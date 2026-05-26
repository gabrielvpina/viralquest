import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from viralquest.biodata import (
    BlastxResult,
    HostViralHit,
    HostViralRecord,
    NucSequence,
    SalmonEntry,
    SalmonQuantReport,
)
from viralquest.salmon_quant import (
    CombinedFastaWriter,
    ConservedFastaLoader,
    QuantsfParser,
    SalmonQuantPipeline,
    TranscriptomeViralAligner,
    _KINGDOM_FILES,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _viral_seq(seq_id: str = "vq1", sequence: str = "ATGC" * 300) -> NucSequence:
    seq = NucSequence(id=seq_id, sequence=sequence)
    seq.is_viral = True
    seq.blastx_nr_hits.append(BlastxResult(
        query_id=seq_id, subject_id="NR_001",
        subject_title=f"polymerase [Influenza A virus]",
        pct_identity=92.0, aln_length=300, mismatches=8,
        gap_opens=0, query_start=1, query_end=300,
        subject_start=1, subject_end=300,
        e_value=1e-90, bit_score=450.0,
    ))
    return seq


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id} description\n{seq}\n")


def _write_quant_sf(path: Path, rows: list[dict]) -> None:
    with open(path, "w") as fh:
        fh.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for r in rows:
            fh.write(
                f"{r['name']}\t{r['length']}\t{r['eff_length']}"
                f"\t{r['tpm']}\t{r['num_reads']}\n"
            )


# ---------------------------------------------------------------------------
# ConservedFastaLoader
# ---------------------------------------------------------------------------

class TestConservedFastaLoader:
    def test_load_all_returns_dict(self, tmp_path):
        fasta = tmp_path / "mammals.fasta"
        fasta.write_text(">NM_001.1 Homo sapiens ACTB\nATGCATGCATGC\n")
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"mammals": fasta}):
            result = ConservedFastaLoader.load_all()
        assert "MAMMALS" in result

    def test_prefix_applied(self, tmp_path):
        fasta = tmp_path / "mammals.fasta"
        fasta.write_text(">NM_001.1 Homo sapiens ACTB\nATGCATGC\n")
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"mammals": fasta}):
            result = ConservedFastaLoader.load_all()
        prefixed_id = result["MAMMALS"][0][0]
        assert prefixed_id == "VQ_CONS_MAMMALS_NM_001.1"

    def test_description_preserved(self, tmp_path):
        fasta = tmp_path / "plants.fasta"
        fasta.write_text(">XM_001.2 Arabidopsis thaliana actin\nATGC\n")
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"plants": fasta}):
            result = ConservedFastaLoader.load_all()
        _, desc, _ = result["PLANTS"][0]
        assert desc == "Arabidopsis thaliana actin"

    def test_sequence_preserved(self, tmp_path):
        fasta = tmp_path / "fish.fasta"
        fasta.write_text(">NM_999.1 Danio rerio\nATGCATGCATGC\n")
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"fish": fasta}):
            result = ConservedFastaLoader.load_all()
        _, _, seq = result["FISH"][0]
        assert seq == "ATGCATGCATGC"

    def test_missing_file_skipped(self, tmp_path):
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"mammals": tmp_path / "does_not_exist.fasta"}):
            result = ConservedFastaLoader.load_all()
        assert "MAMMALS" not in result

    def test_multiple_seqs_parsed(self, tmp_path):
        fasta = tmp_path / "arthropods.fasta"
        fasta.write_text(
            ">PP001.1 Apis mellifera actin\nATGC\n"
            ">PP002.1 Drosophila melanogaster actin\nGCTA\n"
        )
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"arthropods": fasta}):
            result = ConservedFastaLoader.load_all()
        assert len(result["ARTHROPODS"]) == 2

    def test_accession_with_underscore_prefixed_correctly(self, tmp_path):
        # XM_017004835.2 contains an underscore — prefix must not corrupt it
        fasta = tmp_path / "mammals.fasta"
        fasta.write_text(">XM_017004835.2 Homo sapiens POTEF\nATGC\n")
        with patch.dict("viralquest.salmon_quant._KINGDOM_FILES",
                        {"mammals": fasta}):
            result = ConservedFastaLoader.load_all()
        prefixed_id = result["MAMMALS"][0][0]
        assert prefixed_id == "VQ_CONS_MAMMALS_XM_017004835.2"


# ---------------------------------------------------------------------------
# CombinedFastaWriter
# ---------------------------------------------------------------------------

class TestCombinedFastaWriter:
    def _conserved(self) -> dict[str, list[tuple[str, str, str]]]:
        return {
            "MAMMALS": [("VQ_CONS_MAMMALS_NM_001.1", "ACTB", "AAAA")],
            "PLANTS":  [("VQ_CONS_PLANTS_XM_001.1",  "actin", "CCCC")],
        }

    def test_viral_prefix_written(self, tmp_path):
        seq = _viral_seq("myseq", "ATGC" * 10)
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("tx1", "GGGG")])
        CombinedFastaWriter.write(out, tmp_path / "tx.fa", [seq], {})
        text = out.read_text()
        assert ">VQ_VIRAL_myseq" in text

    def test_conserved_prefix_written(self, tmp_path):
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("tx1", "GGGG")])
        CombinedFastaWriter.write(out, tmp_path / "tx.fa", [], self._conserved())
        text = out.read_text()
        assert ">VQ_CONS_MAMMALS_NM_001.1" in text
        assert ">VQ_CONS_PLANTS_XM_001.1" in text

    def test_host_transcriptome_written_unchanged(self, tmp_path):
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("ENST001", "TTTTGGGG")])
        CombinedFastaWriter.write(out, tmp_path / "tx.fa", [], {})
        text = out.read_text()
        assert ">ENST001" in text
        assert "TTTTGGGG" in text

    def test_returns_viral_ids_set(self, tmp_path):
        seq = _viral_seq("seq1")
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("tx1", "AAAA")])
        viral_ids, _ = CombinedFastaWriter.write(out, tmp_path / "tx.fa", [seq], {})
        assert "VQ_VIRAL_seq1" in viral_ids

    def test_returns_conserved_map(self, tmp_path):
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("tx1", "AAAA")])
        _, conserved_map = CombinedFastaWriter.write(
            out, tmp_path / "tx.fa", [], self._conserved()
        )
        assert conserved_map["VQ_CONS_MAMMALS_NM_001.1"] == "MAMMALS"
        assert conserved_map["VQ_CONS_PLANTS_XM_001.1"]  == "PLANTS"

    def test_sequence_order_viral_then_conserved_then_host(self, tmp_path):
        seq = _viral_seq("v1")
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("host1", "GGGG")])
        conserved = {"MAMMALS": [("VQ_CONS_MAMMALS_NM_001.1", "ACTB", "AAAA")]}
        CombinedFastaWriter.write(out, tmp_path / "tx.fa", [seq], conserved)
        text = out.read_text()
        pos_viral     = text.index("VQ_VIRAL_")
        pos_conserved = text.index("VQ_CONS_")
        pos_host      = text.index(">host1")
        assert pos_viral < pos_conserved < pos_host

    def test_empty_viral_and_conserved(self, tmp_path):
        out = tmp_path / "combined.fa"
        _write_fasta(tmp_path / "tx.fa", [("host1", "GGGG")])
        viral_ids, conserved_map = CombinedFastaWriter.write(
            out, tmp_path / "tx.fa", [], {}
        )
        assert viral_ids      == set()
        assert conserved_map  == {}
        assert ">host1" in out.read_text()


# ---------------------------------------------------------------------------
# QuantsfParser
# ---------------------------------------------------------------------------

class TestQuantsfParser:
    def _sf(self, tmp_path: Path, rows: list[dict]) -> Path:
        p = tmp_path / "quant.sf"
        _write_quant_sf(p, rows)
        return p

    def _row(self, name, tpm=100.0, num_reads=50.0) -> dict:
        return dict(name=name, length=1000, eff_length=950.5, tpm=tpm, num_reads=num_reads)

    def test_viral_tagged(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("VQ_VIRAL_seq1", tpm=500.0)])
        entries = QuantsfParser.parse(sf, {})
        assert entries[0].seq_type == "viral"
        assert entries[0].kingdom  == ""

    def test_conserved_tagged(self, tmp_path):
        cmap = {"VQ_CONS_MAMMALS_NM_001.1": "MAMMALS"}
        sf   = self._sf(tmp_path, [self._row("VQ_CONS_MAMMALS_NM_001.1")])
        entries = QuantsfParser.parse(sf, cmap)
        assert entries[0].seq_type == "conserved"
        assert entries[0].kingdom  == "MAMMALS"

    def test_host_tagged(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("ENST00000001234")])
        entries = QuantsfParser.parse(sf, {})
        assert entries[0].seq_type == "host"
        assert entries[0].kingdom  == ""

    def test_tpm_and_num_reads_parsed(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("VQ_VIRAL_seq1", tpm=123.45, num_reads=67.8)])
        e  = QuantsfParser.parse(sf, {})[0]
        assert e.tpm       == pytest.approx(123.45)
        assert e.num_reads == pytest.approx(67.8)

    def test_multiple_entries(self, tmp_path):
        rows = [
            self._row("VQ_VIRAL_seq1"),
            self._row("VQ_CONS_MAMMALS_NM_001.1"),
            self._row("ENST001"),
        ]
        sf      = self._sf(tmp_path, rows)
        entries = QuantsfParser.parse(sf, {"VQ_CONS_MAMMALS_NM_001.1": "MAMMALS"})
        types   = {e.seq_type for e in entries}
        assert types == {"viral", "conserved", "host"}

    def test_nonexistent_file_returns_empty(self, tmp_path):
        entries = QuantsfParser.parse(tmp_path / "missing.sf", {})
        assert entries == []

    def test_conserved_kingdom_fallback_from_name(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("VQ_CONS_ARTHROPODS_PP001.1")])
        # no entry in conserved_map → fallback to parsing name
        entries = QuantsfParser.parse(sf, {})
        assert entries[0].kingdom == "ARTHROPODS"

    def test_returns_salmon_entry_instances(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("VQ_VIRAL_seq1")])
        entries = QuantsfParser.parse(sf, {})
        assert all(isinstance(e, SalmonEntry) for e in entries)

    def test_header_row_skipped(self, tmp_path):
        sf = self._sf(tmp_path, [self._row("VQ_VIRAL_seq1")])
        # there should be exactly 1 entry, not 2 (header must not parse as data)
        assert len(QuantsfParser.parse(sf, {})) == 1


# ---------------------------------------------------------------------------
# TranscriptomeViralAligner
# ---------------------------------------------------------------------------

class TestTranscriptomeViralAligner:
    def _aligner(self) -> TranscriptomeViralAligner:
        return TranscriptomeViralAligner()

    def test_no_viral_seqs_returns_empty(self, tmp_path):
        aligner = self._aligner()
        hits = aligner.align(tmp_path / "tx.fa", [], tmp_path)
        assert hits == []

    def test_empty_tsv_returns_empty(self, tmp_path):
        aligner = self._aligner()
        (tmp_path / "host_vs_viral.tsv").write_text("")
        hits = aligner._parse_tsv(tmp_path / "host_vs_viral.tsv")
        assert hits == []

    def test_parses_valid_tsv_row(self, tmp_path):
        tsv = tmp_path / "host_vs_viral.tsv"
        tsv.write_text("tx1\tseq1\t85.0\t500\t600\t80\t1e-20\t150.0\n")
        hits = self._aligner()._parse_tsv(tsv)
        assert len(hits) == 1
        h = hits[0]
        assert h.host_transcript_id == "tx1"
        assert h.viral_seq_id       == "seq1"
        assert h.pident             == 85.0
        assert h.qcovhsp            == 80
        assert h.bit_score          == 150.0

    def test_malformed_row_skipped(self, tmp_path):
        tsv = tmp_path / "host_vs_viral.tsv"
        tsv.write_text("only\ttwo\tcolumns\n" + "tx1\tseq1\t85.0\t500\t600\t80\t1e-20\t150.0\n")
        hits = self._aligner()._parse_tsv(tsv)
        assert len(hits) == 1

    def test_qcov_filter_applied(self, tmp_path):
        tsv = tmp_path / "host_vs_viral.tsv"
        # qcovhsp = 5, below default min_qcov = 20
        tsv.write_text("tx1\tseq1\t85.0\t50\t1000\t5\t1e-20\t150.0\n")
        aligner = TranscriptomeViralAligner(min_qcov=20)
        hits    = aligner._parse_tsv(tsv)
        hits    = [h for h in hits if h.qcovhsp >= aligner.min_qcov]
        assert hits == []

    def test_write_viral_fasta(self, tmp_path):
        seqs = [_viral_seq("v1", "ATGC"), _viral_seq("v2", "GGGG")]
        out  = tmp_path / "viral.fa"
        TranscriptomeViralAligner._write_fasta(seqs, out)
        text = out.read_text()
        assert ">v1\nATGC" in text
        assert ">v2\nGGGG" in text

    def test_returns_host_viral_hit_instances(self, tmp_path):
        tsv = tmp_path / "host_vs_viral.tsv"
        tsv.write_text("tx1\tseq1\t85.0\t500\t600\t80\t1e-20\t150.0\n")
        hits = self._aligner()._parse_tsv(tsv)
        assert all(isinstance(h, HostViralHit) for h in hits)


# ---------------------------------------------------------------------------
# SalmonQuantPipeline._assemble  (unit test of the pure logic)
# ---------------------------------------------------------------------------

class TestAssemble:
    def _entry(self, name: str, seq_type: str, kingdom: str = "", tpm: float = 100.0) -> SalmonEntry:
        return SalmonEntry(
            name=name, length=1000, eff_length=950.0,
            tpm=tpm, num_reads=50.0, seq_type=seq_type, kingdom=kingdom,
        )

    def _hit(self, host_id: str, viral_id: str, bit_score: float = 200.0) -> HostViralHit:
        return HostViralHit(
            host_transcript_id=host_id, viral_seq_id=viral_id,
            pident=85.0, qcovhsp=70, evalue=1e-20, bit_score=bit_score,
        )

    def test_viral_quant_separated(self):
        entries = [
            self._entry("VQ_VIRAL_seq1", "viral"),
            self._entry("ENST001",        "host"),
        ]
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, [])
        assert len(report.viral_quant) == 1
        assert report.viral_quant[0].name == "VQ_VIRAL_seq1"

    def test_conserved_quant_separated(self):
        entries = [
            self._entry("VQ_CONS_MAMMALS_NM_001.1", "conserved", "MAMMALS"),
            self._entry("VQ_CONS_PLANTS_XM_001.1",  "conserved", "PLANTS"),
        ]
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, [])
        assert len(report.conserved_quant) == 2

    def test_conserved_sorted_by_kingdom_then_tpm(self):
        entries = [
            self._entry("VQ_CONS_PLANTS_XM_001.1",  "conserved", "PLANTS",  tpm=50.0),
            self._entry("VQ_CONS_MAMMALS_NM_002.1",  "conserved", "MAMMALS", tpm=200.0),
            self._entry("VQ_CONS_MAMMALS_NM_001.1",  "conserved", "MAMMALS", tpm=500.0),
        ]
        report  = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, [])
        kingdoms = [e.kingdom for e in report.conserved_quant]
        assert kingdoms[:2] == ["MAMMALS", "MAMMALS"]
        assert kingdoms[2]  == "PLANTS"
        # within MAMMALS: higher TPM first
        mammals = [e for e in report.conserved_quant if e.kingdom == "MAMMALS"]
        assert mammals[0].tpm > mammals[1].tpm

    def test_host_viral_record_built(self):
        entries = [self._entry("ENST001", "host")]
        hits    = [self._hit("ENST001", "VQ_VIRAL_seq1")]
        report  = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, hits)
        assert len(report.host_viral_hits) == 1
        rec = report.host_viral_hits[0]
        assert rec.transcript.name == "ENST001"
        assert len(rec.blast_hits) == 1

    def test_host_viral_hits_sorted_by_bitscore(self):
        entries = [
            self._entry("ENST001", "host"),
            self._entry("ENST002", "host"),
        ]
        hits = [self._hit("ENST001", "v1", bit_score=100.0),
                self._hit("ENST002", "v1", bit_score=500.0)]
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, hits)
        assert report.host_viral_hits[0].transcript.name == "ENST002"

    def test_blast_hits_within_record_sorted_by_bitscore(self):
        entries = [self._entry("ENST001", "host")]
        hits    = [
            self._hit("ENST001", "v1", bit_score=100.0),
            self._hit("ENST001", "v2", bit_score=400.0),
        ]
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, hits)
        rec    = report.host_viral_hits[0]
        assert rec.blast_hits[0].bit_score == 400.0

    def test_host_not_in_quant_gets_zero_tpm(self):
        hits   = [self._hit("ENST_NOT_IN_QUANT", "v1")]
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, [], hits)
        assert report.host_viral_hits[0].transcript.tpm       == 0.0
        assert report.host_viral_hits[0].transcript.num_reads == 0.0

    def test_total_reads_summed(self):
        entries = [
            self._entry("VQ_VIRAL_seq1", "viral", tpm=100.0),
            self._entry("ENST001",        "host",  tpm=50.0),
        ]
        # num_reads is 50.0 for both entries (default in _entry)
        report = SalmonQuantPipeline._assemble(["reads.fq"], 45.0, entries, [])
        assert report.total_reads == 100

    def test_mapping_rate_preserved(self):
        report = SalmonQuantPipeline._assemble(["r.fq"], 62.5, [], [])
        assert report.mapping_rate == 62.5

    def test_reads_list_preserved(self):
        report = SalmonQuantPipeline._assemble(["r1.fq", "r2.fq"], 0.0, [], [])
        assert report.reads == ["r1.fq", "r2.fq"]

    def test_returns_salmon_quant_report(self):
        report = SalmonQuantPipeline._assemble(["reads.fq"], 0.0, [], [])
        assert isinstance(report, SalmonQuantReport)


# ---------------------------------------------------------------------------
# Exporter integration — salmon_quant section
# ---------------------------------------------------------------------------

class TestExporterSalmonSection:
    def _report(self) -> SalmonQuantReport:
        viral_entry = SalmonEntry(
            name="VQ_VIRAL_seq1", length=1200, eff_length=1150.0,
            tpm=250.0, num_reads=80.0, seq_type="viral", kingdom="",
        )
        conserved_entry = SalmonEntry(
            name="VQ_CONS_MAMMALS_NM_001.1", length=1000, eff_length=950.0,
            tpm=8400.0, num_reads=2000.0, seq_type="conserved", kingdom="MAMMALS",
        )
        host_entry = SalmonEntry(
            name="ENST001", length=900, eff_length=850.0,
            tpm=30.0, num_reads=10.0, seq_type="host", kingdom="",
        )
        hit = HostViralHit(
            host_transcript_id="ENST001", viral_seq_id="seq1",
            pident=82.0, qcovhsp=65, evalue=1e-15, bit_score=180.0,
        )
        record = HostViralRecord(transcript=host_entry, blast_hits=[hit])
        return SalmonQuantReport(
            reads=["sample.fastq"], mapping_rate=52.3, total_reads=5000,
            viral_quant=[viral_entry],
            conserved_quant=[conserved_entry],
            host_viral_hits=[record],
        )

    def _make_seq(self) -> NucSequence:
        seq = NucSequence(id="seq1", sequence="ATGC" * 300)
        seq.is_viral = True
        seq.blastx_nr_hits.append(BlastxResult(
            query_id="seq1", subject_id="NR_001",
            subject_title="polymerase [Influenza A virus]",
            pct_identity=92.0, aln_length=300, mismatches=8,
            gap_opens=0, query_start=1, query_end=300,
            subject_start=1, subject_end=300,
            e_value=1e-90, bit_score=450.0,
        ))
        return seq

    def test_salmon_quant_key_present_when_report_given(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [], salmon_report=self._report())
        assert "salmon_quant" in report

    def test_salmon_quant_key_absent_when_no_report(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [])
        assert "salmon_quant" not in report

    def test_salmon_quant_structure(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [], salmon_report=self._report())
        sq     = report["salmon_quant"]
        assert "reads"           in sq
        assert "mapping_rate"    in sq
        assert "total_reads"     in sq
        assert "viral_quant"     in sq
        assert "conserved_quant" in sq
        assert "host_viral_hits" in sq

    def test_salmon_quant_json_serialisable(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [], salmon_report=self._report())
        # must not raise
        json.dumps(report)

    def test_salmon_quant_written_to_file(self, tmp_path):
        from viralquest.exporter import ReportExporter
        out    = tmp_path / "report.json"
        exp    = ReportExporter()
        exp.export([self._make_seq()], [], output_path=out, salmon_report=self._report())
        data   = json.loads(out.read_text())
        assert "salmon_quant" in data
        assert data["salmon_quant"]["mapping_rate"] == pytest.approx(52.3)

    def test_viral_quant_entries_present(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [], salmon_report=self._report())
        vq     = report["salmon_quant"]["viral_quant"]
        assert len(vq) == 1
        assert vq[0]["tpm"] == pytest.approx(250.0)

    def test_host_viral_hits_entries_present(self):
        from viralquest.exporter import ReportExporter
        exp    = ReportExporter()
        report = exp.export([self._make_seq()], [], salmon_report=self._report())
        hits   = report["salmon_quant"]["host_viral_hits"]
        assert len(hits) == 1
        assert hits[0]["transcript"]["name"] == "ENST001"
        assert hits[0]["blast_hits"][0]["pident"] == pytest.approx(82.0)
