import csv
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from viralquest.biodata import BlastnResult, BlastxResult, NucSequence
from viralquest.blastn import (
    BlastnMode,
    BlastnOutputParser,
    BlastnResultAttacher,
    BlastnRunner,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_seq(
    seq_id: str = "seq1",
    sequence: str = "ATGC" * 300,
    is_viral: bool = False,
    with_nr_hit: bool = False,
) -> NucSequence:
    seq = NucSequence(id=seq_id, sequence=sequence)
    seq.is_viral = is_viral
    if with_nr_hit:
        seq.blastx_nr_hits.append(
            BlastxResult(
                query_id=seq_id, subject_id="NR_001",
                subject_title=f"polymerase [Influenza A virus]",
                pct_identity=92.0, aln_length=300, mismatches=8,
                gap_opens=0, query_start=1, query_end=300,
                subject_start=1, subject_end=300,
                e_value=1e-90, bit_score=450.0,
            )
        )
    return seq


def _tsv_row(
    qseqid: str  = "seq1",
    sseqid: str  = "NC_002017",
    stitle: str  = "Influenza A virus segment 1 RNA polymerase PB2",
    pident: str  = "92.5",
    length: str  = "1200",
    qlen: str    = "1200",
    slen: str    = "13600",
    qcovhsp: str = "80",
    evalue: str  = "1e-60",
    bitscore: str = "200.0",
) -> str:
    return "\t".join([qseqid, sseqid, stitle, pident, length, qlen, slen,
                      qcovhsp, evalue, bitscore])


def _write_tsv(tmp_path: Path, rows: list[str]) -> Path:
    p = tmp_path / "out.tsv"
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# BlastnMode
# ---------------------------------------------------------------------------

class TestBlastnMode:
    def test_values(self):
        assert BlastnMode.LOCAL.value  == "local"
        assert BlastnMode.ONLINE.value == "online"


# ---------------------------------------------------------------------------
# BlastnRunner  — init
# ---------------------------------------------------------------------------

class TestBlastnRunnerInit:
    def test_local_without_db_path_raises(self):
        with pytest.raises(ValueError, match="db_path"):
            BlastnRunner(mode=BlastnMode.LOCAL)

    def test_local_with_db_path_succeeds(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path=str(tmp_path / "db"))
        assert runner.db_path == str(tmp_path / "db")

    def test_online_without_db_path_succeeds(self):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        assert runner.mode == BlastnMode.ONLINE

    def test_default_params(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/db")
        assert runner.blastn_bin      == "blastn"
        assert runner.threads         == 2
        assert runner.e_value         == 1e-5
        assert runner.max_target_seqs == 5
        assert runner.batch_size      == 200
        assert runner.request_delay   == 0.4

    def test_outdir_created(self, tmp_path):
        outdir = tmp_path / "blast_out"
        BlastnRunner(mode=BlastnMode.LOCAL, db_path="/db", outdir=str(outdir))
        assert outdir.exists()


# ---------------------------------------------------------------------------
# BlastnRunner.run  — filtering logic
# ---------------------------------------------------------------------------

class TestBlastnRunnerRun:
    def _runner(self, tmp_path) -> BlastnRunner:
        return BlastnRunner(
            mode=BlastnMode.LOCAL,
            db_path="/fake/db",
            outdir=str(tmp_path),
        )

    def test_skips_non_viral(self, tmp_path):
        seq = make_seq(is_viral=False, with_nr_hit=True)
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local", return_value=[]) as mock:
            runner.run([seq])
        mock.assert_not_called()

    def test_skips_viral_without_nr_hits(self, tmp_path):
        seq = make_seq(is_viral=True, with_nr_hit=False)
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local", return_value=[]) as mock:
            runner.run([seq])
        mock.assert_not_called()

    def test_passes_viral_with_nr_hits(self, tmp_path):
        seq = make_seq(is_viral=True, with_nr_hit=True)
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local", return_value=[]) as mock:
            runner.run([seq])
        mock.assert_called_once()
        called_with = mock.call_args[0][0]
        assert len(called_with) == 1
        assert called_with[0].id == "seq1"

    def test_empty_input_returns_empty(self, tmp_path):
        runner = self._runner(tmp_path)
        assert runner.run([]) == []

    def test_mixed_seqs_only_qualified_passed(self, tmp_path):
        viral_confirmed = make_seq("v1", is_viral=True, with_nr_hit=True)
        viral_no_nr     = make_seq("v2", is_viral=True,  with_nr_hit=False)
        non_viral       = make_seq("v3", is_viral=False, with_nr_hit=True)
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local", return_value=[]) as mock:
            runner.run([viral_confirmed, viral_no_nr, non_viral])
        passed = mock.call_args[0][0]
        assert len(passed) == 1
        assert passed[0].id == "v1"

    def test_routes_online_mode(self, tmp_path):
        seq = make_seq(is_viral=True, with_nr_hit=True)
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path))
        with patch.object(runner, "_run_online", return_value=[]) as mock:
            runner.run([seq])
        mock.assert_called_once()


# ---------------------------------------------------------------------------
# BlastnRunner._chunk  (static)
# ---------------------------------------------------------------------------

class TestChunk:
    def _runner(self) -> BlastnRunner:
        return BlastnRunner(mode=BlastnMode.ONLINE)

    def test_exact_multiple(self):
        seqs = [object()] * 6
        chunks = BlastnRunner._chunk(seqs, 3)
        assert len(chunks) == 2
        assert all(len(c) == 3 for c in chunks)

    def test_remainder(self):
        seqs = [object()] * 7
        chunks = BlastnRunner._chunk(seqs, 3)
        assert len(chunks) == 3
        assert len(chunks[-1]) == 1

    def test_single_chunk(self):
        seqs = [object()] * 3
        chunks = BlastnRunner._chunk(seqs, 10)
        assert len(chunks) == 1

    def test_empty(self):
        assert BlastnRunner._chunk([], 5) == []


# ---------------------------------------------------------------------------
# BlastnRunner._write_batch_fasta
# ---------------------------------------------------------------------------

class TestWriteBatchFasta:
    def test_fasta_format(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        seq1 = make_seq("seq1", "ATGCATGC")
        seq2 = make_seq("seq2", "TTTTCCCC")
        out = tmp_path / "batch.fa"
        runner._write_batch_fasta([seq1, seq2], out)
        text = out.read_text()
        assert ">seq1\nATGCATGC\n" in text
        assert ">seq2\nTTTTCCCC\n" in text

    def test_empty_list_creates_empty_file(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        out = tmp_path / "empty.fa"
        runner._write_batch_fasta([], out)
        assert out.read_text() == ""


# ---------------------------------------------------------------------------
# BlastnOutputParser.parse_tsv
# ---------------------------------------------------------------------------

class TestParseTsv:
    def test_parses_valid_row(self, tmp_path):
        tsv = _write_tsv(tmp_path, [_tsv_row()])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 1
        h = hits[0]
        assert h.qseqid    == "seq1"
        assert h.pident    == 92.5
        assert h.qlen      == 1200
        assert h.slen      == 13600
        assert h.qcovhsp   == 80
        assert h.evalue    == pytest.approx(1e-60)
        assert h.bit_score == 200.0
        assert "Influenza" in h.stitle

    def test_parses_multiple_rows(self, tmp_path):
        rows = [_tsv_row("seq1"), _tsv_row("seq2")]
        tsv  = _write_tsv(tmp_path, rows)
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 2
        assert {h.qseqid for h in hits} == {"seq1", "seq2"}

    def test_nonexistent_file_returns_empty(self, tmp_path):
        hits = BlastnOutputParser.parse_tsv(tmp_path / "missing.tsv")
        assert hits == []

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "empty.tsv"
        p.write_text("")
        assert BlastnOutputParser.parse_tsv(p) == []

    def test_malformed_row_skipped(self, tmp_path):
        rows = ["only\ttwo\tcolumns", _tsv_row()]
        tsv  = _write_tsv(tmp_path, rows)
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 1

    def test_qcovhsp_float_truncated_to_int(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row(qcovhsp="79.9")])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert hits[0].qcovhsp == 79

    def test_scientific_evalue(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row(evalue="1.2e-100")])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert hits[0].evalue == pytest.approx(1.2e-100)

    def test_returns_blastn_result_instances(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row()])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert all(isinstance(h, BlastnResult) for h in hits)


# ---------------------------------------------------------------------------
# BlastnOutputParser.parse_xml_record
# ---------------------------------------------------------------------------

def _make_hsp(bits=200.0, identities=90, align_length=100,
              query_start=1, query_end=100, expect=1e-50):
    hsp = MagicMock()
    hsp.bits         = bits
    hsp.identities   = identities
    hsp.align_length = align_length
    hsp.query_start  = query_start
    hsp.query_end    = query_end
    hsp.expect       = expect
    return hsp


def _make_alignment(title="NC_002017 Influenza PB2", hit_def="Influenza PB2",
                    length=13600, hsps=None):
    aln = MagicMock()
    aln.title   = title
    aln.hit_def = hit_def
    aln.length  = length
    aln.hsps    = hsps if hsps is not None else [_make_hsp()]
    return aln


def _make_blast_record(query_length=1200, alignments=None):
    rec = MagicMock()
    rec.query_length = query_length
    rec.alignments   = alignments if alignments is not None else [_make_alignment()]
    return rec


class TestParseXmlRecord:
    def test_basic_parse(self):
        rec  = _make_blast_record()
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert len(hits) == 1
        h = hits[0]
        assert h.qseqid    == "seq1"
        assert h.qlen      == 1200
        assert h.slen      == 13600
        assert h.bit_score == 200.0
        assert h.evalue    == pytest.approx(1e-50)
        assert h.stitle    == "Influenza PB2"

    def test_pident_computed_from_hsp(self):
        hsp = _make_hsp(identities=95, align_length=100)
        rec = _make_blast_record(alignments=[_make_alignment(hsps=[hsp])])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].pident == 95.0

    def test_qcovhsp_computed(self):
        # query_start=1, query_end=600, qlen=1200 → span=600, cov=50%
        hsp = _make_hsp(query_start=1, query_end=600)
        rec = _make_blast_record(query_length=1200, alignments=[_make_alignment(hsps=[hsp])])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].qcovhsp == 50

    def test_best_hsp_selected_by_bitscore(self):
        hsp_low  = _make_hsp(bits=100.0, identities=60)
        hsp_high = _make_hsp(bits=400.0, identities=95)
        aln = _make_alignment(hsps=[hsp_low, hsp_high])
        hits = BlastnOutputParser.parse_xml_record(_make_blast_record(alignments=[aln]), "seq1")
        assert hits[0].bit_score == 400.0
        assert hits[0].pident    == 95.0

    def test_multiple_alignments(self):
        rec = _make_blast_record(alignments=[_make_alignment(), _make_alignment(hit_def="Other virus")])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert len(hits) == 2

    def test_empty_alignments_returns_empty(self):
        rec = _make_blast_record(alignments=[])
        assert BlastnOutputParser.parse_xml_record(rec, "seq1") == []

    def test_alignment_without_hsps_skipped(self):
        aln = _make_alignment(hsps=[])
        rec = _make_blast_record(alignments=[aln])
        assert BlastnOutputParser.parse_xml_record(rec, "seq1") == []

    def test_hit_def_none_falls_back_to_title(self):
        hsp = _make_hsp()
        aln = _make_alignment(title=">NC_002017 Fallback title", hit_def=None, hsps=[hsp])
        aln.hit_def = None
        rec = _make_blast_record(alignments=[aln])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert "Fallback title" in hits[0].stitle

    def test_zero_query_length_gives_zero_coverage(self):
        rec = _make_blast_record(query_length=0)
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].qcovhsp == 0

    def test_returns_blastn_result_instances(self):
        rec  = _make_blast_record()
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert all(isinstance(h, BlastnResult) for h in hits)


# ---------------------------------------------------------------------------
# BlastnOutputParser._accession_from_alignment  (online Hit_id fallback)
# ---------------------------------------------------------------------------

def _aln(hit_id=None, accession="__unset__", hsps=None):
    """Alignment stub. accession is omitted from the object when '__unset__',
    mirroring Biopython alignments parsed from online qblast XML (no Hit_accession)."""
    from types import SimpleNamespace
    ns = SimpleNamespace(title="t", hit_def="Some virus", length=1000,
                         hsps=hsps if hsps is not None else [_make_hsp()])
    if hit_id is not None:
        ns.hit_id = hit_id
    if accession != "__unset__":
        ns.accession = accession
    return ns


class TestAccessionResolution:
    _resolve = staticmethod(BlastnOutputParser._accession_from_alignment)

    def test_explicit_hit_accession_preferred(self):
        # Local BLAST+/XML v2 sets alignment.accession — it wins over Hit_id.
        aln = _aln(hit_id="gi|123|ref|X99999.9|", accession="NC_045512")
        assert self._resolve(aln) == "NC_045512"

    def test_ref_tag_in_pipe_hit_id(self):
        aln = _aln(hit_id="gi|1798174254|ref|NC_045512.2|")
        assert self._resolve(aln) == "NC_045512.2"

    def test_gb_tag_in_pipe_hit_id(self):
        aln = _aln(hit_id="gi|555|gb|MN908947.3|")
        assert self._resolve(aln) == "MN908947.3"

    def test_plain_accession_hit_id_no_pipes(self):
        aln = _aln(hit_id="NC_045512.2")
        assert self._resolve(aln) == "NC_045512.2"

    def test_no_tag_falls_back_to_last_token(self):
        # Pipe-delimited but no recognised db tag → last token.
        aln = _aln(hit_id="12345|NC_045512.2")
        assert self._resolve(aln) == "NC_045512.2"

    def test_no_accession_and_no_hit_id_returns_none(self):
        aln = _aln(hit_id=None)
        assert self._resolve(aln) is None

    def test_parse_xml_record_uses_hit_id_fallback(self):
        # Online path: no Hit_accession → accession pulled from Hit_id.
        aln  = _aln(hit_id="gi|1798174254|ref|NC_045512.2|")
        rec  = _make_blast_record(alignments=[aln])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].accession == "NC_045512.2"


# ---------------------------------------------------------------------------
# BlastnResultAttacher.attach
# ---------------------------------------------------------------------------

class TestBlastnResultAttacher:
    def _hit(self, qseqid: str, bit_score: float) -> BlastnResult:
        return BlastnResult(
            qseqid=qseqid, qlen=1200, slen=13600,
            qcovhsp=80, pident=92.5, evalue=1e-50,
            bit_score=bit_score,
            stitle="Influenza A",
        )

    def test_hits_attached_to_matching_seq(self):
        seq  = make_seq("seq1")
        hit  = self._hit("seq1", 200.0)
        BlastnResultAttacher.attach([hit], [seq])
        assert len(seq.blastn_hits) == 1

    def test_hits_sorted_by_bitscore_descending(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", b) for b in [100.0, 300.0, 200.0]]
        BlastnResultAttacher.attach(hits, [seq])
        scores = [h.bit_score for h in seq.blastn_hits]
        assert scores == sorted(scores, reverse=True)

    def test_max_hits_capped(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i)) for i in range(10)]
        BlastnResultAttacher.attach(hits, [seq], max_hits_per_query=3)
        assert len(seq.blastn_hits) == 3

    def test_top_hits_kept_after_capping(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i * 10)) for i in range(5)]
        BlastnResultAttacher.attach(hits, [seq], max_hits_per_query=2)
        scores = [h.bit_score for h in seq.blastn_hits]
        assert scores[0] == 40.0
        assert scores[1] == 30.0

    def test_unknown_query_id_ignored(self):
        seq = make_seq("seq1")
        hit = self._hit("seq_unknown", 200.0)
        attached = BlastnResultAttacher.attach([hit], [seq])
        assert len(seq.blastn_hits) == 0
        assert attached == 0

    def test_returns_total_attached_count(self):
        seq1 = make_seq("seq1")
        seq2 = make_seq("seq2")
        hits = [self._hit("seq1", 100.0), self._hit("seq2", 200.0), self._hit("seq2", 300.0)]
        count = BlastnResultAttacher.attach(hits, [seq1, seq2])
        assert count == 3

    def test_empty_hits_returns_zero(self):
        seq = make_seq("seq1")
        assert BlastnResultAttacher.attach([], [seq]) == 0
        assert seq.blastn_hits == []

    def test_empty_seqs_returns_zero(self):
        hit = self._hit("seq1", 200.0)
        assert BlastnResultAttacher.attach([hit], []) == 0

    def test_multiple_seqs_independent(self):
        seq1 = make_seq("seq1")
        seq2 = make_seq("seq2")
        hits = [self._hit("seq1", 100.0), self._hit("seq2", 200.0)]
        BlastnResultAttacher.attach(hits, [seq1, seq2])
        assert len(seq1.blastn_hits) == 1
        assert len(seq2.blastn_hits) == 1
        assert seq1.blastn_hits[0].bit_score == 100.0
        assert seq2.blastn_hits[0].bit_score == 200.0

    def test_default_max_is_five(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i)) for i in range(10)]
        BlastnResultAttacher.attach(hits, [seq])
        assert len(seq.blastn_hits) == 5
