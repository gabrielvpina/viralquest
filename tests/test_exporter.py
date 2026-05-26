import json
import uuid
from pathlib import Path

import pytest

from viralquest.biodata import (
    BlastnResult,
    BlastxResult,
    ClusterMember,
    InputFasta,
    LlmOutput,
    NucSequence,
    Orf,
    Taxonomy,
    ViralCluster,
    ViralFamilyInfo,
)
from viralquest.exporter import ReportExporter, _to_serializable


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _taxonomy(family: str = "Orthomyxoviridae") -> Taxonomy:
    return Taxonomy(
        tax_id=11520,
        scientific_name="Influenza A virus",
        no_rank=None,
        clade=None,
        kingdom=None,
        phylum=None,
        class_=None,
        order="Articulavirales",
        family=family,
        subfamily=None,
        genus="Alphainfluenzavirus",
        species="Influenza A virus",
        genome="ssRNA(-)",
    )


def _llm_output(seq_id: str = "seq1") -> LlmOutput:
    return LlmOutput(
        seq_id=seq_id, model="test-model", mode="low",
        vq_score=85, classification="viral-known",
        analysis="Confirmed influenza.", blastn_species="Influenza A virus",
    )


def _viral_family_info() -> ViralFamilyInfo:
    return ViralFamilyInfo(
        source="ViralZone", type="Family", name="Orthomyxoviridae",
        info_high="Full description.", info_low="Short description.",
    )


def _blastx_nr_hit(seq_id: str = "seq1") -> BlastxResult:
    return BlastxResult(
        query_id=seq_id, subject_id="NR_001",
        subject_title="polymerase [Influenza A virus]",
        pct_identity=92.0, aln_length=300, mismatches=8,
        gap_opens=0, query_start=1, query_end=300,
        subject_start=1, subject_end=300,
        e_value=1e-90, bit_score=450.0,
    )


def _blastn_hit(seq_id: str = "seq1") -> BlastnResult:
    return BlastnResult(
        qseqid=seq_id, qlen=1200, slen=13600,
        qcovhsp=80, pident=92.5, evalue=1e-60,
        bit_score=200.0, stitle="Influenza A virus PB2",
    )


def _make_confirmed_seq(seq_id: str = "seq1") -> NucSequence:
    seq = NucSequence(id=seq_id, sequence="ATGC" * 300)
    seq.is_viral  = True
    seq.taxonomy  = _taxonomy()
    seq.llm_output = _llm_output(seq_id)
    seq.blastx_nr_hits.append(_blastx_nr_hit(seq_id))
    seq.blastn_hits.append(_blastn_hit(seq_id))
    return seq


def _make_cluster(rep_id: str = "seq1") -> ViralCluster:
    member = ClusterMember(
        seq_id=rep_id, length=1200, is_representative=True,
        identity=100.0, query_coverage=100.0, aln_start=1, aln_end=1200,
    )
    cluster = ViralCluster(
        cluster_id="VQ_CLU_001", species="Influenza A virus",
        representative_id=rep_id, members=[member],
    )
    return cluster


# ---------------------------------------------------------------------------
# _to_serializable
# ---------------------------------------------------------------------------

class TestToSerializable:
    def test_uuid_becomes_string(self):
        u = uuid.uuid4()
        assert _to_serializable(u) == str(u)

    def test_path_becomes_string(self):
        p = Path("/some/path/file.json")
        assert _to_serializable(p) == "/some/path/file.json"

    def test_list_recursed(self):
        u = uuid.uuid4()
        result = _to_serializable([u, 42])
        assert result == [str(u), 42]

    def test_dict_recursed(self):
        u = uuid.uuid4()
        result = _to_serializable({"key": u, "num": 1})
        assert result == {"key": str(u), "num": 1}

    def test_dataclass_converted_to_dict(self):
        hit = _blastn_hit()
        result = _to_serializable(hit)
        assert isinstance(result, dict)
        assert result["qseqid"]    == "seq1"
        assert result["pident"]    == 92.5
        assert result["bit_score"] == 200.0

    def test_dataclass_fields_only_no_properties(self):
        cluster = _make_cluster()
        result = _to_serializable(cluster)
        # ViralCluster.size is a @property — it must NOT appear here
        assert "size" not in result

    def test_none_preserved(self):
        assert _to_serializable(None) is None

    def test_primitives_passthrough(self):
        assert _to_serializable(42)      == 42
        assert _to_serializable(3.14)    == 3.14
        assert _to_serializable("hello") == "hello"
        assert _to_serializable(True)    is True

    def test_nested_dataclass(self):
        seq  = NucSequence(id="seq1", sequence="ATGC")
        hit  = _blastn_hit()
        seq.blastn_hits.append(hit)
        result = _to_serializable(seq)
        assert isinstance(result, dict)
        assert isinstance(result["blastn_hits"], list)
        assert result["blastn_hits"][0]["qseqid"] == "seq1"

    def test_type_object_not_converted(self):
        # dataclasses.is_dataclass(NucSequence) is True, but it's a class — must skip
        result = _to_serializable(NucSequence)
        assert result is NucSequence


# ---------------------------------------------------------------------------
# ReportExporter — meta
# ---------------------------------------------------------------------------

class TestBuildMeta:
    def test_meta_has_required_keys(self):
        exp    = ReportExporter()
        report = exp.export([], [])
        meta   = report["meta"]
        assert "viralquest_version" in meta
        assert "timestamp"          in meta

    def test_meta_version_embedded(self):
        exp    = ReportExporter()
        report = exp.export([], [], version="1.2.3")
        assert report["meta"]["viralquest_version"] == "1.2.3"

    def test_meta_no_input_file_when_none(self):
        exp    = ReportExporter()
        report = exp.export([], [])
        assert "input_file" not in report["meta"]

    def test_meta_input_file_when_provided(self):
        fasta = InputFasta(name="sample.fa", num_seqs=10, path="/data/sample.fa", size=4096)
        exp   = ReportExporter()
        report = exp.export([], [], input_fasta=fasta)
        assert report["meta"]["input_file"]["name"]     == "sample.fa"
        assert report["meta"]["input_file"]["num_seqs"] == 10
        assert report["meta"]["input_file"]["size_bytes"] == 4096


# ---------------------------------------------------------------------------
# ReportExporter — filtering
# ---------------------------------------------------------------------------

class TestExporterFiltering:
    def test_confirmed_seq_included(self):
        seq    = _make_confirmed_seq()
        exp    = ReportExporter()
        report = exp.export([seq], [])
        assert len(report["sequences"]) == 1
        assert report["sequences"][0]["id"] == "seq1"

    def test_non_viral_excluded(self):
        seq = NucSequence(id="s1", sequence="ATGC" * 100)
        # is_viral=False by default
        exp    = ReportExporter()
        report = exp.export([seq], [])
        assert len(report["sequences"]) == 0

    def test_viral_without_nr_hits_excluded(self):
        seq = NucSequence(id="s1", sequence="ATGC" * 100)
        seq.is_viral = True   # no blastx_nr_hits → excluded
        exp = ReportExporter()
        assert len(exp.export([seq], [])["sequences"]) == 0

    def test_force_exports_all(self):
        non_viral = NucSequence(id="s1", sequence="ATGC" * 100)
        confirmed = _make_confirmed_seq("s2")
        exp = ReportExporter(force=True)
        report = exp.export([non_viral, confirmed], [])
        assert len(report["sequences"]) == 2

    def test_force_exports_all_clusters(self):
        non_viral = NucSequence(id="s1", sequence="ATGC" * 100)
        cluster   = _make_cluster("s1")
        exp = ReportExporter(force=True)
        report = exp.export([non_viral], [cluster])
        assert len(report["clusters"]) == 1

    def test_cluster_filtered_to_confirmed_representatives(self):
        confirmed = _make_confirmed_seq("s1")
        dropped   = NucSequence(id="s2", sequence="TTTT" * 100)   # not viral
        c_confirmed = _make_cluster("s1")
        c_dropped   = _make_cluster("s2")
        exp    = ReportExporter()
        report = exp.export([confirmed, dropped], [c_confirmed, c_dropped])
        assert len(report["sequences"]) == 1
        assert len(report["clusters"])  == 1
        assert report["clusters"][0]["representative_id"] == "s1"


# ---------------------------------------------------------------------------
# ReportExporter — sequence fields
# ---------------------------------------------------------------------------

class TestSeqFields:
    def _report_seq(self, **kwargs) -> dict:
        seq = _make_confirmed_seq()
        exp = ReportExporter(**kwargs)
        return exp.export([seq], [])["sequences"][0]

    def test_basic_fields_present(self):
        s = self._report_seq()
        for key in ("id", "uid", "sequence", "length", "gc_content", "n_count",
                    "is_viral", "cluster_id", "orfs", "blastx_hits",
                    "blastx_nr_hits", "blastn_hits", "taxonomy"):
            assert key in s, f"missing key: {key}"

    def test_uid_is_string(self):
        s = self._report_seq()
        assert isinstance(s["uid"], str)
        # must parse as UUID
        uuid.UUID(s["uid"])

    def test_viral_family_info_excluded(self):
        seq = _make_confirmed_seq()
        seq.viral_family_info = _viral_family_info()
        exp = ReportExporter()
        s   = exp.export([seq], [])["sequences"][0]
        assert "viral_family_info" not in s

    def test_llm_output_included_by_default(self):
        s = self._report_seq(include_llm=True)
        assert "llm_output" in s
        assert s["llm_output"] is not None

    def test_llm_output_excluded_when_flag_false(self):
        s = self._report_seq(include_llm=False)
        assert "llm_output" not in s

    def test_llm_output_none_when_absent(self):
        seq = _make_confirmed_seq()
        seq.llm_output = None
        exp = ReportExporter(include_llm=True)
        s   = exp.export([seq], [])["sequences"][0]
        assert s["llm_output"] is None

    def test_blastn_hits_serialised(self):
        s = self._report_seq()
        assert len(s["blastn_hits"]) == 1
        h = s["blastn_hits"][0]
        assert h["stitle"] == "Influenza A virus PB2"
        assert h["pident"] == 92.5

    def test_blastx_nr_hits_serialised(self):
        s = self._report_seq()
        assert len(s["blastx_nr_hits"]) == 1

    def test_sequence_preserved(self):
        s = self._report_seq()
        assert s["sequence"] == "ATGC" * 300


# ---------------------------------------------------------------------------
# ReportExporter — taxonomy
# ---------------------------------------------------------------------------

class TestTaxonomySerialisation:
    def test_class_renamed(self):
        seq = _make_confirmed_seq()
        tax = _taxonomy()
        # set a non-None class_ to verify rename
        seq.taxonomy = Taxonomy(
            tax_id=1, scientific_name="Virus X", no_rank=None, clade=None,
            kingdom=None, phylum=None, class_="Insthoviricetes", order=None,
            family="Orthomyxoviridae", subfamily=None, genus=None,
            species=None, genome="ssRNA(-)",
        )
        exp = ReportExporter()
        s   = exp.export([seq], [])["sequences"][0]
        tax_d = s["taxonomy"]
        assert "class"  in tax_d
        assert "class_" not in tax_d
        assert tax_d["class"] == "Insthoviricetes"

    def test_taxonomy_none_preserved(self):
        seq = _make_confirmed_seq()
        seq.taxonomy = None
        exp = ReportExporter()
        s   = exp.export([seq], [])["sequences"][0]
        assert s["taxonomy"] is None

    def test_all_taxonomy_fields_present(self):
        seq = _make_confirmed_seq()
        exp = ReportExporter()
        s   = exp.export([seq], [])["sequences"][0]
        for key in ("tax_id", "scientific_name", "family", "genus", "order",
                    "species", "genome", "class"):
            assert key in s["taxonomy"], f"missing taxonomy key: {key}"


# ---------------------------------------------------------------------------
# ReportExporter — clusters
# ---------------------------------------------------------------------------

class TestClusterSerialisation:
    def _cluster_dict(self) -> dict:
        confirmed = _make_confirmed_seq("s1")
        cluster   = _make_cluster("s1")
        cluster.members.append(
            ClusterMember(seq_id="s2", length=800, is_representative=False,
                          identity=88.0, query_coverage=75.0, aln_start=50, aln_end=850)
        )
        exp = ReportExporter()
        return exp.export([confirmed], [cluster])["clusters"][0]

    def test_cluster_keys_present(self):
        c = self._cluster_dict()
        for key in ("cluster_id", "species", "representative_id", "size", "members"):
            assert key in c, f"missing key: {key}"

    def test_cluster_size_is_member_count(self):
        c = self._cluster_dict()
        assert c["size"] == len(c["members"])
        assert c["size"] == 2

    def test_cluster_members_serialised(self):
        c = self._cluster_dict()
        ids = {m["seq_id"] for m in c["members"]}
        assert "s1" in ids

    def test_representative_flagged_in_members(self):
        c    = self._cluster_dict()
        reps = [m for m in c["members"] if m["is_representative"]]
        assert len(reps) == 1
        assert reps[0]["seq_id"] == "s1"

    def test_cluster_size_is_int(self):
        c = self._cluster_dict()
        assert isinstance(c["size"], int)


# ---------------------------------------------------------------------------
# ReportExporter — file writing
# ---------------------------------------------------------------------------

class TestFileWriting:
    def test_file_created(self, tmp_path):
        out = tmp_path / "report.json"
        exp = ReportExporter()
        exp.export([], [], output_path=out)
        assert out.exists()

    def test_written_json_is_valid(self, tmp_path):
        out  = tmp_path / "report.json"
        exp  = ReportExporter()
        exp.export([], [], output_path=out)
        data = json.loads(out.read_text())
        assert "meta"      in data
        assert "sequences" in data
        assert "clusters"  in data

    def test_nested_parents_created(self, tmp_path):
        out = tmp_path / "a" / "b" / "report.json"
        ReportExporter().export([], [], output_path=out)
        assert out.exists()

    def test_returned_dict_matches_file_content(self, tmp_path):
        out     = tmp_path / "report.json"
        seq     = _make_confirmed_seq()
        cluster = _make_cluster("seq1")
        exp     = ReportExporter()
        returned = exp.export([seq], [cluster], output_path=out)
        written  = json.loads(out.read_text())
        assert returned["meta"]["viralquest_version"] == written["meta"]["viralquest_version"]
        assert len(returned["sequences"]) == len(written["sequences"])
        assert len(returned["clusters"])  == len(written["clusters"])

    def test_no_file_when_output_path_none(self, tmp_path):
        exp = ReportExporter()
        exp.export([], [], output_path=None)
        assert list(tmp_path.iterdir()) == []


# ---------------------------------------------------------------------------
# ReportExporter — JSON round-trip
# ---------------------------------------------------------------------------

class TestJsonRoundTrip:
    def test_full_report_json_serialisable(self):
        seq     = _make_confirmed_seq()
        cluster = _make_cluster("seq1")
        fasta   = InputFasta(name="test.fa", num_seqs=1, path="/data/test.fa", size=2048)
        exp     = ReportExporter(include_llm=True)
        report  = exp.export([seq], [cluster], input_fasta=fasta, version="0.9.0")
        # must not raise
        serialised = json.dumps(report)
        restored   = json.loads(serialised)
        assert restored["sequences"][0]["id"] == "seq1"
        assert restored["clusters"][0]["size"] == 1

    def test_top_level_keys(self):
        report = ReportExporter().export([], [])
        assert set(report.keys()) == {"meta", "sequences", "clusters"}
