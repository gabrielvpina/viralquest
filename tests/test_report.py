import json
from pathlib import Path
from unittest.mock import patch

import pytest

from viralquest.report_builder import _render_template
from viralquest.report_clusters import (
    build_general_clusters,
    write_artifacts,
)
from viralquest.report_data import build_report_data
from viralquest.report_loader import (
    GID_SEP,
    discover_result_dirs,
    load_samples,
)

STUB_D3 = "/* d3-stub */"


# ── Fixtures ────────────────────────────────────────────────────────────────

def _seq(seq_id, family="Picornaviridae", length=600, nr=True, tpm=None, heur=None):
    s = {
        "id": seq_id,
        "length": length,
        "is_viral": True,
        "taxonomy": {"family": family},
        "blastx_hits": [{"species": "Some refseq virus"}],
        "blastx_nr_hits": [{"species": "Some NR virus"}] if nr else [],
        "blastn_hits": [{"x": 1}],
        "orfs": [],
        "heuristic_output": {"vq_score": heur, "classification": "viral-known"} if heur is not None else None,
    }
    return s


def _report(seqs, *, salmon=False, llm=False, cap3=False, nr=True):
    rep = {
        "meta": {"viralquest_version": "3.0.0", "input_file": {"name": "contigs.fasta"}},
        "pipeline_stats": {
            "blast": {
                "total_input": 1000,
                "total_viral_flagged": 50,
                "total_confirmed": len(seqs),
                "nr_unique_seqs": len(seqs) if nr else 0,
                "blastn_unique_seqs": len(seqs),
            },
            "llm": {"present": llm},
        },
        "sequences": seqs,
        "clusters": [],
    }
    if cap3:
        rep["pipeline_stats"]["cap3"] = {"used": True, "contigs": 5, "singlets": 9}
    if salmon:
        rep["salmon_quant"] = {
            "mapping_rate": 0.8,
            "viral_quant": [
                {"name": f"VQ_VIRAL_{s['id']}", "tpm": 10.0 + i}
                for i, s in enumerate(seqs)
            ],
        }
    return rep


def _write_sample(root: Path, name: str, report: dict) -> Path:
    d = root / name
    d.mkdir(parents=True)
    (d / f"{name}_viralquest.json").write_text(json.dumps(report))
    return d


@pytest.fixture
def results_root(tmp_path):
    _write_sample(tmp_path, "sampleA",
                  _report([_seq("k141_1"), _seq("k141_2")], salmon=True))
    _write_sample(tmp_path, "sampleB",
                  _report([_seq("k141_1", family="Flaviviridae", nr=False)], nr=False))
    _write_sample(tmp_path, "sampleC",
                  _report([_seq("k141_1", heur=70.0)], llm=True, cap3=True))
    return tmp_path


# ── Discovery ───────────────────────────────────────────────────────────────

def test_discover_finds_child_dirs(results_root):
    found = discover_result_dirs(results_root)
    names = sorted(d.name for d, _ in found)
    assert names == ["sampleA", "sampleB", "sampleC"]


def test_discover_root_is_single_result(tmp_path):
    _write_sample(tmp_path, "solo", _report([_seq("k141_1")]))
    found = discover_result_dirs(tmp_path / "solo")
    assert len(found) == 1


def test_load_samples_empty_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_samples(tmp_path)


# ── Loading / stamping / namespacing ───────────────────────────────────────

def test_sample_attribution_and_gid(results_root):
    samples = load_samples(results_root)
    by_name = {s.sample: s for s in samples}
    seq = by_name["sampleA"].sequences[0]
    assert seq["sample"] == "sampleA"
    assert seq["gid"] == f"sampleA{GID_SEP}k141_1"
    # Colliding raw ids across samples become distinct gids.
    gids = {s.sequences[0]["gid"] for s in samples}
    assert gids == {"sampleA::k141_1", "sampleB::k141_1", "sampleC::k141_1"}


def test_salmon_tpm_stamped(results_root):
    samples = load_samples(results_root)
    a = next(s for s in samples if s.sample == "sampleA")
    assert a.sequences[0]["tpm"] == 10.0
    b = next(s for s in samples if s.sample == "sampleB")
    assert "tpm" not in b.sequences[0]  # no salmon → not stamped


# ── Aggregation ─────────────────────────────────────────────────────────────

def test_build_report_data_steps_and_summaries(results_root):
    samples = load_samples(results_root)
    data = build_report_data(samples, input_root=str(results_root), version="3.0.0")

    assert data["meta"]["n_samples"] == 3
    assert data["has_clusters"] is False
    assert len(data["sequences"]) == 4  # 2 + 1 + 1

    steps = {s["sample"]: s["steps"] for s in data["samples"]}
    assert steps["sampleA"] == {"nr": True, "blastn": True, "salmon": True, "llm": False, "cap3": False}
    assert steps["sampleB"]["nr"] is False and steps["sampleB"]["salmon"] is False
    assert steps["sampleC"]["llm"] is True and steps["sampleC"]["cap3"] is True

    c = next(s for s in data["samples"] if s["sample"] == "sampleC")
    assert c["heuristic_scores"] == [70.0]


# ── Rendering ───────────────────────────────────────────────────────────────

# ── Cross-sample clustering (BLASTn mocked) ─────────────────────────────────

def _gseq(gid, sample, species="Virus A", db="nr", length=600, tpm=None):
    """A stamped sequence dict as produced by report_loader."""
    s = {
        "gid": gid, "sample": sample, "id": gid.split(GID_SEP)[1],
        "sequence": "ACGT" * (length // 4), "length": length,
        "blastx_nr_hits": [{"species": species}] if db == "nr" else [],
        "blastx_hits": [{"species": species}] if db == "refseq" else [],
        "blastn_hits": [],
    }
    if tpm is not None:
        s["tpm"] = tpm
    return s


def _hit(pident=99.0, qcov=95.0):
    return {"pident": pident, "qcov": qcov, "bitscore": 500.0,
            "qstart": 1, "qend": 580, "sstart": 1, "send": 580}


def test_cluster_spans_two_samples():
    seqs = [_gseq("sA::s1", "sA", length=600), _gseq("sB::s1", "sB", length=590)]
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit()}):
        clusters = build_general_clusters(seqs)
    assert len(clusters) == 1
    c = clusters[0]
    assert c.gid == "VQR_CLU_0001"
    assert c.samples == ["sA", "sB"]
    assert c.representative == "sA::s1"  # longest is representative
    assert any(m.is_representative for m in c.members)


def test_no_cluster_when_single_sample():
    seqs = [_gseq("sA::s1", "sA", length=600), _gseq("sA::s2", "sA", length=590)]
    with patch("viralquest.report_clusters._run_blastn", return_value={"sA::s2": _hit()}):
        clusters = build_general_clusters(seqs)
    assert clusters == []  # needs >=2 distinct samples


def test_member_below_floor_excluded():
    seqs = [_gseq("sA::s1", "sA", length=600), _gseq("sB::s1", "sB", length=590)]
    # coverage below the 70 floor → member dropped → no 2-sample cluster
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit(qcov=50.0)}):
        clusters = build_general_clusters(seqs)
    assert clusters == []


def test_homogeneity_flag():
    same = [_gseq("sA::s1", "sA", species="Virus A"), _gseq("sB::s1", "sB", species="Virus A")]
    diff = [_gseq("sA::s1", "sA", species="Virus A"), _gseq("sB::s1", "sB", species="Virus B")]
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit()}):
        assert build_general_clusters(same)[0].homogeneous is True
        assert build_general_clusters(diff)[0].homogeneous is False


def test_write_artifacts(tmp_path):
    seqs = [_gseq("sA::s1", "sA", tpm=12.0), _gseq("sB::s1", "sB", tpm=None)]
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit()}):
        clusters = build_general_clusters(seqs)
    write_artifacts(clusters, {s["gid"]: s for s in seqs}, tmp_path)
    assert (tmp_path / "clusters" / "VQR_CLU_0001.fasta").exists()
    assert (tmp_path / "quantification" / "VQR_CLU_0001.tsv").exists()
    table = (tmp_path / "blastn" / "general_blastn.tsv").read_text()
    assert "sA::s1" in table and "sB::s1" in table


def test_cluster_to_dict_shape():
    seqs = [_gseq("sA::s1", "sA"), _gseq("sB::s1", "sB")]
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit()}):
        c = build_general_clusters(seqs)[0].to_dict()
    assert c["cluster_id"] == c["gid"]
    assert c["size"] == 2
    for m in c["members"]:
        assert m["query_coverage"] == m["coverage"]
        assert {"aln_start", "aln_end", "query_start"} <= set(m)


def test_render_template_no_unresolved_placeholders(results_root):
    samples = load_samples(results_root)
    data = build_report_data(samples, version="3.0.0")
    html = _render_template(data, d3_js=STUB_D3)
    import re
    assert not re.findall(r"\{\{VQ_\w+\}\}", html)
    assert "{{TITLE}}" not in html
    assert "vqInitOverview" in html
    assert STUB_D3 in html
