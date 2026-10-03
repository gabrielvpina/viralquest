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


def test_virome_tab_in_consolidated_report(results_root):
    data = build_report_data(load_samples(results_root))
    html = _render_template(data, d3_js=STUB_D3)
    assert 'data-section="virome"' in html and 'id="section-virome"' in html
    assert "window.vqInitVirusesReport = vqInitVirusesReport" in html   # sample-filter wrapper
    assert "window.vqInitViruses = vqInitViruses" in html               # shared Virome tab
    assert "vqInitVirusesReport(VQ_REPORT.sequences" in html
    # overview helpers the Virome tab needs are exposed by the report's export.js
    for name in ("statChip:", "fmtNum:", "pct:", "cardGroup:", "redrawOnResize:"):
        assert name in html, name


def test_virome_tab_has_a_single_source():
    """The consolidated report inlines the per-sample report's section_viruses.js."""
    from pathlib import Path
    import viralquest.report_builder as rb
    assert (rb._SAMPLE_COMPONENTS / "section_viruses.js").is_file()
    assert not (Path(rb._COMPONENTS) / "section_viruses.js").exists()


def test_export_helpers_identical_in_both_reports():
    """Overview helpers are copied into components-report/export.js — keep them in sync."""
    from pathlib import Path
    import re
    import viralquest.report_builder as rb

    def block(path):
        text = Path(path).read_text(encoding="utf-8")
        return text[text.index("// ── Overview tabs (Run · Virome)"):text.index("// ── Expose on window.VQ")]

    assert block(rb._SAMPLE_COMPONENTS / "export.js") == block(rb._COMPONENTS / "export.js")


# ── Pipeline workflow per sample ────────────────────────────────────────────

_WORKFLOW = {
    "started_at": "2026-10-02T17:46:51+00:00",
    "finished_at": "2026-10-02T18:28:51+00:00",
    "total_seconds": 2519.6,
    "options": {"input": "contigs.fasta", "threads": 4, "cap3": False, "nr_db": None,
                "blastn": "online · nt", "reads": [], "llm": None, "force": False},
    "steps": [
        {"key": "parse", "label": "Parse FASTA", "status": "done", "seconds": 0.1,
         "details": {"Input sequences": 12}, "message": None},
        {"key": "nr", "label": "Diamond BLASTx NR", "status": "skipped", "seconds": None,
         "details": {}, "message": "--nr-db not set"},
        {"key": "blastn", "label": "BLASTn online", "status": "partial", "seconds": 90.5,
         "details": {"Failed": 3}, "message": "3 sequences failed"},
    ],
}


def test_workflow_summarized_per_sample(tmp_path):
    rep = _report([_seq("k141_1")])
    rep["pipeline_stats"]["workflow"] = _WORKFLOW
    _write_sample(tmp_path, "withWf", rep)
    _write_sample(tmp_path, "legacy", _report([_seq("k141_1")]))

    data = build_report_data(load_samples(tmp_path))
    by = {s["sample"]: s for s in data["samples"]}

    wf = by["withWf"]["workflow"]
    assert wf["version"] == "3.0.0"
    assert wf["total_seconds"] == 2519.6
    assert wf["options"]["blastn"] == "online · nt"
    assert [s["key"] for s in wf["steps"]] == ["parse", "nr", "blastn"]
    assert wf["steps"][2] == {"key": "blastn", "label": "BLASTn online", "status": "partial",
                              "seconds": 90.5, "details": {"Failed": 3},
                              "message": "3 sequences failed"}
    # Runs older than the workflow record keep the detected-steps fallback.
    assert by["legacy"]["workflow"] is None
    assert by["legacy"]["steps"]["blastn"] is True


def test_overview_renders_workflow_card(results_root):
    html = _render_template(build_report_data(load_samples(results_root)), d3_js=STUB_D3)
    assert "Pipeline Workflow per Sample" in html
    assert "_renderWorkflowMatrix(samples)" in html
    assert ".ov-wf-matrix" in html and ".vq-wf-pill--error" in html   # card CSS shipped


def test_salmon_summary_has_viral_reads(tmp_path):
    rep = _report([_seq("k141_1"), _seq("k141_2")], salmon=True)
    rep["salmon_quant"].update(total_reads=1000)
    for i, e in enumerate(rep["salmon_quant"]["viral_quant"]):
        e["num_reads"] = 10.5 + i
    _write_sample(tmp_path, "s1", rep)
    salmon = build_report_data(load_samples(tmp_path))["samples"][0]["salmon"]
    assert salmon["viral_reads"] == 22.0
    assert salmon["total_reads"] == 1000


def test_overview_cards(results_root):
    html = _render_template(build_report_data(load_samples(results_root)), d3_js=STUB_D3)
    for card in ("ov-funnel-card", "ov-bubbles-card", "ov-salmon-card", "Putative Novel",
                 "_renderLengthDensity(samples)", "_renderScoreBoxes('ov-heur'"):
        assert card in html, card
    # Removed from the overview.
    for gone in ("Viral Family Diversity per Sample", "Cross-sample Clusters per Sample",
                 "Family Distribution per Sample"):
        assert gone not in html, gone
    # The overview reuses the Virome tab's novelty rule rather than a copy.
    assert "window.vqNovelty = {" in html


def test_detected_viruses_table_is_per_contig_with_csv(results_root):
    """Both reports inline the same section_viruses.js: per-contig table + CSV export."""
    html = _render_template(build_report_data(load_samples(results_root)), d3_js=STUB_D3)
    assert "function _contigRows(viral)" in html
    assert 'id="stats-viruses-csv"' in html and "_virusTableCsv(contigs, tableState)" in html
    for col in ("'BLASTx id %'", "'BLASTx cov %'", "'BLASTn id %'", "'BLASTn cov %'",
                "'BLASTx best hit'", "'BLASTn best hit'", "'Novelty'"):
        assert col in html, col
    assert "Contigs by novelty" not in html


def test_clusters_tab_is_cluster_centred(results_root):
    """General Clusters: selector table, single-cluster detail, presence matrix, sample similarity."""
    html = _render_template(build_report_data(load_samples(results_root)), d3_js=STUB_D3)
    for marker in ('id="clu-table"', 'id="clu-detail"', 'id="clu-matrix"',
                   'id="clu-csv"', "VQ.upgma(sim.M)", "function _clustersCsv()"):
        assert marker in html, marker
    # The old one-card-per-cluster list and ego network are gone.
    assert "function _drawEgo" not in html and 'id="clu-cards"' not in html
    assert 'id="clu-sim"' not in html             # sample-similarity card removed


# ── RNA quantification ──────────────────────────────────────────────────────

def test_salmon_reads_stamped_and_library_size(tmp_path):
    rep = _report([_seq("k141_1"), _seq("k141_2")], salmon=True)
    rep["salmon_quant"].update(mapping_rate=80.0, total_reads=8000)
    rep["salmon_quant"]["viral_quant"][0]["num_reads"] = 120.5
    _write_sample(tmp_path, "s1", rep)
    sample = load_samples(tmp_path)[0]
    by_id = {q["id"]: q for q in sample.sequences}
    assert by_id["k141_1"]["reads"] == 120.5 and by_id["k141_1"]["tpm"] == 10.0
    assert by_id["k141_2"]["reads"] is None          # Salmon ran, no NumReads for it
    salmon = build_report_data([sample])["samples"][0]["salmon"]
    assert salmon["input_reads"] == 10000             # 8000 mapped / 80 %


def test_cluster_members_carry_reads():
    seqs = [_gseq("sA::s1", "sA", tpm=12.0), _gseq("sB::s1", "sB", tpm=3.0)]
    seqs[0]["reads"], seqs[1]["reads"] = 40.0, 7.5
    with patch("viralquest.report_clusters._run_blastn", return_value={"sB::s1": _hit()}):
        members = build_general_clusters(seqs)[0].to_dict()["members"]
    assert {m["gid"]: m["reads"] for m in members} == {"sA::s1": 40.0, "sB::s1": 7.5}


def test_quant_tab_is_target_centred(results_root):
    html = _render_template(build_report_data(load_samples(results_root)), d3_js=STUB_D3)
    # One target picker at the top; every card below follows it.
    for marker in ('id="qt-group"', 'id="qt-select"', 'id="qt-search"', 'id="qt-metric"',
                   'id="qt-load"', 'id="qt-share"', 'id="qt-hk"', 'id="qt-parts"', 'id="qt-contigs"',
                   'id="qt-co"', 'id="qt-load-csv"', 'id="qt-contigs-csv"', "vq:cluster-select",
                   'id="qt-hm"', 'id="qt-hm-csv"', 'id="qt-hm-level"', 'id="qt-hm-filter"',
                   'id="qt-hm-metric"', "function _heatmapRows()", 'id="qt-hm-pal"', "const _HM_PALETTES"):
        assert marker in html, marker
    # Whole-virome views and the old pickers are gone from this tab.
    for gone in ('id="qt-mx"', 'id="qt-comp"', 'id="qt-div"', 'id="qt-list"', 'id="hm-clu-list"', 'id="qt-kg"'):
        assert gone not in html, gone
    # The host kingdom matrix (not target-specific) lives in the Overview again.
    assert 'id="ov-kingdoms-card"' in html
    assert "upgma:          vqUpgma" in html


# ── Named dataset exports (JSON + SQLite) ───────────────────────────────────

def test_dataset_name_validation():
    from viralquest.report_db import validate_dataset_name
    for ok in ("PRJNA123456", "Aedes_aegypti_2024", "study-1.v2"):
        assert validate_dataset_name(ok) == ok
    for bad in ("", "bad name", "-starts-with-dash", "a/b", "x" * 65, None):
        with pytest.raises(ValueError):
            validate_dataset_name(bad)


def test_cli_requires_a_dataset_name(results_root, tmp_path):
    from viralquest.cli_report import _build_parser
    with pytest.raises(SystemExit):
        _build_parser().parse_args(["-in", str(results_root), "-out", str(tmp_path)])
    with pytest.raises(SystemExit):
        _build_parser().parse_args(["-in", str(results_root), "-out", str(tmp_path), "-n", "bad name"])
    args = _build_parser().parse_args(["-in", str(results_root), "-out", str(tmp_path), "-n", "PRJ1"])
    assert args.name == "PRJ1"


def _dataset_files(results_root, tmp_path, name):
    from viralquest.report_db import output_paths, write_json, write_sqlite
    data = build_report_data(load_samples(results_root), version="3.0.0", dataset=name)
    j, db = output_paths(tmp_path, name)
    return write_json(data, j), write_sqlite(data, db), data


def test_json_export_carries_name_and_schema(results_root, tmp_path):
    j, _, data = _dataset_files(results_root, tmp_path, "PRJ1")
    assert j.name == "PRJ1.viralquest.json"
    doc = json.loads(j.read_text())
    assert doc["dataset"] == "PRJ1" and doc["meta"]["dataset"] == "PRJ1"
    assert doc["schema"] == "viralquest-report" and doc["schema_version"] >= 1
    assert len(doc["sequences"]) == len(data["sequences"]) == 4


def test_sqlite_export_tables(results_root, tmp_path):
    import sqlite3
    _, db, _ = _dataset_files(results_root, tmp_path, "PRJ1")
    con = sqlite3.connect(db)
    assert con.execute("SELECT dataset, n_samples, n_sequences FROM dataset_info").fetchall() == [("PRJ1", 3, 4)]
    assert con.execute("SELECT count(*) FROM samples WHERE dataset='PRJ1'").fetchone()[0] == 3
    assert con.execute("SELECT count(DISTINCT dataset) FROM sequences").fetchone()[0] == 1
    fam = dict(con.execute("SELECT gid, family FROM sequences").fetchall())
    assert fam["sampleB::k141_1"] == "Flaviviridae"
    assert con.execute("SELECT count(*) FROM blastx_hits WHERE source='nr'").fetchone()[0] == 3
    # Salmon reads → RPM is only set where the sample ran Salmon.
    assert con.execute("SELECT count(*) FROM sequences WHERE tpm IS NOT NULL").fetchone()[0] == 2
    assert con.execute("PRAGMA user_version").fetchone()[0] >= 1


def test_two_datasets_merge_without_clashes(results_root, tmp_path):
    """Same sample / contig ids in two datasets: the dataset column keeps them apart."""
    import sqlite3
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    _, db_a, _ = _dataset_files(results_root, tmp_path / "a", "STUDY_A")
    _, db_b, _ = _dataset_files(results_root, tmp_path / "b", "STUDY_B")
    con = sqlite3.connect(db_a)
    con.execute("ATTACH ? AS b", (str(db_b),))
    for table in ("dataset_info", "samples", "sequences", "blastx_hits", "blastn_hits", "orfs", "domains"):
        con.execute(f"INSERT INTO {table} SELECT * FROM b.{table}")
    rows = con.execute("SELECT dataset, count(*) FROM sequences GROUP BY dataset ORDER BY dataset").fetchall()
    assert rows == [("STUDY_A", 4), ("STUDY_B", 4)]
