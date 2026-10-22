"""Tests for the run's output layout: OutputOrganizer, novelty FASTAs, read tables, SQLite."""
import json
import sqlite3
from pathlib import Path

from viralquest.biodata import HeuristicScore, NucSequence
from viralquest.output import OutputOrganizer
from viralquest.report_db import write_sample_sqlite

EXAMPLE = Path(__file__).parent / "example.json"


def _seq(sid, novelty=None):
    s = NucSequence(id=sid, sequence="ACGT" * 5)
    if novelty:
        s.heuristic_output = HeuristicScore(seq_id=sid, vq_score=80, classification="viral-known",
                                            analysis="", novelty=novelty)
    return s


def test_layout_directories(tmp_path):
    org = OutputOrganizer(tmp_path, "S1")
    assert org.diamond_dir == tmp_path / "features" / "diamond"
    assert org.hmm_dir == tmp_path / "features" / "hmm"
    assert org.blastn_dir == tmp_path / "features" / "blastn"
    assert org.synteny_dir.parent == org.features_dir
    assert org.coverage_dir == tmp_path / "reads" / "coverage"
    assert org.salmon_dir == tmp_path / "reads" / "salmon"
    assert org.report_path("sqlite") == tmp_path / "S1_viralquest.sqlite"
    # reads/ only appears once a read step writes into it
    assert org.fastas_dir.is_dir() and not org.reads_dir.exists()


def test_main_fasta_lives_in_fastas(tmp_path):
    org = OutputOrganizer(tmp_path, "S1")
    path = org.save_viral_contigs([_seq("a"), _seq("b")])
    assert path == tmp_path / "fastas" / "S1_viral_contigs.fasta"
    assert path.read_text().count(">") == 2


def test_one_fasta_per_novelty_tier(tmp_path):
    org = OutputOrganizer(tmp_path, "S1")
    seqs = [_seq("a", "known"), _seq("b", "divergent"), _seq("c", "known"), _seq("d")]
    out = org.save_novelty_fastas(seqs)
    assert list(out) == ["known", "divergent", "unscored"]
    assert out["known"].name == "S1_known.fasta"
    assert out["known"].read_text() == ">a\nACGTACGTACGTACGTACGT\n>c\nACGTACGTACGTACGTACGT\n"
    assert not (tmp_path / "fastas" / "S1_variant.fasta").exists()


def test_read_tables_only_with_read_data(tmp_path):
    org = OutputOrganizer(tmp_path, "S1")
    assert org.save_read_tables({"sequences": [{"id": "a"}]}) == []
    report = {
        "sequences": [{"id": "a", "coverage": {"length": 10, "mean_depth": 5.0, "breadth_1x": 1.0}}],
        "salmon_quant": {"viral_quant": [{"name": "VQ_VIRAL_a", "tpm": 1.0, "num_reads": 2.0,
                                          "seq_type": "viral"}]},
    }
    paths = org.save_read_tables(report)
    assert [p.name for p in paths] == ["coverage_summary.tsv", "quantification.tsv"]
    lines = (tmp_path / "reads" / "quantification.tsv").read_text().splitlines()
    assert lines[0].startswith("Target\tType") and lines[1].startswith("VQ_VIRAL_a\tviral")


def test_sample_sqlite_from_example(tmp_path):
    report = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    db = write_sample_sqlite(report, tmp_path / "ex_viralquest.sqlite")
    con = sqlite3.connect(db)
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    assert tables == {"sequences", "blastx_hits", "blastn_hits", "blastn_taxonomy",
                      "orfs", "domains", "clusters", "cluster_members", "coverage",
                      "seq_quality", "synteny_hits", "quantification", "host_viral_hits"}
    cols = {r[1] for r in con.execute("PRAGMA table_info(sequences)")}
    assert "seq_id" in cols and not {"dataset", "sample", "gid"} & cols
    assert con.execute("SELECT COUNT(*) FROM sequences").fetchone()[0] == len(report["sequences"])
    n_orfs = len({(s["id"], o["name"]) for s in report["sequences"] for o in s.get("orfs") or []})
    assert con.execute("SELECT COUNT(*) FROM orfs").fetchone()[0] == n_orfs
    assert con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0] == len(report["clusters"])
    con.close()


def test_sample_sqlite_quantification_and_rpm(tmp_path):
    report = {
        "sequences": [{"id": "a", "length": 100, "sequence": "A" * 100}],
        "clusters": [],
        "salmon_quant": {"total_reads": 800, "mapping_rate": 80.0,
                         "viral_quant": [{"name": "VQ_VIRAL_a", "tpm": 12.0, "num_reads": 10.0}]},
    }
    con = sqlite3.connect(write_sample_sqlite(report, tmp_path / "x.sqlite"))
    tpm, reads, rpm = con.execute("SELECT tpm, reads, rpm FROM sequences").fetchone()
    con.close()
    assert (tpm, reads) == (12.0, 10.0)
    assert rpm == 10.0 / 1000 * 1e6          # 1000 input reads = 800 mapped / 80 %
