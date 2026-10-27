"""Tests for viralquest.report_synteny — the multi-sample general synteny."""
import json
import shutil
from pathlib import Path

from viralquest.biodata import BlastnResult, NucSequence, Orf
from viralquest.exporter import ReportExporter
from viralquest.report_loader import SampleReport, _stamp
from viralquest.report_synteny import _groups, _similar_pairs, build_general_synteny
from viralquest.synteny import GenBankFetcher, SyntenyAnalyzer, parse_genbank

MOCK = Path(__file__).parent / "mock" / "synteny"
DENGUE, HIV = MOCK / "NC_001477.gb", MOCK / "NC_001802.gb"


def _offline(params):
    raise RuntimeError("offline")


class _Aligner:
    """Perfect HSPs of every query on the first subject; *proteome* rows for the
    reference ↔ reference search; contig ↔ contig: every query on every other."""

    def __init__(self, proteome=()):
        self.proteome, self.tags = list(proteome), []

    def align(self, queries, subjects, tag=""):
        self.tags.append(tag)
        if tag == "proteome":
            return self.proteome
        if tag == "contigs":
            return [_row(a, b, len(queries[a])) for a in queries for b in queries if a != b]
        sid = next(iter(subjects))
        return [_row(q, sid, len(aa)) for q, aa in queries.items()]


def _row(q, s, n, pident=99.0):
    return {"qseqid": q, "sseqid": s, "pident": pident, "length": n, "qstart": 1, "qend": n,
            "sstart": 1, "send": n, "evalue": 1e-50, "bitscore": 200.0, "qlen": n, "slen": 3392}


def _sample(tmp_path, name, acc, title, gb):
    """A loaded sample with one sequence whose per-sample synteny used *acc*."""
    s = NucSequence(id="c1", sequence="A" * 300)
    s.is_viral = True
    s.blastn_hits = [BlastnResult(qseqid="c1", qlen=300, slen=10000, qcovhsp=90, pident=97.0,
                                  evalue=0.0, bit_score=500.0, stitle=title, accession=acc)]
    s.orfs = [Orf(start_codon="ATG", stop_codon="TAA", start_position=0, stop_position=300,
                  strand="+", frame=1, length_aa=100, length_nt=300, bigger_than_50=True,
                  aa_sequence="M" * 100, nuc_sequence="", orf_type="complete", name="c1_orf")]
    d = tmp_path / name
    refs_dir = d / "features" / "synteny" / "refs"
    refs_dir.mkdir(parents=True)
    shutil.copyfile(gb, refs_dir / gb.name)
    refs = SyntenyAnalyzer(d / "features" / "synteny", fetcher=GenBankFetcher(refs_dir, http=_offline),
                           aligner=_Aligner()).run([s])
    rep = json.loads(json.dumps(ReportExporter(force=True).export([s], [], synteny_refs=refs)))
    _stamp(rep, name)
    return SampleReport(sample=name, directory=d, json_path=d / f"{name}_viralquest.json", report=rep)


def test_same_reference_in_two_samples_is_one_shared_group(tmp_path):
    a = _sample(tmp_path, "A", "NC_001477", "Dengue virus 1", DENGUE)
    b = _sample(tmp_path, "B", "NC_001477", "Dengue virus 1", DENGUE)
    out_dir = tmp_path / "out"
    gs = build_general_synteny([a, b], out_dir,
                               fetcher=GenBankFetcher(out_dir / "refs", http=_offline),
                               aligner=_Aligner())
    assert gs["clusters"] == [{"reference": "NC_001477.1", "members": ["A::c1", "B::c1"],
                               "samples": ["A", "B"],
                               "merged": [{"accession": "NC_001477.1", "organism": "dengue virus type 1",
                                           "length": 10735, "n_members": 2}]}]
    assert (out_dir / "refs" / "NC_001477.gb").is_file()          # copied from a sample
    syn = a.sequences[0]["general_synteny"]
    assert syn["reference"] == "NC_001477.1" and syn["hits"]
    assert [l["other_seq"] for l in syn["links"]] == ["B::c1"]     # cross-sample link
    assert a.sequences[0]["synteny"]["reference"] == "NC_001477.1"  # per-sample untouched


def test_similar_references_are_pooled_on_one_anchor(tmp_path):
    a = _sample(tmp_path, "A", "NC_001477", "Dengue virus 1", DENGUE)
    b = _sample(tmp_path, "B", "NC_001477", "Dengue virus 1", DENGUE)
    c = _sample(tmp_path, "C", "NC_001802", "Human immunodeficiency virus 1", HIV)
    n_den = len(parse_genbank(DENGUE).proteins[0])
    proteome = [_row("NC_001802|0", "NC_001477|0", 300, pident=85.0)]  # HIV gag ~ a Dengue protein
    out_dir = tmp_path / "out"
    gs = build_general_synteny([a, b, c], out_dir,
                               fetcher=GenBankFetcher(out_dir / "refs", http=_offline),
                               aligner=_Aligner(proteome), min_coverage=1)
    assert n_den > 300
    [cl] = gs["clusters"]
    assert cl["reference"] == "NC_001477.1"                      # two members beat one
    assert cl["samples"] == ["A", "B", "C"]
    assert [r["accession"] for r in cl["merged"]] == ["NC_001477.1", "NC_001802.1"]
    assert c.sequences[0]["general_synteny"]["reference"] == "NC_001477.1"


def test_dissimilar_references_stay_apart(tmp_path):
    a = _sample(tmp_path, "A", "NC_001477", "Dengue virus 1", DENGUE)
    c = _sample(tmp_path, "C", "NC_001802", "Human immunodeficiency virus 1", HIV)
    proteome = [_row("NC_001802|0", "NC_001477|0", 300, pident=40.0)]
    out_dir = tmp_path / "out"
    gs = build_general_synteny([a, c], out_dir,
                               fetcher=GenBankFetcher(out_dir / "refs", http=_offline),
                               aligner=_Aligner(proteome), min_coverage=1)
    assert sorted(c_["reference"] for c_ in gs["clusters"]) == ["NC_001477.1", "NC_001802.1"]


def test_no_synteny_references_returns_none(tmp_path):
    s = SampleReport(sample="A", directory=tmp_path, json_path=tmp_path / "a.json",
                     report={"sequences": [{"id": "c1", "gid": "A::c1"}]})
    assert build_general_synteny([s], tmp_path / "out") is None


def test_similar_pairs_needs_identity_and_coverage():
    class _R:
        def __init__(self, n):
            self.proteins = ["M" * n]
    refs = {"X": _R(100), "Y": _R(1000)}
    rows = [_row("X|0", "Y|0", 80, pident=90.0)]
    al = type("A", (), {"align": lambda self, q, s, tag="": rows})()
    assert _similar_pairs(refs, al, 70, 50) == [("X", "Y")]      # 80 % of X's proteome
    assert _similar_pairs(refs, al, 95, 50) == []
    assert _similar_pairs(refs, al, 70, 90) == []


def test_groups_are_connected_components():
    gs = _groups(["a", "b", "c", "d"], [("a", "b"), ("b", "c")])
    assert sorted(sorted(g) for g in gs) == [["a", "b", "c"], ["d"]]


def test_html_tab_and_workbook_sheet(tmp_path):
    from viralquest.report_builder import _render_template, _synteny_css
    from viralquest.report_data import build_report_data
    from viralquest.seq_table import _general_synteny
    a = _sample(tmp_path, "A", "NC_001477", "Dengue virus 1", DENGUE)
    b = _sample(tmp_path, "B", "NC_001477", "Dengue virus 1", DENGUE)
    gs = build_general_synteny([a, b], tmp_path / "out",
                               fetcher=GenBankFetcher(tmp_path / "out" / "refs", http=_offline),
                               aligner=_Aligner())
    data = build_report_data([a, b], general_synteny=gs)
    html = _render_template(data, d3_js="")
    assert 'id="tab-synteny"' in html and "{{VQ_SYNTENY}}" not in html
    assert "function vqInitSynteny" in html and ".vq-syn-member" in _synteny_css()
    sheet = _general_synteny(data)
    assert sheet.columns[0] == "Sample"
    assert [(r["Sample"], r["Sequence ID"], r["Reference"]) for r in sheet.rows] == \
        [("A", "c1", "NC_001477.1"), ("B", "c1", "NC_001477.1")]
