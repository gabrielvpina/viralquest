"""Tests for viralquest.synteny — reference download/parsing, blastp mapping, clusters."""
import shutil
from pathlib import Path

import pytest
from Bio.Seq import Seq

from viralquest.biodata import BlastnResult, NucSequence, Orf, SyntenyHit
from viralquest.exporter import ReportExporter
from viralquest.synteny import (
    GenBankFetcher,
    SyntenyAligner,
    SyntenyAnalyzer,
    aligned_segments,
    best_viral_hit,
    cds_aa_to_genome,
    orf_aa_to_contig,
    parse_genbank,
    placement,
)

MOCK = Path(__file__).parent / "mock" / "synteny"
DENGUE, HIV = MOCK / "NC_001477.gb", MOCK / "NC_001802.gb"


def _offline(params):
    raise RuntimeError("offline")


def _hit(acc="NC_001477", title="Dengue virus 1, complete genome", bits=500.0, slen=10735):
    return BlastnResult(qseqid="c", qlen=1000, slen=slen, qcovhsp=90, pident=97.0,
                        evalue=0.0, bit_score=bits, stitle=title, accession=acc)


def _orf(start, stop, strand="+", aa="M" * 50, name=None):
    return Orf(start_codon="ATG", stop_codon="TAA", start_position=start, stop_position=stop,
               strand=strand, frame=1, length_aa=len(aa), length_nt=stop - start,
               bigger_than_50=True, aa_sequence=aa, nuc_sequence="", orf_type="complete",
               name=name or f"orf_{start}_{stop}{strand}")


# ── Reference selection ─────────────────────────────────────────────────────

class TestBestViralHit:

    def test_picks_highest_bit_score_viral_title(self):
        s = NucSequence(id="c", sequence="A" * 100)
        s.blastn_hits = [_hit("H1", "Homo sapiens chromosome 3", bits=900),
                         _hit("V1", bits=300), _hit("V2", bits=400)]
        assert best_viral_hit(s).accession == "V2"

    def test_none_without_viral_hit(self):
        s = NucSequence(id="c", sequence="A" * 100)
        s.blastn_hits = [_hit("H1", "Aedes aegypti mRNA")]
        assert best_viral_hit(s) is None

    def test_oversized_subject_is_ignored(self):
        s = NucSequence(id="c", sequence="A" * 100)
        s.blastn_hits = [_hit("V1", slen=5_000_000)]
        assert best_viral_hit(s) is None


# ── GenBank fetch / cache / parse ───────────────────────────────────────────

class TestGenBankFetcher:

    def test_splits_multi_record_response_into_cache(self, tmp_path):
        text = DENGUE.read_text() + HIV.read_text()
        calls = []
        f = GenBankFetcher(tmp_path, http=lambda p: calls.append(p) or text, sleep=lambda s: None)
        paths = f.fetch(["NC_001477.1", "NC_001802"])
        assert set(paths) == {"NC_001477", "NC_001802"}
        assert calls[0]["rettype"] == "gbwithparts" and calls[0]["db"] == "nuccore"
        assert parse_genbank(paths["NC_001802"]).genome.accession == "NC_001802.1"

    def test_cached_records_are_not_fetched_again(self, tmp_path):
        shutil.copy(DENGUE, tmp_path / "NC_001477.gb")
        calls = []
        f = GenBankFetcher(tmp_path, http=lambda p: calls.append(p) or "")
        assert "NC_001477" in f.fetch(["NC_001477"])
        assert calls == []

    def test_download_failure_returns_nothing(self, tmp_path):
        assert GenBankFetcher(tmp_path, http=_offline).fetch(["X1"]) == {}


class TestParseGenbank:

    def test_polyprotein_with_mat_peptides(self):
        ref = parse_genbank(DENGUE)
        g = ref.genome
        assert g.organism.lower().startswith("dengue")
        assert len(g.cds) == 1 and g.cds[0].product == "polyprotein"
        assert g.peptides and all(p.parent == 0 for p in g.peptides)
        assert len(ref.proteins[0]) == g.cds[0].length_aa

    def test_joined_and_minus_strand_cds(self):
        ref = parse_genbank(HIV)
        by_name = {c.product: (i, c) for i, c in enumerate(ref.genome.cds)}
        tat_i, _ = by_name["Tat"]
        assert len(ref.parts[tat_i]) == 2                    # spliced
        assert any(c.strand == "-" for c in ref.genome.cds)  # antisense Asp

    def test_duplicate_mat_peptides_are_dropped(self):
        peps = parse_genbank(HIV).genome.peptides
        keys = [(p.start, p.end, p.strand) for p in peps]
        assert len(keys) == len(set(keys))


# ── Coordinate mapping ──────────────────────────────────────────────────────

def test_cds_aa_to_genome_plus_and_minus():
    assert cds_aa_to_genome([(100, 400, 1)], 1, 1, 10) == [100, 130]
    assert cds_aa_to_genome([(100, 400, -1)], 1, 1, 10) == [370, 400]


def test_orf_aa_to_contig_plus_and_minus():
    assert orf_aa_to_contig(10, 310, "+", 2, 5) == [13, 25]
    assert orf_aa_to_contig(10, 310, "-", 2, 5) == [295, 307]


def test_aligned_segments_split_across_exons():
    parts = [(0, 30, 1), (100, 130, 1)]        # 10 aa + 10 aa
    segs = aligned_segments(parts, 1, 6, 15, [500, 530], "+")
    assert segs == [[15, 30, 500, 515], [100, 115, 515, 530]]


def test_placement_detects_reverse_complement():
    hits = [SyntenyHit(orf_name="o", cds_index=0, pident=99, evalue=0, bit_score=100,
                       q_start=1, q_end=10, s_start=1, s_end=10, q_cov=100, s_cov=10,
                       orf_nt=[0, 30], ref_nt=[970, 1000])]
    flipped, offset = placement(hits, {"o": "-"}, ["+"])
    assert flipped and offset == 1000


# ── End to end (DIAMOND mocked) ─────────────────────────────────────────────

class _FakeAligner:
    """Returns one perfect HSP per ORF against the first CDS of its reference."""

    def __init__(self, rows=None):
        self.rows, self.queries, self.subjects = rows, None, None

    def align(self, queries, subjects, tag=""):
        if tag == "contigs":                      # contig ↔ contig search
            self.contig_queries = queries
            return []
        self.queries, self.subjects = queries, subjects
        if self.rows is not None:
            return self.rows
        if not subjects:
            return []
        sid = next(iter(subjects))
        return [{"qseqid": q, "sseqid": sid, "pident": 99.0, "length": len(aa),
                 "qstart": 1, "qend": len(aa), "sstart": 1, "send": len(aa),
                 "evalue": 1e-50, "bitscore": 200.0, "qlen": len(aa), "slen": 3392}
                for q, aa in queries.items()]


def _dengue_seq(rc=False):
    from Bio import SeqIO
    den = str(SeqIO.read(str(DENGUE), "genbank").seq)
    sq = den[94:394]
    s = NucSequence(id="c1", sequence=str(Seq(sq).reverse_complement()) if rc else sq)
    s.blastn_hits = [_hit()]
    s.orfs = [_orf(0, 300, "-" if rc else "+", aa="M" * 100, name="c1_orf")]
    return s


class TestSyntenyAnalyzer:

    def test_sets_synteny_and_returns_reference(self, tmp_path):
        s = _dengue_seq()
        refs = SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline),
                               aligner=_FakeAligner()).run([s])
        assert list(refs) == ["NC_001477.1"]
        assert s.synteny.reference == "NC_001477.1"
        assert len(s.synteny.hits) == 1
        assert s.synteny.hits[0].ref_nt == [94, 394]
        assert not s.synteny.flipped
        assert (tmp_path / "synteny.tsv").read_text().count("\n") == 2

    def test_reverse_complement_contig_is_flipped(self, tmp_path):
        s = _dengue_seq(rc=True)
        SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline),
                        aligner=_FakeAligner()).run([s])
        assert s.synteny.flipped

    def test_hits_to_another_reference_are_ignored(self, tmp_path):
        s = _dengue_seq()
        rows = [{"qseqid": "q0", "sseqid": "NC_001802|0", "pident": 90.0, "length": 50,
                 "qstart": 1, "qend": 50, "sstart": 1, "send": 50, "evalue": 1e-9,
                 "bitscore": 80.0, "qlen": 100, "slen": 1435}]
        SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline),
                        aligner=_FakeAligner(rows)).run([s])
        assert s.synteny is not None and s.synteny.hits == []

    def test_no_reference_available_leaves_sequence_untouched(self, tmp_path):
        s = _dengue_seq()
        s.blastn_hits = [_hit("XX_000001", "Some virus")]
        refs = SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(tmp_path / "refs", http=_offline),
                               aligner=_FakeAligner()).run([s])
        assert refs == {} and s.synteny is None

    def test_exporter_builds_clusters(self, tmp_path):
        a, b = _dengue_seq(), _dengue_seq()
        b.id = "c2"
        b.blastn_hits = [_hit()]
        for s in (a, b):
            s.is_viral = True
        refs = SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline),
                               aligner=_FakeAligner()).run([a, b])
        rep = ReportExporter(force=True).export([a, b], [], synteny_refs=refs)
        assert rep["synteny"]["clusters"] == [{"reference": "NC_001477.1", "members": ["c1", "c2"]}]
        assert rep["synteny"]["references"]["NC_001477.1"]["cds"][0]["product"] == "polyprotein"
        assert rep["sequences"][0]["synteny"]["reference"] == "NC_001477.1"

    def test_exporter_without_synteny_has_no_key(self):
        rep = ReportExporter(force=True).export([_dengue_seq()], [])
        assert "synteny" not in rep


def test_aligner_failure_returns_no_rows(tmp_path):
    class _P:
        returncode, stderr = 1, "boom"
    al = SyntenyAligner(tmp_path, run=lambda cmd, **kw: _P())
    assert al.align({"q0": "MK"}, {"r|0": "MK"}) == []


@pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond not installed")
def test_real_diamond_blastp(tmp_path):
    s = _dengue_seq()
    from viralquest.orfs import OrfAnalyzer
    s.orfs = []
    OrfAnalyzer(150).find_n_save_orfs(s)
    SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline)).run([s])
    assert s.synteny.hits and s.synteny.hits[0].pident > 95


@pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond not installed")
def test_real_diamond_links_overlapping_contigs(tmp_path):
    """Two contigs covering the same reference stretch get contig ↔ contig links."""
    from Bio import SeqIO
    from viralquest.orfs import OrfAnalyzer
    den = str(SeqIO.read(str(DENGUE), "genbank").seq)
    seqs = []
    for sid, (a, b) in {"a": (1000, 4000), "b": (2500, 5500)}.items():
        s = NucSequence(id=sid, sequence=den[a:b])
        s.blastn_hits = [_hit()]
        OrfAnalyzer(150).find_n_save_orfs(s)
        seqs.append(s)
    SyntenyAnalyzer(tmp_path, fetcher=GenBankFetcher(MOCK, http=_offline)).run(seqs)
    links = seqs[0].synteny.links
    assert links and all(l.other_seq == "b" for l in links)
    assert seqs[1].synteny.links == []          # stored once per pair
    assert max(l.pident for l in links) > 99
