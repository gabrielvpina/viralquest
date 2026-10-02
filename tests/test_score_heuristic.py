"""Tests for score_heuristic: classification and novelty tiers."""
import pytest

from viralquest.biodata import (
    NOVELTY_FLAGS,
    NOVELTY_TIERS,
    BlastnResult,
    BlastxResult,
    NucSequence,
)
from viralquest.score_heuristic import HeuristicScorer, _novelty


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bx(identity: float, coverage: float = 90.0, title: str = "polymerase [Some virus]") -> BlastxResult:
    h = BlastxResult(
        query_id="s1", subject_id="YP_000001.1", subject_title=title,
        pct_identity=identity, aln_length=300, mismatches=0, gap_opens=0,
        query_start=1, query_end=900, subject_start=1, subject_end=300,
        e_value=1e-50, bit_score=300.0,
    )
    h.query_coverage = coverage
    return h


def bn(identity: float, coverage: int = 90, title: str = "Some virus segment 1, complete") -> BlastnResult:
    return BlastnResult(
        qseqid="s1", qlen=1000, slen=1000, qcovhsp=coverage, pident=identity,
        evalue=1e-100, bit_score=500.0, stitle=title,
    )


def seq_with(blastx=None, blastn=None) -> NucSequence:
    s = NucSequence(id="s1", sequence="ATGC" * 300)
    s.is_viral = True
    if blastx is not None:
        s.blastx_nr_hits.append(blastx)
    if blastn is not None:
        s.blastn_hits.append(blastn)
    return s


# ---------------------------------------------------------------------------
# _novelty — tier boundaries
# ---------------------------------------------------------------------------

class TestNoveltyTiers:
    @pytest.mark.parametrize("nt, cov, tier", [
        (99.0, 90, "known"),
        (95.0, 70, "known"),          # both thresholds inclusive
        (94.9, 90, "variant"),        # just below species-level nt
        (85.0, 70, "variant"),
    ])
    def test_nucleotide_tiers(self, nt, cov, tier):
        assert _novelty("viral-known", bn(nt, cov), None)[0] == tier

    @pytest.mark.parametrize("aa, tier", [
        (99.0, "variant"),            # aa alone never reaches "known"
        (90.0, "variant"),
        (89.9, "novel-species"),
        (70.0, "novel-species"),
        (69.9, "divergent"),
        (40.0, "divergent"),
        (39.9, "highly-divergent"),
    ])
    def test_amino_acid_tiers(self, aa, tier):
        assert _novelty("viral-unknown", None, bx(aa))[0] == tier

    def test_known_needs_nucleotide_coverage(self):
        # 99% nt over 40% of the contig is not a species-level call.
        assert _novelty("viral-unknown", bn(99.0, 40), bx(50.0))[0] == "divergent"

    def test_high_aa_low_coverage_is_flagged_variant(self):
        # Long contigs of large DNA viruses: one protein covers little of the contig.
        assert _novelty("viral-unknown", None, bx(96.0, coverage=12.0)) == ("variant", ["low-coverage"])

    def test_hmm_only(self):
        assert _novelty("viral-unknown", None, None) == ("highly-divergent", ["hmm-only"])

    def test_nucleotide_only(self):
        assert _novelty("viral-unknown", bn(78.0, 50), None) == ("novel-species", ["nt-only"])
        assert _novelty("viral-unknown", bn(65.0, 50), None) == ("divergent", ["nt-only"])
        tier, flags = _novelty("viral-unknown", bn(97.0, 30), None)
        assert tier == "novel-species" and flags == ["nt-only", "low-coverage"]

    def test_non_viral_classification_wins(self):
        assert _novelty("non-viral", bn(99.0, 99), bx(99.0)) == ("non-viral", [])

    def test_tiers_and_flags_are_declared(self):
        cases = [("viral-unknown", None, bx(a)) for a in (95, 80, 50, 20)] + [("viral-unknown", None, None)]
        for args in cases:
            tier, flags = _novelty(*args)
            assert tier in NOVELTY_TIERS
            assert set(flags) <= set(NOVELTY_FLAGS)


# ---------------------------------------------------------------------------
# HeuristicScorer — novelty on the output
# ---------------------------------------------------------------------------

class TestScorerNovelty:
    def test_novelty_attached_to_output(self):
        out = HeuristicScorer().score([seq_with(blastx=bx(78.0))])[0]
        assert out.classification == "viral-unknown"
        assert out.novelty == "novel-species" and out.novelty_flags == []

    def test_known_classification_with_variant_novelty(self):
        # >=90% aa is species-level for the classification (viral-known) but,
        # without a close nucleotide hit, only a "variant" for novelty.
        out = HeuristicScorer().score([seq_with(blastx=bx(97.0))])[0]
        assert out.classification == "viral-known" and out.novelty == "variant"

    def test_blastn_known(self):
        out = HeuristicScorer().score([seq_with(blastx=bx(98.0), blastn=bn(99.0))])[0]
        assert out.classification == "viral-known" and out.novelty == "known"

    def test_non_viral_blastn_host_hit(self):
        host = bn(99.0, 95, title="Homo sapiens chromosome 3, complete sequence")
        out = HeuristicScorer().score([seq_with(blastx=bx(60.0), blastn=host)])[0]
        assert out.classification == "non-viral" and out.novelty == "non-viral"

    def test_analysis_mentions_novelty(self):
        out = HeuristicScorer().score([seq_with(blastx=bx(96.0, coverage=10.0))])[0]
        assert "Novelty: variant (low-coverage)." in out.analysis
