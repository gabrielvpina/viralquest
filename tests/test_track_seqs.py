import subprocess
from unittest.mock import patch

import pytest

from viralquest.biodata import BlastxResult, ClusterMember, NucSequence, ViralCluster
from viralquest.track_seqs import (
    BlastnAligner,
    ClusterBuilder,
    RepresentativeSelector,
    SequenceTracker,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_seq(seq_id: str, sequence: str, species: str = "") -> NucSequence:
    seq = NucSequence(id=seq_id, sequence=sequence)
    if species:
        hit = BlastxResult(
            query_id=seq_id, subject_id="NC_001",
            subject_title=f"polymerase [{species}]",
            pct_identity=95.0, aln_length=200, mismatches=5,
            gap_opens=0, query_start=1, query_end=200,
            subject_start=1, subject_end=200,
            e_value=1e-80, bit_score=300.0,
        )
        seq.blastx_hits.append(hit)
    return seq


def _hit(pident: float = 92.5, qcov: float = 85.0,
         sstart: int = 10, send: int = 350,
         qstart: int = 1,  qend: int = 300) -> dict:
    """A best-hit record as returned by BlastnAligner._parse_blastn."""
    return {"pident": pident, "qcov": qcov, "qstart": qstart, "qend": qend,
            "sstart": sstart, "send": send, "bitscore": 200.0}


def _align_all(rep, members):
    """_run_blastn stand-in: every member aligns well to the representative."""
    return {m.id: _hit() for m in members}


def _tabline(
    qid: str   = "q1",   sid: str  = "s1",
    pident: str = "95.0", aln: str = "200",
    qlen: str  = "200",  slen: str = "500",
    qstart: str = "1",   qend: str = "200",
    sstart: str = "10",  send: str = "209",
    evalue: str = "1e-50", bitscore: str = "300.0",
) -> str:
    return "\t".join([qid, sid, pident, aln, qlen, slen,
                      qstart, qend, sstart, send, evalue, bitscore])


# ---------------------------------------------------------------------------
# ClusterBuilder
# ---------------------------------------------------------------------------

class TestClusterBuilder:
    def test_groups_by_species(self):
        seqs = [
            make_seq("s1", "ATGC" * 100, "Influenza A virus"),
            make_seq("s2", "ATGC" * 80,  "Influenza A virus"),
            make_seq("s3", "ATGC" * 90,  "Tobacco mosaic virus"),
        ]
        groups = ClusterBuilder.build(seqs)
        assert set(groups) == {"Influenza A virus", "Tobacco mosaic virus"}
        assert len(groups["Influenza A virus"]) == 2
        assert len(groups["Tobacco mosaic virus"]) == 1

    def test_skips_seqs_without_species(self):
        seqs = [
            make_seq("s1", "ATGC" * 100, "Virus A"),
            make_seq("s2", "ATGC" * 80),   # no blastx hit → no species
        ]
        groups = ClusterBuilder.build(seqs)
        assert list(groups.keys()) == ["Virus A"]

    def test_sorted_by_species_name(self):
        seqs = [
            make_seq("s1", "ATGC" * 100, "Zymovirus"),
            make_seq("s2", "ATGC" * 80,  "Alphavirus"),
        ]
        groups = ClusterBuilder.build(seqs)
        assert list(groups.keys()) == ["Alphavirus", "Zymovirus"]

    def test_empty_input(self):
        assert ClusterBuilder.build([]) == {}

    def test_all_skipped_returns_empty(self):
        seqs = [make_seq("s1", "ATGC"), make_seq("s2", "TTTT")]
        assert ClusterBuilder.build(seqs) == {}

    def test_nr_hits_preferred_for_species(self):
        seq = make_seq("s1", "ATGC" * 100)
        seq.blastx_hits.append(BlastxResult(
            query_id="s1", subject_id="r1",
            subject_title="protein [Refseq Virus]",
            pct_identity=90.0, aln_length=100, mismatches=10,
            gap_opens=0, query_start=1, query_end=100,
            subject_start=1, subject_end=100,
            e_value=1e-40, bit_score=150.0,
        ))
        seq.blastx_nr_hits.append(BlastxResult(
            query_id="s1", subject_id="n1",
            subject_title="protein [NR Virus]",
            pct_identity=95.0, aln_length=100, mismatches=5,
            gap_opens=0, query_start=1, query_end=100,
            subject_start=1, subject_end=100,
            e_value=1e-60, bit_score=300.0,
        ))
        groups = ClusterBuilder.build([seq])
        assert "NR Virus" in groups


# ---------------------------------------------------------------------------
# RepresentativeSelector
# ---------------------------------------------------------------------------

class TestRepresentativeSelector:
    def test_selects_longest(self):
        seqs = [
            make_seq("short",  "ATGC" * 10),
            make_seq("long",   "ATGC" * 100),
            make_seq("medium", "ATGC" * 50),
        ]
        assert RepresentativeSelector.select(seqs).id == "long"

    def test_single_seq_is_representative(self):
        seq = make_seq("only", "ATGCATGC")
        assert RepresentativeSelector.select([seq]).id == "only"

    def test_first_longest_wins_on_tie(self):
        s1 = make_seq("s1", "ATGC" * 100)
        s2 = make_seq("s2", "ATGC" * 100)
        assert RepresentativeSelector.select([s1, s2]).id == "s1"


# ---------------------------------------------------------------------------
# BlastnAligner._parse_blastn  (static — tested directly)
# ---------------------------------------------------------------------------

class TestParseBlastn:
    def test_basic_parse(self):
        result = BlastnAligner._parse_blastn(_tabline())
        assert "q1" in result
        hit = result["q1"]
        assert hit["pident"] == 95.0
        assert hit["qcov"]   == 100.0   # aln_len=200, qlen=200

    def test_best_bitscore_kept(self):
        lines = "\n".join([
            _tabline(bitscore="100.0"),
            _tabline(bitscore="300.0"),
            _tabline(bitscore="200.0"),
        ])
        result = BlastnAligner._parse_blastn(lines)
        assert result["q1"]["bitscore"] == 300.0

    def test_minus_strand_normalised(self):
        result = BlastnAligner._parse_blastn(_tabline(sstart="209", send="10"))
        assert result["q1"]["sstart"] == 10
        assert result["q1"]["send"]   == 209

    def test_plus_strand_unchanged(self):
        result = BlastnAligner._parse_blastn(_tabline(sstart="10", send="209"))
        assert result["q1"]["sstart"] == 10
        assert result["q1"]["send"]   == 209

    def test_empty_input_returns_empty(self):
        assert BlastnAligner._parse_blastn("") == {}

    def test_whitespace_only_input(self):
        assert BlastnAligner._parse_blastn("   \n  \n") == {}

    def test_short_line_skipped(self):
        assert BlastnAligner._parse_blastn("q1\ts1\t95.0") == {}

    def test_multiple_queries(self):
        lines = "\n".join([
            _tabline(qid="q1", bitscore="200.0"),
            _tabline(qid="q2", bitscore="500.0"),
        ])
        result = BlastnAligner._parse_blastn(lines)
        assert "q1" in result and "q2" in result

    def test_qcov_calculation(self):
        # qstart=1, qend=100, qlen=200 → aln_len=abs(100-1)+1=100, qcov=50.0
        result = BlastnAligner._parse_blastn(
            _tabline(qstart="1", qend="100", qlen="200")
        )
        assert result["q1"]["qcov"] == 50.0

    def test_pident_rounded(self):
        result = BlastnAligner._parse_blastn(_tabline(pident="95.1234"))
        assert result["q1"]["pident"] == 95.12


# ---------------------------------------------------------------------------
# BlastnAligner.align  (mocks _run_blastn)
# ---------------------------------------------------------------------------

class TestBlastnAlignerAlign:
    def _aligner(self) -> BlastnAligner:
        return BlastnAligner()

    def test_representative_is_first_and_marked(self):
        rep = make_seq("rep", "ATGC" * 100)
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", return_value={}):
            members = aligner.align(rep, [rep])
        assert len(members) == 1
        assert members[0].is_representative
        assert members[0].seq_id == "rep"

    def test_representative_identity_and_coverage(self):
        rep = make_seq("rep", "ATGC" * 100)
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", return_value={}):
            members = aligner.align(rep, [rep])
        assert members[0].identity       == 100.0
        assert members[0].query_coverage == 100.0

    def test_representative_aln_spans_full_length(self):
        rep = make_seq("rep", "ATGC" * 50)   # 200 nt
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", return_value={}):
            members = aligner.align(rep, [rep])
        assert members[0].aln_start == 1
        assert members[0].aln_end   == 200

    def test_member_with_hit_gets_alignment_data(self):
        rep = make_seq("rep",  "ATGC" * 100)
        mem = make_seq("mem1", "ATGC" * 80)
        hit_map = {"mem1": _hit(pident=92.5, qcov=85.0, sstart=10, send=350,
                                qstart=5, qend=300)}
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", return_value=hit_map):
            members = aligner.align(rep, [rep, mem])
        m = members[1]
        assert m.seq_id         == "mem1"
        assert m.identity       == 92.5
        assert m.query_coverage == 85.0
        assert m.aln_start      == 10
        assert m.aln_end        == 350
        assert m.query_start    == 5
        assert m.query_end      == 300
        assert not m.is_representative

    def test_member_without_hit_is_dropped(self):
        rep = make_seq("rep",     "ATGC" * 100)
        mem = make_seq("nomatch", "TTTT" * 80)
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", return_value={}):
            members = aligner.align(rep, [rep, mem])
        assert [m.seq_id for m in members] == ["rep"]

    def test_member_below_min_coverage_is_dropped(self):
        rep = make_seq("rep", "ATGC" * 100)
        mem = make_seq("low", "ATGC" * 80)
        aligner = BlastnAligner(min_coverage=50.0)
        with patch.object(aligner, "_run_blastn", return_value={"low": _hit(qcov=30.0)}):
            members = aligner.align(rep, [rep, mem])
        assert [m.seq_id for m in members] == ["rep"]

    def test_single_seq_skips_blastn(self):
        rep = make_seq("rep", "ATGC" * 100)
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn") as mock_run:
            aligner.align(rep, [rep])
        mock_run.assert_not_called()

    def test_result_length_equals_aligned_members_count(self):
        rep  = make_seq("rep", "ATGC" * 100)
        mem1 = make_seq("m1",  "ATGC" * 80)
        mem2 = make_seq("m2",  "ATGC" * 60)
        aligner = self._aligner()
        with patch.object(aligner, "_run_blastn", side_effect=_align_all):
            members = aligner.align(rep, [rep, mem1, mem2])
        assert len(members) == 3


# ---------------------------------------------------------------------------
# BlastnAligner._run_blastn  (error paths)
# ---------------------------------------------------------------------------

class TestBlastnAlignerRunBlastn:
    def test_binary_not_found_returns_empty(self):
        aligner = BlastnAligner(blastn_bin="no_such_binary_xyz_123")
        rep = make_seq("rep", "ATGCATGC")
        mem = make_seq("mem", "ATGCATGC")
        result = aligner._run_blastn(rep, [mem])
        assert result == {}

    def test_subprocess_error_returns_empty(self):
        aligner = BlastnAligner()
        rep = make_seq("rep", "ATGCATGC")
        mem = make_seq("mem", "ATGCATGC")
        with patch("viralquest.track_seqs.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                1, "blastn", stderr="BLAST error"
            )
            result = aligner._run_blastn(rep, [mem])
        assert result == {}


# ---------------------------------------------------------------------------
# SequenceTracker.track  (aligner mocked)
# ---------------------------------------------------------------------------

class TestSequenceTracker:
    """A cluster needs >= 2 sequences of the same species that align to the
    representative; _run_blastn is mocked so every member aligns."""

    def _seqs(self):
        s1 = make_seq("s1", "ATGC" * 100, "Influenza A virus")     # 400 nt
        s2 = make_seq("s2", "ATGC" * 80,  "Influenza A virus")     # 320 nt
        s3 = make_seq("s3", "ATGC" * 90,  "Tobacco mosaic virus")
        s4 = make_seq("s4", "ATGC" * 70,  "Tobacco mosaic virus")
        return s1, s2, s3, s4

    def _track(self, seqs, run_blastn=_align_all, **kw):
        tracker = SequenceTracker(**kw)
        with patch.object(tracker._aligner, "_run_blastn", side_effect=run_blastn):
            return tracker.track(list(seqs))

    def test_returns_list_of_viral_clusters(self):
        clusters = self._track(self._seqs())
        assert len(clusters) == 2
        assert all(isinstance(c, ViralCluster) for c in clusters)

    def test_cluster_id_format(self):
        for c in self._track(self._seqs()):
            assert c.cluster_id.startswith("VQ_CLU_")

    def test_cluster_ids_are_unique(self):
        ids = [c.cluster_id for c in self._track(self._seqs())]
        assert len(ids) == len(set(ids))

    def test_cluster_id_assigned_to_nuc_seqs(self):
        seqs = self._seqs()
        self._track(seqs)
        assert all(s.cluster_id is not None for s in seqs)

    def test_same_species_gets_same_cluster_id(self):
        s1, s2, s3, s4 = self._seqs()
        self._track([s1, s2, s3, s4])
        assert s1.cluster_id == s2.cluster_id
        assert s3.cluster_id == s4.cluster_id
        assert s1.cluster_id != s3.cluster_id

    def test_representative_is_longest_in_cluster(self):
        s1, s2, _, _ = self._seqs()
        clusters = self._track([s2, s1])
        # s1 = 400 nt, s2 = 320 nt
        assert clusters[0].representative_id == "s1"

    def test_cluster_species_matches_blastx(self):
        species_set = {c.species for c in self._track(self._seqs())}
        assert species_set == {"Influenza A virus", "Tobacco mosaic virus"}

    def test_cluster_size_property(self):
        clusters = self._track(self._seqs())
        flu = next(c for c in clusters if c.species == "Influenza A virus")
        assert flu.size == 2

    def test_single_sequence_species_not_clustered(self):
        s1, s2, s3, _ = self._seqs()      # s3 is the only TMV sequence
        clusters = self._track([s1, s2, s3])
        assert [c.species for c in clusters] == ["Influenza A virus"]
        assert s3.cluster_id is None

    def test_no_aligned_member_means_no_cluster(self):
        s1, s2, _, _ = self._seqs()
        clusters = self._track([s1, s2], run_blastn=lambda rep, members: {})
        assert clusters == []
        assert s1.cluster_id is None and s2.cluster_id is None

    def test_unaligned_member_left_out_of_cluster(self):
        s1, s2, _, _ = self._seqs()
        s5 = make_seq("s5", "ATGC" * 60, "Influenza A virus")
        only_s2 = lambda rep, members: {"s2": _hit()}
        clusters = self._track([s1, s2, s5], run_blastn=only_s2)
        assert {m.seq_id for m in clusters[0].members} == {"s1", "s2"}
        assert s5.cluster_id is None

    def test_seqs_without_species_excluded_from_clusters(self):
        no_hit = make_seq("no_hit", "ATGC" * 50)
        a      = make_seq("a",      "ATGC" * 50, "Virus X")
        b      = make_seq("b",      "ATGC" * 40, "Virus X")
        clusters = self._track([no_hit, a, b])
        assert len(clusters) == 1
        assert no_hit.cluster_id is None

    def test_empty_input_returns_empty(self):
        assert SequenceTracker().track([]) == []

    def test_clusters_ordered_alphabetically_by_species(self):
        seqs = [
            make_seq("z1", "ATGC" * 100, "Zymovirus X"),
            make_seq("z2", "ATGC" * 90,  "Zymovirus X"),
            make_seq("a1", "ATGC" * 80,  "Alphavirus Y"),
            make_seq("a2", "ATGC" * 70,  "Alphavirus Y"),
        ]
        clusters = self._track(seqs)
        assert [c.species for c in clusters] == ["Alphavirus Y", "Zymovirus X"]

    def test_custom_min_identity_passed_to_aligner(self):
        tracker = SequenceTracker(min_identity=90.0)
        assert tracker._aligner.min_identity == 90.0

    def test_custom_blastn_bin_passed_to_aligner(self):
        tracker = SequenceTracker(blastn_bin="/usr/local/bin/blastn")
        assert tracker._aligner.blastn_bin == "/usr/local/bin/blastn"
