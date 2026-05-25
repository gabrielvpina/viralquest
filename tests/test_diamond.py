
"""
tests/test_diamond.py
 
pytest test suite for viralquest/diamond.py.
 
Run from the project root:
    pytest tests/test_diamond.py -v
 
Covers:
  - is_viral_subject()               — regex matching
  - DiamondRunner._chunk()           — sequence batching
  - DiamondRunner._write_batch_fasta() — FASTA output
  - DiamondOutputParser.parse()      — TSV parsing (valid, empty, malformed)
  - DiamondResultAttacher.attach()   — hit attachment, is_viral flag,
                                       single-HSP coverage, merged coverage,
                                       NR viral filtering
  - DiamondFilterer.filter()         — quality thresholds
  - merge_query_intervals()          — interval merging edge cases
  - DiamondRunner.run_all()          — empty input guard
"""
 
import csv
from pathlib import Path
 
import pytest
 
from viralquest.biodata import BlastxResult, NucSequence, merge_query_intervals
from viralquest.diamond import (
    DiamondFilterer,
    DiamondOutputParser,
    DiamondPhase,
    DiamondResultAttacher,
    DiamondRunner,
    is_viral_subject,
)


# ===========================================================================
# Fixtures
# ===========================================================================

def make_hit(
    query_id="seq1",
    subject_id="prot1",
    subject_title="hypothetical protein [Tobacco mosaic virus]",
    pct_identity=95.0,
    aln_length=100,
    mismatches=5,
    gap_opens=0,
    query_start=1,
    query_end=300,
    subject_start=1,
    subject_end=100,
    e_value=1e-10,
    bit_score=200.0,
) -> BlastxResult:
    """Factory for BlastxResult with sensible defaults."""
    return BlastxResult(
        query_id=query_id,
        subject_id=subject_id,
        subject_title=subject_title,
        pct_identity=pct_identity,
        aln_length=aln_length,
        mismatches=mismatches,
        gap_opens=gap_opens,
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        e_value=e_value,
        bit_score=bit_score,
    )


def make_seq(seq_id="seq1", length=1000) -> NucSequence:
    return NucSequence(id=seq_id, sequence="A" * length)


@pytest.fixture
def tmp_tsv(tmp_path):
    """Returns a helper that writes TSV rows and gives back the path."""
    def _write(rows: list[list]) -> Path:
        p = tmp_path / "hits.tsv"
        with open(p, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerows(rows)
        return p
    return _write


# ===========================================================================
# is_viral_subject
# ===========================================================================

class TestIsViralSubject:

    @pytest.mark.parametrize("title", [
        "polyprotein [Tobacco mosaic virus]",
        "capsid protein [Tomato bushy stunt virus]",
        "hypothetical protein [Escherichia phage T4]",
        "coat protein [bacteriophage lambda]",
        "RNA-dependent RNA polymerase [Influenza A viral RNA]",
        "glycoprotein [Retroviridae sp.]",
        "movement protein [Virales unclassified]",
        "membrane protein [Viridae family member]",
    ])
    def test_viral_titles_match(self, title):
        assert is_viral_subject(title) is True

    @pytest.mark.parametrize("title", [
        "hypothetical protein [Homo sapiens]",
        "DNA polymerase [Escherichia coli]",
        "ribosomal protein [Arabidopsis thaliana]",
        "",
        "NaN",
        "no brackets here virus",           # 'virus' outside brackets — no match
        "uncharacterized protein [Mus musculus]",
    ])
    def test_non_viral_titles_no_match(self, title):
        assert is_viral_subject(title) is False

    def test_case_insensitive(self):
        assert is_viral_subject("protein [TOBACCO MOSAIC VIRUS]") is True
        assert is_viral_subject("protein [Tobacco Mosaic Virus]") is True

    def test_virus_outside_brackets_not_matched(self):
        # 'virus' must be inside square brackets to match
        assert is_viral_subject("virus-like particle [Homo sapiens]") is False


# ===========================================================================
# DiamondRunner._chunk
# ===========================================================================

class TestDiamondRunnerChunk:

    def _runner(self, tmp_path):
        return DiamondRunner(db_path="fake.dmnd", outdir=str(tmp_path))

    def test_exact_batch(self, tmp_path):
        seqs   = [make_seq(f"s{i}") for i in range(6)]
        chunks = DiamondRunner._chunk(seqs, 2)
        assert len(chunks) == 3
        assert all(len(c) == 2 for c in chunks)

    def test_partial_last_batch(self, tmp_path):
        seqs   = [make_seq(f"s{i}") for i in range(7)]
        chunks = DiamondRunner._chunk(seqs, 3)
        assert len(chunks) == 3
        assert len(chunks[-1]) == 1

    def test_single_sequence(self, tmp_path):
        seqs   = [make_seq("s0")]
        chunks = DiamondRunner._chunk(seqs, 500)
        assert len(chunks) == 1
        assert len(chunks[0]) == 1

    def test_empty_list(self, tmp_path):
        chunks = DiamondRunner._chunk([], 100)
        assert chunks == []

    def test_chunk_larger_than_list(self, tmp_path):
        seqs   = [make_seq(f"s{i}") for i in range(3)]
        chunks = DiamondRunner._chunk(seqs, 1000)
        assert len(chunks) == 1
        assert len(chunks[0]) == 3


# ===========================================================================
# DiamondRunner._write_batch_fasta
# ===========================================================================

class TestWriteBatchFasta:

    def test_fasta_format(self, tmp_path):
        runner = DiamondRunner(db_path="x", outdir=str(tmp_path))
        seqs   = [NucSequence(id="seq1", sequence="ATCG"),
                  NucSequence(id="seq2", sequence="GGCC")]
        fasta  = tmp_path / "out.fa"
        runner._write_batch_fasta(seqs, fasta)

        text = fasta.read_text()
        assert ">seq1\nATCG\n" in text
        assert ">seq2\nGGCC\n" in text

    def test_file_created(self, tmp_path):
        runner = DiamondRunner(db_path="x", outdir=str(tmp_path))
        fasta  = tmp_path / "out.fa"
        runner._write_batch_fasta([make_seq()], fasta)
        assert fasta.exists()

    def test_empty_batch_writes_empty_file(self, tmp_path):
        runner = DiamondRunner(db_path="x", outdir=str(tmp_path))
        fasta  = tmp_path / "out.fa"
        runner._write_batch_fasta([], fasta)
        assert fasta.read_text() == ""


# ===========================================================================
# DiamondOutputParser.parse
# ===========================================================================

class TestDiamondOutputParser:

    def _row(self, query_id="seq1", subject_id="prot1",
             stitle="protein [Virus X]", pident=90.0, length=100,
             mismatch=5, gapopen=0, qstart=1, qend=300,
             sstart=1, send=100, evalue=1e-10, bitscore=200.0):
        return [query_id, subject_id, stitle, pident, length,
                mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore]

    def test_parses_valid_row(self, tmp_tsv):
        path = tmp_tsv([self._row()])
        hits = DiamondOutputParser.parse(path)
        assert len(hits) == 1
        h = hits[0]
        assert h.query_id == "seq1"
        assert h.pct_identity == 90.0
        assert h.aln_length == 100
        assert h.e_value == 1e-10
        assert h.bit_score == 200.0

    def test_parses_multiple_rows(self, tmp_tsv):
        rows = [self._row(query_id=f"seq{i}") for i in range(5)]
        hits = DiamondOutputParser.parse(tmp_tsv(rows))
        assert len(hits) == 5

    def test_skips_short_rows(self, tmp_tsv):
        short_row = ["seq1", "prot1", "title"]          # only 3 columns
        full_row  = self._row()
        hits = DiamondOutputParser.parse(tmp_tsv([short_row, full_row]))
        assert len(hits) == 1

    def test_skips_malformed_numeric_fields(self, tmp_tsv):
        bad_row = self._row()
        bad_row[3] = "not_a_float"                      # pident field
        hits = DiamondOutputParser.parse(tmp_tsv([bad_row]))
        assert hits == []

    def test_empty_file_returns_empty_list(self, tmp_path):
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        assert DiamondOutputParser.parse(empty) == []

    def test_nonexistent_file_returns_empty_list(self, tmp_path):
        assert DiamondOutputParser.parse(tmp_path / "ghost.tsv") == []

    def test_subject_title_with_tabs_in_title_handled(self, tmp_tsv):
        # stitle can contain spaces but not tabs — ensure normal titles parse
        row = self._row(stitle="RNA polymerase [Tobacco mosaic virus] partial")
        hits = DiamondOutputParser.parse(tmp_tsv([row]))
        assert hits[0].subject_title == "RNA polymerase [Tobacco mosaic virus] partial"


# ===========================================================================
# DiamondResultAttacher  — REFSEQ_FILTER phase
# ===========================================================================

class TestAttacherRefseqPhase:

    def test_hit_attached_and_is_viral_set(self):
        seq  = make_seq("seq1", length=1000)
        hit  = make_hit(query_id="seq1", query_start=1, query_end=300, bit_score=200)
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.REFSEQ_FILTER)

        assert seq.is_viral is True
        assert len(seq.blastx_hits) == 1

    def test_unknown_query_id_skipped(self):
        seq = make_seq("seq1")
        hit = make_hit(query_id="unknown")
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.REFSEQ_FILTER)
        assert seq.blastx_hits == []
        assert seq.is_viral is False

    def test_returns_correct_attached_count(self):
        seqs = [make_seq(f"seq{i}", 500) for i in range(3)]
        hits = [make_hit(query_id=f"seq{i}", query_start=1, query_end=100) for i in range(3)]
        n = DiamondResultAttacher.attach(hits, seqs, DiamondPhase.REFSEQ_FILTER)
        assert n == 3

    def test_single_hsp_coverage_computed(self):
        seq = make_seq("seq1", length=1000)
        hit = make_hit(query_id="seq1", query_start=1, query_end=500)
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.REFSEQ_FILTER)
        # best hit gets merged coverage = (500-1)/1000 * 100
        best = seq.blastx_hits[0]
        assert best.query_coverage == pytest.approx(49.9, abs=0.1)

    def test_multiple_hits_same_query_all_attached(self):
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=1,   query_end=200, bit_score=100),
            make_hit(query_id="seq1", query_start=300,  query_end=500, bit_score=150),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        assert len(seq.blastx_hits) == 2


# ===========================================================================
# DiamondResultAttacher  — NR_CHARACTERIZE phase
# ===========================================================================

class TestAttacherNrPhase:

    def test_viral_subject_attached_to_nr_hits(self):
        seq = make_seq("seq1", length=1000)
        hit = make_hit(
            query_id="seq1",
            subject_title="polyprotein [Tobacco mosaic virus]",
            query_start=1, query_end=300,
        )
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.NR_CHARACTERIZE)
        assert len(seq.blastx_nr_hits) == 1

    def test_non_viral_subject_not_attached(self):
        seq = make_seq("seq1", length=1000)
        hit = make_hit(
            query_id="seq1",
            subject_title="ribosomal protein [Homo sapiens]",
            query_start=1, query_end=300,
        )
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.NR_CHARACTERIZE)
        assert seq.blastx_nr_hits == []

    def test_nr_phase_does_not_set_is_viral(self):
        seq = make_seq("seq1", length=1000)
        hit = make_hit(
            query_id="seq1",
            subject_title="capsid [Influenza A virus]",
            query_start=1, query_end=300,
        )
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.NR_CHARACTERIZE)
        # is_viral is only set by REFSEQ_FILTER
        assert seq.is_viral is False

    def test_mixed_viral_nonviral_only_viral_kept(self):
        seq = make_seq("seq1", length=2000)
        hits = [
            make_hit(query_id="seq1", subject_title="protein [Virus X]",
                     query_start=1, query_end=500, bit_score=300),
            make_hit(query_id="seq1", subject_title="protein [Homo sapiens]",
                     query_start=600, query_end=900, bit_score=100),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.NR_CHARACTERIZE)
        assert len(seq.blastx_nr_hits) == 1
        assert "Virus" in seq.blastx_nr_hits[0].subject_title


# ===========================================================================
# Merged coverage  (the core multi-HSP logic)
# ===========================================================================

class TestMergedCoverage:

    def test_non_overlapping_hsps_sum_correctly(self):
        """[10-50] + [200-300] = 40 + 100 = 140 bases → 14% of 1000"""
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=10,  query_end=50,  bit_score=100),
            make_hit(query_id="seq1", query_start=200, query_end=300, bit_score=200),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        best = max(seq.blastx_hits, key=lambda h: h.bit_score)
        assert best.query_coverage == pytest.approx(14.0, abs=0.1)

    def test_overlapping_hsps_not_double_counted(self):
        """[10-80] overlaps [30-120] → merged [10-120] = 110 bases → 11% of 1000"""
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=10, query_end=80,  bit_score=100),
            make_hit(query_id="seq1", query_start=30, query_end=120, bit_score=50),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        best = max(seq.blastx_hits, key=lambda h: h.bit_score)
        assert best.query_coverage == pytest.approx(11.0, abs=0.1)

    def test_single_hsp_coverage_equals_interval_span(self):
        seq = make_seq("seq1", length=500)
        hit = make_hit(query_id="seq1", query_start=0, query_end=250, bit_score=100)
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.REFSEQ_FILTER)
        assert seq.blastx_hits[0].query_coverage == pytest.approx(50.0, abs=0.1)

    def test_minus_strand_hit_normalised(self):
        """Minus strand: query_end < query_start — abs() must handle this."""
        seq = make_seq("seq1", length=1000)
        hit = make_hit(query_id="seq1", query_start=500, query_end=100, bit_score=100)
        DiamondResultAttacher.attach([hit], [seq], DiamondPhase.REFSEQ_FILTER)
        best = seq.blastx_hits[0]
        assert best.query_coverage == pytest.approx(40.0, abs=0.1)

    def test_best_hit_receives_merged_coverage(self):
        """Highest bit_score hit should carry the merged value."""
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=1,   query_end=200, bit_score=50),
            make_hit(query_id="seq1", query_start=400, query_end=600, bit_score=300),  # best
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        best = max(seq.blastx_hits, key=lambda h: h.bit_score)
        # merged = (200-1) + (600-400) = 199 + 200 = 399 → 39.9%
        assert best.query_coverage == pytest.approx(39.9, abs=0.1)

    def test_adjacent_intervals_merged(self):
        """Touching intervals [0-100] and [100-200] → single span of 200."""
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=0,   query_end=100, bit_score=100),
            make_hit(query_id="seq1", query_start=100, query_end=200, bit_score=50),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        best = max(seq.blastx_hits, key=lambda h: h.bit_score)
        assert best.query_coverage == pytest.approx(20.0, abs=0.1)

    def test_three_overlapping_intervals(self):
        """[0-100], [50-150], [120-200] → merged [0-200] = 200 bases → 20%"""
        seq  = make_seq("seq1", length=1000)
        hits = [
            make_hit(query_id="seq1", query_start=0,   query_end=100, bit_score=100),
            make_hit(query_id="seq1", query_start=50,  query_end=150, bit_score=80),
            make_hit(query_id="seq1", query_start=120, query_end=200, bit_score=60),
        ]
        DiamondResultAttacher.attach(hits, [seq], DiamondPhase.REFSEQ_FILTER)
        best = max(seq.blastx_hits, key=lambda h: h.bit_score)
        assert best.query_coverage == pytest.approx(20.0, abs=0.1)


# ===========================================================================
# merge_query_intervals  (unit tests on the helper itself)
# ===========================================================================

class TestMergeQueryIntervals:

    def _hit(self, qs, qe, bs=100):
        return make_hit(query_start=qs, query_end=qe, bit_score=bs)

    def test_empty_list(self):
        assert merge_query_intervals([]) == 0

    def test_single_interval(self):
        assert merge_query_intervals([self._hit(0, 100)]) == 100

    def test_two_non_overlapping(self):
        hits = [self._hit(0, 100), self._hit(200, 400)]
        assert merge_query_intervals(hits) == 300

    def test_two_overlapping(self):
        hits = [self._hit(0, 100), self._hit(50, 150)]
        assert merge_query_intervals(hits) == 150

    def test_one_contained_in_other(self):
        hits = [self._hit(0, 200), self._hit(50, 100)]
        assert merge_query_intervals(hits) == 200

    def test_identical_intervals(self):
        hits = [self._hit(10, 50), self._hit(10, 50)]
        assert merge_query_intervals(hits) == 40

    def test_minus_strand_normalised(self):
        # query_end < query_start on minus strand
        hits = [self._hit(300, 100)]   # span = 200
        assert merge_query_intervals(hits) == 200

    def test_three_intervals_partial_overlap(self):
        # [0-50], [30-80], [200-300] → [0-80] + [200-300] = 80 + 100 = 180
        hits = [self._hit(0, 50), self._hit(30, 80), self._hit(200, 300)]
        assert merge_query_intervals(hits) == 180


# ===========================================================================
# DiamondFilterer
# ===========================================================================

class TestDiamondFilterer:

    def _filtered_hit(self, pct_identity=90.0, e_value=1e-10, coverage=80.0):
        h = make_hit(pct_identity=pct_identity, e_value=e_value)
        h.query_coverage = coverage
        return h

    def test_all_pass_with_permissive_thresholds(self):
        hits = [self._filtered_hit() for _ in range(5)]
        out  = DiamondFilterer.filter(hits, min_identity=0, max_evalue=1, min_coverage=0)
        assert len(out) == 5

    def test_identity_threshold(self):
        hits = [
            self._filtered_hit(pct_identity=95.0),
            self._filtered_hit(pct_identity=25.0),   # below threshold
        ]
        out = DiamondFilterer.filter(hits, min_identity=30.0)
        assert len(out) == 1
        assert out[0].pct_identity == 95.0

    def test_evalue_threshold(self):
        hits = [
            self._filtered_hit(e_value=1e-20),
            self._filtered_hit(e_value=1e-3),         # above threshold
        ]
        out = DiamondFilterer.filter(hits, max_evalue=1e-5)
        assert len(out) == 1

    def test_coverage_threshold(self):
        hits = [
            self._filtered_hit(coverage=80.0),
            self._filtered_hit(coverage=10.0),        # below threshold
        ]
        out = DiamondFilterer.filter(hits, min_coverage=50.0)
        assert len(out) == 1

    def test_all_fail_returns_empty(self):
        hits = [self._filtered_hit(pct_identity=5.0) for _ in range(3)]
        out  = DiamondFilterer.filter(hits, min_identity=90.0)
        assert out == []

    def test_empty_input_returns_empty(self):
        assert DiamondFilterer.filter([]) == []

    def test_combined_thresholds(self):
        hits = [
            self._filtered_hit(pct_identity=95.0, e_value=1e-20, coverage=80.0),  # pass
            self._filtered_hit(pct_identity=20.0, e_value=1e-20, coverage=80.0),  # fail identity
            self._filtered_hit(pct_identity=95.0, e_value=1e-1,  coverage=80.0),  # fail evalue
            self._filtered_hit(pct_identity=95.0, e_value=1e-20, coverage=10.0),  # fail coverage
        ]
        out = DiamondFilterer.filter(hits, min_identity=30.0, max_evalue=1e-5, min_coverage=50.0)
        assert len(out) == 1


# ===========================================================================
# DiamondRunner.run_all  — empty input guard (no subprocess needed)
# ===========================================================================

class TestDiamondRunnerRunAll:

    def test_run_all_empty_sequences_returns_empty_list(self, tmp_path):
        runner = DiamondRunner(db_path="fake.dmnd", outdir=str(tmp_path))
        result = runner.run_all([], DiamondPhase.REFSEQ_FILTER)
        assert result == []

    def test_outdir_created_on_init(self, tmp_path):
        new_dir = tmp_path / "subdir" / "deep"
        DiamondRunner(db_path="x", outdir=str(new_dir))
        assert new_dir.exists()