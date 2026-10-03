from viralquest.biodata import NucSequence
from viralquest.orfs import OrfAnalyzer

def test_none_for_seqs_without_orfs():

    short_seq = NucSequence(id="virus_1", sequence="ATGCGTTAA")
    analyzer = OrfAnalyzer(min_len_nt=150)

    analyzer.find_n_save_orfs(short_seq)
    biggest_orf = analyzer.get_longest_orf(short_seq)

    assert len(short_seq.orfs) == 0
    assert biggest_orf is None


def test_partial_orfs_at_contig_ends_are_found():
    """A gene fragment with no start and no stop inside the contig (common in
    assemblies) is still called, labelled by orfipy as a partial ORF."""
    frag = NucSequence(id="frag", sequence="GCT" * 200)        # 600 nt, no ATG, no stop
    OrfAnalyzer(min_len_nt=150).find_n_save_orfs(frag)
    assert frag.orfs, "partial ORFs must be called"
    assert any("partial" in o.orf_type or o.orf_type == "no_start_no_stop" for o in frag.orfs)
    assert all(o.length_nt % 3 == 0 and o.aa_sequence for o in frag.orfs)


def test_complete_orf_still_called():
    s = NucSequence(id="c", sequence="ATG" + "GCT" * 60 + "TAA")
    OrfAnalyzer(min_len_nt=150).find_n_save_orfs(s)
    assert any(o.orf_type == "complete" for o in s.orfs)
