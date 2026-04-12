from viralquest.biodata import NucSequence
from viralquest.orfs import OrfAnalyzer

def test_none_for_seqs_without_orfs():

    short_seq = NucSequence(id="virus_1", sequence="ATGCGTTAA")
    analyzer = OrfAnalyzer(min_len_nt=150)

    analyzer.find_n_save_orfs(short_seq)
    biggest_orf = analyzer.get_longest_orf(short_seq)

    assert len(short_seq.orfs) == 0
    assert biggest_orf is None
