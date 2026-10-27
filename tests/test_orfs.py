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


# ── --reload of reports older than partial-ORF calling ──────────────────────

def _old_orf(start, stop, strand, domains=()):
    return {"start_codon": "ATG", "stop_codon": "TAA", "start_position": start,
            "stop_position": stop, "strand": strand, "frame": 1, "length_aa": (stop - start) // 3,
            "length_nt": stop - start, "bigger_than_50": True, "aa_sequence": "", "nuc_sequence": "",
            "orf_type": "complete", "name": f"c_ORF_{start}_{stop}{strand}", "domains": list(domains)}


def test_needs_orf_recall_by_version():
    from viralquest.orfs import needs_orf_recall
    assert needs_orf_recall({"meta": {"viralquest_version": "3.0.1"}})
    assert not needs_orf_recall({"meta": {"viralquest_version": "3.0.9"}})
    assert not needs_orf_recall({"meta": {"viralquest_version": "3.0.1",
                                          "orfs_recalled_with": "3.0.14"}})


def test_recall_adds_partial_orfs_and_keeps_domains():
    from viralquest.orfs import recall_orfs
    # ATG + 100 codons + TAA after 60 open nt: the old caller saw only the complete ORF;
    # the 3'-partial frame (no ATG) and the reverse strand fragment are new.
    seq = {"id": "c", "sequence": "GCT" * 20 + "ATG" + "GCT" * 100 + "TAA"}
    dom = {"database": "Pfam", "target": "X", "start": 10, "stop": 50}
    seq["orfs"] = [_old_orf(60, 363, "+", [dom])]
    assert recall_orfs(seq, OrfAnalyzer(min_len_nt=150)) == 2
    kept = [o for o in seq["orfs"] if (o["start_position"], o["strand"]) == (60, "+")]
    assert [(d["start"], d["stop"]) for d in kept[0]["domains"]] == [(10, 50)]
    assert {o["orf_type"] for o in seq["orfs"]} == {"complete", "3-prime-partial", "5-prime-partial"}


def test_recall_keeps_old_orfs_it_cannot_match():
    from viralquest.orfs import recall_orfs
    seq = {"id": "c", "sequence": "A" * 30, "orfs": [_old_orf(0, 30, "+")]}
    assert recall_orfs(seq, OrfAnalyzer(min_len_nt=150)) == 0
    assert [o["name"] for o in seq["orfs"]] == ["c_ORF_0_30+"]
