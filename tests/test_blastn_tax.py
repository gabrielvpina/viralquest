import pytest
"""Tests for blastn_tax: taxonomy of the best BLASTn hit from the virus name in its title."""
from viralquest.biodata import BlastnResult, NucSequence
from viralquest.blastn_tax import BlastnTaxonomyResolver, is_viral_title, name_candidates


def _rec(taxid, name, species=None, family="Rhabdoviridae", genus=None, genome=None):
    return {"TaxId": taxid, "ScientificName": name, "Species": species, "Family": family,
            "Genus": genus, "No_Rank": None, "Clade": None, "Kingdom": None, "Phylum": None,
            "Class": None, "Order": None, "Subfamily": None, "Genome": genome}


RECORDS = [
    _rec(152177, "Dichorhavirus orchidaceae", genus="Dichorhavirus", genome="ssRNA(-)"),
    _rec(1213, "Dichorhavirus", genus="Dichorhavirus"),
    _rec(1922884, "Hubei diptera virus 3", species="Beidivirus muscae", family="Phenuiviridae"),
]


def _resolver():
    sci = {r["ScientificName"].lower(): r for r in RECORDS}
    sp = {r["Species"].lower(): r for r in RECORDS if r["Species"]}
    return BlastnTaxonomyResolver(sci, sp)


def _hit(title, acc="AB000001", bits=100.0):
    return BlastnResult(qseqid="s1", qlen=1000, slen=1000, qcovhsp=90, pident=97.0, evalue=1e-50,
                        bit_score=bits, stitle=title, accession=acc)


def _seq(*hits):
    s = NucSequence(id="s1", sequence="ACGT" * 50)
    s.blastn_hits.extend(hits)
    return s


class TestNameCandidates:
    def test_windows_around_the_anchor_longest_first(self):
        c = name_candidates("Hubei diptera virus 3 strain HB1 segment L, complete sequence")
        assert c.index("Hubei diptera virus 3") < c.index("diptera virus 3") < c.index("virus 3")
        assert "virus" not in [x.lower() for x in c]          # a bare generic word is never a candidate

    def test_prefixes_and_brackets(self):
        c = name_candidates("MAG: Some protein [Hubei diptera virus 3]")
        assert c[0] == "Hubei diptera virus 3"
        assert not any(x.startswith("MAG") for x in c)

    def test_flank_is_bounded(self):
        longest = name_candidates("a b c d e f g Dichorhavirus h i j k l m n", max_flank=5)[0].split()
        assert len(longest) == 11                                # 5 + anchor + 5

    def test_rank_suffixes_are_anchors(self):
        assert "Rhabdoviridae" in name_candidates("unclassified Rhabdoviridae sequence")


class TestResolver:
    def test_longest_name_wins(self):
        s = _seq(_hit("Dichorhavirus orchidaceae isolate S9 segment RNA1"))
        _resolver().annotate([s])
        bt = s.blastn_taxonomy
        assert (bt.status, bt.method, bt.matched_name, bt.taxid) == (
            "resolved", "title", "Dichorhavirus orchidaceae", 152177)
        assert bt.taxonomy.genome == "ssRNA(-)" and bt.accession == "AB000001"

    def test_species_name_index(self):
        bt = _resolver().resolve_hit(_hit("Beidivirus muscae isolate X"))
        assert bt.matched_name == "Beidivirus muscae" and bt.taxonomy.family == "Phenuiviridae"

    def test_genus_only(self):
        bt = _resolver().resolve_hit(_hit("Orchid fleck dichorhavirus isolate Smilax"))
        assert bt.matched_name.lower() == "dichorhavirus"

    def test_unresolved(self):
        s = _seq(_hit("PREDICTED: Homo sapiens actin, mRNA"))
        counts = _resolver().annotate([s])
        assert s.blastn_taxonomy.status == "unresolved" and s.blastn_taxonomy.taxonomy is None
        assert counts == {"resolved": 0, "unresolved": 1}

    def test_best_hit_by_bit_score(self):
        s = _seq(_hit("PREDICTED: Homo sapiens actin", bits=50),
                 _hit("Hubei diptera virus 3 segment S", acc="MK1", bits=400))
        _resolver().annotate([s])
        assert s.blastn_taxonomy.accession == "MK1"

    def test_no_hits_no_annotation(self):
        s = NucSequence(id="x", sequence="ACGT")
        _resolver().annotate([s])
        assert s.blastn_taxonomy is None


def test_exported_shape():
    from viralquest.exporter import ReportExporter
    d = ReportExporter._blastn_taxonomy_to_dict(_resolver().resolve_hit(_hit("Hubei diptera virus 3")))
    assert d["status"] == "resolved" and d["method"] == "title" and d["taxonomy"]["family"] == "Phenuiviridae"
    assert "class" in d["taxonomy"] and "class_" not in d["taxonomy"]
    un = ReportExporter._blastn_taxonomy_to_dict(_resolver().resolve_hit(_hit("PREDICTED: Homo sapiens")))
    assert un["status"] == "unresolved" and un["taxonomy"] is None


# ── is_viral_title ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("title, expected", [
    ("Dengue virus 2 isolate X, complete genome", True),
    ("MAG: Hubei diptera virus 3 strain Y", True),
    ("Escherichia phage T4, complete genome", True),
    ("polyprotein [Zika virus]", True),
    ("Homo sapiens chromosome 1, GRCh38.p14", False),
    ("Aedes aegypti mRNA for actin", False),
    ("", False),
])
def test_is_viral_title(title, expected):
    assert is_viral_title(title) is expected
