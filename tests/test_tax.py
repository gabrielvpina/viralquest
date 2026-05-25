import json
import lzma

import pytest

from viralquest.biodata import BlastxResult, NucSequence, Taxonomy, ViralFamilyInfo
from viralquest.tax import (
    TaxonomyAnnotator,
    TaxonomyLoader,
    TaxonomyMatcher,
    ViralFamilyAnnotator,
    ViralFamilyLoader,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_blastx(subject_title: str, bit_score: float = 200.0) -> BlastxResult:
    return BlastxResult(
        query_id="q1",
        subject_id="s1",
        subject_title=subject_title,
        pct_identity=95.0,
        aln_length=100,
        mismatches=5,
        gap_opens=0,
        query_start=1,
        query_end=100,
        subject_start=1,
        subject_end=100,
        e_value=1e-50,
        bit_score=bit_score,
    )


def make_nuc_seq(seq_id: str = "seq1") -> NucSequence:
    return NucSequence(id=seq_id, sequence="ATGCGT")


def make_taxonomy(
    family: str | None = "Orthomyxoviridae",
    genus:  str | None = "Alphainfluenzavirus",
    order:  str | None = "Articulavirales",
) -> Taxonomy:
    return Taxonomy(
        tax_id=1,
        scientific_name="Influenza A virus",
        no_rank=None,
        clade="Riboviria",
        kingdom="Orthornavirae",
        phylum="Negarnaviricota",
        class_="Insthoviricetes",
        order=order,
        family=family,
        subfamily=None,
        genus=genus,
        species="Alphainfluenzavirus influenzae",
        genome="ssRNA(-)",
    )


def write_tax_json(tmp_path, records: list[dict]) -> str:
    path = tmp_path / "viralTax.json.xz"
    with lzma.open(path, "wt", encoding="utf-8") as fh:
        json.dump(records, fh)
    return str(path)


def write_family_json(tmp_path, records: list[dict], filename: str) -> str:
    path = tmp_path / filename
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(records, fh)
    return str(path)


# ---------------------------------------------------------------------------
# Sample data
# ---------------------------------------------------------------------------

SAMPLE_TAX_RECORDS = [
    {
        "TaxId": 1,
        "ScientificName": "Influenza A virus",
        "No_Rank": None,
        "Clade": "Riboviria",
        "Kingdom": "Orthornavirae",
        "Phylum": "Negarnaviricota",
        "Class": "Insthoviricetes",
        "Order": "Articulavirales",
        "Family": "Orthomyxoviridae",
        "Subfamily": None,
        "Genus": "Alphainfluenzavirus",
        "Species": "Alphainfluenzavirus influenzae",
        "Genome": "ssRNA(-)",
    },
    {
        "TaxId": 2,
        "ScientificName": "Tobacco mosaic virus",
        "No_Rank": "unclassified viruses",
        "Clade": "Riboviria",
        "Kingdom": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": "Virgaviridae",
        "Subfamily": None,
        "Genus": "Tobamovirus",
        "Species": "Tobamovirus virginianum",
        "Genome": "ssRNA(+)",
    },
    {
        "TaxId": 3,
        "ScientificName": "Alphatorquevirus zetaone",
        "No_Rank": None,
        "Clade": None,
        "Kingdom": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": None,
        "Subfamily": None,
        "Genus": None,
        "Species": None,
        "Genome": None,
    },
    {
        "TaxId": 4,
        "ScientificName": "Hepatitis B virus genotype A",
        "No_Rank": None,
        "Clade": "Riboviria",
        "Kingdom": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": "Hepadnaviridae",
        "Subfamily": None,
        "Genus": "Orthohepadnavirus",
        "Species": "Hepatitis B virus",
        "Genome": "dsDNA-RT",
    },
    {
        "TaxId": 5,
        "ScientificName": "Hepatitis B virus genotype B",
        "No_Rank": None,
        "Clade": "Riboviria",
        "Kingdom": None,
        "Phylum": None,
        "Class": None,
        "Order": None,
        "Family": "Hepadnaviridae",
        "Subfamily": None,
        "Genus": "Orthohepadnavirus",
        "Species": "Hepatitis B virus",   # same Species as TaxId=4 — TaxId=4 wins
        "Genome": "dsDNA-RT",
    },
]

# Shared family-info records used across ViralFamily tests
HIGH_RECORDS = [
    {"source": "ICTV",      "type": "Family", "name": "Orthomyxoviridae", "info": "Orthomyxo HIGH text"},
    {"source": "ICTV",      "type": "Family", "name": "Virgaviridae",     "info": "Virgaviridae HIGH text"},
    {"source": "ICTV",      "type": "Genus",  "name": "Tobamovirus",      "info": "Tobamovirus HIGH text"},
    {"source": "ViralZone", "type": "Order",  "name": "Articulavirales",  "info": "Articulavirales HIGH text"},
]

LOW_RECORDS = [
    {"source": "ICTV",      "type": "Family", "name": "Orthomyxoviridae", "info": "Orthomyxo LOW text"},
    {"source": "ICTV",      "type": "Family", "name": "Virgaviridae",     "info": "Virgaviridae LOW text"},
    {"source": "ICTV",      "type": "Genus",  "name": "Tobamovirus",      "info": "Tobamovirus LOW text"},
    {"source": "ViralZone", "type": "Order",  "name": "Articulavirales",  "info": "Articulavirales LOW text"},
]


# ---------------------------------------------------------------------------
# BlastxResult.species extraction  (biodata, tested here for proximity)
# ---------------------------------------------------------------------------

class TestBlastxSpeciesExtraction:
    def test_standard_bracket_format(self):
        hit = make_blastx("RNA-dependent RNA polymerase [Influenza A virus]")
        assert hit.species == "Influenza A virus"

    def test_bracket_with_strain_inside(self):
        hit = make_blastx("coat protein [Tobacco mosaic virus (strain U1)]")
        assert hit.species == "Tobacco mosaic virus (strain U1)"

    def test_no_brackets_returns_empty(self):
        hit = make_blastx("hypothetical protein")
        assert hit.species == ""

    def test_trailing_whitespace_ignored(self):
        hit = make_blastx("polymerase [Virus X]   ")
        assert hit.species == "Virus X"

    def test_only_last_bracket_used(self):
        hit = make_blastx("protein [region] functional [Real Virus]")
        assert hit.species == "Real Virus"

    def test_empty_subject_title(self):
        hit = make_blastx("")
        assert hit.species == ""


# ---------------------------------------------------------------------------
# TaxonomyLoader
# ---------------------------------------------------------------------------

class TestTaxonomyLoader:
    def test_sci_idx_keyed_by_lowercase(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        sci, _ = TaxonomyLoader.load(path)
        assert "influenza a virus" in sci
        assert "tobacco mosaic virus" in sci

    def test_sp_idx_keyed_by_lowercase_species_field(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        assert "alphainfluenzavirus influenzae" in sp
        assert "tobamovirus virginianum" in sp

    def test_records_without_species_field_excluded_from_sp_idx(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        assert None not in sp

    def test_duplicate_species_field_first_record_wins(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        assert sp["hepatitis b virus"]["TaxId"] == 4

    def test_missing_file_returns_empty_dicts(self):
        sci, sp = TaxonomyLoader.load("/nonexistent/path.json.xz")
        assert sci == {}
        assert sp == {}

    def test_empty_records_list(self, tmp_path):
        path = write_tax_json(tmp_path, [])
        sci, sp = TaxonomyLoader.load(path)
        assert sci == {}
        assert sp == {}


# ---------------------------------------------------------------------------
# TaxonomyMatcher
# ---------------------------------------------------------------------------

class TestTaxonomyMatcher:
    @pytest.fixture
    def matcher(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        sci, sp = TaxonomyLoader.load(path)
        return TaxonomyMatcher(sci, sp)

    def test_exact_scientific_name_match(self, matcher):
        tax = matcher.match("Influenza A virus")
        assert tax is not None
        assert tax.tax_id == 1
        assert tax.family == "Orthomyxoviridae"

    def test_match_is_case_insensitive(self, matcher):
        assert matcher.match("INFLUENZA A VIRUS") is not None
        assert matcher.match("influenza a virus") is not None

    def test_species_field_fallback(self, matcher):
        tax = matcher.match("Tobamovirus virginianum")
        assert tax is not None
        assert tax.tax_id == 2

    def test_unknown_species_returns_none(self, matcher):
        assert matcher.match("Unknown virus XYZ") is None

    def test_empty_string_returns_none(self, matcher):
        assert matcher.match("") is None

    def test_taxonomy_fields_populated_correctly(self, matcher):
        tax = matcher.match("Influenza A virus")
        assert isinstance(tax, Taxonomy)
        assert tax.scientific_name == "Influenza A virus"
        assert tax.clade == "Riboviria"
        assert tax.kingdom == "Orthornavirae"
        assert tax.phylum == "Negarnaviricota"
        assert tax.class_ == "Insthoviricetes"
        assert tax.order == "Articulavirales"
        assert tax.genus == "Alphainfluenzavirus"
        assert tax.species == "Alphainfluenzavirus influenzae"
        assert tax.genome == "ssRNA(-)"

    def test_null_fields_are_none(self, matcher):
        tax = matcher.match("Alphatorquevirus zetaone")
        assert tax is not None
        assert tax.family is None
        assert tax.genus is None
        assert tax.species is None

    def test_scientific_name_takes_priority_over_species_field(self, matcher):
        tax_sci = matcher.match("Hepatitis B virus genotype A")
        tax_sp  = matcher.match("Hepatitis B virus")
        assert tax_sci.tax_id == 4
        assert tax_sp.tax_id == 4


# ---------------------------------------------------------------------------
# TaxonomyAnnotator
# ---------------------------------------------------------------------------

class TestTaxonomyAnnotator:
    @pytest.fixture
    def annotator(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_TAX_RECORDS)
        return TaxonomyAnnotator(str(path))

    def test_annotates_sequence_from_best_blastx(self, annotator):
        seq = make_nuc_seq()
        seq.blastx_hits.append(make_blastx("RNA polymerase [Influenza A virus]"))
        count = annotator.annotate([seq])
        assert count == 1
        assert seq.taxonomy is not None
        assert seq.taxonomy.tax_id == 1

    def test_uses_best_blastx_not_first(self, annotator):
        seq = make_nuc_seq()
        seq.blastx_hits.append(make_blastx("protein [Unknown virus XYZ]", bit_score=100.0))
        seq.blastx_hits.append(make_blastx("coat protein [Tobacco mosaic virus]", bit_score=500.0))
        annotator.annotate([seq])
        assert seq.taxonomy.tax_id == 2

    def test_no_blastx_hits_not_annotated(self, annotator):
        seq = make_nuc_seq()
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.taxonomy is None

    def test_unresolvable_species_not_annotated(self, annotator):
        seq = make_nuc_seq()
        seq.blastx_hits.append(make_blastx("protein [Unknown virus XYZ]"))
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.taxonomy is None

    def test_blastx_title_without_brackets_not_annotated(self, annotator):
        seq = make_nuc_seq()
        seq.blastx_hits.append(make_blastx("hypothetical protein"))
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.taxonomy is None

    def test_empty_sequence_list(self, annotator):
        assert annotator.annotate([]) == 0

    def test_multiple_sequences_partial_resolution(self, annotator):
        seq1 = make_nuc_seq("seq1")
        seq1.blastx_hits.append(make_blastx("RdRp [Influenza A virus]"))
        seq2 = make_nuc_seq("seq2")
        seq2.blastx_hits.append(make_blastx("protein [Unknown virus XYZ]"))
        seq3 = make_nuc_seq("seq3")

        count = annotator.annotate([seq1, seq2, seq3])
        assert count == 1
        assert seq1.taxonomy is not None
        assert seq2.taxonomy is None
        assert seq3.taxonomy is None

    def test_prefers_nr_hits_over_refseq_hits(self, annotator):
        seq = make_nuc_seq()
        seq.blastx_hits.append(make_blastx("protein [Influenza A virus]", bit_score=100.0))
        seq.blastx_nr_hits.append(make_blastx("coat protein [Tobacco mosaic virus]", bit_score=50.0))
        annotator.annotate([seq])
        assert seq.taxonomy.tax_id == 2


# ---------------------------------------------------------------------------
# ViralFamilyLoader
# ---------------------------------------------------------------------------

class TestViralFamilyLoader:
    @pytest.fixture
    def paths(self, tmp_path):
        high = write_family_json(tmp_path, HIGH_RECORDS, "high.json")
        low  = write_family_json(tmp_path, LOW_RECORDS,  "low.json")
        return high, low

    def test_load_count_matches_union_of_keys(self, paths):
        records = ViralFamilyLoader.load(*paths)
        assert len(records) == 4

    def test_info_high_and_info_low_populated(self, paths):
        records = ViralFamilyLoader.load(*paths)
        ortho = next(r for r in records if r.name == "Orthomyxoviridae")
        assert ortho.info_high == "Orthomyxo HIGH text"
        assert ortho.info_low  == "Orthomyxo LOW text"

    def test_returns_viral_family_info_instances(self, paths):
        records = ViralFamilyLoader.load(*paths)
        assert all(isinstance(r, ViralFamilyInfo) for r in records)

    def test_source_type_name_fields_populated(self, paths):
        records = ViralFamilyLoader.load(*paths)
        ortho = next(r for r in records if r.name == "Orthomyxoviridae")
        assert ortho.source == "ICTV"
        assert ortho.type   == "Family"

    def test_record_only_in_high_gets_empty_info_low(self, tmp_path):
        high = write_family_json(tmp_path, [
            {"source": "ICTV", "type": "Family", "name": "OnlyHigh", "info": "HIGH only"},
        ], "high.json")
        low  = write_family_json(tmp_path, [], "low.json")
        records = ViralFamilyLoader.load(high, low)
        assert len(records) == 1
        assert records[0].info_high == "HIGH only"
        assert records[0].info_low  == ""

    def test_record_only_in_low_gets_empty_info_high(self, tmp_path):
        high = write_family_json(tmp_path, [], "high.json")
        low  = write_family_json(tmp_path, [
            {"source": "ICTV", "type": "Family", "name": "OnlyLow", "info": "LOW only"},
        ], "low.json")
        records = ViralFamilyLoader.load(high, low)
        assert len(records) == 1
        assert records[0].info_high == ""
        assert records[0].info_low  == "LOW only"

    def test_missing_high_file_loads_from_low_only(self, tmp_path):
        low = write_family_json(tmp_path, LOW_RECORDS, "low.json")
        records = ViralFamilyLoader.load("/nonexistent/high.json", low)
        assert len(records) == len(LOW_RECORDS)
        assert all(r.info_high == "" for r in records)

    def test_missing_low_file_loads_from_high_only(self, tmp_path):
        high = write_family_json(tmp_path, HIGH_RECORDS, "high.json")
        records = ViralFamilyLoader.load(high, "/nonexistent/low.json")
        assert len(records) == len(HIGH_RECORDS)
        assert all(r.info_low == "" for r in records)

    def test_both_files_missing_returns_empty_list(self):
        records = ViralFamilyLoader.load("/no/high.json", "/no/low.json")
        assert records == []

    def test_records_missing_type_or_name_are_skipped(self, tmp_path):
        high = write_family_json(tmp_path, [
            {"source": "ICTV", "type": "Family", "name": "GoodRecord", "info": "ok"},
            {"source": "ICTV", "type": "",       "name": "NoType",      "info": "bad"},
            {"source": "ICTV", "type": "Family", "name": "",            "info": "bad"},
        ], "high.json")
        low = write_family_json(tmp_path, [], "low.json")
        records = ViralFamilyLoader.load(high, low)
        assert len(records) == 1
        assert records[0].name == "GoodRecord"

    def test_join_merges_same_type_and_name(self, tmp_path):
        high = write_family_json(tmp_path, [
            {"source": "ICTV", "type": "Family", "name": "SharedFamily", "info": "HIGH"},
        ], "high.json")
        low = write_family_json(tmp_path, [
            {"source": "ICTV", "type": "Family", "name": "SharedFamily", "info": "LOW"},
        ], "low.json")
        records = ViralFamilyLoader.load(high, low)
        assert len(records) == 1
        assert records[0].info_high == "HIGH"
        assert records[0].info_low  == "LOW"


# ---------------------------------------------------------------------------
# ViralFamilyAnnotator
# ---------------------------------------------------------------------------

class TestViralFamilyAnnotator:
    @pytest.fixture
    def annotator(self, tmp_path):
        high = write_family_json(tmp_path, HIGH_RECORDS, "high.json")
        low  = write_family_json(tmp_path, LOW_RECORDS,  "low.json")
        return ViralFamilyAnnotator(high, low)

    def test_annotates_via_family(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="Orthomyxoviridae", genus=None, order=None)
        count = annotator.annotate([seq])
        assert count == 1
        assert seq.viral_family_info is not None
        assert seq.viral_family_info.name == "Orthomyxoviridae"
        assert seq.viral_family_info.type == "Family"

    def test_annotates_via_genus_fallback(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="NoSuchFamily", genus="Tobamovirus", order=None)
        count = annotator.annotate([seq])
        assert count == 1
        assert seq.viral_family_info.name == "Tobamovirus"
        assert seq.viral_family_info.type == "Genus"

    def test_annotates_via_order_fallback(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family=None, genus=None, order="Articulavirales")
        count = annotator.annotate([seq])
        assert count == 1
        assert seq.viral_family_info.name == "Articulavirales"
        assert seq.viral_family_info.type == "Order"

    def test_family_takes_priority_over_genus(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="Orthomyxoviridae", genus="Tobamovirus", order=None)
        annotator.annotate([seq])
        assert seq.viral_family_info.name == "Orthomyxoviridae"

    def test_genus_takes_priority_over_order(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family=None, genus="Tobamovirus", order="Articulavirales")
        annotator.annotate([seq])
        assert seq.viral_family_info.name == "Tobamovirus"

    def test_info_high_and_low_set_on_annotated_seq(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="Orthomyxoviridae", genus=None, order=None)
        annotator.annotate([seq])
        assert seq.viral_family_info.info_high == "Orthomyxo HIGH text"
        assert seq.viral_family_info.info_low  == "Orthomyxo LOW text"

    def test_no_taxonomy_not_annotated(self, annotator):
        seq = make_nuc_seq()
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.viral_family_info is None

    def test_taxonomy_with_no_matching_entry_not_annotated(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="UnknownFamily", genus="UnknownGenus", order="UnknownOrder")
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.viral_family_info is None

    def test_lookup_is_case_insensitive(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family="ORTHOMYXOVIRIDAE", genus=None, order=None)
        count = annotator.annotate([seq])
        assert count == 1
        assert seq.viral_family_info.name == "Orthomyxoviridae"

    def test_empty_sequence_list(self, annotator):
        assert annotator.annotate([]) == 0

    def test_all_taxonomy_fields_none_not_annotated(self, annotator):
        seq = make_nuc_seq()
        seq.taxonomy = make_taxonomy(family=None, genus=None, order=None)
        count = annotator.annotate([seq])
        assert count == 0
        assert seq.viral_family_info is None

    def test_multiple_sequences_partial_annotation(self, annotator):
        seq1 = make_nuc_seq("seq1")
        seq1.taxonomy = make_taxonomy(family="Orthomyxoviridae")

        seq2 = make_nuc_seq("seq2")
        seq2.taxonomy = make_taxonomy(family="UnknownFamily", genus="UnknownGenus", order="UnknownOrder")

        seq3 = make_nuc_seq("seq3")  # no taxonomy

        count = annotator.annotate([seq1, seq2, seq3])
        assert count == 1
        assert seq1.viral_family_info is not None
        assert seq2.viral_family_info is None
        assert seq3.viral_family_info is None
