import json
import lzma

import pytest

from viralquest.biodata import BlastxResult, NucSequence, Taxonomy
from viralquest.tax import TaxonomyAnnotator, TaxonomyLoader, TaxonomyMatcher


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


def write_tax_json(tmp_path, records: list[dict]) -> str:
    path = tmp_path / "viralTax.json.xz"
    with lzma.open(path, "wt", encoding="utf-8") as fh:
        json.dump(records, fh)
    return str(path)


SAMPLE_RECORDS = [
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


# ---------------------------------------------------------------------------
# BlastxResult.species extraction  (biodata, but tested here for proximity)
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
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
        sci, _ = TaxonomyLoader.load(path)
        assert "influenza a virus" in sci
        assert "tobacco mosaic virus" in sci

    def test_sp_idx_keyed_by_lowercase_species_field(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        assert "alphainfluenzavirus influenzae" in sp
        assert "tobamovirus virginianum" in sp

    def test_records_without_species_field_excluded_from_sp_idx(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        # TaxId=3 has Species=None — must not appear as a key
        assert None not in sp

    def test_duplicate_species_field_first_record_wins(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
        _, sp = TaxonomyLoader.load(path)
        # Both TaxId=4 and TaxId=5 share Species="Hepatitis B virus"
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
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
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
        # "Tobamovirus virginianum" is a Species value, not a ScientificName
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
        # TaxId=3 has all classification fields as null
        tax = matcher.match("Alphatorquevirus zetaone")
        assert tax is not None
        assert tax.family is None
        assert tax.genus is None
        assert tax.species is None

    def test_scientific_name_takes_priority_over_species_field(self, matcher):
        # "Hepatitis B virus" matches both ScientificName (no) and Species (yes)
        # but "Hepatitis B virus genotype A" is a ScientificName — test priority
        tax_sci = matcher.match("Hepatitis B virus genotype A")
        tax_sp  = matcher.match("Hepatitis B virus")
        assert tax_sci.tax_id == 4
        assert tax_sp.tax_id == 4   # Species fallback → first record with that Species


# ---------------------------------------------------------------------------
# TaxonomyAnnotator
# ---------------------------------------------------------------------------

class TestTaxonomyAnnotator:
    @pytest.fixture
    def annotator(self, tmp_path):
        path = write_tax_json(tmp_path, SAMPLE_RECORDS)
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
        # best_blastx picks highest bit_score → Tobacco mosaic virus
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

        seq3 = make_nuc_seq("seq3")   # no hits

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
        # best_blastx prefers blastx_nr_hits pool when non-empty
        assert seq.taxonomy.tax_id == 2
