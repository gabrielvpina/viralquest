"""NCBI synonym table: build script, loader and use in BLASTx / BLASTn taxonomy."""
import importlib.util
import io
import json
import lzma
import tarfile
from pathlib import Path

from viralquest.biodata import BlastnResult, BlastxResult, NucSequence
from viralquest.blastn_tax import BlastnTaxonomyResolver
from viralquest.tax import SYNONYMS_FILE, TaxonomyAnnotator, TaxonomyLoader

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("build_viral_synonyms", ROOT / "scripts" / "build_viral_synonyms.py")
build_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_mod)


def _rec(taxid, name, species=None, family="Rhabdoviridae", genus=None):
    return {"TaxId": taxid, "ScientificName": name, "Species": species, "Family": family,
            "Genus": genus, "No_Rank": None, "Clade": None, "Kingdom": None, "Phylum": None,
            "Class": None, "Order": None, "Subfamily": None, "Genome": None}


RECORDS = [
    _rec(152177, "Dichorhavirus orchidaceae", species="Dichorhavirus orchidaceae", genus="Dichorhavirus"),
    _rec(1213, "Dichorhavirus", genus="Dichorhavirus"),
    _rec(500, "Alpha virus A", species="Alphavirus prima", family="Fooviridae"),
    _rec(501, "Alpha virus B", species="Alphavirus prima", family="Fooviridae"),
    _rec(600, "Beta virus", species="Betavirus secunda", family="Fooviridae"),
    _rec(700, "Old merged virus", species="Gammavirus tertia", family="Fooviridae"),
]

NAMES = [
    (152177, "Orchid fleck virus", "equivalent name"),
    (152177, "Orchid fleck dichorhavirus", "equivalent name"),
    (152177, "OFV", "acronym"),                          # acronyms are left out
    (500, "shared alpha virus", "synonym"),            # same species as 501 → kept
    (501, "shared alpha virus", "synonym"),
    (500, "clash virus name", "synonym"),              # different species → ambiguous, dropped
    (600, "clash virus name", "synonym"),
    (600, "Dichorhavirus", "synonym"),                 # already a primary name → skipped
    (600, "xy", "synonym"),                            # too short
    (9999, "not in bundle virus", "synonym"),          # taxid not in the viral bundle
    (7000, "renamed old virus", "synonym"),            # current id of merged 700
]


def _dmp(rows):
    return ("".join("\t|\t".join(map(str, r)) + "\t|\n" for r in rows)).encode()


def _write_inputs(tmp_path):
    tax = tmp_path / "viralTax.json.xz"
    with lzma.open(tax, "wt", encoding="utf-8") as fh:
        json.dump(RECORDS, fh)
    dump = tmp_path / "taxdump.tar.gz"
    with tarfile.open(dump, "w:gz") as tar:
        for name, rows in (("names.dmp", [(t, n, "", c) for t, n, c in NAMES]),
                           ("merged.dmp", [(700, 7000)])):
            data = _dmp(rows)
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return tax, dump


def _build(tmp_path):
    tax, dump = _write_inputs(tmp_path)
    data = build_mod.build(dump, tax)
    out = tmp_path / SYNONYMS_FILE
    with lzma.open(out, "wt", encoding="utf-8") as fh:
        json.dump(data, fh)
    return tax, data


# ── build script ─────────────────────────────────────────────────────────────

class TestBuild:
    def test_kept_and_dropped(self, tmp_path):
        _, data = _build(tmp_path)
        n = data["names"]
        assert n["orchid fleck virus"] == 152177 and n["orchid fleck dichorhavirus"] == 152177
        assert "ofv" not in n                               # acronym class
        assert n["shared alpha virus"] == 500               # same species → lowest taxid
        assert "clash virus name" not in n and data["meta"]["dropped_ambiguous"] == 1
        assert "dichorhavirus" not in n                     # primary name wins
        assert "xy" not in n and "not in bundle virus" not in n
        assert n["renamed old virus"] == 700                # followed through merged.dmp

    def test_meta(self, tmp_path):
        _, data = _build(tmp_path)
        assert data["meta"]["names"] == len(data["names"])
        assert "acronym" not in data["meta"]["name_classes"]


# ── runtime ─────────────────────────────────────────────────────────────────

class TestLoader:
    def test_synonyms_loaded_next_to_taxonomy(self, tmp_path):
        tax, _ = _build(tmp_path)
        *_, syn = TaxonomyLoader.load_indices(str(tax))
        assert syn["orchid fleck virus"]["TaxId"] == 152177

    def test_missing_table_is_optional(self, tmp_path):
        tax, _ = _write_inputs(tmp_path)
        sci, sp, tid, syn = TaxonomyLoader.load_indices(str(tax))
        assert syn == {} and sci                             # primary names still work

    def test_legacy_load_ignores_synonyms(self, tmp_path):
        tax, _ = _build(tmp_path)
        assert len(TaxonomyLoader.load(str(tax))) == 2


def _bn(title):
    return BlastnResult(qseqid="s1", qlen=1000, slen=1000, qcovhsp=90, pident=97.0, evalue=1e-50,
                        bit_score=100.0, stitle=title, accession="NC_009609")


class TestBlastnWithSynonyms:
    def _resolver(self, tmp_path):
        tax, _ = _build(tmp_path)
        return BlastnTaxonomyResolver(*TaxonomyLoader.load_indices(str(tax)))

    def test_former_name_resolves_to_current_species(self, tmp_path):
        bt = self._resolver(tmp_path).resolve_hit(_bn("NC_009609.1 Orchid fleck virus genomic RNA, segment RNA 2"))
        assert bt.status == "resolved" and bt.matched_name == "Orchid fleck virus"
        assert bt.synonym_of == "Dichorhavirus orchidaceae" and bt.taxonomy.species == "Dichorhavirus orchidaceae"

    def test_longer_synonym_beats_shorter_primary_genus(self, tmp_path):
        bt = self._resolver(tmp_path).resolve_hit(_bn("Orchid fleck dichorhavirus isolate Smilax segment RNA 2"))
        assert bt.matched_name == "Orchid fleck dichorhavirus" and bt.taxonomy.species == "Dichorhavirus orchidaceae"

    def test_primary_name_has_no_synonym_note(self, tmp_path):
        bt = self._resolver(tmp_path).resolve_hit(_bn("Dichorhavirus orchidaceae isolate S9"))
        assert bt.synonym_of is None


def test_blastx_taxonomy_uses_synonyms(tmp_path):
    tax, _ = _build(tmp_path)
    seq = NucSequence(id="s1", sequence="ACGT")
    hit = BlastxResult(query_id="s1", subject_id="YP_1", subject_title="polymerase [Orchid fleck virus]",
                       pct_identity=90.0, aln_length=100, mismatches=0, gap_opens=0, query_start=1,
                       query_end=300, subject_start=1, subject_end=100, e_value=1e-30, bit_score=200.0)
    hit.species = "Orchid fleck virus"
    seq.blastx_hits.append(hit)
    TaxonomyAnnotator(str(tax)).annotate([seq])
    assert seq.taxonomy is not None and seq.taxonomy.tax_id == 152177


def test_synonym_table_is_optional_for_the_cli():
    from pathlib import Path as P
    from viralquest.cli import _check_databases

    class Console:
        def print(self, *a, **k):
            pass

    assert _check_databases(Console(), {"viral_syn": P("/nonexistent/viralSynonyms.json.xz")}) is True
    assert _check_databases(Console(), {"viral_tax": P("/nonexistent/viralTax.json.xz")}) is False


def test_bundled_table_resolves_known_renames():
    """The shipped data/viralSynonyms.json.xz carries the renames seen in real runs."""
    path = ROOT / "data" / SYNONYMS_FILE
    with lzma.open(path, "rt", encoding="utf-8") as fh:
        names = json.load(fh)["names"]
    assert names["orchid fleck virus"] == 152177
    assert "cotesia congregata virus" in names
