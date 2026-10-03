import json
import lzma
from pathlib import Path

from loguru import logger

from viralquest.biodata import NucSequence, Taxonomy, ViralFamilyInfo


class ViralFamilyLoader:
    """
    Loads the highToken and lowToken JSON files and merges them into
    ViralFamilyInfo objects, joined on the unique (type, name) key.
    Records present in only one file receive an empty string for the
    missing side.
    """

    @staticmethod
    def _read(path: str) -> dict[tuple[str, str], dict]:
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            logger.error(f"Cannot read {path}: {exc}")
            return {}
        return {(r["type"], r["name"]): r for r in data if r.get("type") and r.get("name")}

    @classmethod
    def load(cls, high_path: str, low_path: str) -> list[ViralFamilyInfo]:
        high = cls._read(high_path)
        low  = cls._read(low_path)
        all_keys = high.keys() | low.keys()
        records: list[ViralFamilyInfo] = []
        for key in all_keys:
            h = high.get(key, {})
            l = low.get(key, {})
            base = h or l
            records.append(ViralFamilyInfo(
                source    = base.get("source", ""),
                type      = base.get("type",   ""),
                name      = base.get("name",   ""),
                info_high = h.get("info", ""),
                info_low  = l.get("info", ""),
            ))
        logger.success(
            f"ViralFamilyLoader: {len(records)} entries merged "
            f"({len(high)} high-token, {len(low)} low-token)."
        )
        return records


class ViralFamilyAnnotator:
    """
    Attaches ViralFamilyInfo to NucSequence objects.

    Lookup priority (most-specific first) uses the taxonomy already
    resolved on each sequence:
        1. taxonomy.family  → type "Family"
        2. taxonomy.genus   → type "Genus"
        3. taxonomy.order   → type "Order"
    """

    def __init__(self, high_path: str, low_path: str):
        records = ViralFamilyLoader.load(high_path, low_path)
        self._idx: dict[tuple[str, str], ViralFamilyInfo] = {
            (r.type.lower(), r.name.lower()): r for r in records
        }

    def _lookup(self, type_: str, name: str | None) -> ViralFamilyInfo | None:
        if not name:
            return None
        return self._idx.get((type_, name.lower()))

    def annotate(self, nuc_seqs: list[NucSequence]) -> int:
        """
        Sets seq.viral_family_info for each sequence that has a resolved
        taxonomy. Returns the count of sequences annotated.
        """
        annotated = 0
        for seq in nuc_seqs:
            tax = seq.taxonomy
            if not tax:
                continue
            info = (
                self._lookup("family", tax.family)
                or self._lookup("genus",  tax.genus)
                or self._lookup("order",  tax.order)
            )
            if info:
                seq.viral_family_info = info
                annotated += 1
            else:
                logger.debug(
                    f"No ViralFamilyInfo for seq '{seq.id}' "
                    f"(family={tax.family}, genus={tax.genus}, order={tax.order})"
                )
        logger.success(f"ViralFamilyInfo annotated {annotated}/{len(nuc_seqs)} sequences.")
        return annotated


SYNONYMS_FILE = "viralSynonyms.json.xz"


class TaxonomyLoader:
    """
    Loads viralTax.json.xz and builds the lookup indices:
      - by ScientificName  (primary match)
      - by Species field   (fallback; first record per species name wins)
      - by TaxId
      - by NCBI synonym    (old / equivalent / common names, from
                            viralSynonyms.json.xz built by
                            scripts/build_viral_synonyms.py; optional)
    """

    @staticmethod
    def load(json_xz_path: str) -> tuple[dict[str, dict], dict[str, dict]]:
        sci_idx, sp_idx, _, _ = TaxonomyLoader.load_indices(json_xz_path, synonyms_path=False)
        return sci_idx, sp_idx

    @staticmethod
    def load_indices(json_xz_path: str, synonyms_path: str | None | bool = None
                     ) -> tuple[dict[str, dict], dict[str, dict], dict[int, dict], dict[str, dict]]:
        """
        (scientific-name, species, taxid, synonym) indices.  ``synonyms_path``
        defaults to viralSynonyms.json.xz next to the taxonomy file; ``False``
        skips synonyms.  A missing synonym file only disables that index.
        """
        try:
            with lzma.open(json_xz_path, "rt", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            logger.error(f"Cannot read {json_xz_path}: {exc}")
            return {}, {}, {}, {}

        sci_idx: dict[str, dict] = {}
        sp_idx: dict[str, dict] = {}
        taxid_idx: dict[int, dict] = {}

        for record in data:
            tid = record.get("TaxId")
            if tid:
                taxid_idx[int(tid)] = record
            sci = record.get("ScientificName")
            if sci:
                sci_idx[sci.lower()] = record
            sp = record.get("Species")
            if sp:
                key = sp.lower()
                if key not in sp_idx:
                    sp_idx[key] = record

        syn_idx: dict[str, dict] = {}
        if synonyms_path is not False:
            syn_path = Path(synonyms_path) if synonyms_path else Path(json_xz_path).with_name(SYNONYMS_FILE)
            syn_idx = TaxonomyLoader._load_synonyms(syn_path, taxid_idx)

        logger.success(
            f"TaxonomyLoader: {len(sci_idx)} scientific names, "
            f"{len(sp_idx)} species names, {len(syn_idx)} NCBI synonyms indexed from {json_xz_path}."
        )
        return sci_idx, sp_idx, taxid_idx, syn_idx

    @staticmethod
    def _load_synonyms(path: Path, taxid_idx: dict[int, dict]) -> dict[str, dict]:
        if not path.exists():
            logger.debug(f"No synonym table at {path} — matching on primary names only.")
            return {}
        try:
            with lzma.open(path, "rt", encoding="utf-8") as fh:
                names = json.load(fh).get("names", {})
        except Exception as exc:
            logger.warning(f"Cannot read synonym table {path}: {exc}")
            return {}
        return {name: taxid_idx[tid] for name, tid in names.items() if tid in taxid_idx}


class TaxonomyMatcher:
    """Matches a species string against the loaded taxonomy indices
    (scientific name, then species, then NCBI synonym)."""

    def __init__(self, sci_idx: dict[str, dict], sp_idx: dict[str, dict],
                 syn_idx: dict[str, dict] | None = None):
        self._sci = sci_idx
        self._sp  = sp_idx
        self._syn = syn_idx or {}

    def lookup(self, name: str) -> tuple[dict | None, bool]:
        """(record, via_synonym) for a name, or (None, False)."""
        if not name:
            return None, False
        key = name.lower()
        rec = self._sci.get(key) or self._sp.get(key)
        if rec:
            return rec, False
        rec = self._syn.get(key)
        return (rec, True) if rec else (None, False)

    def match(self, species: str) -> Taxonomy | None:
        record, _ = self.lookup(species)
        return record_to_taxonomy(record) if record else None


def record_to_taxonomy(record: dict) -> Taxonomy:
    """Build a Taxonomy from one viralTax.json record."""
    return Taxonomy(
        tax_id        = record.get("TaxId", 0),
        scientific_name = record.get("ScientificName", ""),
        no_rank       = record.get("No_Rank"),
        clade         = record.get("Clade"),
        kingdom       = record.get("Kingdom"),
        phylum        = record.get("Phylum"),
        class_        = record.get("Class"),
        order         = record.get("Order"),
        family        = record.get("Family"),
        subfamily     = record.get("Subfamily"),
        genus         = record.get("Genus"),
        species       = record.get("Species"),
        genome        = record.get("Genome"),
    )


class TaxonomyAnnotator:
    """
    Annotates NucSequence objects with taxonomy data derived from
    the species extracted from their best BLASTx hit's subject_title.
    """

    def __init__(self, json_xz_path: str, synonyms_path: str | None = None):
        sci_idx, sp_idx, taxid_idx, syn_idx = TaxonomyLoader.load_indices(json_xz_path, synonyms_path)
        self._matcher   = TaxonomyMatcher(sci_idx, sp_idx, syn_idx)
        # Shared with BlastnTaxonomyResolver so the files are decompressed once.
        self.indices    = (sci_idx, sp_idx, taxid_idx, syn_idx)

    def annotate(self, nuc_seqs: list[NucSequence]) -> int:
        """
        Sets seq.taxonomy for each sequence that has a resolvable species.
        Returns the count of sequences that received a taxonomy annotation.
        """
        resolved = 0
        for seq in nuc_seqs:
            best = seq.best_blastx
            if not best or not best.species:
                continue
            tax = self._matcher.match(best.species)
            if tax:
                seq.taxonomy = tax
                resolved += 1
            else:
                logger.debug(f"No taxonomy match for '{best.species}' (seq: {seq.id})")
        logger.success(f"Taxonomy resolved for {resolved}/{len(nuc_seqs)} sequences.")
        return resolved
