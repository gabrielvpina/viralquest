import json
import lzma

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


class TaxonomyLoader:
    """
    Loads viralTax.json.xz and builds two lookup indices:
      - by ScientificName  (primary match)
      - by Species field   (fallback; first record per species name wins)
    """

    @staticmethod
    def load(json_xz_path: str) -> tuple[dict[str, dict], dict[str, dict]]:
        try:
            with lzma.open(json_xz_path, "rt", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            logger.error(f"Cannot read {json_xz_path}: {exc}")
            return {}, {}

        sci_idx: dict[str, dict] = {}
        sp_idx: dict[str, dict] = {}

        for record in data:
            sci = record.get("ScientificName")
            if sci:
                sci_idx[sci.lower()] = record
            sp = record.get("Species")
            if sp:
                key = sp.lower()
                if key not in sp_idx:
                    sp_idx[key] = record

        logger.success(
            f"TaxonomyLoader: {len(sci_idx)} scientific names, "
            f"{len(sp_idx)} species names indexed from {json_xz_path}."
        )
        return sci_idx, sp_idx


class TaxonomyMatcher:
    """Matches a species string against the loaded taxonomy indices."""

    def __init__(self, sci_idx: dict[str, dict], sp_idx: dict[str, dict]):
        self._sci = sci_idx
        self._sp  = sp_idx

    def match(self, species: str) -> Taxonomy | None:
        if not species:
            return None
        key    = species.lower()
        record = self._sci.get(key) or self._sp.get(key)
        if record is None:
            return None
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

    def __init__(self, json_xz_path: str):
        sci_idx, sp_idx = TaxonomyLoader.load(json_xz_path)
        self._matcher   = TaxonomyMatcher(sci_idx, sp_idx)

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
