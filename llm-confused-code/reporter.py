"""
reporter.py — ViralQuest JSON output builder.

Three additions over the previous version:
  1. BLASTx fields now correctly map to the right BlastxResult attributes.
  2. Domain overlap filtering per-database before serialisation.
  3. Taxonomy join from viralTax.json.xz matched on ScientificName → Species.
"""

import json
import lzma
import re
from pathlib import Path

from loguru import logger
from .biodata import NucSequence, BlastxResult, BlastnResult, HmmDomain, Orf


# ---------------------------------------------------------------------------
# Taxonomy loader
# ---------------------------------------------------------------------------

class ViralTaxonomyDB:
    """
    Loads viralTax.json.xz and provides fast lookup by organism name.

    Lookup order for a given subject title:
      1. Exact match on ScientificName
      2. Exact match on Species
      3. Substring search of subject title against ScientificName values
         (handles titles like "gi|...|gb|...| Tomato bushy stunt virus, ...")
      4. Returns empty taxonomy dict if nothing matches.
    """

    _EMPTY: dict = {
        "TaxId":         None,
        "ScientificName": None,
        "NoRank":        None,
        "Clade":         None,
        "Kingdom":       None,
        "Phylum":        None,
        "Class":         None,
        "Order":         None,
        "Family":        None,
        "Subfamily":     None,
        "Genus":         None,
        "Species":       None,
        "Genome":        None,
    }

    def __init__(self, tax_path: str | Path):
        self._by_scientific: dict[str, dict] = {}
        self._by_species:    dict[str, dict] = {}
        self._all:           list[dict]      = []
        self._load(Path(tax_path))

    def _load(self, path: Path) -> None:
        if not path.exists():
            logger.warning(f"viralTax not found at {path} — taxonomy will be empty.")
            return
        try:
            opener = lzma.open if path.suffix == ".xz" else open
            with opener(path, "rb") as fh:
                records: list[dict] = json.load(fh)
        except Exception as exc:
            logger.error(f"Failed to load viralTax: {exc}")
            return

        for rec in records:
            sci = rec.get("ScientificName") or ""
            spe = rec.get("Species") or ""
            if sci:
                self._by_scientific[sci.lower()] = rec
            if spe and spe != sci:
                self._by_species[spe.lower()] = rec
            self._all.append(rec)

        logger.success(f"viralTax loaded — {len(self._all)} records.")

    def lookup(self, subject_title: str) -> dict:
        """Return taxonomy dict for the best match, or _EMPTY."""
        if not subject_title or subject_title == "NaN":
            return dict(self._EMPTY)

        title_lower = subject_title.lower()

        # 1. exact ScientificName
        hit = self._by_scientific.get(title_lower)
        if hit:
            return hit

        # 2. exact Species
        hit = self._by_species.get(title_lower)
        if hit:
            return hit

        # 3. substring scan — subject titles often contain the organism name
        #    e.g. "gi|...|gb|MK204389.1| Sharp-tailed sandpiper Picornavirus..."
        for sci_name, rec in self._by_scientific.items():
            if sci_name and sci_name in title_lower:
                return rec

        # 4. Species substring scan
        for spe_name, rec in self._by_species.items():
            if spe_name and spe_name in title_lower:
                return rec

        return dict(self._EMPTY)


# ---------------------------------------------------------------------------
# Domain overlap filter
# ---------------------------------------------------------------------------

class DomainOverlapFilter:
    """
    For each ORF, removes overlapping domains **within the same database**,
    keeping the one with the highest score.

    Two domains overlap when one starts before the other ends:
        domain_a.start < domain_b.stop  AND  domain_b.start < domain_a.stop
    """

    @staticmethod
    def filter_orf(orf: Orf) -> None:
        """Mutates orf.domains in-place — replaces with filtered list."""
        if not orf.domains:
            return

        # group by database
        db_groups: dict[str, list[HmmDomain]] = {}
        for d in orf.domains:
            db_groups.setdefault(d.database, []).append(d)

        kept: list[HmmDomain] = []
        for db, domains in db_groups.items():
            # sort by start position, then resolve overlaps greedily
            sorted_domains = sorted(domains, key=lambda d: (d.start, -d.score))
            resolved: list[HmmDomain] = []
            for current in sorted_domains:
                if not resolved:
                    resolved.append(current)
                    continue
                last = resolved[-1]
                # check overlap: current starts before last ends
                if current.start < last.stop:
                    # keep the higher-scoring one
                    if current.score > last.score:
                        resolved[-1] = current
                    # else discard current (last stays)
                else:
                    resolved.append(current)
            kept.extend(resolved)

        removed = len(orf.domains) - len(kept)
        if removed:
            logger.debug(f"  {orf.name}: removed {removed} overlapping domain(s).")
        orf.domains = kept

    @classmethod
    def filter_all(cls, nuc_seqs: list[NucSequence]) -> int:
        """Filters all ORFs in all sequences. Returns total domains removed."""
        before = sum(len(orf.domains) for seq in nuc_seqs for orf in seq.orfs)
        for seq in nuc_seqs:
            for orf in seq.orfs:
                cls.filter_orf(orf)
        after = sum(len(orf.domains) for seq in nuc_seqs for orf in seq.orfs)
        removed = before - after
        logger.info(f"Domain overlap filter: {removed} redundant domain(s) removed.")
        return removed


# ---------------------------------------------------------------------------
# Reporter
# ---------------------------------------------------------------------------

_NAN = "NaN"


class ViralQuestReporter:
    """
    Builds and writes the three-section ViralQuest JSON:
        Viral_Hits  — one entry per viral sequence (BLAST + taxonomy)
        HMM_hits    — one entry per viral sequence (best ORF domains)
        ORF_Data    — all ORFs per viral sequence
    """

    def __init__(self, sample_name: str, tax_db: ViralTaxonomyDB | None = None):
        self.sample_name = sample_name
        self.tax_db      = tax_db   # optional; taxonomy fields are NaN if absent

    # ------------------------------------------------------------------
    # Viral_Hits builder
    # ------------------------------------------------------------------

    def _build_viral_hit(self, seq: NucSequence) -> dict:
        bx = seq.best_blastx   # best NR hit, or best refseq hit as fallback
        bn = seq.best_blastn

        # ---- BLASTx fields ----
        if bx:
            blastx = {
                "BLASTx_Qlength":       str(seq.length),
                "BLASTx_Slength":       str(bx.aln_length),   # alignment length
                "BLASTx_Cover":         f"{bx.query_coverage:.1f}",
                "BLASTx_Ident":         f"{bx.pct_identity:.1f}",
                "BLASTx_evalue":        str(bx.e_value),
                "BLASTx_Subject_Title": bx.subject_title,
                "BLASTx_Subject_ID":    bx.subject_id,
            }
        else:
            blastx = {
                "BLASTx_Qlength":       str(seq.length),
                "BLASTx_Slength":       _NAN,
                "BLASTx_Cover":         _NAN,
                "BLASTx_Ident":         _NAN,
                "BLASTx_evalue":        _NAN,
                "BLASTx_Subject_Title": _NAN,
                "BLASTx_Subject_ID":    _NAN,
            }

        # ---- BLASTn fields ----
        if bn:
            blastn = {
                "BLASTn_Qlength":       str(bn.query_length),
                "BLASTn_Slength":       str(bn.subject_length),
                "BLASTn_Cover":         str(bn.query_coverage),
                "BLASTn_Ident":         f"{bn.pct_identity:.1f}",
                "BLASTn_evalue":        str(bn.e_value),
                "BLASTn_Subject_Title": bn.subject_title,
            }
        else:
            blastn = {
                "BLASTn_Qlength":       _NAN,
                "BLASTn_Slength":       _NAN,
                "BLASTn_Cover":         _NAN,
                "BLASTn_Ident":         _NAN,
                "BLASTn_evalue":        _NAN,
                "BLASTn_Subject_Title": _NAN,
            }

        # ---- Taxonomy ----
        # Try BLASTn title first (nucleotide match is more specific),
        # fall back to BLASTx title.
        tax_source = (bn.subject_title if bn else None) or (bx.subject_title if bx else None)
        tax = self.tax_db.lookup(tax_source) if self.tax_db else ViralTaxonomyDB._EMPTY.copy()

        entry = {
            "Sample_name": self.sample_name,
            "QueryID":     seq.id,
            **blastx,
            **blastn,
            # taxonomy block
            "TaxId":          tax.get("TaxId"),
            "ScientificName": tax.get("ScientificName"),
            "NoRank":         tax.get("NoRank"),
            "Clade":          tax.get("Clade"),
            "Kingdom":        tax.get("Kingdom"),
            "Phylum":         tax.get("Phylum"),
            "Class":          tax.get("Class"),
            "Order":          tax.get("Order"),
            "Family":         tax.get("Family"),
            "Subfamily":      tax.get("Subfamily"),
            "Genus":          tax.get("Genus"),
            "Species":        tax.get("Species"),
            "Genome":         tax.get("Genome"),
            # sequence
            "GC_content":     f"{seq.gc_content:.2f}",
            "N_count":        seq.n_count,
            "FullSeq":        seq.sequence,
        }
        return entry

    # ------------------------------------------------------------------
    # HMM_hits builder
    # ------------------------------------------------------------------

    def _build_hmm_hit(self, seq: NucSequence) -> dict | None:
        """
        One entry per viral sequence.
        Uses the ORF with the most domain hits as representative.
        Domains are already overlap-filtered at this point.
        """
        orfs_with_domains = [o for o in seq.orfs if o.domains]
        if not orfs_with_domains:
            return None

        best_orf = max(orfs_with_domains, key=lambda o: len(o.domains))

        entry: dict = {
            "Sample_name": self.sample_name,
            "Query_name":  best_orf.name,
            "QueryID":     seq.id,
        }

        db_order = ["RVDB", "Vfam", "EggNOG", "Pfam"]
        for db in db_order:
            db_domains = [d for d in best_orf.domains if d.database == db]
            for idx, dom in enumerate(db_domains):
                # first hit: no suffix; subsequent hits: _2, _3, …
                suffix = "" if idx == 0 else f"_{idx + 1}"
                prefix = f"{db}{suffix}"
                entry[f"{prefix}_TargetID"]    = dom.target
                entry[f"{prefix}_Description"] = dom.description
                entry[f"{prefix}_Score"]       = str(dom.score)
                entry[f"{prefix}_Start"]       = str(dom.start)
                entry[f"{prefix}_End"]         = str(dom.stop)
                entry[f"{prefix}_length"]      = str(dom.length)
                if db == "Pfam":
                    entry[f"{prefix}_Type"]    = dom.type
                    entry[f"{prefix}_Details"] = dom.details

        entry["FullSequence"] = best_orf.aa_sequence
        return entry

    # ------------------------------------------------------------------
    # ORF_Data builder
    # ------------------------------------------------------------------

    def _build_orf_data(self, seq: NucSequence) -> list[dict]:
        return [
            {
                "Query_name":  orf.name,
                "QueryID":     seq.id,
                "start":       orf.start_position,
                "end":         orf.stop_position,
                "strand":      orf.strand,
                "sequence":    orf.aa_sequence,
                "type":        orf.orf_type,
                "length":      orf.length_nt,
                "frame":       str(orf.frame),
                "start_codon": orf.start_codon,
                "stop_codon":  orf.stop_codon if orf.stop_codon else "NA",
                # include per-ORF domain summary for convenience
                "domains": [
                    {
                        "database":    d.database,
                        "target":      d.target,
                        "description": d.description,
                        "score":       d.score,
                        "e_value":     d.e_value,
                        "start":       d.start,
                        "stop":        d.stop,
                        "length":      d.length,
                    }
                    for d in orf.domains
                ],
            }
            for orf in seq.orfs
        ]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self, nuc_seqs: list[NucSequence]) -> dict:
        """
        Runs domain overlap filtering, then builds the full JSON dict.
        Only sequences with is_viral=True are included.
        """
        viral_seqs = [s for s in nuc_seqs if s.is_viral]
        logger.info(
            f"Reporter: {len(viral_seqs)} viral / {len(nuc_seqs)} total sequences."
        )

        # overlap filtering happens here — mutates orf.domains in-place
        DomainOverlapFilter.filter_all(viral_seqs)

        viral_hits: list[dict] = []
        hmm_hits:   list[dict] = []
        orf_data:   dict       = {}

        for seq in viral_seqs:
            viral_hits.append(self._build_viral_hit(seq))

            hmm_entry = self._build_hmm_hit(seq)
            if hmm_entry:
                hmm_hits.append(hmm_entry)

            orf_data[seq.id] = self._build_orf_data(seq)

        return {
            "Viral_Hits": viral_hits,
            "HMM_hits":   hmm_hits,
            "ORF_Data":   orf_data,
        }

    def save(self, nuc_seqs: list[NucSequence], output_path: Path) -> Path:
        """Builds, overlap-filters, and writes the JSON. Returns output path."""
        data = self.build(nuc_seqs)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=4, ensure_ascii=False)
        logger.success(
            f"JSON saved → {output_path}  "
            f"({len(data['Viral_Hits'])} viral hits, "
            f"{len(data['HMM_hits'])} HMM entries)"
        )
        return output_path