"""
blastn_tax.py — Taxonomy of the best BLASTn hit of each sequence (offline).

Runs after BLASTn (local or online), on the finished hits.  BLASTn subject
titles are free text ("Orchid fleck dichorhavirus isolate Smilax segment RNA 2,
complete sequence", "MAG: Hubei diptera virus 3 strain X ..."), so the species
is read out of the title: every token that looks like a virus name (…virus,
…viridae, …virales, …viricetes, …phage, …) anchors a window that grows up to
``MAX_FLANK`` words to the left and to the right; candidate phrases are tried
longest first against the bundled viral taxonomy (viralTax.json.xz), so the
most complete name found wins; former / equivalent names from NCBI
(viralSynonyms.json.xz) are tried when no primary name matches.  No network access.

Hits whose title names no known viral taxon (e.g. a host mRNA) are recorded
as ``unresolved``.
"""

from __future__ import annotations

import re

from loguru import logger

from viralquest.biodata import BlastnResult, BlastnTaxonomy, NucSequence
from viralquest.tax import TaxonomyMatcher, record_to_taxonomy

MAX_FLANK = 5


# ── Virus names in a free-text title ────────────────────────────────────────

# A token that names a virus or a viral taxon (rank suffixes from ICTV).
_ANCHOR_RE = re.compile(
    r"(?:vir(?:us|uses|idae|inae|ales|icetes|icota|iae|ia|ae|oid|oids|ophyta)|phages?|satellites?)$",
    re.IGNORECASE,
)
# Words that alone say nothing about which virus it is.
_GENERIC = {"virus", "viruses", "phage", "phages", "satellite", "satellites", "viroid", "viroids"}
# Database / annotation prefixes in GenBank definitions.
_PREFIX_RE = re.compile(r"^(?:(?:UNVERIFIED|MAG|TPA|TSA|PREDICTED|TPA_asm|TPA_exp|TPA_inf)\s*:\s*)+", re.IGNORECASE)
_EDGE = ".,;:()[]{}\"'"


def _tokens(title: str) -> list[str]:
    title = _PREFIX_RE.sub("", title or "")
    return [w.strip(_EDGE) for w in title.split() if w.strip(_EDGE)]


def is_viral_title(title: str) -> bool:
    """True if a free-text hit title carries a virus-like token (…virus,
    …viridae, …phage, …).  No taxonomy lookup — a cheap pre-filter."""
    return any(_ANCHOR_RE.search(w) for w in _tokens(title))


def name_candidates(title: str, max_flank: int = MAX_FLANK) -> list[str]:
    """
    Candidate virus names from a free-text hit title, longest first.

    A bracketed organism (``... [Some virus]``) is tried first.  Then, around
    every virus-like token, windows of 0..max_flank words on each side.
    """
    out: list[str] = []
    seen: set[str] = set()

    def add(c: str) -> None:
        k = c.lower()
        if c and k not in seen and k not in _GENERIC:
            seen.add(k)
            out.append(c)

    m = re.search(r"\[([^\]]+)\]\s*$", title or "")
    if m:
        add(m.group(1).strip())

    toks = _tokens(title)
    spans: set[tuple[int, int]] = set()
    for i, w in enumerate(toks):
        if not _ANCHOR_RE.search(w):
            continue
        for left in range(max_flank + 1):
            for right in range(max_flank + 1):
                a, b = i - left, i + right
                if a >= 0 and b < len(toks):
                    spans.add((a, b))
    # Longest phrase first; among equals, the one starting further left.
    for a, b in sorted(spans, key=lambda s: (-(s[1] - s[0]), s[0])):
        add(" ".join(toks[a:b + 1]))
    return out


class TitleMatcher:
    """Virus name in the hit title → bundled viral taxonomy (primary names,
    then NCBI synonyms such as former species names)."""

    def __init__(self, sci_idx: dict[str, dict], sp_idx: dict[str, dict],
                 syn_idx: dict[str, dict] | None = None, max_flank: int = MAX_FLANK):
        self._matcher = TaxonomyMatcher(sci_idx, sp_idx, syn_idx)
        self._flank = max_flank

    def match(self, title: str) -> tuple[dict | None, str | None, bool]:
        """(record, matched name, via_synonym).  Candidates are tried longest
        first, so a full synonym ("Orchid fleck dichorhavirus" → the species)
        beats a shorter primary name inside it (the genus "Dichorhavirus");
        for one candidate the primary name is preferred over a synonym."""
        for cand in name_candidates(title, self._flank):
            rec, syn = self._matcher.lookup(cand)
            if rec:
                return rec, cand, syn
        return None, None, False


def _best_hit(seq: NucSequence) -> BlastnResult | None:
    return max(seq.blastn_hits, key=lambda h: h.bit_score) if seq.blastn_hits else None


class BlastnTaxonomyResolver:
    """
    Sets ``seq.blastn_taxonomy`` for every sequence with a BLASTn hit, from
    the virus name in its best hit's title.  The indices are the ones
    TaxonomyAnnotator loaded (``.indices``), so the bundle is read once.
    """

    def __init__(self, sci_idx: dict[str, dict], sp_idx: dict[str, dict],
                 taxid_idx: dict[int, dict] | None = None, syn_idx: dict[str, dict] | None = None,
                 max_flank: int = MAX_FLANK):
        self._title = TitleMatcher(sci_idx, sp_idx, syn_idx, max_flank)

    def resolve_hit(self, hit: BlastnResult) -> BlastnTaxonomy:
        rec, name, via_syn = self._title.match(hit.stitle)
        if rec:
            return BlastnTaxonomy(status="resolved", method="title", accession=hit.accession,
                                  subject=hit.stitle, taxid=rec.get("TaxId"), matched_name=name,
                                  synonym_of=rec.get("ScientificName") if via_syn else None,
                                  taxonomy=record_to_taxonomy(rec))
        return BlastnTaxonomy(status="unresolved", method=None, accession=hit.accession, subject=hit.stitle)

    def annotate(self, nuc_seqs: list[NucSequence]) -> dict[str, int]:
        counts = {"resolved": 0, "unresolved": 0}
        for seq in nuc_seqs:
            hit = _best_hit(seq)
            if hit is None:
                continue
            seq.blastn_taxonomy = self.resolve_hit(hit)
            counts[seq.blastn_taxonomy.status] += 1
        logger.success(
            f"BLASTn hit taxonomy: {counts['resolved']} resolved from the hit title, "
            f"{counts['unresolved']} unresolved."
        )
        return counts
