"""
score_heuristic.py — deterministic, rule-based viral sequence scoring (no LLM).

Produces the same metric as score_ai.SequenceScorer (a 0-100 ``vq_score`` plus a
classification) but computes it from explicit rules over the alignment / HMM
evidence already attached to each ``NucSequence``. It needs no API key, no NR
database, and runs instantly, so the pipeline scores every confirmed sequence
by default; the LLM scorer remains an optional layer on top.

Scoring model
-------------
Three components, each on a 0-100 scale, combined by a weighted average:

    blastn 30%   — best *viral* BLASTn hit (identity x coverage)
    blastx 30%   — best BLASTx hit          (identity x coverage)
    hmm    40%   — FILTER-database HMM evidence (best bit score + hit multiplicity)

Several rules shape the raw average, and they are kept deliberately separate:

  * Renormalisation — a component that is *absent* (e.g. BLASTn was never run)
    is dropped and its weight redistributed over the present components, so the
    score is never artificially capped. Absent is not the same as zero.

  * Low-coverage penalty — applied *inside* each BLAST component. A hit covering
    only a tiny fraction of the query is damped, so high identity over a short
    footprint can no longer carry the component.

  * Single-evidence penalty — applied *after* the weighted sum. A score that
    rests on a single component (e.g. HMM only) is damped, because the other two
    evidence channels are missing.

  * False-positive penalty — applied *after* the weighted sum. If the dominant
    BLASTn hit is to a clearly NON-viral subject at high identity and coverage,
    the whole score is damped multiplicatively. This can override otherwise
    strong BLASTx / HMM signal, flagging a likely host / contaminant contig.

Only the false-positive penalty affects the *classification*; the low-coverage
and single-evidence penalties lower the numeric score alone.

Because RefSeq BLASTx hits are viral by construction and NR BLASTx hits are
pre-filtered to viral subjects upstream (diamond.py), the viral-subject check
and the penalty effectively act on BLASTn — which is exactly where a
non-standard or general user database can introduce false positives.

All thresholds, weights, and the penalty magnitude are module-level constants
grouped under "CALIBRATION" below, so they can be tuned against reference data.
"""

import re

from loguru import logger

from viralquest.biodata import (
    BlastnResult,
    BlastxResult,
    HeuristicScore,
    NucSequence,
    ScoreComponents,
)

# ===========================================================================
# CALIBRATION — all tunable parameters live here
# ===========================================================================

# Base component weights (rubric: blastn 30 / blastx 30 / hmm 40). Renormalised
# over whichever components are present for a given sequence.
WEIGHTS: dict[str, float] = {"blastn": 0.30, "blastx": 0.30, "hmm": 0.40}

# Within an alignment component, how identity and coverage are mixed.
ID_FRAC:  float = 0.6
COV_FRAC: float = 0.4

# Low-coverage penalty (BLASTx and BLASTn). A hit covering less than
# LOWCOV_COV_MIN of the query is weak evidence no matter how high its identity
# is — high identity over a tiny footprint says little about the whole contig.
# When coverage is below the threshold the alignment component is multiplied by
# LOWCOV_PENALTY. This is baked straight into the component value (so it shows up
# in components.blastn / components.blastx), NOT in components.penalty.
LOWCOV_COV_MIN: float = 30.0
LOWCOV_PENALTY: float = 0.5

# Single-evidence penalty. A score resting on a single present component (e.g.
# HMM only, with no BLASTn/BLASTx corroboration) is damped by
# SINGLE_EVIDENCE_PENALTY because the other two evidence channels are missing.
# Unlike the false-positive penalty this only lowers the score — it does NOT
# change the classification.
SINGLE_EVIDENCE_PENALTY: float = 0.7

# Identity ramps (value <= lo scores 0, >= hi scores 100, linear between).
# BLASTx is protein-level so meaningful identity is lower than for nucleotide.
BX_ID_LO, BX_ID_HI = 20.0, 90.0     # BLASTx percent identity
BN_ID_LO, BN_ID_HI = 70.0, 98.0     # BLASTn percent identity
# Shared query-coverage ramp for both BLAST components.
COV_LO,   COV_HI   = 10.0, 80.0

# HMM component.
FILTER_DBS = {"RVDB", "Vfam", "EggNOG"}     # viral-confirmation databases
HMM_SCORE_LO, HMM_SCORE_HI = 50.0, 150.0    # best FILTER-domain bit score ramp
HMM_COUNT_LO, HMM_COUNT_HI = 1.0,  6.0      # FILTER-hit multiplicity ramp
HMM_QUAL_W,   HMM_QUANT_W  = 0.7,  0.3      # qualitative vs quantitative mix

# Classification: "viral-known" demands near-species-level evidence.
KNOWN_ID, KNOWN_COV = 90.0, 70.0

# False-positive penalty: a dominant BLASTn hit to a non-viral subject at this
# strength damps the final score by NON_VIRAL_PENALTY (heavy, can override).
FP_ID, FP_COV       = 90.0, 80.0
NON_VIRAL_PENALTY   = 0.2

# Subject-title viral-term matcher. Searches anywhere in the title (not only
# inside NCBI brackets) so non-standard BLASTn databases are still classified.
_VIRAL_RE = re.compile(
    r"vir(?:us|al|idae|ales|inae|oid|ion|aceae)|phage|bacteriophage",
    re.IGNORECASE,
)

# Tokens that mark the end of a species name in a BLASTn subject title.
_SPECIES_CUT_RE = re.compile(
    r"\b(?:segment|isolate|strain|gene|genes|complete|partial|polyprotein|"
    r"chromosome|clone|mrna|cds|protein|nonstructural|structural|rna|dna|"
    r"genome|sequence|assembly|scaffold|contig)\b",
    re.IGNORECASE,
)

_VALID_CLASSIFICATIONS = {"viral-known", "viral-unknown", "non-viral"}


# ===========================================================================
# Small numeric / text helpers
# ===========================================================================

def _ramp(x: float, lo: float, hi: float) -> float:
    """Linear ramp: 0 at/below lo, 100 at/above hi, interpolated between."""
    if hi <= lo:
        return 100.0 if x >= hi else 0.0
    if x <= lo:
        return 0.0
    if x >= hi:
        return 100.0
    return (x - lo) / (hi - lo) * 100.0


def _alignment_component(
    pident: float, qcov: float,
    id_lo: float, id_hi: float, cov_lo: float, cov_hi: float,
) -> float:
    """
    0-100 from a single alignment's percent identity and query coverage.

    A hit whose query coverage is below LOWCOV_COV_MIN is damped by
    LOWCOV_PENALTY: high identity over a tiny footprint is weak evidence about
    the whole contig, so identity alone must not carry the component.
    """
    comp = ID_FRAC * _ramp(pident, id_lo, id_hi) + COV_FRAC * _ramp(qcov, cov_lo, cov_hi)
    if qcov < LOWCOV_COV_MIN:
        comp *= LOWCOV_PENALTY
    return comp


def _looks_viral(title: str) -> bool:
    """True if the subject title contains a viral / phage term anywhere."""
    return bool(title) and bool(_VIRAL_RE.search(title))


def _extract_species(title: str) -> str:
    """
    Light, deterministic species-name extraction from a BLASTn subject title.
    Drops a leading accession token, cuts at the first descriptor keyword and at
    the first comma / parenthesis. Best-effort only — kept minimal on purpose.
    """
    if not title:
        return ""
    t = title.strip()
    # strip a leading "ACCESSION " token (e.g. "NC_001498.1 Measles virus ...")
    parts = t.split(None, 1)
    if len(parts) == 2 and re.fullmatch(r"[A-Za-z]{1,4}_?\d+(?:\.\d+)?", parts[0]):
        t = parts[1]
    t = re.split(r"[(,]", t, maxsplit=1)[0]
    m = _SPECIES_CUT_RE.search(t)
    if m:
        t = t[: m.start()]
    return t.strip()


# ===========================================================================
# Per-sequence evidence selection
# ===========================================================================

def _best_blastn_overall(seq: NucSequence) -> BlastnResult | None:
    """Strongest BLASTn hit by bit score, viral or not — used for FP detection."""
    return max(seq.blastn_hits, key=lambda h: h.bit_score, default=None)


def _best_viral_blastn(seq: NucSequence) -> BlastnResult | None:
    """Strongest BLASTn hit whose subject looks viral — used for positive signal."""
    viral = [h for h in seq.blastn_hits if _looks_viral(h.stitle)]
    return max(viral, key=lambda h: h.bit_score, default=None)


def _hmm_component(seq: NucSequence) -> float | None:
    """
    FILTER-database HMM evidence as 0-100, or None if no FILTER domain is present.
    Combines the best domain bit score (qualitative) with the raw hit
    multiplicity across FILTER databases (quantitative).
    """
    filter_domains = [d for d in seq.hmm_domains if d.database in FILTER_DBS]
    raw            = seq.hmm_raw_counts
    filter_count   = sum(raw.get(db, 0) for db in FILTER_DBS)
    if not filter_domains and filter_count == 0:
        return None
    best_score = max((d.score for d in filter_domains), default=0.0)
    qual = _ramp(best_score, HMM_SCORE_LO, HMM_SCORE_HI)
    quant = _ramp(filter_count, HMM_COUNT_LO, HMM_COUNT_HI)
    return HMM_QUAL_W * qual + HMM_QUANT_W * quant


# ===========================================================================
# Analysis text
# ===========================================================================

def _build_analysis(
    seq: NucSequence,
    components: dict[str, float],
    applied: dict[str, float],
    fp_penalty: float,
    single_penalty: float,
    fp_hit: BlastnResult | None,
    viral_bn: BlastnResult | None,
    bx: BlastxResult | None,
) -> str:
    def _lowcov_note(qcov: float) -> str:
        return (f" — low-coverage penalty x{LOWCOV_PENALTY:g}"
                if qcov < LOWCOV_COV_MIN else "")

    parts: list[str] = []
    if "blastx" in components and bx is not None:
        parts.append(
            f"BLASTx best hit {bx.pct_identity:.1f}% identity / "
            f"{bx.query_coverage:.1f}% coverage (component {components['blastx']:.0f}"
            f"{_lowcov_note(bx.query_coverage)})."
        )
    if "blastn" in components and viral_bn is not None:
        parts.append(
            f"Viral BLASTn best hit {viral_bn.pident:.1f}% identity / "
            f"{viral_bn.qcovhsp:.0f}% coverage (component {components['blastn']:.0f}"
            f"{_lowcov_note(viral_bn.qcovhsp)})."
        )
    if "hmm" in components:
        counts = ", ".join(
            f"{db}:{n}" for db, n in seq.hmm_raw_counts.items() if db in FILTER_DBS and n
        ) or "none above threshold"
        parts.append(
            f"HMM FILTER evidence (component {components['hmm']:.0f}); raw hits {counts}."
        )
    if not parts:
        parts.append("No viral alignment or HMM evidence available.")
    if single_penalty < 1.0 and components:
        only = next(iter(components)).upper()
        parts.append(
            f"Single-evidence penalty x{single_penalty:g} applied: score rests on the "
            f"{only} component alone, with no corroboration from the other two channels."
        )
    if fp_penalty < 1.0 and fp_hit is not None:
        parts.append(
            f"False-positive penalty x{fp_penalty:g} applied: dominant BLASTn hit is "
            f"non-viral at {fp_hit.pident:.1f}% identity / {fp_hit.qcovhsp:.0f}% coverage "
            f"('{fp_hit.stitle[:80]}') — possible false positive."
        )
    if applied:
        w = ", ".join(f"{k} {applied[k]*100:.0f}%" for k in applied)
        parts.append(f"Weights applied: {w}.")
    return " ".join(parts)


# ===========================================================================
# Classification
# ===========================================================================

def _classify(
    present: dict[str, float],
    viral_bn: BlastnResult | None,
    bx: BlastxResult | None,
    fp_penalty: float,
) -> str:
    # A triggered false-positive penalty overrides everything. Only the
    # false-positive penalty reaches here — the low-coverage and single-evidence
    # penalties lower the numeric score but never change the classification.
    if fp_penalty < 1.0:
        return "non-viral"
    # No viral evidence of any kind.
    if not present:
        return "non-viral"
    # Species-level evidence in either BLAST channel → known.
    if viral_bn is not None and viral_bn.pident >= KNOWN_ID and viral_bn.qcovhsp >= KNOWN_COV:
        return "viral-known"
    if bx is not None and bx.pct_identity >= KNOWN_ID and bx.query_coverage >= KNOWN_COV:
        return "viral-known"
    return "viral-unknown"


# ===========================================================================
# Scorer
# ===========================================================================

class HeuristicScorer:
    """Scores NucSequence objects deterministically, attaching a HeuristicScore."""

    def score(self, nuc_seqs: list[NucSequence]) -> list[HeuristicScore]:
        """
        Scores each sequence, sets seq.heuristic_output, and returns the
        HeuristicScore list in the same order.
        """
        outputs: list[HeuristicScore] = []
        for i, seq in enumerate(nuc_seqs):
            result = self._score_one(seq)
            seq.heuristic_output = result
            outputs.append(result)
            logger.info(
                f"[{i+1}/{len(nuc_seqs)}] {seq.id}: "
                f"vq_score={result.vq_score}, class={result.classification}"
            )
        return outputs

    def _score_one(self, seq: NucSequence) -> HeuristicScore:
        bx       = seq.best_blastx                       # viral by construction
        viral_bn = _best_viral_blastn(seq)               # positive BLASTn signal
        fp_hit   = _best_blastn_overall(seq)             # dominant hit for FP check

        # --- components (None = absent → excluded from the weighted average) ---
        present: dict[str, float] = {}
        if bx is not None:
            present["blastx"] = _alignment_component(
                bx.pct_identity, bx.query_coverage, BX_ID_LO, BX_ID_HI, COV_LO, COV_HI
            )
        if viral_bn is not None:
            present["blastn"] = _alignment_component(
                viral_bn.pident, viral_bn.qcovhsp, BN_ID_LO, BN_ID_HI, COV_LO, COV_HI
            )
        hmm_comp = _hmm_component(seq)
        if hmm_comp is not None:
            present["hmm"] = hmm_comp

        # --- weighted average with renormalisation over present components ---
        if present:
            total_w = sum(WEIGHTS[k] for k in present)
            applied = {k: WEIGHTS[k] / total_w for k in present}
            raw_score = sum(present[k] * applied[k] for k in present)
        else:
            applied = {}
            raw_score = 0.0

        # --- penalties applied AFTER the weighted sum (multiplicative) ---
        # (the low-coverage penalty is already baked into the BLAST components.)
        #
        # false-positive: dominant BLASTn hit is strongly non-viral → host /
        # contaminant. This one also forces the "non-viral" classification.
        fp_penalty = (
            NON_VIRAL_PENALTY
            if (fp_hit is not None and not _looks_viral(fp_hit.stitle)
                and fp_hit.pident >= FP_ID and fp_hit.qcovhsp >= FP_COV)
            else 1.0
        )
        # single-evidence: score rests on a single component (no corroboration).
        single_penalty = SINGLE_EVIDENCE_PENALTY if len(present) == 1 else 1.0
        # components.penalty stores the combined post-sum dampers (1.0 = none).
        penalty = fp_penalty * single_penalty

        vq_score = int(max(0, min(100, round(raw_score * penalty))))
        classification = _classify(present, viral_bn, bx, fp_penalty)
        blastn_species = _extract_species(viral_bn.stitle) if viral_bn else ""
        analysis = _build_analysis(
            seq, present, applied, fp_penalty, single_penalty, fp_hit, viral_bn, bx
        )

        return HeuristicScore(
            seq_id=seq.id,
            vq_score=vq_score,
            classification=classification,
            blastn_species=blastn_species,
            analysis=analysis,
            components=ScoreComponents(
                blastn=round(present["blastn"], 2) if "blastn" in present else None,
                blastx=round(present["blastx"], 2) if "blastx" in present else None,
                hmm=round(present["hmm"], 2) if "hmm" in present else None,
                weights={k: round(v, 4) for k, v in applied.items()},
                penalty=penalty,
            ),
        )
