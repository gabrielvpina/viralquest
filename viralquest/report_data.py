"""
report_data.py — Assemble the single ``VQ_REPORT`` structure injected into the
consolidated HTML, from a list of loaded :class:`SampleReport` objects.

Shape produced by :func:`build_report_data`::

    {
      "meta":      {generated, viralquest_version, n_samples, input_root},
      "samples":   [ <per-sample overview summary, incl. salmon rollup>, ... ],
      "sequences": [ <every stamped sequence, for the viewer>, ... ],
      "clusters":  [ <cross-sample clusters>, ... ],   # filled in stage 2
      "has_clusters": bool,
    }

The per-sample *overview* summaries are intentionally compact (counts and
small arrays) so Section 1 can render many samples without re-scanning the
full sequence list client-side.  The full ``sequences`` array is still shipped
because the Sequence Viewer (Section 3) needs the raw records.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from statistics import median

from viralquest.report_loader import SampleReport


# ── Per-sample step detection ──────────────────────────────────────────────

def _detect_steps(s: SampleReport) -> dict[str, bool]:
    """Which optional pipeline stages ran for this sample (for the uniformity matrix)."""
    ps    = s.pipeline_stats
    blast = ps.get("blast", {}) or {}
    seqs  = s.sequences

    nr = bool(blast.get("nr_unique_seqs")) or any(seq.get("blastx_nr_hits") for seq in seqs)
    blastn = bool(blast.get("blastn_unique_seqs")) or any(seq.get("blastn_hits") for seq in seqs)
    salmon = s.salmon_quant is not None
    llm = bool((ps.get("llm", {}) or {}).get("present")) or any(seq.get("llm_output") for seq in seqs)
    cap3 = bool((ps.get("cap3", {}) or {}).get("used"))

    return {"nr": nr, "blastn": blastn, "salmon": salmon, "llm": llm, "cap3": cap3}


def _best_species(seq: dict) -> str | None:
    """Best BLASTx species: NR first hit, fallback RefSeq first hit."""
    for key in ("blastx_nr_hits", "blastx_hits"):
        hits = seq.get(key) or []
        if hits:
            return hits[0].get("species") or None
    return None


def _summarize_salmon(s: SampleReport) -> dict | None:
    """
    Kingdom-level rollup of a sample's Salmon quantification.

    Deliberately a rollup and never the rows: ``conserved_quant`` carries ~11k
    bundled housekeeping genes per sample and only the handful belonging to the
    real host ever collects reads, so shipping the table would bloat the HTML
    for data the report only ever aggregates.

    ``detected`` (genes with TPM > 0) is the load-bearing number, not a measure
    of position — a median over a kingdom is 0.0 for every kingdom, and a median
    over just the detected genes lets a single stray transcript outrank the
    actual host.  Totals and detection counts are what separate the two.
    """
    sq = s.salmon_quant
    if not sq:
        return None

    def _tpm(entry: dict) -> float:
        return entry.get("tpm") or 0.0

    conserved: dict[str, dict] = {}
    for entry in sq.get("conserved_quant", []) or []:
        kingdom = entry.get("kingdom") or "Unknown"
        row = conserved.setdefault(kingdom, {"n": 0, "detected": 0, "tpm_sum": 0.0})
        row["n"] += 1
        tpm = _tpm(entry)
        row["tpm_sum"] += tpm
        if tpm > 0:
            row["detected"] += 1
    for row in conserved.values():
        row["tpm_sum"] = round(row["tpm_sum"], 3)

    ref_hk    = [_tpm(e) for e in (sq.get("ref_hk_quant") or [])]
    ref_found = sorted(v for v in ref_hk if v > 0)

    return {
        "pathway":       sq.get("pathway"),
        "mapping_rate":  sq.get("mapping_rate"),
        "total_reads":   sq.get("total_reads"),
        "viral_tpm_sum": round(sum(_tpm(e) for e in (sq.get("viral_quant") or [])), 3),
        "conserved":     conserved,
        "ref_hk": {
            "n":        len(ref_hk),
            "detected": len(ref_found),
            "tpm_sum":  round(sum(ref_hk), 3),
            "median":   round(median(ref_found), 3) if ref_found else None,
        },
    }


def _summarize_sample(s: SampleReport) -> dict:
    """Compact overview record for one sample."""
    ps    = s.pipeline_stats
    blast = ps.get("blast", {}) or {}
    seqs  = s.sequences

    families = Counter(
        (seq.get("taxonomy") or {}).get("family") or "Unclassified"
        for seq in seqs
    )
    lengths = [seq.get("length") for seq in seqs if seq.get("length")]

    heur_scores = [
        seq["heuristic_output"]["vq_score"]
        for seq in seqs
        if (seq.get("heuristic_output") or {}).get("vq_score") is not None
    ]
    llm_scores = [
        seq["llm_output"]["vq_score"]
        for seq in seqs
        if (seq.get("llm_output") or {}).get("vq_score") is not None
        and seq["llm_output"].get("classification") != "api-error"
    ]

    return {
        "sample":       s.sample,
        "n_input":      blast.get("total_input"),
        "n_viral":      blast.get("total_viral_flagged"),
        "n_confirmed":  blast.get("total_confirmed", len(seqs)),
        "families":     dict(families),
        "lengths":      lengths,
        "heuristic_scores": heur_scores,
        "llm_scores":   llm_scores,
        "n_clusters":   0,   # filled in stage 2 (cross-sample clusters per sample)
        "steps":        _detect_steps(s),
        "salmon":       _summarize_salmon(s),
    }


# ── Public entry point ─────────────────────────────────────────────────────

def build_report_data(
    samples:    list[SampleReport],
    input_root: str = "",
    version:    str = "unknown",
    clusters:   list[dict] | None = None,
) -> dict:
    """Build the ``VQ_REPORT`` dict from loaded samples (+ optional clusters)."""
    clusters = clusters or []

    sequences: list[dict] = []
    for s in samples:
        sequences.extend(s.sequences)

    summaries = [_summarize_sample(s) for s in samples]

    # Per-sample cross-sample cluster counts (stage 2 populates `clusters`).
    if clusters:
        per_sample = Counter()
        for c in clusters:
            for sample in c.get("samples", []) or []:
                per_sample[sample] += 1
        for summary in summaries:
            summary["n_clusters"] = per_sample.get(summary["sample"], 0)

    return {
        "meta": {
            "generated":          datetime.now(timezone.utc).isoformat(),
            "viralquest_version": version,
            "n_samples":          len(samples),
            "input_root":         str(input_root),
        },
        "samples":      summaries,
        "sequences":    sequences,
        "clusters":     clusters,
        "has_clusters": bool(clusters),
    }
