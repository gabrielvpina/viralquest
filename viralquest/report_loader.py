"""
report_loader.py — Discover and load multiple ViralQuest result directories
into a single in-memory structure for the consolidated ``viralquest-report``.

A "result directory" is any sub-directory of the input root that contains a
``*_viralquest.json`` file (the per-sample report written by ReportExporter).

Two correctness rules are enforced at the load boundary — everything
downstream (overview, cross-sample clusters, viewer filter, quantification)
depends on them:

  1. Sample attribution comes from the result *directory name*, not from the
     per-sequence ``sample_name`` field (which is empty in the real schema)
     nor from ``meta.input_file.name`` (often identical across samples, e.g.
     ``contigs_megahit.fasta``).
  2. Sequence IDs are namespaced as ``{sample}::{id}`` (exposed as ``gid``)
     because assemblers emit colliding IDs (``k141_*``) in every sample.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

# Separator for the global (cross-sample) unique id.
GID_SEP = "::"


@dataclass
class SampleReport:
    """One loaded per-sample ViralQuest JSON, with sample attribution stamped."""

    sample:    str
    directory: Path
    json_path: Path
    report:    dict

    @property
    def sequences(self) -> list[dict]:
        return self.report.get("sequences", []) or []

    @property
    def clusters(self) -> list[dict]:
        return self.report.get("clusters", []) or []

    @property
    def salmon_quant(self) -> dict | None:
        return self.report.get("salmon_quant")

    @property
    def pipeline_stats(self) -> dict:
        return self.report.get("pipeline_stats", {}) or {}


# ── Discovery ──────────────────────────────────────────────────────────────

def _find_report_json(directory: Path) -> Path | None:
    """Return the single ``*_viralquest.json`` in *directory*, or None."""
    matches = sorted(directory.glob("*_viralquest.json"))
    if not matches:
        return None
    if len(matches) > 1:
        logger.warning(
            f"{directory.name}: {len(matches)} '*_viralquest.json' files found; "
            f"using '{matches[0].name}'."
        )
    return matches[0]


def discover_result_dirs(root: Path) -> list[tuple[Path, Path]]:
    """
    Find ViralQuest result directories under *root*.

    Returns a list of ``(directory, json_path)`` tuples, sorted by directory
    name.  Both *root* itself and its immediate sub-directories are inspected,
    so passing either a single result dir or a parent of many works.
    """
    root = Path(root)
    if not root.is_dir():
        raise NotADirectoryError(f"Input is not a directory: {root}")

    found: list[tuple[Path, Path]] = []

    # The root itself may be a single result directory.
    self_json = _find_report_json(root)
    if self_json is not None:
        found.append((root, self_json))

    for child in sorted(p for p in root.iterdir() if p.is_dir()):
        jp = _find_report_json(child)
        if jp is not None:
            found.append((child, jp))

    return found


# ── Loading + stamping ─────────────────────────────────────────────────────

def _salmon_tpm_map(report: dict) -> dict[str, float]:
    """
    Map original sequence id → TPM from ``salmon_quant.viral_quant``.

    Entry names look like ``VQ_VIRAL_{seq_id}``; the prefix is stripped so the
    key matches the sequence ``id``.
    """
    sq = report.get("salmon_quant")
    if not sq:
        return {}
    out: dict[str, float] = {}
    for entry in sq.get("viral_quant", []) or []:
        name = entry.get("name", "")
        seq_id = name[len("VQ_VIRAL_"):] if name.startswith("VQ_VIRAL_") else name
        tpm = entry.get("tpm")
        if seq_id and tpm is not None:
            out[seq_id] = tpm
    return out


def _stamp(report: dict, sample: str) -> None:
    """
    Mutate *report* in place: stamp ``sample``/``gid`` (and ``tpm`` when salmon
    ran) onto every sequence.  The original ``id`` is preserved for display.
    """
    tpm_map = _salmon_tpm_map(report)
    for seq in report.get("sequences", []) or []:
        sid = seq.get("id", "")
        seq["sample"] = sample
        seq["gid"]    = f"{sample}{GID_SEP}{sid}"
        if tpm_map:
            seq["tpm"] = tpm_map.get(sid)


def load_samples(root: str | Path) -> list[SampleReport]:
    """
    Discover, load, and stamp every ViralQuest result directory under *root*.

    Raises FileNotFoundError when no result directories are found.
    """
    root = Path(root)
    discovered = discover_result_dirs(root)
    if not discovered:
        raise FileNotFoundError(
            f"No ViralQuest results found under '{root}'. "
            f"Expected sub-directories containing a '*_viralquest.json' file."
        )

    samples: list[SampleReport] = []
    seen_names: dict[str, int] = {}

    for directory, json_path in discovered:
        try:
            report = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            logger.error(f"Skipping '{json_path}': {exc}")
            continue

        sample = directory.name
        # Guard against duplicate directory names (e.g. root + child collision).
        if sample in seen_names:
            seen_names[sample] += 1
            sample = f"{sample}_{seen_names[sample]}"
        else:
            seen_names[directory.name] = 0

        _stamp(report, sample)
        samples.append(SampleReport(
            sample=sample,
            directory=directory,
            json_path=json_path,
            report=report,
        ))
        logger.info(
            f"Loaded sample '{sample}' "
            f"({len(report.get('sequences', []) or [])} sequences) "
            f"from '{json_path.name}'."
        )

    logger.success(f"Loaded {len(samples)} sample(s).")
    return samples
