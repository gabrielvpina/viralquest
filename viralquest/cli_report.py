"""
cli_report.py — ``viralquest-report`` console entry point.

Consolidates several ViralQuest result directories into a single
self-contained HTML report.

    viralquest-report -in my_results -out final_report

The output directory contains the HTML plus (from stage 2 onward) the
``clusters/``, ``quantification/`` and ``blastn/`` artefact folders.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from loguru import logger

from viralquest.report_builder import write_report
from viralquest.report_clusters import (
    FLOOR_COVERAGE,
    FLOOR_IDENTITY,
    build_general_clusters,
    write_artifacts,
)
from viralquest.report_data import build_report_data
from viralquest.report_loader import load_samples

__version__ = "3.0.1"

HTML_NAME = "viralquest_report.html"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viralquest-report",
        description="Consolidate multiple ViralQuest result directories into a single HTML report.",
    )
    req = parser.add_argument_group("required")
    req.add_argument(
        "-in", "--input", dest="input", type=str, required=True,
        help="Directory containing ViralQuest result directories (each with a *_viralquest.json).",
    )
    req.add_argument(
        "-out", "--outdir", dest="outdir", type=str, required=True,
        help="Output directory for the consolidated report.",
    )
    clu = parser.add_argument_group("cross-sample clustering")
    clu.add_argument(
        "--min-identity", dest="min_identity", type=float, default=FLOOR_IDENTITY,
        help=f"BLASTn identity floor for cross-sample clusters (default {FLOOR_IDENTITY}).",
    )
    clu.add_argument(
        "--min-coverage", dest="min_coverage", type=float, default=FLOOR_COVERAGE,
        help=f"BLASTn coverage floor for cross-sample clusters (default {FLOOR_COVERAGE}).",
    )
    clu.add_argument(
        "--blastn-bin", dest="blastn_bin", type=str, default=None,
        help="blastn binary (default: the pixi installation, then PATH).",
    )
    clu.add_argument(
        "--no-clusters", dest="no_clusters", action="store_true",
        help="Skip cross-sample clustering (Overview + Viewer only).",
    )
    parser.add_argument(
        "--version", action="version", version=f"viralquest-report {__version__}",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    # Same binary resolution as the main pipeline: the pixi environment is put
    # on PATH so blastn (and anything else shelled out to) comes from there.
    from viralquest.setup_env import activate_pixi_env, resolve_tool

    activate_pixi_env()
    if not args.no_clusters and resolve_tool("blastn", args.blastn_bin) is None:
        logger.error(
            "blastn not found (pixi environment or PATH). "
            "Run 'viralquest-setup' to install it, pass --blastn-bin with an "
            "explicit path, or use --no-clusters to skip cross-sample clustering."
        )
        return 1

    in_root = Path(args.input)
    out_dir = Path(args.outdir)

    try:
        samples = load_samples(in_root)
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(str(exc))
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)

    # Merge sequences once; used both for clustering and as the viewer source.
    all_sequences: list[dict] = []
    for s in samples:
        all_sequences.extend(s.sequences)

    cluster_dicts: list[dict] = []
    if not args.no_clusters:
        try:
            clusters = build_general_clusters(
                all_sequences,
                blastn_bin=args.blastn_bin,
                min_identity=args.min_identity,
                min_coverage=args.min_coverage,
            )
        except RuntimeError as exc:
            logger.error(f"{exc}  Continuing without cross-sample clusters.")
            clusters = []

        if clusters:
            seq_by_gid = {s["gid"]: s for s in all_sequences}
            write_artifacts(clusters, seq_by_gid, out_dir)
            cluster_dicts = [c.to_dict() for c in clusters]

    report_data = build_report_data(
        samples,
        input_root=str(in_root),
        version=__version__,
        clusters=cluster_dicts,
    )

    html_path = out_dir / HTML_NAME
    write_report(report_data, html_path)
    logger.success(
        f"Consolidated report → '{html_path}'  "
        f"({report_data['meta']['n_samples']} samples)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
