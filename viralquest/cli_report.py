"""
cli_report.py — ``viralquest-report`` console entry point.

Consolidates several ViralQuest result directories into one named dataset:
a self-contained HTML report plus merge-ready data exports.

    viralquest-report -in my_results -out final_report -n PRJNA123456

The output directory mirrors a single run's layout:

    features/blastn/general_blastn.tsv   cross-sample cluster BLASTn (when clustering ran)
    features/synteny/                    general synteny: pooled references + DIAMOND blastp
                                         (see report_synteny.py)
    fastas/                              every viral sequence + one FASTA per novelty tier;
                                         clusters/<VQR_CLU_####>.fasta
    reads/                               coverage_summary.tsv · quantification.tsv ·
                                         host_viral_hits.tsv (Sample column); clusters/<id>.tsv
    viralquest_report.html
    <name>.viralquest.json / .sqlite     (see report_db.py)
    <name>.viralquest.xlsx               every per-sample sheet + Samples + General clusters
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
from viralquest.report_db import output_paths, validate_dataset_name, write_json, write_sqlite
from viralquest.report_loader import load_samples
from viralquest.report_synteny import MIN_COVERAGE as SYN_COVERAGE
from viralquest.report_synteny import MIN_IDENTITY as SYN_IDENTITY
from viralquest.report_synteny import build_general_synteny
from viralquest.seq_table import write_multi_tables


def write_multi_fastas(report_data: dict, fasta_dir: Path, name: str) -> dict[str, Path]:
    """
    fastas/<name>_viral_contigs.fasta with every sequence of the dataset and one
    <name>_<tier>.fasta per heuristic novelty tier present; headers are the
    dataset-wide ids (sample::seq_id).
    """
    from viralquest.biodata import NOVELTY_TIERS

    seqs = [q for q in report_data.get("sequences") or [] if q.get("sequence")]
    fasta_dir.mkdir(parents=True, exist_ok=True)
    groups: dict[str, list[dict]] = {"viral_contigs": seqs}
    for q in seqs:
        tier = (q.get("heuristic_output") or {}).get("novelty") or "unscored"
        groups.setdefault(tier, []).append(q)
    order = ["viral_contigs"] + [t for t in NOVELTY_TIERS if t in groups] + \
            sorted(set(groups) - set(NOVELTY_TIERS) - {"viral_contigs"})
    out = {}
    for key in order:
        path = fasta_dir / f"{name}_{key}.fasta"
        with open(path, "w", encoding="utf-8") as fh:
            for q in groups[key]:
                fh.write(f">{q.get('gid') or q.get('id')}\n{q['sequence']}\n")
        out[key] = path
    logger.info(f"FASTAs → {fasta_dir}/ ({len(seqs)} sequences, {len(out) - 1} novelty tier(s)).")
    return out


def write_multi_read_tables(samples, reads_dir: Path) -> list[Path]:
    """reads/coverage_summary.tsv, quantification.tsv, host_viral_hits.tsv over every
    sample (Sample column); each only when some sample has that data."""
    from viralquest.seq_table import merged_sample_sheets, write_sheet_tsv

    files = {"Genome coverage": "coverage_summary.tsv", "Quantification": "quantification.tsv",
             "Host-viral hits": "host_viral_hits.tsv"}
    written = []
    for sh in merged_sample_sheets(samples):
        if sh.name in files and sh.rows:
            reads_dir.mkdir(parents=True, exist_ok=True)
            written.append(write_sheet_tsv(sh, reads_dir / files[sh.name]))
    return written

__version__ = "3.0.15"

HTML_NAME = "viralquest_report.html"


def _dataset_name(value: str) -> str:
    try:
        return validate_dataset_name(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


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
    req.add_argument(
        "-n", "--name", dest="name", type=_dataset_name, required=True, metavar="NAME",
        help="Name of this dataset (e.g. a BioProject, a host species, a study): 1-64 letters, "
             "digits, '.', '_' or '-'. Stored in the JSON / SQLite exports and in every SQLite row, "
             "so datasets from different viralquest-report runs can be merged later.",
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
    syn = parser.add_argument_group("general synteny")
    syn.add_argument(
        "--no-synteny", dest="no_synteny", action="store_true",
        help="Skip the general synteny (the samples' reference genomes pooled and grouped "
             "by similarity, every sequence aligned to its group's reference).",
    )
    syn.add_argument(
        "--synteny-identity", dest="synteny_identity", type=float, default=SYN_IDENTITY,
        metavar="PCT",
        help=f"aa identity for two reference genomes to share a group (default {SYN_IDENTITY:g}).",
    )
    syn.add_argument(
        "--synteny-coverage", dest="synteny_coverage", type=float, default=SYN_COVERAGE,
        metavar="PCT",
        help=f"Proteome coverage for two reference genomes to share a group (default {SYN_COVERAGE:g}).",
    )
    syn.add_argument(
        "-cpu", "--cpu", dest="cpu", type=int, default=2, metavar="N",
        help="DIAMOND threads for the general synteny (default 2).",
    )
    syn.add_argument(
        "--email", dest="email", type=str, default=None,
        help="E-mail sent to NCBI when a reference GenBank record must be downloaded.",
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

    from viralquest.orfs import needs_orf_recall
    old = [s.sample for s in samples if needs_orf_recall(s.report)]
    if old:
        logger.warning(
            f"{len(old)} sample(s) ran ViralQuest < 3.0.9, whose ORFs are complete only "
            f"({', '.join(old[:5])}{' …' if len(old) > 5 else ''}): gene fragments have no ORF, "
            f"so the synteny misses them. Run 'viralquest --reload <sample dir>' on each first."
        )

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

    general_synteny = None
    if not args.no_synteny:
        if resolve_tool("diamond") is None:
            logger.warning("diamond not found (pixi environment or PATH) — general synteny skipped.")
        else:
            general_synteny = build_general_synteny(
                samples, out_dir / "features" / "synteny", threads=args.cpu, email=args.email,
                min_identity=args.synteny_identity, min_coverage=args.synteny_coverage,
            )

    report_data = build_report_data(
        samples,
        input_root=str(in_root),
        version=__version__,
        clusters=cluster_dicts,
        dataset=args.name,
        general_synteny=general_synteny,
    )

    html_path = out_dir / HTML_NAME
    write_report(report_data, html_path)
    json_path, db_path = output_paths(out_dir, args.name)
    write_json(report_data, json_path)
    write_sqlite(report_data, db_path, sample_reports=samples)
    xlsx_path = out_dir / f"{args.name}.viralquest.xlsx"
    write_multi_tables(samples, report_data, xlsx_path)
    write_multi_fastas(report_data, out_dir / "fastas", args.name)
    write_multi_read_tables(samples, out_dir / "reads")
    logger.success(
        f"Dataset '{args.name}' ({report_data['meta']['n_samples']} samples) → "
        f"'{html_path}', '{json_path.name}', '{db_path.name}', '{xlsx_path.name}'."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
