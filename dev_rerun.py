#!/usr/bin/env python3
"""
dev_rerun.py — Reload precomputed ViralQuest results and run only the final
pipeline steps without re-running Diamond, HMMsearch, or BLASTn.

Steps skipped (loaded from precomputed/):
  diamond/refseq.tsv   — RefSeq BLASTx hits
  diamond/nr.tsv       — NR BLASTx hits (optional)
  hmm/RVDB.tsv etc.    — HMM domain hits
  blastn/blastn.tsv    — BLASTn hits (optional)

Steps executed:
  Parse FASTA → ORF finding → attach cached results →
  Taxonomy → Clustering → LLM (optional) → JSON/HTML export

Usage:
  python dev_rerun.py -in data/test.fasta -pre testFast -out dev_out
  python dev_rerun.py -in data/test.fasta -pre testFast -out dev_out \\
      --model-type ollama --model-name gemma4:e2b --llm-tokens low
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT))

from loguru import logger

from viralquest.cli import _DB, _HMM_FILTER, _INDEX_DIR, __version__
from viralquest.biodata import BlastnResult
from viralquest.diamond import (DiamondOutputParser, DiamondResultAttacher,
                                 DiamondPhase)
from viralquest.exporter import ReportExporter
from viralquest.hmm import (HmmMetadataLoader, HmmResultAttacher,
                              HmmSequencePreparer, HmmViralFlagSetter)
from viralquest.html_report import write_report
from viralquest.blastn import BlastnResultAttacher
from viralquest.output import OutputOrganizer
from viralquest.orfs import OrfAnalyzer
from viralquest.parser import FastaParser
from viralquest.tax import TaxonomyAnnotator, ViralFamilyAnnotator
from viralquest.track_seqs import SequenceTracker


# ---------------------------------------------------------------------------
# TSV loaders for OutputOrganizer-format files
# ---------------------------------------------------------------------------

def _load_blastn_tsv(path: Path) -> list[BlastnResult]:
    """
    Parse the 8-column OutputOrganizer BLASTn TSV (has header row):
      qseqid  qlen  slen  qcovhsp  pident  evalue  bit_score  stitle
    """
    hits: list[BlastnResult] = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                hits.append(BlastnResult(
                    qseqid    = row["qseqid"],
                    qlen      = int(row["qlen"]),
                    slen      = int(row["slen"]),
                    qcovhsp   = int(float(row["qcovhsp"])),
                    pident    = float(row["pident"]),
                    evalue    = float(row["evalue"]),
                    bit_score = float(row["bit_score"]),
                    stitle    = row["stitle"],
                ))
            except (KeyError, ValueError) as exc:
                logger.warning(f"Skipping malformed BLASTn row: {exc}")
    return hits


def _load_hmm_tsv(path: Path) -> list[tuple]:
    """
    Parse the 6-column OutputOrganizer HMM TSV (has header row):
      hmm_profile  orf_id  score  i_evalue  env_from  env_to
    Returns raw hit tuples compatible with HmmResultAttacher.attach().
    """
    hits: list[tuple] = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                hits.append((
                    row["hmm_profile"],
                    row["orf_id"],
                    float(row["score"]),
                    float(row["i_evalue"]),
                    int(row["env_from"]),
                    int(row["env_to"]),
                ))
            except (KeyError, ValueError) as exc:
                logger.warning(f"Skipping malformed HMM row: {exc}")
    return hits


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dev_rerun",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-in",  "--input",       dest="input",       required=True,
                   help="Original input FASTA file.")
    p.add_argument("-pre", "--precomputed", dest="precomputed", required=True,
                   help="Directory containing precomputed results (diamond/, hmm/, blastn/).")
    p.add_argument("-out", "--outdir",      dest="outdir",      required=True,
                   help="Output directory for final reports.")
    p.add_argument("-cpu", "--cpu",         dest="cpu",         type=int, default=2,
                   help="CPU threads (default: 2).")
    p.add_argument("--min-identity",        dest="min_identity", type=float, default=70.0,
                   help="Min %% identity for clustering (default: 70.0).")
    p.add_argument("--force",               action="store_true",
                   help="Export all sequences, not only confirmed viral ones.")
    # Salmon (optional)
    p.add_argument("--transcriptome", dest="transcriptome", type=str, metavar="HOST.fasta",
                   help="Host transcriptome FASTA for Salmon quantification.")
    p.add_argument("--reads",         dest="reads", nargs="+", type=str, metavar="READS.fastq",
                   help="FASTQ file(s): one = single-end, two = paired-end.")
    # LLM (optional)
    p.add_argument("--model-type",  dest="model_type",
                   choices=["ollama", "openai", "anthropic", "google"])
    p.add_argument("--model-name",  dest="model_name",  type=str)
    p.add_argument("--llm-tokens",  dest="llm_tokens",  choices=["high", "low"])
    p.add_argument("--api-key",     dest="api_key",     type=str)
    return p


def main() -> None:
    args   = _build_parser().parse_args()

    if bool(getattr(args, "transcriptome", None)) ^ bool(getattr(args, "reads", None)):
        logger.error("--transcriptome and --reads must be provided together.")
        sys.exit(1)

    pre    = Path(args.precomputed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem      = Path(args.input).stem
    organizer = OutputOrganizer(outdir, stem)
    log_id    = organizer.add_log_sink()

    # ── 1. Parse FASTA ────────────────────────────────────────────────────────
    logger.info("Parsing FASTA...")
    fp = FastaParser(args.input)
    fp.read_input_file()
    seqs        = fp.sequences
    input_fasta = fp.input_fasta
    logger.success(f"Parsed {len(seqs)} sequences.")

    # ── 2. ORF finding ────────────────────────────────────────────────────────
    logger.info("Finding ORFs...")
    analyzer = OrfAnalyzer()
    for seq in seqs:
        analyzer.find_n_save_orfs(seq)
    logger.success("ORF finding done.")

    # build orf_map once — shared by all HMM attach calls
    _, orf_map = HmmSequencePreparer.prepare(seqs)

    # ── 3. Diamond RefSeq hits ────────────────────────────────────────────────
    refseq_tsv = pre / "diamond" / "refseq.tsv"
    logger.info(f"Loading Diamond RefSeq → {refseq_tsv}")
    hits = DiamondOutputParser.parse(refseq_tsv)
    DiamondResultAttacher.attach(hits, seqs, DiamondPhase.REFSEQ_FILTER)
    logger.success(f"RefSeq: {len(hits)} hits attached.")

    # ── 4. HMM viral confirmation (RVDB / Vfam / EggNOG) ─────────────────────
    for _db_key, json_path, db_name in _HMM_FILTER:
        tsv_path = pre / "hmm" / f"{db_name}.tsv"
        logger.info(f"Loading {db_name} HMM → {tsv_path}")
        meta     = HmmMetadataLoader.load(str(json_path), db_name)
        hmm_hits = _load_hmm_tsv(tsv_path)
        HmmResultAttacher.attach(hmm_hits, orf_map, meta, db_name)
        logger.success(f"{db_name}: {len(hmm_hits)} hits attached.")
    HmmViralFlagSetter.flag(seqs)

    # ── 5. Diamond NR hits (optional) ─────────────────────────────────────────
    nr_tsv = pre / "diamond" / "nr.tsv"
    if nr_tsv.exists() and nr_tsv.stat().st_size > 0:
        logger.info(f"Loading Diamond NR → {nr_tsv}")
        nr_hits = DiamondOutputParser.parse(nr_tsv)
        DiamondResultAttacher.attach(nr_hits, seqs, DiamondPhase.NR_CHARACTERIZE)
        logger.success(f"NR: {len(nr_hits)} hits attached.")

    # ── 6. Pfam hits (optional) ───────────────────────────────────────────────
    pfam_tsv = pre / "hmm" / "Pfam.tsv"
    if pfam_tsv.exists() and pfam_tsv.stat().st_size > 0:
        logger.info(f"Loading Pfam HMM → {pfam_tsv}")
        pfam_meta = HmmMetadataLoader.load(str(_INDEX_DIR / "Pfam-index.json"), "Pfam")
        pfam_hits = _load_hmm_tsv(pfam_tsv)
        HmmResultAttacher.attach(pfam_hits, orf_map, pfam_meta, "Pfam")
        logger.success(f"Pfam: {len(pfam_hits)} hits attached.")

    # ── 7. BLASTn hits (optional) ─────────────────────────────────────────────
    blastn_tsv = pre / "blastn" / "blastn.tsv"
    if blastn_tsv.exists() and blastn_tsv.stat().st_size > 0:
        logger.info(f"Loading BLASTn → {blastn_tsv}")
        blastn_hits = _load_blastn_tsv(blastn_tsv)
        BlastnResultAttacher.attach(blastn_hits, seqs)
        logger.success(f"BLASTn: {len(blastn_hits)} hits attached.")

    # ── 8. Taxonomy annotation ────────────────────────────────────────────────
    logger.info("Annotating taxonomy...")
    TaxonomyAnnotator(str(_DB["viral_tax"])).annotate(seqs)
    ViralFamilyAnnotator(str(_DB["fam_high"]), str(_DB["fam_low"])).annotate(seqs)

    # ── 9. Clustering ─────────────────────────────────────────────────────────
    logger.info("Clustering sequences...")
    tracker      = SequenceTracker(min_identity=args.min_identity)
    viral_for_cl = [s for s in seqs if s.is_viral]
    clusters     = tracker.track(viral_for_cl)

    # ── 10. LLM scoring (optional) ────────────────────────────────────────────
    if args.model_type and args.model_name and args.llm_tokens:
        from viralquest.score_ai import LlmMode, SequenceScorer
        mode       = LlmMode.HIGH if args.llm_tokens == "high" else LlmMode.LOW
        viral_seqs = [s for s in seqs if s.blastx_nr_hits]
        logger.info(
            f"LLM scoring ({args.llm_tokens}-token) — "
            f"{args.model_type}/{args.model_name} — {len(viral_seqs)} sequences"
        )
        SequenceScorer(
            model_type = args.model_type,
            model_name = args.model_name,
            mode       = mode,
            api_key    = args.api_key,
        ).score(viral_seqs)

    # ── 11. Salmon quantification (optional) ─────────────────────────────────
    salmon_report = None
    if args.transcriptome and args.reads:
        from viralquest.salmon_quant import SalmonQuantPipeline
        logger.info("Running Salmon quantification...")
        try:
            salmon_report = SalmonQuantPipeline().run(
                user_transcriptome = Path(args.transcriptome),
                reads              = args.reads,
                viral_seqs         = seqs,
                outdir             = outdir / "salmon",
            )
        except Exception as exc:
            logger.error(f"Salmon quantification failed — skipping: {exc}")

    # ── 12. Export viral FASTA + JSON + HTML ──────────────────────────────────
    logger.info("Exporting reports...")
    viral_fasta = organizer.save_viral_contigs(seqs)
    json_path   = outdir / f"{stem}_viralquest.json"
    report      = ReportExporter(force=args.force).export(
        nuc_seqs      = seqs,
        clusters      = clusters,
        input_fasta   = input_fasta,
        output_path   = json_path,
        version       = __version__,
        salmon_report = salmon_report,
    )
    html_path = outdir / f"{stem}_viralquest.html"
    write_report(report, html_path)

    viral_count = sum(1 for s in seqs if s.is_viral)
    logger.success(
        f"Done — {viral_count} viral sequences, {len(clusters)} clusters.\n"
        f"  JSON  → {json_path}\n"
        f"  HTML  → {html_path}\n"
        f"  FASTA → {viral_fasta}"
    )
    OutputOrganizer.remove_log_sink(log_id)


if __name__ == "__main__":
    main()
