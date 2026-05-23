# import time
# from viralquest.parser import FastaParser
# from viralquest.orfs import OrfAnalyzer

# Parser = FastaParser("data/test.fasta")
# Parser.read_input_file()

# orf_anlyzer = OrfAnalyzer(min_len_nt=150)


# for seq in Parser.sequences:

#     orf_anlyzer.find_n_save_orfs(seq)
#     biggest_orf = orf_anlyzer.get_longest_orf(seq)

#     print(f"\n# ===================================")
#     print("-> SEQ INFO")
#     print(f"ID: {seq.id}")
#     print(f"Length: {seq.length} bp")
#     print(f"N Bases: {seq.n_count}")
#     print(f"GC Content: {seq.gc_content:.2f}%")
#     print(f"UUID: {seq.uid}\n")
#     print(f"Seq: {seq.sequence}\n")

#     if biggest_orf is not None:
#         print(f"-> ORF INFO:")
#         print(f"ORF name: {biggest_orf.name}")
#         print(f"ORF type: {biggest_orf.orf_type}")
#         print(f"Length: {biggest_orf.length_aa} AA")
#         print(f"Frame: {biggest_orf.frame} | Strand: {biggest_orf.strand}")
#         print(f"Sequence (AA): {biggest_orf.aa_sequence}")
#         print(f"# ===================================\n")

#     else:
#         print(f"-> ORF INFO")
#         print(f"No valid ORFs :(\n")

# print(f"\n\nTotal sequences processed: {len(Parser.sequences)}\n")



"""
main.py — ViralQuest pipeline entry point.

Demonstrates both sequential and parallel execution of:
  - FASTA parsing
  - ORF finding
  - HMM domain search  (pyhmmer, parallel across HMM databases)
  - Diamond BLASTx     (batched, parallel)
  - BLASTn local       (batched, parallel)
  - BLASTn online      (sequential — NCBI rate-limit policy)
"""

import math
import os
import sys
from pathlib import Path

from loguru import logger

# --- viralquest modules ---
from viralquest.parser  import FastaParser
from viralquest.orfs    import OrfAnalyzer
from viralquest.hmm     import (
    HmmMetadataLoader,
    HmmSequencePreparer,
    HmmSearcher,
    HmmResultAttacher,
)
from viralquest.diamond import (
    DiamondRunner,
    DiamondOutputParser,
    DiamondResultAttacher,
    DiamondFilterer,
)
from viralquest.blastn  import (
    BlastnRunner,
    BlastnOnlineRunner,
    BlastnOutputParser,
    BlastnResultAttacher,
)


# ===========================================================================
# Configuration — edit these paths before running
# ===========================================================================

FASTA_INPUT   = "data/test.fasta"
EMPTY_FASTA   = "data/empty.fasta"

# HMM databases: list of (hmm_file, db_name, json_metadata)
HMM_DATABASES = [
    ("data/hmm-dbs/Pfam-A.hmm",   "Pfam",   "data/hmm-index/Pfam-index.json"),
    ("data/hmm-dbs/U-RVDBv29.0-prot.hmm",   "RVDB",   "data/hmm-index/RVDB-index.json"),
    ("data/hmm-dbs/Vfam-228.hmm",   "Vfam",   "data/hmm-index/Vfam-index.json"),
    ("data/hmm-dbs/EggNOG-4.5.hmm", "EggNOG", "data/hmm-index/eggNOG-4.5-index.json"),
]

DIAMOND_DB    = "/run/media/gabriel/DATA_01/viralDB.dmnd"
BLASTN_DB     = "/run/media/gabriel/DATA_01/viral_blastn"
NCBI_EMAIL    = "gvprodrigues.ppggbm@uesc.br"          # required for online BLAST

OUTPUT_DIR    = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)

CPU_COUNT     = 4 # os.cpu_count() or 4


# ===========================================================================
# Helpers
# ===========================================================================

def _optimal_workers(threads_per_job: int) -> int:
    """How many parallel jobs before we saturate the CPU."""
    return max(1, math.ceil(CPU_COUNT / threads_per_job))


# ===========================================================================
# Pipeline steps
# ===========================================================================

def step_parse(fasta_path: str):
    """Parse the input FASTA and return a list of NucSequence objects."""
    logger.info(f"[1/5] Parsing {fasta_path}")
    parser = FastaParser(fasta_path)
    parser.read_input_file()
    logger.success(f"Parsed {len(parser.sequences)} sequences.")
    return parser.sequences


def step_find_orfs(nuc_seqs):
    """Find ORFs in all sequences (sequential — fast enough in pure Python)."""
    logger.info("[2/5] Finding ORFs")
    analyzer = OrfAnalyzer(min_len_nt=150)
    for seq in nuc_seqs:
        analyzer.find_n_save_orfs(seq)
    total = sum(len(s.orfs) for s in nuc_seqs)
    logger.success(f"Found {total} ORFs across {len(nuc_seqs)} sequences.")


def step_hmm_search(nuc_seqs, parallel: bool = True):
    """
    Run HMM domain search against all configured databases.

    parallel=True  → all HMM files dispatched concurrently (ProcessPoolExecutor)
    parallel=False → sequential, one database at a time
    """
    logger.info(f"[3/5] HMM search (parallel={parallel})")

    preparer = HmmSequencePreparer()
    seq_block, orf_map = preparer.prepare(nuc_seqs)

    if seq_block is None:
        logger.warning("No valid ORF amino-acid sequences — skipping HMM search.")
        return

    # Load all metadata upfront
    all_metadata = {
        hmm_path: HmmMetadataLoader.load(json_path, db_name)
        for hmm_path, db_name, json_path in HMM_DATABASES
    }

    searcher = HmmSearcher(cpus=0, score_threshold=50.0)
    hmm_paths = [p for p, _, _ in HMM_DATABASES]

    if parallel:
        # All databases searched concurrently; each search is CPU-bound
        # so we cap workers to avoid over-subscription.
        all_hits = searcher.search_parallel(
            seq_block,
            hmm_paths,
            max_workers=min(len(hmm_paths), _optimal_workers(threads_per_job=1)),
        )
        for hmm_path, db_name, _ in HMM_DATABASES:
            hits = all_hits.get(hmm_path, [])
            attached = HmmResultAttacher.attach(hits, orf_map, all_metadata[hmm_path], db_name)
            logger.info(f"  {db_name}: {attached} domains attached")

    else:
        for hmm_path, db_name, _ in HMM_DATABASES:
            hits = searcher.search(seq_block, hmm_path)
            attached = HmmResultAttacher.attach(hits, orf_map, all_metadata[hmm_path], db_name)
            logger.info(f"  {db_name}: {attached} domains attached")

    total = sum(len(orf.domains) for seq in nuc_seqs for orf in seq.orfs)
    logger.success(f"HMM search complete — {total} total domains.")


def step_diamond(nuc_seqs, parallel: bool = True):
    """
    Run Diamond BLASTx on all sequences.

    parallel=True  → batches of 1000 sequences dispatched concurrently
    parallel=False → batches run sequentially (useful for debugging)
    """
    logger.info(f"[4/5] Diamond BLASTx (parallel={parallel})")

    runner = DiamondRunner(
        db_path=DIAMOND_DB,
        threads=2,
        e_value=1e-5,
        max_target_seqs=5,
        outdir=str(OUTPUT_DIR / "diamond_tmp"),
        batch_size=1000,
    )

    workers = _optimal_workers(threads_per_job=runner.threads) if parallel else 1
    tsv_paths = runner.run_all(nuc_seqs, max_workers=workers)

    seq_map = {seq.id: seq for seq in nuc_seqs}

    for tsv in tsv_paths:
        hits = DiamondOutputParser.parse(tsv)
        filtered = DiamondFilterer.filter(hits, min_identity=30.0, max_evalue=1e-5)

        for hit in filtered:
            seq = seq_map.get(hit.query_id)
            if seq:
                DiamondResultAttacher.attach([hit], seq)

    total = sum(len(seq.blast_hits) for seq in nuc_seqs)
    logger.success(f"Diamond complete — {total} hits attached.")


def step_blastn_local(nuc_seqs, parallel: bool = True):
    """
    Run local BLASTn on all sequences.

    parallel=True  → batches run concurrently via ThreadPoolExecutor
    parallel=False → one batch at a time
    """
    logger.info(f"[5/5] BLASTn local (parallel={parallel})")

    runner = BlastnRunner(
        database=BLASTN_DB,
        threads=2,
        e_value=1e-5,
        max_target_seqs=1,
        outdir=str(OUTPUT_DIR / "blastn_tmp"),
        batch_size=1000,
    )

    workers = _optimal_workers(threads_per_job=runner.threads) if parallel else 1
    tsv_paths = runner.run_all(nuc_seqs, max_workers=workers)

    all_hits = []
    for tsv in tsv_paths:
        all_hits.extend(BlastnOutputParser.parse(tsv))

    BlastnResultAttacher.attach(all_hits, nuc_seqs)
    logger.success(f"BLASTn local complete — {len(all_hits)} hits.")


def step_blastn_online(nuc_seqs):
    """
    Fallback: query NCBI qblast sequentially.
    Use only when a local BLASTn database is unavailable.
    """
    logger.info("[5/5] BLASTn online (NCBI qblast — sequential)")

    runner = BlastnOnlineRunner(
        database="nt",
        email=NCBI_EMAIL,
        hitlist_size=1,
        sleep_interval=1.0,
        outdir=str(OUTPUT_DIR / "blastn_tmp"),
    )

    tsv_path = runner.run_all(nuc_seqs)
    hits = BlastnOutputParser.parse(tsv_path)
    BlastnResultAttacher.attach(hits, nuc_seqs)
    logger.success(f"BLASTn online complete — {len(hits)} hits.")


# ===========================================================================
# Report helper (minimal — replace with a proper reporter class later)
# ===========================================================================

def print_summary(nuc_seqs):
    logger.info("=" * 60)
    logger.info("Pipeline summary")
    logger.info("=" * 60)
    for seq in nuc_seqs:
        orfs_with_domains = [o for o in seq.orfs if o.domains]
        logger.info(
            f"  {seq.id} | len={seq.length} nt | gc={seq.gc_content:.1f}% | "
            f"orfs={len(seq.orfs)} | orfs_with_domains={len(orfs_with_domains)} | "
            f"blast_hits={len(seq.blast_hits)}"
        )


# ===========================================================================
# Entry points
# ===========================================================================

def run_sequential(fasta_path: str, use_online_blast: bool = False):
    """Full pipeline, all steps sequential."""
    logger.info("Running pipeline — SEQUENTIAL mode")

    nuc_seqs = step_parse(fasta_path)
    if not nuc_seqs:
        logger.warning("No sequences to process.")
        return

    step_find_orfs(nuc_seqs)
    step_hmm_search(nuc_seqs, parallel=False)
    step_diamond(nuc_seqs, parallel=False)

    if use_online_blast:
        step_blastn_online(nuc_seqs)
    else:
        step_blastn_local(nuc_seqs, parallel=False)

    print_summary(nuc_seqs)
    return nuc_seqs


def run_parallel(fasta_path: str, use_online_blast: bool = False):
    """Full pipeline with parallelism enabled at every step that supports it."""
    logger.info("Running pipeline — PARALLEL mode")

    nuc_seqs = step_parse(fasta_path)
    if not nuc_seqs:
        logger.warning("No sequences to process.")
        return

    step_find_orfs(nuc_seqs)            # sequential — fast enough
    step_hmm_search(nuc_seqs, parallel=True)
    step_diamond(nuc_seqs, parallel=True)

    if use_online_blast:
        # NCBI enforces rate limits — always sequential regardless of mode
        step_blastn_online(nuc_seqs)
    else:
        step_blastn_local(nuc_seqs, parallel=True)

    print_summary(nuc_seqs)
    return nuc_seqs


# ===========================================================================
# CLI
# ===========================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="ViralQuest — viral sequence analysis pipeline"
    )
    parser.add_argument(
        "fasta", nargs="?", default=FASTA_INPUT,
        help="Input FASTA file (default: %(default)s)"
    )
    parser.add_argument(
        "--sequential", action="store_true",
        help="Disable parallelism (useful for debugging)"
    )
    parser.add_argument(
        "--online-blast", action="store_true",
        help="Use NCBI qblast instead of a local BLASTn database"
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR"],
        help="Loguru log level (default: INFO)"
    )

    args = parser.parse_args()

    # configure loguru
    logger.remove()
    logger.add(sys.stderr, level=args.log_level)
    logger.add(
        OUTPUT_DIR / "viralquest.log",
        level="DEBUG",
        rotation="10 MB",
        retention=3,
    )

    if args.sequential:
        run_sequential(args.fasta, use_online_blast=args.online_blast)
    else:
        run_parallel(args.fasta, use_online_blast=args.online_blast)