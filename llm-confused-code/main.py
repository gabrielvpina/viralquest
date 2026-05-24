"""
main.py — ViralQuest pipeline  (sequential)

Pipeline order:
  ┌─ FILTER PHASE ──────────────────────────────────────────────────────────┐
  │ 1. Parse FASTA                                                           │
  │ 2. Diamond BLASTx → RefSeq viral proteins  → flag is_viral candidates   │
  │ 3. ORF finding on ALL sequences                                          │
  │ 4. HMM search (RVDB + Vfam + EggNOG) → flag additional is_viral seqs   │
  │    → union of BLASTx hits + HMM hits = confirmed viral candidates        │
  └──────────────────────────────────────────────────────────────────────────┘
  ┌─ CHARACTERISATION PHASE (viral candidates only) ─────────────────────────┐
  │ 5. Diamond BLASTx → NR global db → keep only viral subject titles        │
  │ 6. BLASTn → NT database (local binary) or NCBI qblast (online fallback)  │
  │ 7. HMM search → Pfam only  (functional domain annotation)                │
  │ 8. Build and save JSON report                                             │
  └──────────────────────────────────────────────────────────────────────────┘

Set CPUS at the top. All tools respect it.
"""

import sys
from pathlib import Path

from loguru import logger

from viralquest.parser   import FastaParser
from viralquest.orfs     import OrfAnalyzer
from llm-confused-code.hmm      import (
    HmmMetadataLoader,
    HmmSequencePreparer,
    HmmSearcher,
    HmmResultAttacher,
    HmmViralFlagSetter,
    HmmRole,
    HMM_ROLES,
)
from llm-confused-code.diamond  import (
    DiamondRunner,
    DiamondOutputParser,
    DiamondResultAttacher,
    DiamondFilterer,
    DiamondPhase,
)
from llm-confused-code.blastn   import (
    BlastnRunner,
    BlastnOnlineRunner,
    BlastnOutputParser,
    BlastnResultAttacher,
)
from llm-confused-code.reporter import ViralQuestReporter, ViralTaxonomyDB


# ===========================================================================
# Configuration
# ===========================================================================

CPUS        = 4
SAMPLE_NAME = "MySample"
FASTA_INPUT = "data/test.fasta"

# Diamond databases
DIAMOND_REFSEQ_DB = "/run/media/gabriel/DATA_01/viralDB.dmnd"
DIAMOND_NR_DB     = "/run/media/gabriel/DATA_01/nr/nr.dmnd"

# BLASTn
USE_ONLINE_BLAST = False
BLASTN_DB        = "/run/media/gabriel/DATA_01/viral_blastn/viral_genomic.fna"
NCBI_EMAIL       = "gvprodrigues.ppggbm@uesc.br"

# HMM databases — order matters for logging but not for logic
HMM_DATABASES = [
    ("data/hmm-dbs/U-RVDBv29.0-prot.hmm", "RVDB",   "data/hmm-index/RVDB-index.json"),
    ("data/hmm-dbs/Vfam-228.hmm",          "Vfam",   "data/hmm-index/Vfam-index.json"),
    ("data/hmm-dbs/EggNOG-4.5.hmm",        "EggNOG", "data/hmm-index/eggNOG-4.5-index.json"),
    ("data/hmm-dbs/Pfam-A.hmm",            "Pfam",   "data/hmm-index/Pfam-index.json"),
]

# Quality thresholds
BLASTX_MIN_IDENTITY = 30.0
BLASTX_MAX_EVALUE   = 1e-5
BLASTX_MIN_COVERAGE = 50.0
HMM_SCORE_THRESHOLD = 50.0

OUTPUT_DIR    = Path("results")
OUTPUT_DIR.mkdir(exist_ok=True)

# Taxonomy database (lzma-compressed JSON)
VIRAL_TAX_PATH = Path("data/viralTax.json.xz")


# ===========================================================================
# Shared HMM helper
# ===========================================================================

def _run_hmm_databases(
    nuc_seqs: list,
    role_filter: HmmRole,
    label: str,
) -> None:
    """
    Runs all HMM databases whose HmmRole matches `role_filter`.
    Prepares sequences once, then searches each matching database sequentially.
    """
    target_dbs = [
        (path, name, json_path)
        for path, name, json_path in HMM_DATABASES
        if HMM_ROLES.get(name) == role_filter
    ]
    if not target_dbs:
        return

    preparer = HmmSequencePreparer()
    seq_block, orf_map = preparer.prepare(nuc_seqs)
    if seq_block is None:
        logger.warning(f"HMM [{label}]: no valid ORF sequences — skipping.")
        return

    searcher = HmmSearcher(cpus=CPUS, score_threshold=HMM_SCORE_THRESHOLD)

    for hmm_path, db_name, json_path in target_dbs:
        logger.info(f"  HMM [{label}] searching {db_name} …")
        metadata = HmmMetadataLoader.load(json_path, db_name)
        hits     = searcher.search(seq_block, hmm_path)
        attached = HmmResultAttacher.attach(hits, orf_map, metadata, db_name)
        logger.info(f"  {db_name}: {attached} domains attached")


# ===========================================================================
# ─── FILTER PHASE ───────────────────────────────────────────────────────────
# ===========================================================================

def step_1_parse(fasta_path: str, sample_name: str) -> list:
    logger.info(f"[1/7] Parsing {fasta_path}")
    parser = FastaParser(fasta_path)
    parser.read_input_file()
    for seq in parser.sequences:
        seq.sample_name = sample_name
    logger.success(f"  {len(parser.sequences)} sequences parsed.")
    return parser.sequences


def step_2_blastx_refseq(nuc_seqs: list) -> None:
    """
    BLASTx against RefSeq viral proteins.
    Sets is_viral=True on any sequence with a passing hit.
    The RefSeq db is already viral-only so no subject-title filtering needed.
    """
    logger.info(f"[2/7] BLASTx filter — RefSeq viral proteins")

    runner = DiamondRunner(
        db_path=DIAMOND_REFSEQ_DB,
        threads=CPUS,
        e_value=BLASTX_MAX_EVALUE,
        max_target_seqs=5,
        outdir=str(OUTPUT_DIR / "diamond_refseq_tmp"),
        batch_size=1000,
    )
    tsv_paths = runner.run_all(nuc_seqs, phase=DiamondPhase.REFSEQ_FILTER, max_workers=1)

    all_hits = []
    for tsv in tsv_paths:
        raw = DiamondOutputParser.parse(tsv)
        all_hits.extend(DiamondFilterer.filter(
            raw,
            min_identity=BLASTX_MIN_IDENTITY,
            max_evalue=BLASTX_MAX_EVALUE,
            min_coverage=BLASTX_MIN_COVERAGE,
        ))

    DiamondResultAttacher.attach(all_hits, nuc_seqs, DiamondPhase.REFSEQ_FILTER)
    flagged = sum(1 for s in nuc_seqs if s.is_viral)
    logger.success(f"  BLASTx refseq: {flagged} viral candidates flagged.")


def step_3_find_orfs(nuc_seqs: list) -> None:
    """ORF finding on ALL sequences (needed for HMM filter step)."""
    logger.info(f"[3/7] ORF finding — {len(nuc_seqs)} sequences")
    analyzer = OrfAnalyzer(min_len_nt=150)
    for seq in nuc_seqs:
        analyzer.find_n_save_orfs(seq)
    total = sum(len(s.orfs) for s in nuc_seqs)
    logger.success(f"  {total} ORFs found.")


def step_4_hmm_filter(nuc_seqs: list) -> None:
    """
    HMM filter using RVDB, Vfam, EggNOG.
    Sequences not already flagged by BLASTx but with a hit here are also
    marked is_viral=True (union logic).
    """
    logger.info(f"[4/7] HMM filter — RVDB / Vfam / EggNOG")
    _run_hmm_databases(nuc_seqs, HmmRole.FILTER, "filter")
    newly = HmmViralFlagSetter.flag(nuc_seqs)
    total = sum(1 for s in nuc_seqs if s.is_viral)
    logger.success(
        f"  HMM filter done. {newly} additional sequences flagged. "
        f"Total viral candidates: {total}."
    )


# ===========================================================================
# ─── CHARACTERISATION PHASE ─────────────────────────────────────────────────
# ===========================================================================

def step_5_blastx_nr(viral_seqs: list) -> None:
    """
    BLASTx against NR (global, all organisms).
    Only hits whose subject title matches viral keywords are kept
    (handled inside DiamondResultAttacher for NR_CHARACTERIZE phase).
    """
    logger.info(f"[5/7] BLASTx characterisation — NR global db ({len(viral_seqs)} seqs)")

    runner = DiamondRunner(
        db_path=DIAMOND_NR_DB,
        threads=CPUS,
        e_value=BLASTX_MAX_EVALUE,
        max_target_seqs=5,        # more targets — we filter by title afterwards
        outdir=str(OUTPUT_DIR / "diamond_nr_tmp"),
        batch_size=500,
    )
    tsv_paths = runner.run_all(viral_seqs, phase=DiamondPhase.NR_CHARACTERIZE, max_workers=1)

    all_hits = []
    for tsv in tsv_paths:
        raw = DiamondOutputParser.parse(tsv)
        all_hits.extend(DiamondFilterer.filter(
            raw,
            min_identity=BLASTX_MIN_IDENTITY,
            max_evalue=BLASTX_MAX_EVALUE,
            min_coverage=BLASTX_MIN_COVERAGE,
        ))

    # attacher applies is_viral_subject() internally for NR phase
    DiamondResultAttacher.attach(all_hits, viral_seqs, DiamondPhase.NR_CHARACTERIZE)
    attached = sum(len(s.blastx_nr_hits) for s in viral_seqs)
    logger.success(f"  BLASTx NR done — {attached} viral NR hits attached.")


def step_6_blastn(viral_seqs: list) -> None:
    """BLASTn on viral candidates — local binary or NCBI qblast."""
    if USE_ONLINE_BLAST:
        _blastn_online(viral_seqs)
    else:
        _blastn_local(viral_seqs)


def _blastn_local(viral_seqs: list) -> None:
    logger.info(f"[6/7] BLASTn local — {len(viral_seqs)} sequences")
    runner = BlastnRunner(
        database=BLASTN_DB,
        threads=CPUS,
        e_value=1e-5,
        max_target_seqs=1,
        outdir=str(OUTPUT_DIR / "blastn_tmp"),
        batch_size=500,
    )
    tsv_paths = runner.run_all(viral_seqs, max_workers=1)
    all_hits  = []
    for tsv in tsv_paths:
        all_hits.extend(BlastnOutputParser.parse(tsv))
    BlastnResultAttacher.attach(all_hits, viral_seqs)
    logger.success(f"  BLASTn local done — {len(all_hits)} hits.")


def _blastn_online(viral_seqs: list) -> None:
    logger.info(f"[6/7] BLASTn online (NCBI qblast) — {len(viral_seqs)} sequences")
    runner = BlastnOnlineRunner(
        database="nt",
        email=NCBI_EMAIL,
        hitlist_size=1,
        sleep_interval=1.5,
        outdir=str(OUTPUT_DIR / "blastn_tmp"),
    )
    tsv_path = runner.run_all(viral_seqs)
    hits = BlastnOutputParser.parse(tsv_path)
    BlastnResultAttacher.attach(hits, viral_seqs)
    logger.success(f"  BLASTn online done — {len(hits)} hits.")


def step_7_pfam(viral_seqs: list) -> None:
    """Pfam domain annotation on confirmed viral sequences."""
    logger.info(f"[7/7] Pfam annotation — {len(viral_seqs)} sequences")
    _run_hmm_databases(viral_seqs, HmmRole.CHARACTERIZE, "Pfam")
    total = sum(
        len([d for d in orf.domains if d.database == "Pfam"])
        for seq in viral_seqs
        for orf in seq.orfs
    )
    logger.success(f"  Pfam done — {total} domains annotated.")


# ===========================================================================
# Full pipeline
# ===========================================================================

def run_pipeline(fasta_path: str = FASTA_INPUT, sample_name: str = SAMPLE_NAME):
    logger.info(f"ViralQuest — sample: {sample_name} | CPUs: {CPUS}")

    # ── Filter phase ──────────────────────────────────────────────────────
    nuc_seqs = step_1_parse(fasta_path, sample_name)
    if not nuc_seqs:
        logger.warning("No sequences found. Exiting.")
        return

    step_2_blastx_refseq(nuc_seqs)
    step_3_find_orfs(nuc_seqs)
    step_4_hmm_filter(nuc_seqs)

    viral_seqs = [s for s in nuc_seqs if s.is_viral]
    if not viral_seqs:
        logger.warning("No viral candidates after filter phase. Exiting.")
        return

    logger.info(
        f"Filter phase complete — {len(viral_seqs)}/{len(nuc_seqs)} "
        f"sequences confirmed as viral candidates."
    )

    # ── Characterisation phase ────────────────────────────────────────────
    step_5_blastx_nr(viral_seqs)
    step_6_blastn(viral_seqs)
    step_7_pfam(viral_seqs)

    # ── Report ────────────────────────────────────────────────────────────
    logger.info("Loading taxonomy database …")
    tax_db   = ViralTaxonomyDB(VIRAL_TAX_PATH)
    reporter = ViralQuestReporter(sample_name=sample_name, tax_db=tax_db)
    out_json  = OUTPUT_DIR / f"{sample_name}_viralquest.json"
    reporter.save(nuc_seqs, out_json)   # reporter filters is_viral internally

    logger.success(f"Pipeline complete → {out_json}")
    return nuc_seqs


# ===========================================================================
# CLI
# ===========================================================================

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="ViralQuest — viral sequence analysis pipeline")
    ap.add_argument("fasta",    nargs="?",  default=FASTA_INPUT)
    ap.add_argument("--sample", default=SAMPLE_NAME)
    ap.add_argument("--online-blast", action="store_true")
    ap.add_argument(
        "--log-level", default="INFO",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR"],
    )
    args = ap.parse_args()

    if args.online_blast:
        USE_ONLINE_BLAST = True

    logger.remove()
    logger.add(sys.stderr, level=args.log_level)
    logger.add(OUTPUT_DIR / "viralquest.log", level="DEBUG", rotation="10 MB", retention=3)

    run_pipeline(fasta_path=args.fasta, sample_name=args.sample)