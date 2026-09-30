#!/usr/bin/env python3
"""
cli.py — ViralQuest command-line interface.

Auto-detected databases (no flags needed):
  data/viralDB.dmnd          — RefSeq viral Diamond filter DB
  data/hmm-dbs/*.hmm         — RVDB, Vfam, EggNOG, Pfam profiles

Run  python download_dbs.py  first if databases are missing.

Output modes
------------
Default  : loguru writes directly to stderr — full timestamped log stream.
--live   : ASCII banner + scrolling log box + step progress bar via Rich Live.
"""

from __future__ import annotations

import platform
import shutil
import sys
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path

# ── Version ───────────────────────────────────────────────────────────────────

__version__ = "3.0.2"
__author__  = "Gabriel Rodrigues"
__link__ = "https://github.com/gabrielvpina/viralquest"

# ── Database path resolution ──────────────────────────────────────────────────

_DEFAULT_DATA = Path(__file__).parent.parent / "data"


def _resolve_db_paths(ext_db_dir: Path | None = None) -> tuple[dict, Path, list]:
    """Return (db, index_dir, hmm_filter).

    ext_db_dir — flat directory supplied via --db-dir that contains the five
    binary database files (*.hmm + viralDB.dmnd) with no subdirectories.
    When None the default nested data/ layout is used.

    Index files, viralTax.json.xz, and viral-family-info/ are always read
    from the bundled data/ directory regardless of ext_db_dir.
    """
    data      = _DEFAULT_DATA
    index_dir = data / "hmm-index"
    fam_dir   = data / "viral-family-info"

    if ext_db_dir is not None:
        hmm_base  = ext_db_dir
        dmnd_base = ext_db_dir
    else:
        hmm_base  = data / "hmm-dbs"
        dmnd_base = data / "viral-db"

    db = {
        "viral_dmnd": dmnd_base / "viralDB.dmnd",
        "rvdb":       hmm_base  / "U-RVDBv29.0-prot.hmm",
        "vfam":       hmm_base  / "Vfam-228.hmm",
        "eggnog":     hmm_base  / "EggNOG-4.5.hmm",
        "pfam":       hmm_base  / "Pfam-A.hmm",
        "viral_tax":  data      / "viralTax.json.xz",
        "fam_high":   fam_dir   / "viral_info_highToken.json",
        "fam_low":    fam_dir   / "viral_info_lowToken.json",
    }

    hmm_filter = [
        ("rvdb",   index_dir / "RVDB-index.json",       "RVDB"),
        ("vfam",   index_dir / "Vfam-index.json",       "Vfam"),
        ("eggnog", index_dir / "eggNOG-4.5-index.json", "EggNOG"),
    ]

    return db, index_dir, hmm_filter

# ── Banner ────────────────────────────────────────────────────────────────────

_BANNER_RAW = r"""
██╗   ██╗██╗██████╗  █████╗ ██╗      ██████╗ ██╗   ██╗███████╗███████╗████████╗
██║   ██║██║██╔══██╗██╔══██╗██║     ██╔═══██╗██║   ██║██╔════╝██╔════╝╚══██╔══╝
██║   ██║██║██████╔╝███████║██║     ██║   ██║██║   ██║█████╗  ███████╗   ██║
╚██╗ ██╔╝██║██╔══██╗██╔══██║██║     ██║▄▄ ██║██║   ██║██╔══╝  ╚════██║   ██║
 ╚████╔╝ ██║██║  ██║██║  ██║███████╗╚██████╔╝╚██████╔╝███████╗███████║   ██║
  ╚═══╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚══▀▀═╝  ╚═════╝ ╚══════╝╚══════╝   ╚═╝
"""

# ── Argument parser ───────────────────────────────────────────────────────────

def _build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        prog="viralquest",
        description=_BANNER_RAW,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-h", "--help",
        action="store_true", default=False,
        help="Show this help message and exit.")

    # Required ─────────────────────────────────────────────────────────────────
    req = parser.add_argument_group("required")
    req.add_argument("-in", "--input", dest="input", type=str,
        metavar="INPUT.fasta",
        help="Input nucleotide FASTA file.")
    req.add_argument("-out", "--outdir", dest="outdir", type=str,
        metavar="OUTPUT_DIR",
        help="Directory for all output files (created if absent).")

    # Assembly ─────────────────────────────────────────────────────────────────
    asm = parser.add_argument_group("assembly")
    asm.add_argument("--cap3", action="store_true",
        help="Run CAP3 before analysis to assemble overlapping reads into contigs.")

    # Databases (optional overrides) ───────────────────────────────────────────
    dbs = parser.add_argument_group("databases (optional overrides)")
    dbs.add_argument("-nr", "--nr-db", dest="nr_db", type=str,
        metavar="NR.dmnd",
        help="Path to a Diamond-format NCBI NR database for protein-level "
             "characterisation of confirmed viral sequences. "
             "This is a large (>100 GB) database and is optional.")
    dbs.add_argument("--nr-block-size", dest="nr_block_size", type=float,
        default=None, metavar="N",
        help="Diamond --block-size for NR search: GB of RAM per thread pass "
             "(default: Diamond's built-in default ~2.0). "
             "Increase to reduce database passes and speed up the search "
             "(e.g. 6 uses ~45 GB, 12 uses ~90 GB).")
    dbs.add_argument("--nr-index-chunks", dest="nr_index_chunks", type=int,
        default=None, metavar="N",
        help="Diamond --index-chunks for NR search: number of seed-index chunks "
             "(default: 4). Lower values (e.g. 1 or 2) mean fewer disk passes "
             "but higher RAM usage.")
    dbs.add_argument("--nr-tmpdir", dest="nr_tmpdir", type=str,
        default=None, metavar="DIR",
        help="Directory for Diamond temporary files during the NR search. "
             "Pointing this at a fast NVMe drive or RAM disk (/dev/shm) "
             "removes disk I/O as a bottleneck.")

    # BLASTn ───────────────────────────────────────────────────────────────────
    bln = parser.add_argument_group("blastn (choose one)")
    bln_ex = bln.add_mutually_exclusive_group()
    bln_ex.add_argument("-n", "--blastn-local", dest="blastn_local", type=str,
        metavar="BLAST_DB",
        help="Path to a local BLASTn nucleotide database.")
    bln_ex.add_argument("--blastn-online", dest="blastn_online", type=str,
        metavar="EMAIL",
        help="NCBI e-mail for web BLASTn (slower, no local DB required).")
    bln.add_argument("--blastn-online-db", dest="blastn_online_db", type=str,
        default="nt", metavar="DB",
        help="NCBI nucleotide database for web BLASTn (default: nt).")

    # Pipeline tuning ──────────────────────────────────────────────────────────
    tun = parser.add_argument_group("pipeline tuning")
    tun.add_argument("-cpu", "--cpu", dest="cpu", type=int,
        default=2, metavar="N",
        help="CPU threads for Diamond, HMMsearch, BLASTn, and Salmon (default: 2).")
    tun.add_argument("--force", action="store_true",
        help="Export all input sequences, not only confirmed viral ones.")
    tun.add_argument("--min-identity", dest="min_identity", type=float,
        default=90.0, metavar="PCT",
        help="Minimum %% identity for intra-cluster BLASTn alignment (default: 90.0).")
    tun.add_argument("--min-coverage", dest="min_coverage", type=float,
        default=50.0, metavar="PCT",
        help="Minimum query coverage (%%) for a sequence to qualify as a cluster member "
             "(default: 50.0). Clusters with fewer than 2 qualifying members are dropped.")

    # Salmon quantification ────────────────────────────────────────────────────
    sal = parser.add_argument_group("salmon quantification (optional)")
    sal.add_argument("--transcriptome", dest="transcriptome", type=str,
        metavar="HOST.fasta",
        help="Host transcriptome FASTA for Salmon quantification.")
    sal.add_argument("--reads", dest="reads", nargs="+", type=str,
        metavar="READS.fastq",
        help="FASTQ file(s) for Salmon: one = single-end, two = paired-end.")
    sal.add_argument("--read-type", dest="read_type", type=str,
        choices=["sr", "ont", "pb", "hifi"],
        metavar="TYPE",
        help="Sequencing technology of the reads supplied to --reads. "
             "Required when --reads is used. "
             "Choices: sr (Illumina short reads), ont (Oxford Nanopore), "
             "pb (PacBio CLR), hifi (PacBio HiFi/CCS). "
             "NOTE: Salmon quantification runs for short reads (sr) only — its "
             "short-read mapping is invalid for long reads, so for ont/pb/hifi the "
             "Salmon step is skipped and only minimap2 read coverage is produced.")
    sal.add_argument("--hk-genes", dest="hk_genes", type=str,
        default=None, metavar="IDS.txt",
        help="Text file with reference housekeeping gene IDs (one per line) for "
             "normalization. Matched against the --transcriptome headers by "
             "accession (with or without the .version suffix, case-insensitive) "
             "or by gene symbol found in the header. Requires --transcriptome.")
    sal.add_argument("--low-memory", dest="low_memory", action="store_true",
        help="Low-RAM mode for Salmon: builds the index with a smaller k-mer size "
             "(-k 21 instead of 31) and, in de-novo mode, excludes background "
             "contigs shorter than 500 bp from the index (viral and HK-matched "
             "contigs are always included). Both changes directly reduce the SSHash "
             "index footprint. Recommended when salmon index runs out of memory.")
    sal.add_argument("--skip-salmon", dest="skip_salmon", action="store_true",
        help="Skip Salmon quantification but still process the reads: read "
             "coverage (incl. sense/antisense) and sequence-quality signals run "
             "as usual.")

    # Sequence quality ──────────────────────────────────────────────────────────
    # Always runs (reads are not required): dustmasker, self-BLASTn and jellyfish
    # flag structural issues (low complexity, repeats) on the assembled sequences.
    rq = parser.add_argument_group("sequence quality (always on)")
    rq.add_argument("--kmer", dest="kmer", type=int, default=15, metavar="K",
        help="k-mer size for the jellyfish repetitiveness score (default: 15).")
    rq.add_argument("--skip-seq-quality", dest="skip_seq_quality", action="store_true",
        help="Skip the sequence-quality module (dustmasker low-complexity, "
             "self-BLASTn repeats, jellyfish k-mer score and dot plot). It runs "
             "by default on the confirmed viral sequences, with or without --reads.")

    # AI scoring ───────────────────────────────────────────────────────────────
    ai = parser.add_argument_group("AI scoring (optional)")
    ai.add_argument("--model-type", dest="model_type",
        choices=["ollama", "openai", "anthropic", "google"],
        help="AI provider for LLM viral sequence scoring.")
    ai.add_argument("--model-name", dest="model_name", type=str,
        metavar="MODEL",
        help="Model name (e.g. 'qwen3:4b', 'gpt-4o', 'claude-opus-4-7').")
    ai.add_argument("--llm-tokens", dest="llm_tokens",
        choices=["high", "low"],
        help="Token usage mode: 'high' (all hits, full details) or "
             "'low' (best hits only, compact). Required when --model-type is set.")
    ai.add_argument("--api-key", dest="api_key", type=str,
        metavar="KEY",
        help="API key for cloud AI providers (not needed for ollama).")

    # Output mode / misc ──────────────────────────────────────────────────────
    parser.add_argument("--live", action="store_true",
        help="Rich Live display: banner + scrolling log box + step progress bar.")
    parser.add_argument("--db-dir", dest="db_dir", type=str,
        default=None, metavar="DIR",
        help="Flat directory containing the five binary database files "
             "(U-RVDBv29.0-prot.hmm, Vfam-228.hmm, EggNOG-4.5.hmm, Pfam-A.hmm, "
             "viralDB.dmnd). No subdirectories needed. Indices and metadata are "
             "always read from the bundled data/ folder.")

    # Version ──────────────────────────────────────────────────────────────────
    parser.add_argument("-v", "--version",
        action="version", version=f"ViralQuest v{__version__}")

    return parser


# ── Rich help ─────────────────────────────────────────────────────────────────

def _show_rich_help() -> None:
    from rich.console import Console
    from rich.panel import Panel
    from rich.text import Text
    from rich import box

    console = Console()
    console.print(Text(_BANNER_RAW, style="bold cyan"))

    tagline = Text()
    tagline.append("  Viral diversity discovery and characterisation pipeline\n", style="italic")
    tagline.append(f"  v{__version__} · {__author__}\n", style="italic")
    console.print(tagline)

    console.print(Panel(
        "[bold cyan]-in / --input[/]    [dim]FASTA[/]\n"
        "  Input nucleotide FASTA file.\n\n"
        "[bold cyan]-out / --outdir[/]  [dim]DIR[/]\n"
        "  Output directory (created if absent).",
        title="[bold red]REQUIRED[/bold red]",
        border_style="red", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]--cap3[/]\n"
        "  Assemble overlapping sequences with CAP3 before analysis.\n\n"
        "[bold cyan]-cpu / --cpu[/]       [dim]N[/]  (default: 2)\n"
        "  CPU threads used by Diamond, HMMsearch, BLASTn, and Salmon.\n\n"
        "[bold cyan]--min-identity[/]    [dim]%[/]  (default: 90.0)\n"
        "  Minimum %% identity for intra-cluster BLASTn alignment.\n\n"
        "[bold cyan]--min-coverage[/]    [dim]%[/]  (default: 50.0)\n"
        "  Minimum query coverage for a sequence to qualify as a cluster member.\n"
        "  Clusters with fewer than 2 qualifying members are dropped.\n\n"
        "[bold cyan]--force[/]\n"
        "  Export all input sequences even if viral confirmation fails.",
        title="[bold green]PIPELINE OPTIONS[/bold green]",
        border_style="green", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]-nr / --nr-db[/]         [dim]NR.dmnd[/]\n"
        "  Diamond-format NCBI NR database for full protein-level characterisation\n"
        "  of confirmed viral sequences. Optional but recommended for complete\n"
        "  annotation. NR is a large database (>100 GB); build with:\n"
        "  [dim]diamond makedb --in nr.fasta --db nr[/dim]\n\n"
        "[bold cyan]--nr-block-size[/]       [dim]N[/]  (default: Diamond built-in ~2.0)\n"
        "  GB of RAM loaded per database pass. Higher = fewer passes = faster.\n"
        "  Examples: 6 → ~45 GB RAM · 12 → ~90 GB RAM.\n\n"
        "[bold cyan]--nr-index-chunks[/]     [dim]N[/]  (default: 4)\n"
        "  Seed-index chunks. Lower values (1–2) reduce disk passes but use\n"
        "  more RAM. Combine with --nr-block-size for maximum speed.\n\n"
        "[bold cyan]--nr-tmpdir[/]           [dim]DIR[/]\n"
        "  Directory for Diamond temporary files. Use a fast NVMe drive or\n"
        "  RAM disk ([dim]/dev/shm[/dim]) to eliminate I/O as a bottleneck.",
        title="[bold yellow]NR DATABASE (optional)[/bold yellow]",
        border_style="yellow", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]-n / --blastn-local[/]   [dim]BLAST_DB[/]\n"
        "  Path to a local BLAST nucleotide database (e.g. nt).\n\n"
        "[bold cyan]--blastn-online[/]       [dim]EMAIL[/]\n"
        "  NCBI e-mail for web BLASTn — no local database required.\n\n"
        "[bold cyan]--blastn-online-db[/]    [dim]DB[/]  (default: nt)\n"
        "  NCBI database to query when using --blastn-online.\n\n"
        "[bold]Note:[/] --blastn-local and --blastn-online are mutually exclusive.",
        title="[bold yellow]BLASTN (choose one)[/bold yellow]",
        border_style="yellow", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]--reads[/]          [dim]R1.fastq [R2.fastq][/]\n"
        "  FASTQ file(s): one = single-end, two = paired-end.\n\n"
        "[bold cyan]--read-type[/]      [dim]sr | ont | pb | hifi[/]\n"
        "  Sequencing technology. Required when --reads is used.\n"
        "  sr = Illumina short reads, ont = Oxford Nanopore,\n"
        "  pb = PacBio CLR, hifi = PacBio HiFi/CCS.\n"
        "  [yellow]Salmon quantification runs for [bold]sr[/bold] only[/yellow] — its short-read\n"
        "  mapping is invalid for long reads. With ont/pb/hifi the Salmon\n"
        "  step is skipped; only minimap2 read coverage is produced.\n\n"
        "[bold cyan]--transcriptome[/]  [dim]HOST.fasta[/]\n"
        "  Host transcriptome FASTA. When provided, runs the [bold]reference pathway[/bold]:\n"
        "  viral + bundled HK + ref-HK + transcriptome as a combined Salmon index.\n"
        "  Without --transcriptome, runs the [bold]de-novo pathway[/bold]: quantifies\n"
        "  assembled contigs directly, using BLASTn to find HK-matching contigs.\n\n"
        "[bold cyan]--hk-genes[/]       [dim]IDS.txt[/]\n"
        "  Text file with reference HK gene IDs (one per line) for normalization.\n"
        "  Matched against the --transcriptome headers by accession (with or\n"
        "  without the [dim].version[/dim] suffix, case-insensitive) or by gene symbol\n"
        "  present in the header. Requires --transcriptome.\n\n"
        "[bold cyan]--low-memory[/]\n"
        "  Low-RAM mode for [dim]salmon index[/dim]: uses [dim]-k 21[/dim] (instead of 31) and,\n"
        "  in de-novo mode, excludes background contigs shorter than 500 bp\n"
        "  from the index. Viral and HK-matched contigs are always included\n"
        "  regardless of length. Both changes reduce the SSHash index footprint.\n"
        "  Recommended when [dim]salmon index[/dim] runs out of memory.\n\n"
        "[bold cyan]--skip-salmon[/]\n"
        "  Skip Salmon quantification but still process the reads: read coverage\n"
        "  (incl. sense/antisense strands) and sequence-quality signals still run.\n\n"
        "[bold]Note:[/] --reads alone triggers de-novo mode. --transcriptome requires --reads.\n"
        "[bold]Note:[/] For metagenomics or metatranscriptomics samples, omit --transcriptome — "
        "the de-novo pathway is the appropriate mode and no host reference is needed.",
        title="[bold cyan]SALMON QUANTIFICATION (optional)[/bold cyan]",
        border_style="cyan", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "Structural checks on the confirmed viral sequences: [dim]dustmasker[/dim]\n"
        "low-complexity regions, [dim]self-BLASTn[/dim] direct/inverted repeats (with a\n"
        "dot plot) and a [dim]jellyfish[/dim] k-mer repetitiveness score.\n"
        "Runs by default in every mode — [bold]reads are not required[/bold].\n\n"
        "[bold cyan]--kmer[/]              [dim]K[/]\n"
        "  k-mer size for the jellyfish repetitiveness score (default: 15).\n\n"
        "[bold cyan]--skip-seq-quality[/]\n"
        "  Skip this module entirely (saves the self-BLASTn / dot-plot time on\n"
        "  runs with many or very long contigs).",
        title="[bold cyan]SEQUENCE QUALITY[/bold cyan]",
        border_style="cyan", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]--model-type[/]   [dim]ollama | openai | anthropic | google[/]\n"
        "  AI provider for LLM viral sequence scoring.\n\n"
        "[bold cyan]--model-name[/]   [dim]MODEL[/]\n"
        "  Model identifier (e.g. 'qwen3:4b', 'gpt-4o', 'claude-opus-4-7',\n"
        "  'gemini-2.5-pro').\n\n"
        "[bold cyan]--llm-tokens[/]   [dim]high | low[/]  [bold red](required with --model-type)[/]\n"
        "  Token usage mode for LLM prompts:\n"
        "    [bold]high[/] — all hits, Pfam details, full taxonomy, full family description.\n"
        "    [bold]low[/]  — best hits only, no Pfam details, compact taxonomy.\n\n"
        "[bold cyan]--api-key[/]      [dim]KEY[/]\n"
        "  API key for cloud providers (not required for ollama).",
        title="[bold magenta]AI SCORING (optional)[/bold magenta]",
        border_style="magenta", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold cyan]--db-dir[/]  [dim]DIR[/]\n"
        "  Flat directory that contains the five binary database files:\n"
        "  [dim]U-RVDBv29.0-prot.hmm  Vfam-228.hmm  EggNOG-4.5.hmm\n"
        "  Pfam-A.hmm  viralDB.dmnd[/dim]\n"
        "  Files must sit directly in DIR — no subdirectories.\n"
        "  Useful when databases are stored outside the package.\n"
        "  Download them with: [dim]viralquest-download --db-dir DIR[/dim]\n"
        "  Indices and metadata are always read from the bundled data/ folder.\n\n"
        "[bold cyan]--live[/]\n"
        "  Rich Live display: ASCII banner + scrolling log box + step progress bar.\n"
        "  Default (without flag): raw loguru log stream to stderr.\n\n"
        "[bold cyan]-h / --help[/]     Show this help message and exit.\n"
        "[bold cyan]-v / --version[/]  Show version and exit.",
        title="[bold cyan]OTHER[/bold cyan]",
        border_style="blue", width=85, box=box.ROUNDED,
    ))


# ── Validation ────────────────────────────────────────────────────────────────

def _check_databases(console, db: dict) -> bool:
    missing = [name for name, path in db.items() if not path.exists()]
    if not missing:
        return True
    console.print(f"\n[bold red]ERROR:[/bold red] {len(missing)} database file(s) not found:\n")
    for name in missing:
        console.print(f"  [red]✗[/red]  {name}  [dim]{db[name]}[/dim]")
    console.print(
        "\n[bold yellow]Fix (option 1):[/bold yellow]  "
        "[bold cyan]viralquest-download[/bold cyan]"
        "  — downloads into the package data/ folder.\n"
        "[bold yellow]Fix (option 2):[/bold yellow]  "
        "[bold cyan]viralquest --db-dir DIR ...[/bold cyan]"
        "  — point to a flat directory that already contains the database files.\n"
    )
    return False


def _validate_args(args, console) -> None:
    errors = []
    if not args.input:
        errors.append("--input is required.")
    elif not Path(args.input).exists():
        errors.append(f"input file not found: {args.input}")
    if not args.outdir:
        errors.append("--outdir is required.")
    if args.transcriptome and not args.reads:
        errors.append("--transcriptome requires --reads.")
    if args.hk_genes and not args.transcriptome:
        errors.append("--hk-genes requires --transcriptome.")
    if args.reads and len(args.reads) > 2:
        errors.append("--reads accepts at most two files (R1 and R2).")
    if args.reads and not args.read_type:
        errors.append("--read-type is required when --reads is used (choices: sr, ont, pb, hifi).")
    if args.read_type and not args.reads:
        errors.append("--read-type requires --reads.")
    if args.kmer is not None and args.kmer < 2:
        errors.append("--kmer must be >= 2.")
    if args.model_type and not args.model_name:
        errors.append("--model-name is required when --model-type is set.")
    if args.model_type and not args.llm_tokens:
        errors.append("--llm-tokens (high|low) is required when --model-type is set.")
    if args.llm_tokens and not args.model_type:
        errors.append("--llm-tokens requires --model-type to be set.")
    if args.nr_db and not Path(args.nr_db).exists():
        errors.append(f"NR database not found: {args.nr_db}")
    if args.cap3 and platform.system() == "Darwin" and not shutil.which("cap3"):
        errors.append(
            "--cap3 is unavailable: CAP3 has no native macOS build on bioconda "
            "(neither Apple Silicon nor Intel), so it isn't installed by the macOS "
            "pixi environment. Run without --cap3, "
            "or install a cap3 binary on PATH yourself."
        )
    for msg in errors:
        console.print(f"[bold red]ERROR:[/bold red] {msg}")
    if errors:
        sys.exit(1)
    if args.model_type and not args.nr_db:
        console.print(
            "[bold yellow]WARNING:[/bold yellow] --model-type is set but --nr-db was not provided. "
            "LLM scoring requires NR-confirmed sequences and will be skipped."
        )


# ── Live-mode display helpers ─────────────────────────────────────────────────

_BANNER_GRADIENT = ("#8be9fd", "#6fd3f7", "#4fc3f7", "#29b6f6", "#0ea5e9", "#0284c7")

_LEVEL_STYLE = {
    "TRACE":    "dim",
    "DEBUG":    "dim",
    "INFO":     "cyan",
    "SUCCESS":  "green",
    "WARNING":  "yellow",
    "ERROR":    "bold red",
    "CRITICAL": "bold white on red",
}


def _fmt_clock(sec: float) -> str:
    """Elapsed time as m:ss (or h:mm:ss) for the live screen."""
    sec = int(max(sec, 0))
    h, rem = divmod(sec, 3600)
    m, s   = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


def _fmt_step_time(sec: float) -> str:
    if sec < 60:
        return f"{sec:.1f}s"
    return _fmt_clock(sec)


class _LiveState:
    """Mutable pipeline state read by the live screen on every refresh."""

    def __init__(self, steps: list[str]):
        self.steps       = steps
        self.current     = 0                 # index of the running step
        self.timings: dict[int, float] = {}
        self.run_start   = time.time()
        self.step_start  = time.time()
        self.batch: str | None = None        # e.g. "3/12" during Diamond batches
        self.failed: int | None = None       # index of the step that raised
        self.done        = False
        self.end: float | None = None        # set when finished/failed → clock stops
        self.warnings    = 0
        self.errors      = 0

    def finish_step(self, idx: int, elapsed: float) -> None:
        self.timings[idx] = elapsed
        self.current      = idx + 1
        self.step_start   = time.time()
        self.batch        = None
        self.done         = self.current >= len(self.steps)
        if self.done:
            self.end = time.time()

    def fail(self) -> None:
        self.failed = self.current
        self.end    = time.time()

    def elapsed(self) -> float:
        return (self.end or time.time()) - self.run_start


class _LogTail:
    """Last log records that fit the panel height, coloured by level.

    Records are kept as plain data and rendered as Text (never parsed as
    markup), so messages like "diamond [refseq]" display verbatim.
    """

    def __init__(self, records: deque):
        self.records = records

    def __rich_console__(self, console, options):
        from rich.text import Text

        height = options.height or 20
        rows   = list(self.records)[-height:]
        if not rows:
            yield Text("Waiting for output…", style="dim italic")
            return
        for ts, level, message in rows:
            line = Text(no_wrap=True, overflow="ellipsis")
            line.append(f"{ts} ", style="dim")
            line.append(f"{level:<8} ", style=_LEVEL_STYLE.get(level, ""))
            line.append(message.splitlines()[0] if message else "",
                        style="red" if level in ("ERROR", "CRITICAL") else "")
            yield line


class _LiveScreen:
    """Full live layout: header · step checklist | log · progress footer.

    Rebuilt on every Live refresh, so elapsed timers, the spinner and the
    log tail keep moving during long steps (not only when a step finishes).
    """

    def __init__(self, state: _LiveState, records: deque, info: dict):
        from rich.spinner import Spinner

        self.state   = state
        self.records = records
        self.info    = info
        self.spinner = Spinner("dots", style="bold cyan")   # reused → animates

    # ── pieces ────────────────────────────────────────────────────────────
    def _banner(self):
        from rich.text import Text

        text = Text(no_wrap=True, overflow="crop")
        lines = _BANNER_RAW.strip("\n").splitlines()
        for i, line in enumerate(lines):
            text.append(line + ("\n" if i < len(lines) - 1 else ""),
                        style=f"bold {_BANNER_GRADIENT[i % len(_BANNER_GRADIENT)]}")
        return text

    def _info(self):
        from rich.table import Table

        o    = self.info
        grid = Table.grid(padding=(0, 1))
        grid.add_column(style="dim", justify="right", no_wrap=True, min_width=7)
        grid.add_column(no_wrap=True, overflow="ellipsis")
        grid.add_row("version", f"[bold]v{__version__}[/]")
        grid.add_row("input",   f"[bold]{o['input']}[/]")
        grid.add_row("output",  str(o["outdir"]))
        grid.add_row("threads", str(o["threads"]))

        def mod(name, on, value=""):
            return (f"[green]●[/] {name} [dim]{value}[/]".rstrip() if on
                    else f"[dim]○ {name}[/]")
        reads = f"{len(o['reads'])} files · {o['read_type']}" if o["reads"] else ""
        grid.add_row("modules", "  ".join([
            mod("CAP3",   o["cap3"]),
            mod("NR",     o["nr_db"],  o["nr_db"] or ""),
            mod("BLASTn", o["blastn"], o["blastn"] or ""),
        ]))
        grid.add_row("", "  ".join([
            mod("reads", o["reads"], reads),
            mod("LLM",   o["llm"],   o["llm"] or ""),
        ]))
        return grid

    def _header(self, width: int):
        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        if width < 84:   # banner would not fit — compact title line
            title = Text.assemble(("ViralQuest", "bold #29b6f6"),
                                  (f"  v{__version__}  ·  {self.info['input']}", "dim"))
            return Panel(title, border_style="#0284c7", padding=(0, 1)), 3

        if width >= 128:
            row = Table.grid(padding=(0, 3), expand=True)
            row.add_column(no_wrap=True)
            row.add_column(ratio=1, vertical="middle")
            row.add_row(self._banner(), self._info())
            body = row
        else:
            body = self._banner()
        return Panel(body, border_style="#0284c7", padding=(0, 2)), 8

    def _steps_panel(self):
        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        st    = self.state
        now   = st.end or time.time()
        table = Table.grid(padding=(0, 1), expand=True)
        table.add_column(width=2, no_wrap=True)
        table.add_column(ratio=1, no_wrap=True, overflow="ellipsis")
        table.add_column(justify="right", no_wrap=True)

        for i, label in enumerate(st.steps):
            if st.failed == i:
                icon, name, tm = (Text("✗", style="bold red"),
                                  Text(label, style="bold red"),
                                  Text(_fmt_step_time(now - st.step_start), style="red"))
            elif i in st.timings:
                icon, name, tm = (Text("✓", style="bold green"),
                                  Text(label, style="default"),
                                  Text(_fmt_step_time(st.timings[i]), style="dim green"))
            elif i == st.current and not st.done and st.failed is None:
                extra = f"  [{st.batch}]" if st.batch else ""
                icon, name, tm = (self.spinner,
                                  Text(label + extra, style="bold cyan"),
                                  Text(_fmt_step_time(now - st.step_start), style="bold cyan"))
            else:
                icon, name, tm = Text("·", style="dim"), Text(label, style="dim"), Text("")
            table.add_row(icon, name, tm)

        done_n = len(st.timings)
        return Panel(table, title="[bold]Pipeline[/]", title_align="left",
                     subtitle=f"[dim]{done_n}/{len(st.steps)} steps[/]", subtitle_align="right",
                     border_style="#0284c7", padding=(0, 1))

    def _log_panel(self):
        from rich.panel import Panel

        st = self.state
        if st.failed is not None:
            title = "[bold]Log[/]  [red]failed[/]"
        elif st.done:
            title = "[bold]Log[/]  [green]complete[/]"
        else:
            title = f"[bold]Log[/]  [dim]{st.steps[st.current]}[/]"
        return Panel(_LogTail(self.records), title=title, title_align="left",
                     border_style="#0284c7", padding=(0, 1))

    def _footer(self):
        from rich.panel import Panel
        from rich.progress_bar import ProgressBar
        from rich.table import Table

        st    = self.state
        total = len(st.steps)
        done  = len(st.timings)
        if st.failed is not None:
            status, colour = "[bold red]failed[/]", "red"
        elif st.done:
            status, colour = "[bold green]complete[/]", "green"
        else:
            status, colour = "[cyan]running[/]", "cyan"

        counters = []
        if st.warnings:
            counters.append(f"[yellow]▲ {st.warnings} warning{'s' if st.warnings > 1 else ''}[/]")
        if st.errors:
            counters.append(f"[red]✗ {st.errors} error{'s' if st.errors > 1 else ''}[/]")

        row = Table.grid(padding=(0, 2), expand=True)
        row.add_column(no_wrap=True)
        row.add_column(ratio=1)
        row.add_column(no_wrap=True, justify="right")
        row.add_row(
            status,
            ProgressBar(total=total, completed=done, complete_style=colour,
                        finished_style=colour, style="grey23"),
            f"[bold]{done}[/][dim]/{total}[/]   "
            + ("   ".join(counters) + "   " if counters else "")
            + f"[dim]elapsed[/] [bold]{_fmt_clock(st.elapsed())}[/]",
        )
        return Panel(row, border_style="#0284c7", padding=(0, 1))

    # ── layout ────────────────────────────────────────────────────────────
    def __rich_console__(self, console, options):
        from rich.layout import Layout

        width = options.max_width
        header, header_h = self._header(width)
        steps_w = min(max(len(s) for s in self.state.steps) + 14, 72, int(width * 0.55))

        layout = Layout()
        layout.split_column(
            Layout(header,          name="header", size=header_h),
            Layout(name="body",     ratio=1),
            Layout(self._footer(),  name="footer", size=3),
        )
        if width >= 100:
            layout["body"].split_row(
                Layout(self._steps_panel(), name="steps", size=steps_w),
                Layout(self._log_panel(),   name="log",   ratio=1),
            )
        else:   # narrow terminal: stack steps over log
            layout["body"].split_column(
                Layout(self._steps_panel(), name="steps", size=len(self.state.steps) + 2),
                Layout(self._log_panel(),   name="log",   ratio=1),
            )
        yield layout


# ── Pipeline ──────────────────────────────────────────────────────────────────

# Long-read technologies. Salmon quantifies via short-read selective-alignment
# mapping (k-mer index + --validateMappings), which is not valid for long,
# error-prone reads, so the Salmon step is skipped for these. Read coverage
# (minimap2) still runs — it is long-read aware via per-technology presets.
_LONG_READ_TYPES = ("ont", "pb", "hifi")


def _salmon_enabled(args) -> bool:
    """
    Salmon runs only for short reads and only when not explicitly skipped;
    long reads use coverage (minimap2) only.
    """
    return (bool(args.reads)
            and not args.skip_salmon
            and args.read_type not in _LONG_READ_TYPES)


def _build_steps(args) -> list[str]:
    steps = [
        f"Parse FASTA{'  +  CAP3' if args.cap3 else ''}",
        "Find ORFs",
        "Diamond BLASTx  —  RefSeq viral filter",
        "HMMsearch  —  RVDB · Vfam · EggNOG  (viral confirmation)",
    ]
    if args.nr_db:
        steps.append(f"Diamond BLASTx NR  —  {Path(args.nr_db).name}")
    if args.blastn_local:
        steps.append(f"BLASTn local  —  {Path(args.blastn_local).name}")
    elif args.blastn_online:
        steps.append(f"BLASTn online  —  {args.blastn_online_db}")
    steps.append("HMMsearch  —  Pfam  (functional annotation)")
    steps.append("Taxonomy annotation")
    steps.append("Cluster sequences by species")
    if _salmon_enabled(args):
        mode = "reference" if args.transcriptome else "de novo"
        steps.append(f"Salmon quantification  —  {mode}")
    if not args.skip_seq_quality:
        steps.append("Sequence quality  —  dustmask · self-BLAST · jellyfish")
    if args.reads:
        steps.append("Read coverage profiling")
    steps.append("Heuristic scoring  —  rule-based vq_score")
    if args.model_type and args.model_name:
        token_tag = f"[{args.llm_tokens}]" if args.llm_tokens else ""
        steps.append(f"LLM scoring  —  {args.model_type} / {args.model_name}  {token_tag}".rstrip())
    steps.append("Export JSON report")
    steps.append("Build HTML report")
    return steps


def _workflow_options(args) -> dict:
    """User-chosen options shown in the report's workflow card (never secrets)."""
    if args.blastn_local:
        blastn = f"local · {Path(args.blastn_local).name}"
    elif args.blastn_online:
        blastn = f"online · {args.blastn_online_db}"
    else:
        blastn = None
    return {
        "input":            Path(args.input).name,
        "threads":          args.cpu,
        "cap3":             bool(args.cap3),
        "nr_db":            Path(args.nr_db).name if args.nr_db else None,
        "blastn":           blastn,
        "reads":            [Path(r).name for r in (args.reads or [])],
        "read_type":        args.read_type if args.reads else None,
        "transcriptome":    Path(args.transcriptome).name if args.transcriptome else None,
        "skip_salmon":      bool(args.skip_salmon),
        "skip_seq_quality": bool(args.skip_seq_quality),
        "llm":              (f"{args.model_type} / {args.model_name}"
                             if args.model_type and args.model_name else None),
        "force":            bool(args.force),
    }


def _run_pipeline(args):
    """
    Execute every pipeline step.

    Yields tagged event tuples:
      ("step",  step_idx, elapsed_seconds)   — a full pipeline step completed
      ("batch", step_idx, batch_num, total)  — one Diamond batch completed (RefSeq only)
    """
    from loguru import logger

    from .parser      import FastaParser, Cap3Runner
    from .orfs        import OrfAnalyzer
    from .diamond     import (DiamondRunner, DiamondOutputParser,
                              DiamondResultAttacher, DiamondPhase)
    from .hmm         import (HmmMetadataLoader, HmmSequencePreparer,
                              HmmSearcher, HmmResultAttacher, HmmViralFlagSetter)
    from .blastn      import BlastnRunner, BlastnMode, BlastnResultAttacher
    from .tax         import TaxonomyAnnotator, ViralFamilyAnnotator
    from .track_seqs  import SequenceTracker
    from .exporter    import ReportExporter
    from .html_report import write_report
    from .output      import OutputOrganizer

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem       = Path(args.input).stem
    organizer  = OutputOrganizer(outdir, stem)
    step       = 0
    labels     = _build_steps(args)
    run_start  = time.time()
    started_at = datetime.now(timezone.utc).isoformat()
    workflow: list[dict] = []   # one entry per step (incl. skipped) → report's workflow card

    def _record(key: str, elapsed: float, status: str, details, message) -> None:
        workflow.append({
            "key":     key,
            "label":   labels[step] if step < len(labels) else key,
            "status":  status,               # done | partial | error | skipped
            "seconds": round(elapsed, 2),
            "details": details or {},
            "message": message,
        })

    def _tick(start: float, key: str, details: dict | None = None,
              status: str = "done", message: str | None = None):
        nonlocal step
        elapsed = time.time() - start
        _record(key, elapsed, status, details, message)
        yield ("step", step, elapsed)
        step += 1

    def _skip(key: str, label: str, reason: str) -> None:
        workflow.append({"key": key, "label": label, "status": "skipped",
                         "seconds": None, "details": {}, "message": reason})

    # ── 1. Parse / assemble ───────────────────────────────────────────────────
    t  = time.time()
    fp = FastaParser(args.input)
    fp.read_input_file()
    input_fasta = fp.input_fasta

    cap3_info = None   # populated only when --cap3 is used; drives the report's CAP3 card
    blastn_info = None # populated only when BLASTn runs; records failed queries
    if args.cap3:
        cap3_res   = Cap3Runner(args.input, outdir=str(outdir / "cap3")).cap3_runner()
        contigs_p  = FastaParser(str(cap3_res.contigs));  contigs_p.read_input_file()
        singlets_p = FastaParser(str(cap3_res.singlets)); singlets_p.read_input_file()
        seqs = contigs_p.sequences + singlets_p.sequences
        cap3_info = {
            "used":     True,
            "contigs":  len(contigs_p.sequences),
            "singlets": len(singlets_p.sequences),
        }
        assembled_fasta = outdir / "cap3" / "assembled.fasta"
        cap3_res.get_combined_fasta(assembled_fasta)
    else:
        seqs = fp.sequences
        assembled_fasta = Path(args.input)
    parse_details = {"Input sequences": len(fp.sequences)}
    if cap3_info:
        parse_details |= {"CAP3 contigs": cap3_info["contigs"],
                          "CAP3 singlets": cap3_info["singlets"]}
    parse_details["Sequences carried forward"] = len(seqs)
    yield from _tick(t, "parse", parse_details)

    # ── 2. ORF finding ────────────────────────────────────────────────────────
    t        = time.time()
    analyzer = OrfAnalyzer()
    for seq in seqs:
        analyzer.find_n_save_orfs(seq)
    yield from _tick(t, "orfs", {"ORFs found": sum(len(s.orfs) for s in seqs)})

    # ── 3. Diamond RefSeq filter ──────────────────────────────────────────────
    t    = time.time()
    dmnd = DiamondRunner(db_path=str(args._db["viral_dmnd"]), threads=args.cpu,
                         outdir=str(organizer.diamond_dir))
    tsvs:  list[Path] = []
    hits:  list        = []
    for batch_num, n_batches, tsv in dmnd.run_batched(seqs, DiamondPhase.REFSEQ_FILTER):
        tsvs.append(tsv)
        hits.extend(DiamondOutputParser.parse(tsv))
        yield ("batch", step, batch_num, n_batches)
    DiamondResultAttacher.attach(
        hits, seqs, DiamondPhase.REFSEQ_FILTER,
        min_identity = 0.0  if args.nr_db else 50.0,
        min_coverage = 0.0  if args.nr_db else 40.0,
    )
    organizer.finalize_diamond_refseq(tsvs)
    yield from _tick(t, "refseq", {
        "Database":            Path(args._db["viral_dmnd"]).name,
        "Sequences with hits": sum(1 for s in seqs if s.blastx_hits),
        "Flagged viral":       sum(1 for s in seqs if s.is_viral),
        "Filter":              "none (NR confirms)" if args.nr_db else "identity ≥ 50% · coverage ≥ 40%",
    })

    # ── 4. HMM viral confirmation (RVDB + Vfam + EggNOG) ─────────────────────
    t        = time.time()
    searcher = HmmSearcher(cpus=args.cpu)
    seq_block, orf_map = HmmSequencePreparer.prepare(seqs)
    if seq_block is not None:
        for db_key, json_path, db_name in args._hmm_filter:
            meta      = HmmMetadataLoader.load(str(json_path), db_name)
            hmm_hits  = searcher.search(seq_block, str(args._db[db_key]))
            HmmResultAttacher.attach(hmm_hits, orf_map, meta, db_name)
            organizer.save_hmm_table(hmm_hits, db_name)
    HmmViralFlagSetter.flag(seqs)
    yield from _tick(t, "hmm", {
        "Databases":                  " · ".join(db_name for _, _, db_name in args._hmm_filter),
        "Flagged viral (RefSeq+HMM)": sum(1 for s in seqs if s.is_viral),
    })

    # ── 5. Diamond NR characterisation (optional) ─────────────────────────────
    if args.nr_db:
        t       = time.time()
        viral   = [s for s in seqs if s.is_viral]
        dmnd_nr = DiamondRunner(
            db_path      = str(args.nr_db),
            threads      = args.cpu,
            outdir       = str(organizer.diamond_dir),
            block_size   = args.nr_block_size,
            index_chunks = args.nr_index_chunks,
            tmpdir       = args.nr_tmpdir,
        )
        nr_tsv  = dmnd_nr.run_single(viral, DiamondPhase.NR_CHARACTERIZE)
        if nr_tsv:
            hits_nr = DiamondOutputParser.parse(nr_tsv)
            DiamondResultAttacher.attach(hits_nr, seqs, DiamondPhase.NR_CHARACTERIZE)
        organizer.finalize_diamond_nr(nr_tsv)
        yield from _tick(t, "nr", {
            "Database":        Path(args.nr_db).name,
            "Queried":         len(viral),
            "Viral NR hits":   sum(1 for s in seqs if s.blastx_nr_hits),
        })
    else:
        _skip("nr", "Diamond BLASTx NR", "--nr-db not set")

    # ── 6. BLASTn ─────────────────────────────────────────────────────────────
    if args.blastn_local or args.blastn_online:
        t = time.time()
        if args.blastn_local:
            blastn = BlastnRunner(mode=BlastnMode.LOCAL, db_path=args.blastn_local,
                                  threads=args.cpu, outdir=str(organizer.blastn_dir))
        else:
            blastn = BlastnRunner(mode=BlastnMode.ONLINE,
                                  outdir=str(organizer.blastn_dir),
                                  online_db=args.blastn_online_db,
                                  email=args.blastn_online)
        # NR run → only NR-confirmed sequences; no NR → all viral sequences
        blastn_seqs = (
            [s for s in seqs if s.blastx_nr_hits]
            if args.nr_db
            else [s for s in seqs if s.is_viral]
        )
        blastn_hits = blastn.run(blastn_seqs)
        BlastnResultAttacher.attach(blastn_hits, seqs)
        blastn_info = {
            "run":        True,
            "mode":       blastn.mode.value,
            "db":         args.blastn_local or args.blastn_online_db,
            "queried":    len(blastn_seqs),
            "failed":     len(blastn.failed_ids),
            "failed_ids": blastn.failed_ids,
        }
        bn_status, bn_msg = "done", None
        if blastn_seqs and len(blastn.failed_ids) == len(blastn_seqs):
            bn_status = "error"
            bn_msg    = "every query failed — see the log for the cause"
            logger.error(
                "BLASTn: every query failed — the report will carry no BLASTn "
                "hits. See the log above for the cause."
            )
        elif blastn.failed_ids:
            bn_status = "partial"
            bn_msg    = f"{len(blastn.failed_ids)} of {len(blastn_seqs)} queries failed"
            logger.warning(
                f"BLASTn: {len(blastn.failed_ids)}/{len(blastn_seqs)} query(ies) "
                f"failed — listed in pipeline_stats.blastn.failed_ids."
            )
        organizer.save_blastn_table(blastn_hits)
        yield from _tick(t, "blastn", {
            "Mode":      blastn_info["mode"],
            "Database":  (Path(args.blastn_local).name if args.blastn_local
                          else args.blastn_online_db),
            "Queried":   len(blastn_seqs),
            "With hits": sum(1 for s in seqs if s.blastn_hits),
            "Failed":    len(blastn.failed_ids),
        }, status=bn_status, message=bn_msg)
    else:
        _skip("blastn", "BLASTn", "--blastn-local / --blastn-online not set")

    # ── 7. Pfam characterisation (confirmed viral only) ───────────────────────
    t          = time.time()
    viral_seqs = [s for s in seqs if s.is_viral]
    if viral_seqs:
        pfam_block, pfam_orf_map = HmmSequencePreparer.prepare(viral_seqs)
        if pfam_block is not None:
            pfam_meta = HmmMetadataLoader.load(
                str(args._index_dir / "Pfam-index.json"), "Pfam"
            )
            pfam_hits = HmmSearcher(cpus=args.cpu).search(
                pfam_block, str(args._db["pfam"])
            )
            HmmResultAttacher.attach(pfam_hits, pfam_orf_map, pfam_meta, "Pfam")
            organizer.save_hmm_table(pfam_hits, "Pfam")
    yield from _tick(t, "pfam", {
        "Queried":      len(viral_seqs),
        "Pfam domains": sum(1 for s in viral_seqs for o in s.orfs
                            for d in o.domains if d.database == "Pfam"),
    })

    # ── 8. Taxonomy annotation ────────────────────────────────────────────────
    t = time.time()
    TaxonomyAnnotator(str(args._db["viral_tax"])).annotate(seqs)
    ViralFamilyAnnotator(
        str(args._db["fam_high"]), str(args._db["fam_low"])
    ).annotate(seqs)
    yield from _tick(t, "taxonomy", {
        "With taxonomy": sum(1 for s in seqs if s.taxonomy is not None),
    })

    # ── 9. Clustering ─────────────────────────────────────────────────────────
    t            = time.time()
    tracker      = SequenceTracker(min_identity=args.min_identity, min_coverage=args.min_coverage)
    viral_for_cl = [s for s in seqs if s.is_viral]
    clusters     = tracker.track(viral_for_cl)
    yield from _tick(t, "clusters", {
        "Viral sequences": len(viral_for_cl),
        "Clusters":        len(clusters),
    })

    # ── 10. Salmon quantification ─────────────────────────────────────────────
    salmon_report = None
    if not args.reads:
        _skip("salmon", "Salmon quantification", "--reads not set")
    elif not _salmon_enabled(args):
        _skip("salmon", "Salmon quantification",
              "--skip-salmon" if args.skip_salmon
              else f"long-read input (--read-type {args.read_type})")
        if args.skip_salmon:
            logger.warning(
                "Salmon quantification skipped (--skip-salmon). Read coverage "
                "(incl. sense/antisense) and sequence-quality signals still run."
            )
        else:
            logger.warning(
                f"Salmon quantification skipped: --read-type '{args.read_type}' is a long-read "
                "technology, and Salmon's short-read mapping mode is not valid for long reads. "
                "Read coverage (minimap2) still runs and is long-read aware."
            )
    if _salmon_enabled(args):
        t = time.time()
        from .salmon_quant import SalmonQuantPipeline
        salmon_error = None
        try:
            salmon_report = SalmonQuantPipeline(
                threads=args.cpu,
                low_memory=args.low_memory,
                pfam_hmm_path=Path(args._db["pfam"]),
            ).run(
                reads              = args.reads,
                viral_seqs         = seqs,
                outdir             = outdir / "salmon",
                assembled_fasta    = assembled_fasta,
                user_transcriptome = Path(args.transcriptome) if args.transcriptome else None,
                hk_genes_file      = Path(args.hk_genes) if args.hk_genes else None,
            )
        except Exception as exc:
            salmon_error = str(exc) or type(exc).__name__
            logger.error(f"Salmon quantification failed — skipping: {exc}")

        # Populate salmon_tpm / salmon_reads on each viral sequence so the
        # LLM scorer (step 11) can include expression data in its prompt.
        if salmon_report:
            tpm_map = {
                e.name.removeprefix("VQ_VIRAL_"): (e.tpm, e.num_reads)
                for e in salmon_report.viral_quant
            }
            for seq in seqs:
                if seq.id in tpm_map:
                    seq.salmon_tpm, seq.salmon_reads = tpm_map[seq.id]

        yield from _tick(t, "salmon", {
            "Mode":                 "reference" if args.transcriptome else "de novo",
            "Read files":           len(args.reads),
            "Viral quantified":     len(salmon_report.viral_quant) if salmon_report else 0,
        }, status="error" if salmon_error else "done", message=salmon_error)

    # ── 10a. Sequence quality (dustmask · self-BLAST · jellyfish · dot plot) ──
    # Structural signals on the confirmed viral sequences — same set the coverage
    # track and exporter use, so every signal aligns with the report viewer.
    # Reads are not needed: the signals come from the assembled sequences alone,
    # so this runs in every mode unless --skip-seq-quality is given.
    if not args.skip_seq_quality:
        t = time.time()
        from .seq_quality import SequenceQualityPipeline
        from .exporter    import select_confirmed_sequences
        viral_for_sq = select_confirmed_sequences(seqs, force=args.force, nr_run=bool(args.nr_db))
        sq_map = SequenceQualityPipeline(
            threads=args.cpu, kmer_size=args.kmer,
        ).run(viral_seqs=viral_for_sq, outdir=outdir / "seq_quality")
        for seq in seqs:
            if seq.id in sq_map:
                seq.seq_quality = sq_map[seq.id]
        yield from _tick(t, "seq_quality", {
            "Sequences analysed": len(viral_for_sq),
            "With signals":       len(sq_map),
            "k-mer size":         args.kmer if args.kmer is not None else "auto",
        })
    else:
        _skip("seq_quality", "Sequence quality", "--skip-seq-quality")

    # ── 10b. Read coverage (per-base depth track) ─────────────────────────────
    # Runs for any --reads input (short or long); minimap2 is long-read aware,
    # so this is the coverage signal used when Salmon is skipped for long reads.
    # Profile exactly the sequences that will appear in the report, so the
    # coverage set matches the viewer in every mode (--nr-db, no --nr-db, --force).
    if args.reads:
        t = time.time()
        from .coverage import CoveragePipeline
        from .exporter import select_confirmed_sequences
        viral_for_cov = select_confirmed_sequences(seqs, force=args.force, nr_run=bool(args.nr_db))
        cov_map = CoveragePipeline(
            threads=args.cpu, low_memory=args.low_memory,
            read_type=args.read_type,
        ).run(
            reads      = args.reads,
            viral_seqs = viral_for_cov,
            outdir     = outdir / "coverage",
        )
        for seq in seqs:
            seq.coverage = cov_map.get(seq.id)
        yield from _tick(t, "coverage", {
            "Read type":          args.read_type,
            "Sequences profiled": len(viral_for_cov),
            "With coverage":      sum(1 for s in viral_for_cov if s.coverage is not None),
        })
    else:
        _skip("coverage", "Read coverage profiling", "--reads not set")

    # ── 11. Heuristic scoring (always; no LLM, no NR, no API key) ─────────────
    # Scores exactly the sequences the exporter will emit, so the report's
    # mandatory heuristic field is populated for every exported sequence.
    t = time.time()
    from .exporter import select_confirmed_sequences
    from .score_heuristic import HeuristicScorer
    heur_seqs = select_confirmed_sequences(seqs, force=args.force, nr_run=bool(args.nr_db))
    HeuristicScorer().score(heur_seqs)
    yield from _tick(t, "heuristic", {"Scored": len(heur_seqs)})

    # ── 12. LLM scoring (optional) ────────────────────────────────────────────
    # Score the same set the report emits (single source of truth): NR-confirmed
    # in the --nr-db pathway, is_viral in RefSeq/HMM-only runs, all under --force.
    if args.model_type and args.model_name:
        t          = time.time()
        from .exporter import select_confirmed_sequences
        viral_seqs = select_confirmed_sequences(seqs, force=args.force, nr_run=bool(args.nr_db))
        from .score_ai import SequenceScorer, LlmMode
        mode = LlmMode.HIGH if args.llm_tokens == "high" else LlmMode.LOW
        SequenceScorer(model_type=args.model_type,
                       model_name=args.model_name,
                       mode=mode,
                       api_key=args.api_key).score(viral_seqs)
        scored     = [s for s in viral_seqs if s.llm_output]
        api_errors = sum(1 for s in scored if s.llm_output.classification == "api-error")
        llm_status, llm_msg = "done", None
        if scored and api_errors == len(scored):
            llm_status, llm_msg = "error", "every request returned an API error"
        elif api_errors:
            llm_status, llm_msg = "partial", f"{api_errors} of {len(scored)} requests returned an API error"
        yield from _tick(t, "llm", {
            "Model":      f"{args.model_type} / {args.model_name}",
            "Scored":     len(scored),
            "API errors": api_errors,
        }, status=llm_status, message=llm_msg)
    else:
        _skip("llm", "LLM scoring", "--model-type / --model-name not set")

    # ── 13. Export: viral FASTA + JSON ────────────────────────────────────────
    t             = time.time()
    fasta_seqs    = (
        [s for s in seqs if s.blastx_nr_hits]
        if args.nr_db
        else [s for s in seqs if s.is_viral]
    )
    viral_fasta   = organizer.save_viral_contigs(fasta_seqs)
    json_path     = outdir / f"{stem}_viralquest.json"
    exporter      = ReportExporter(force=args.force, nr_run=bool(args.nr_db))
    report        = exporter.export(
        nuc_seqs=seqs,
        clusters=clusters,
        input_fasta=input_fasta,
        output_path=None,           # written below, once the workflow is complete
        version=__version__,
        salmon_report=salmon_report,
        cap3=cap3_info,
        blastn=blastn_info,
    )
    # The export step is recorded before the JSON is written so the workflow
    # it carries is complete (JSON write time itself is not included).
    _record("export", time.time() - t, "done", {
        "Confirmed sequences": len(report["sequences"]),
        "Clusters":            len(report["clusters"]),
        "Viral FASTA":         Path(viral_fasta).name if viral_fasta else "—",
    }, None)
    report["pipeline_stats"]["workflow"] = {
        "started_at":    started_at,
        "finished_at":   datetime.now(timezone.utc).isoformat(),
        "total_seconds": round(time.time() - run_start, 2),
        "options":       _workflow_options(args),
        "steps":         workflow,
    }
    ReportExporter.write(report, json_path)
    yield ("step", step, time.time() - t)
    step += 1

    # ── 14. HTML report ───────────────────────────────────────────────────────
    t         = time.time()
    html_path = outdir / f"{stem}_viralquest.html"
    write_report(report, html_path)
    # Not recorded in the workflow: the JSON (and the workflow it carries) is
    # already written, so only the progress event is emitted here.
    yield ("step", step, time.time() - t)
    step += 1

    # Stash for summary
    args._result_seqs         = seqs
    args._result_clusters     = clusters
    args._result_json         = json_path
    args._result_html         = html_path
    args._result_viral_fasta  = viral_fasta
    args._result_diamond_dir  = organizer.diamond_dir
    args._result_blastn_dir   = organizer.blastn_dir if (args.blastn_local or args.blastn_online) else None
    args._result_hmm_dir      = organizer.hmm_dir


# ── Default mode (plain loguru) ───────────────────────────────────────────────

def _run_default(args) -> None:
    from loguru import logger
    from rich.console import Console
    from .output import OutputOrganizer

    console  = Console(stderr=True)
    steps    = _build_steps(args)
    timings: dict[int, float] = {}

    outdir   = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_id   = OutputOrganizer(outdir, Path(args.input).stem).add_log_sink()

    try:
        for event in _run_pipeline(args):
            if event[0] == "step":
                _, step_idx, elapsed = event
                timings[step_idx] = elapsed
                label = steps[step_idx] if step_idx < len(steps) else "?"
                logger.success(f"  ✓  {label}  ({elapsed:.1f}s)")
            # "batch" events are already logged inside DiamondRunner
    finally:
        OutputOrganizer.remove_log_sink(log_id)

    _print_summary(console, args, timings)


# ── Live mode (Rich Layout) ───────────────────────────────────────────────────

def _run_live(args) -> None:
    from loguru import logger
    from rich.console import Console
    from rich.live import Live
    from .output import OutputOrganizer

    console = Console()
    steps   = _build_steps(args)
    records: deque = deque(maxlen=500)   # (time, level, message) — sized for tall terminals
    state   = _LiveState(steps)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def _sink(message) -> None:
        rec   = message.record
        level = rec["level"].name
        records.append((rec["time"].strftime("%H:%M:%S"), level, rec["message"]))
        if level == "WARNING":
            state.warnings += 1
        elif level in ("ERROR", "CRITICAL"):
            state.errors += 1

    # Redirect loguru to the live screen + log file (remove default stderr sink first)
    logger.remove()
    logger.add(_sink, level="INFO")
    log_id = OutputOrganizer(outdir, Path(args.input).stem).add_log_sink()

    info   = _workflow_options(args) | {"outdir": outdir}
    screen = _LiveScreen(state, records, info)

    try:
        with Live(screen, console=console, refresh_per_second=10, screen=False):
            try:
                for event in _run_pipeline(args):
                    if event[0] == "step":
                        _, step_idx, elapsed = event
                        state.finish_step(step_idx, elapsed)
                    elif event[0] == "batch":
                        _, _step_idx, batch_num, total = event
                        state.batch = f"{batch_num}/{total}"
            except BaseException:
                # Freeze the screen on the failing step before the traceback prints.
                state.fail()
                raise
    finally:
        OutputOrganizer.remove_log_sink(log_id)
        # Restore loguru to stderr after Live exits
        logger.remove()
        logger.add(sys.stderr, colorize=True,
                   format="{time:HH:mm:ss} | <level>{level:<8}</level> | {message}")

    _print_summary(console, args, state.timings)


# ── Summary panel ─────────────────────────────────────────────────────────────

def _print_summary(console, args, timings: dict) -> None:
    from rich.panel import Panel

    seqs         = getattr(args, "_result_seqs",         [])
    clusters     = getattr(args, "_result_clusters",     [])
    json_p       = getattr(args, "_result_json",         Path(args.outdir))
    html_p       = getattr(args, "_result_html",         Path(args.outdir))
    viral_fasta  = getattr(args, "_result_viral_fasta",  None)
    diamond_dir  = getattr(args, "_result_diamond_dir",  None)
    hmm_dir      = getattr(args, "_result_hmm_dir",      None)

    viral_count = sum(1 for s in seqs if s.is_viral)
    total_time  = sum(timings.values())

    blastn_dir = getattr(args, "_result_blastn_dir", None)

    lines = (
        f"[bold green]Confirmed viral sequences:[/]  {viral_count}\n"
        f"[bold green]Clusters:[/]                   {len(clusters)}\n"
        f"[bold green]Viral contigs FASTA:[/]        {viral_fasta}\n"
        f"[bold green]Diamond results:[/]             {diamond_dir}\n"
        + (f"[bold green]BLASTn results:[/]              {blastn_dir}\n" if blastn_dir else "")
        + f"[bold green]HMM tables:[/]                 {hmm_dir}\n"
        f"[bold green]JSON report:[/]                {json_p}\n"
        f"[bold green]HTML report:[/]                {html_p}\n"
        f"[bold green]Log file:[/]                   {Path(args.outdir) / 'viralquest.log'}\n"
        f"[bold green]Total time:[/]                 {total_time:.1f}s"
    )

    console.print()
    console.print(Panel(
        lines,
        title=f"[bold green]ViralQuest — {Path(args.input).name} complete[/bold green]",
        border_style="green", width=85,
    ))


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    try:
        import setproctitle
        setproctitle.setproctitle("viralquest")
    except ImportError:
        pass

    # Check bioinformatics binaries before doing anything else.
    # missing_tools() activates the pixi environment first — the user may have
    # run viralquest-setup but not reloaded their shell PATH.
    from viralquest.setup_env import missing_tools
    absent = missing_tools()

    if absent:
        print(
            f"[ERROR] Missing required tools: {', '.join(absent)}\n"
            "Run:  viralquest-setup\n"
            "This will install pixi and all bioinformatics dependencies automatically."
        )
        sys.exit(1)

    from rich.console import Console
    console = Console(stderr=True)
    parser  = _build_parser()

    if len(sys.argv) == 1:
        console.print(
            "[bold red]ERROR:[/bold red] No arguments provided. "
            "Use [bold cyan]-h[/] or [bold cyan]--help[/] for usage."
        )
        sys.exit(1)

    args = parser.parse_args()

    if args.help:
        _show_rich_help()
        sys.exit(0)

    ext_db_dir = Path(args.db_dir) if args.db_dir else None
    args._db, args._index_dir, args._hmm_filter = _resolve_db_paths(ext_db_dir)

    _validate_args(args, console)
    if not _check_databases(console, args._db):
        sys.exit(1)

    if args.live:
        _run_live(args)
    else:
        _run_default(args)


if __name__ == "__main__":
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
    from viralquest.cli import main
    main()
