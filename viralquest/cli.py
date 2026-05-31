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

import sys
import time
from collections import deque
from pathlib import Path

# ── Version ───────────────────────────────────────────────────────────────────

__version__ = "3.0.0"
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

    # Salmon quantification ────────────────────────────────────────────────────
    sal = parser.add_argument_group("salmon quantification (optional)")
    sal.add_argument("--transcriptome", dest="transcriptome", type=str,
        metavar="HOST.fasta",
        help="Host transcriptome FASTA for Salmon quantification.")
    sal.add_argument("--reads", dest="reads", nargs="+", type=str,
        metavar="READS.fastq",
        help="FASTQ file(s) for Salmon: one = single-end, two = paired-end.")
    sal.add_argument("--hk-genes", dest="hk_genes", type=str,
        default=None, metavar="IDS.txt",
        help="Text file with reference housekeeping gene IDs (one per line) for "
             "normalization. IDs must exactly match the first word of the FASTA "
             "header in --transcriptome (case-sensitive). Requires --transcriptome.")
    sal.add_argument("--low-memory", dest="low_memory", action="store_true",
        help="Low-RAM mode for Salmon: builds the index with a smaller k-mer size "
             "(-k 25 instead of 31) to reduce the SSHash index footprint, and runs "
             "quant with GC-bias correction at reduced memory cost "
             "(--gcBias --reduceGCMemory). Recommended on machines with <16 GB RAM.")

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
        "  Minimum %% identity for intra-cluster alignment.\n\n"
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
        "[bold cyan]--transcriptome[/]  [dim]HOST.fasta[/]\n"
        "  Host transcriptome FASTA. When provided, runs the [bold]reference pathway[/bold]:\n"
        "  viral + bundled HK + ref-HK + transcriptome as a combined Salmon index.\n"
        "  Without --transcriptome, runs the [bold]de-novo pathway[/bold]: quantifies\n"
        "  assembled contigs directly, using BLASTn to find HK-matching contigs.\n\n"
        "[bold cyan]--hk-genes[/]       [dim]IDS.txt[/]\n"
        "  Text file with reference HK gene IDs (one per line) for normalization.\n"
        "  IDs must exactly match the first word of the --transcriptome header\n"
        "  (case-sensitive). Requires --transcriptome.\n\n"
        "[bold cyan]--low-memory[/]\n"
        "  Low-RAM mode: builds the Salmon index with a smaller k-mer size\n"
        "  ([dim]-k 25[/dim] instead of the default 31), reducing the SSHash index\n"
        "  memory footprint, and runs quant with GC-bias correction at reduced\n"
        "  memory cost ([dim]--gcBias --reduceGCMemory[/dim]).\n"
        "  Recommended on machines with <16 GB RAM.\n\n"
        "[bold]Note:[/] --reads alone triggers de-novo mode. --transcriptome requires --reads.",
        title="[bold cyan]SALMON QUANTIFICATION (optional)[/bold cyan]",
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
    if args.model_type and not args.model_name:
        errors.append("--model-name is required when --model-type is set.")
    if args.model_type and not args.llm_tokens:
        errors.append("--llm-tokens (high|low) is required when --model-type is set.")
    if args.llm_tokens and not args.model_type:
        errors.append("--llm-tokens requires --model-type to be set.")
    if args.nr_db and not Path(args.nr_db).exists():
        errors.append(f"NR database not found: {args.nr_db}")
    for msg in errors:
        console.print(f"[bold red]ERROR:[/bold red] {msg}")
    if errors:
        sys.exit(1)


# ── Live-mode display helpers ─────────────────────────────────────────────────

def _build_live_display(log_buf: deque, steps: list[str], current: int, progress):
    from rich.layout import Layout
    from rich.panel import Panel
    from rich.text import Text

    layout = Layout()
    layout.split_column(
        Layout(name="banner",   size=9),
        Layout(name="log",      ratio=1),
        Layout(name="progress", size=3),
    )

    layout["banner"].update(
        Panel(Text(_BANNER_RAW, style="bold cyan", no_wrap=True),
              border_style="blue", padding=(0, 1))
    )

    log_content = "\n".join(log_buf) if log_buf else "[dim]Waiting for output…[/dim]"
    layout["log"].update(
        Panel(log_content,
              title=f"[bold cyan]Log[/bold cyan]  "
                    f"[dim]{steps[current] if current < len(steps) else 'done'}[/dim]",
              border_style="blue",
              subtitle=f"[dim]step {min(current+1, len(steps))}/{len(steps)}[/dim]")
    )

    layout["progress"].update(Panel(progress, border_style="blue", padding=(0, 1)))

    return layout


# ── Pipeline ──────────────────────────────────────────────────────────────────

def _build_steps(args) -> list[str]:
    steps = [
        f"Parse FASTA{'  +  CAP3' if args.cap3 else ''}",
        "Find ORFs",
        "Diamond BLASTx  —  RefSeq viral filter",
        "HMMsearch  —  RVDB · Vfam · EggNOG  (viral confirmation)",
    ]
    if args.blastn_local:
        steps.append(f"BLASTn local  —  {Path(args.blastn_local).name}")
    elif args.blastn_online:
        steps.append(f"BLASTn online  —  {args.blastn_online_db}")
    if args.nr_db:
        steps.append(f"Diamond BLASTx NR  —  {Path(args.nr_db).name}")
    steps.append("HMMsearch  —  Pfam  (functional annotation)")
    steps.append("Taxonomy annotation")
    steps.append("Cluster sequences by species")
    if args.model_type:
        token_tag = f"[{args.llm_tokens}]" if args.llm_tokens else ""
        steps.append(f"LLM scoring  —  {args.model_type} / {args.model_name}  {token_tag}".rstrip())
    if args.reads:
        mode = "reference" if args.transcriptome else "de novo"
        steps.append(f"Salmon quantification  —  {mode}")
    steps.append("Export JSON report")
    steps.append("Build HTML report")
    return steps


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

    def _tick(start: float):
        nonlocal step
        elapsed = time.time() - start
        yield ("step", step, elapsed)
        step += 1

    # ── 1. Parse / assemble ───────────────────────────────────────────────────
    t  = time.time()
    fp = FastaParser(args.input)
    fp.read_input_file()
    input_fasta = fp.input_fasta

    if args.cap3:
        cap3_res   = Cap3Runner(args.input, outdir=str(outdir / "cap3")).cap3_runner()
        contigs_p  = FastaParser(str(cap3_res.contigs));  contigs_p.read_input_file()
        singlets_p = FastaParser(str(cap3_res.singlets)); singlets_p.read_input_file()
        seqs = contigs_p.sequences + singlets_p.sequences
        assembled_fasta = outdir / "cap3" / "assembled.fasta"
        cap3_res.get_combined_fasta(assembled_fasta)
    else:
        seqs = fp.sequences
        assembled_fasta = Path(args.input)
    yield from _tick(t)

    # ── 2. ORF finding ────────────────────────────────────────────────────────
    t        = time.time()
    analyzer = OrfAnalyzer()
    for seq in seqs:
        analyzer.find_n_save_orfs(seq)
    yield from _tick(t)

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
    DiamondResultAttacher.attach(hits, seqs, DiamondPhase.REFSEQ_FILTER)
    organizer.finalize_diamond_refseq(tsvs)
    yield from _tick(t)

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
    yield from _tick(t)

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
        yield from _tick(t)

    # ── 6. BLASTn ─────────────────────────────────────────────────────────────
    if args.blastn_local or args.blastn_online:
        t = time.time()
        if args.blastn_local:
            blastn = BlastnRunner(mode=BlastnMode.LOCAL, db_path=args.blastn_local,
                                  threads=args.cpu, outdir=str(organizer.blastn_dir))
        else:
            blastn = BlastnRunner(mode=BlastnMode.ONLINE,
                                  outdir=str(organizer.blastn_dir))
        # NR run → only NR-confirmed sequences; no NR → all viral sequences
        blastn_seqs = (
            [s for s in seqs if s.blastx_nr_hits]
            if args.nr_db
            else [s for s in seqs if s.is_viral]
        )
        blastn_hits = blastn.run(blastn_seqs)
        BlastnResultAttacher.attach(blastn_hits, seqs)
        organizer.save_blastn_table(blastn_hits)
        yield from _tick(t)

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
    yield from _tick(t)

    # ── 8. Taxonomy annotation ────────────────────────────────────────────────
    t = time.time()
    TaxonomyAnnotator(str(args._db["viral_tax"])).annotate(seqs)
    ViralFamilyAnnotator(
        str(args._db["fam_high"]), str(args._db["fam_low"])
    ).annotate(seqs)
    yield from _tick(t)

    # ── 9. Clustering ─────────────────────────────────────────────────────────
    t            = time.time()
    tracker      = SequenceTracker(min_identity=args.min_identity)
    viral_for_cl = [s for s in seqs if s.is_viral]
    clusters     = tracker.track(viral_for_cl)
    yield from _tick(t)

    # ── 10. LLM scoring ───────────────────────────────────────────────────────
    if args.model_type and args.model_name:
        t          = time.time()
        viral_seqs = [s for s in seqs if s.blastx_nr_hits]
        from .score_ai import SequenceScorer, LlmMode
        mode = LlmMode.HIGH if args.llm_tokens == "high" else LlmMode.LOW
        SequenceScorer(model_type=args.model_type,
                       model_name=args.model_name,
                       mode=mode,
                       api_key=args.api_key).score(viral_seqs)
        yield from _tick(t)

    # ── 11. Salmon quantification ─────────────────────────────────────────────
    salmon_report = None
    if args.reads:
        t = time.time()
        from .salmon_quant import SalmonQuantPipeline
        try:
            salmon_report = SalmonQuantPipeline(threads=args.cpu, low_memory=args.low_memory).run(
                reads              = args.reads,
                viral_seqs         = seqs,
                outdir             = outdir / "salmon",
                assembled_fasta    = assembled_fasta,
                user_transcriptome = Path(args.transcriptome) if args.transcriptome else None,
                hk_genes_file      = Path(args.hk_genes) if args.hk_genes else None,
            )
        except Exception as exc:
            logger.error(f"Salmon quantification failed — skipping: {exc}")
        yield from _tick(t)

    # ── 12. Export: viral FASTA + JSON ────────────────────────────────────────
    t             = time.time()
    viral_fasta   = organizer.save_viral_contigs(seqs)
    json_path     = outdir / f"{stem}_viralquest.json"
    exporter      = ReportExporter(force=args.force)
    report        = exporter.export(
        nuc_seqs=seqs,
        clusters=clusters,
        input_fasta=input_fasta,
        output_path=json_path,
        version=__version__,
        salmon_report=salmon_report,
    )
    yield from _tick(t)

    # ── 13. HTML report ───────────────────────────────────────────────────────
    t         = time.time()
    html_path = outdir / f"{stem}_viralquest.html"
    write_report(report, html_path)
    yield from _tick(t)

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
    from rich.panel import Panel
    from rich.progress import (Progress, BarColumn, TextColumn,
                               MofNCompleteColumn, TimeElapsedColumn)
    from .output import OutputOrganizer

    console  = Console()
    steps    = _build_steps(args)
    log_buf  = deque(maxlen=28)
    timings: dict[int, float] = {}

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Redirect loguru to deque + log file (remove default stderr sink first)
    logger.remove()
    logger.add(
        lambda msg: log_buf.append(msg.rstrip()),
        format="{time:HH:mm:ss} | <level>{level:<8}</level> | {message}",
        colorize=False,
    )
    log_id = OutputOrganizer(outdir, Path(args.input).stem).add_log_sink()

    progress = Progress(
        TextColumn("[bold cyan]{task.description}[/]"),
        BarColumn(bar_width=38, style="cyan", complete_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        expand=False,
    )
    task_id = progress.add_task(steps[0], total=len(steps))

    try:
        with Live(
            _build_live_display(log_buf, steps, 0, progress),
            console=console,
            refresh_per_second=10,
            screen=False,
        ) as live:
            for event in _run_pipeline(args):
                if event[0] == "step":
                    _, step_idx, elapsed = event
                    timings[step_idx] = elapsed
                    progress.advance(task_id)
                    next_step = step_idx + 1
                    if next_step < len(steps):
                        progress.update(task_id, description=steps[next_step])
                    else:
                        progress.update(task_id, description="[green]complete[/]")
                    live.update(_build_live_display(log_buf, steps, next_step, progress))

                elif event[0] == "batch":
                    _, step_idx, batch_num, total = event
                    label = steps[step_idx] if step_idx < len(steps) else "?"
                    progress.update(
                        task_id,
                        description=f"{label}  [dim][{batch_num}/{total}][/dim]",
                    )
                    live.update(_build_live_display(log_buf, steps, step_idx, progress))
    finally:
        OutputOrganizer.remove_log_sink(log_id)

    # Restore loguru to stderr after Live exits
    logger.remove()
    logger.add(sys.stderr, colorize=True,
               format="{time:HH:mm:ss} | <level>{level:<8}</level> | {message}")

    _print_summary(console, args, timings)


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
    # If tools are missing, try to auto-activate the pixi environment first —
    # the user may have run viralquest-setup but not reloaded their shell PATH.
    import os
    from viralquest.setup_env import missing_tools, _find_pixi_toml
    absent = missing_tools()
    if absent:
        pixi_toml = _find_pixi_toml()
        if pixi_toml:
            env_bin = pixi_toml.parent / ".pixi" / "envs" / "default" / "bin"
            if env_bin.exists():
                os.environ["PATH"] = str(env_bin) + ":" + os.environ.get("PATH", "")
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
