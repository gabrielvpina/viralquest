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
__author__  = "gabrielvpina"

# ── Auto-detected database paths ──────────────────────────────────────────────

_ROOT    = Path(__file__).parent.parent
_DATA    = _ROOT / "data"
_HMM_DIR = _DATA / "hmm-dbs"

_FAMILY_INFO_DIR = _DATA / "viral-family-info"

_DB = {
    "viral_dmnd":   _DATA            / "viralDB.dmnd",
    "rvdb":         _HMM_DIR         / "U-RVDBv29.0-prot.hmm",
    "vfam":         _HMM_DIR         / "Vfam-228.hmm",
    "eggnog":       _HMM_DIR         / "EggNOG-4.5.hmm",
    "pfam":         _HMM_DIR         / "Pfam-A.hmm",
    "viral_tax":    _DATA            / "viralTax.json.xz",
    "fam_high":     _FAMILY_INFO_DIR / "viral_info_highToken.json",
    "fam_low":      _FAMILY_INFO_DIR / "viral_info_lowToken.json",
}

_INDEX_DIR = _DATA / "hmm-index"

# (hmm_path_key, json_index_path, db_name_for_metadata)
_HMM_FILTER = [
    ("rvdb",   _INDEX_DIR / "RVDB-index.json",       "RVDB"),
    ("vfam",   _INDEX_DIR / "Vfam-index.json",       "Vfam"),
    ("eggnog", _INDEX_DIR / "eggNOG-4.5-index.json", "EggNOG"),
]

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
        default=70.0, metavar="PCT",
        help="Minimum %% identity for intra-cluster BLASTn alignment (default: 90.0).")

    # Salmon quantification ────────────────────────────────────────────────────
    sal = parser.add_argument_group("salmon quantification (optional)")
    sal.add_argument("--transcriptome", dest="transcriptome", type=str,
        metavar="HOST.fasta",
        help="Host transcriptome FASTA for Salmon quantification.")
    sal.add_argument("--reads", dest="reads", nargs="+", type=str,
        metavar="READS.fastq",
        help="FASTQ file(s) for Salmon: one = single-end, two = paired-end.")

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

    # Output mode ──────────────────────────────────────────────────────────────
    parser.add_argument("--live", action="store_true",
        help="Rich Live display: banner + scrolling log box + step progress bar.")

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
    tagline.append(f"  v{__version__} · {__author__}\n", style="dim")
    console.print(tagline)

    console.print(Panel(
        "[bold white]-in / --input[/]    [dim]FASTA[/]\n"
        "  Input nucleotide FASTA file.\n\n"
        "[bold white]-out / --outdir[/]  [dim]DIR[/]\n"
        "  Output directory (created if absent).",
        title="[bold red]REQUIRED[/bold red]",
        border_style="red", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]--cap3[/]\n"
        "  Assemble overlapping sequences with CAP3 before analysis.\n\n"
        "[bold white]-cpu / --cpu[/]       [dim]N[/]  (default: 2)\n"
        "  CPU threads used by Diamond, HMMsearch, BLASTn, and Salmon.\n\n"
        "[bold white]--min-identity[/]    [dim]%[/]  (default: 70.0)\n"
        "  Minimum %% identity for intra-cluster alignment.\n\n"
        "[bold white]--force[/]\n"
        "  Export all input sequences even if viral confirmation fails.",
        title="[bold green]PIPELINE OPTIONS[/bold green]",
        border_style="green", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]-nr / --nr-db[/]         [dim]NR.dmnd[/]\n"
        "  Diamond-format NCBI NR database for full protein-level characterisation\n"
        "  of confirmed viral sequences. Optional but recommended for complete\n"
        "  annotation. NR is a large database (>100 GB); build with:\n"
        "  [dim]diamond makedb --in nr.fasta --db nr[/dim]\n\n"
        "[bold white]--nr-block-size[/]       [dim]N[/]  (default: Diamond built-in ~2.0)\n"
        "  GB of RAM loaded per database pass. Higher = fewer passes = faster.\n"
        "  Examples: 6 → ~45 GB RAM · 12 → ~90 GB RAM.\n\n"
        "[bold white]--nr-index-chunks[/]     [dim]N[/]  (default: 4)\n"
        "  Seed-index chunks. Lower values (1–2) reduce disk passes but use\n"
        "  more RAM. Combine with --nr-block-size for maximum speed.\n\n"
        "[bold white]--nr-tmpdir[/]           [dim]DIR[/]\n"
        "  Directory for Diamond temporary files. Use a fast NVMe drive or\n"
        "  RAM disk ([dim]/dev/shm[/dim]) to eliminate I/O as a bottleneck.",
        title="[bold yellow]NR DATABASE (optional)[/bold yellow]",
        border_style="yellow", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]-n / --blastn-local[/]   [dim]BLAST_DB[/]\n"
        "  Path to a local BLAST nucleotide database (e.g. nt).\n\n"
        "[bold white]--blastn-online[/]       [dim]EMAIL[/]\n"
        "  NCBI e-mail for web BLASTn — no local database required.\n\n"
        "[bold white]--blastn-online-db[/]    [dim]DB[/]  (default: nt)\n"
        "  NCBI database to query when using --blastn-online.\n\n"
        "[bold yellow]Note:[/] --blastn-local and --blastn-online are mutually exclusive.",
        title="[bold yellow]BLASTN (choose one)[/bold yellow]",
        border_style="yellow", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]--transcriptome[/]  [dim]HOST.fasta[/]\n"
        "  Host transcriptome FASTA. Salmon quantifies viral + housekeeping\n"
        "  gene expression in the context of this assembly.\n\n"
        "[bold white]--reads[/]          [dim]R1.fastq [R2.fastq][/]\n"
        "  FASTQ file(s): one = single-end, two = paired-end.\n\n"
        "[bold yellow]Note:[/] --transcriptome and --reads must be used together.",
        title="[bold cyan]SALMON QUANTIFICATION (optional)[/bold cyan]",
        border_style="cyan", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]--model-type[/]   [dim]ollama | openai | anthropic | google[/]\n"
        "  AI provider for LLM viral sequence scoring.\n\n"
        "[bold white]--model-name[/]   [dim]MODEL[/]\n"
        "  Model identifier (e.g. 'qwen3:4b', 'gpt-4o', 'claude-opus-4-7',\n"
        "  'gemini-2.5-pro').\n\n"
        "[bold white]--llm-tokens[/]   [dim]high | low[/]  [bold red](required with --model-type)[/]\n"
        "  Token usage mode for LLM prompts:\n"
        "    [bold]high[/] — all hits, Pfam details, full taxonomy, full family description.\n"
        "    [bold]low[/]  — best hits only, no Pfam details, compact taxonomy.\n\n"
        "[bold white]--api-key[/]      [dim]KEY[/]\n"
        "  API key for cloud providers (not required for ollama).",
        title="[bold magenta]AI SCORING (optional)[/bold magenta]",
        border_style="magenta", width=85, box=box.ROUNDED,
    ))
    console.print()

    console.print(Panel(
        "[bold white]--live[/]\n"
        "  Rich Live display: ASCII banner + scrolling log box + step progress bar.\n"
        "  Default (without flag): raw loguru log stream to stderr.\n\n"
        "[bold white]-h / --help[/]     Show this help message and exit.\n"
        "[bold white]-v / --version[/]  Show version and exit.",
        title="[bold bright_black]OTHER[/bold bright_black]",
        border_style="bright_black", width=85, box=box.ROUNDED,
    ))


# ── Validation ────────────────────────────────────────────────────────────────

def _check_databases(console) -> bool:
    missing = [name for name, path in _DB.items() if not path.exists()]
    if not missing:
        return True
    console.print(f"\n[bold red]ERROR:[/bold red] {len(missing)} database(s) not found:\n")
    for name in missing:
        console.print(f"  [red]✗[/red]  {name}  [dim]{_DB[name]}[/dim]")
    console.print(
        "\n[bold yellow]Fix:[/bold yellow] run  "
        "[bold white]python download_dbs.py[/bold white]  to fetch all required databases.\n"
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
    if bool(args.transcriptome) ^ bool(args.reads):
        errors.append("--transcriptome and --reads must be provided together.")
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
              border_style="bright_black", padding=(0, 1))
    )

    log_content = "\n".join(log_buf) if log_buf else "[dim]Waiting for output…[/dim]"
    layout["log"].update(
        Panel(log_content,
              title=f"[bold cyan]Log[/bold cyan]  "
                    f"[dim]{steps[current] if current < len(steps) else 'done'}[/dim]",
              border_style="bright_black",
              subtitle=f"[dim]step {min(current+1, len(steps))}/{len(steps)}[/dim]")
    )

    layout["progress"].update(Panel(progress, border_style="bright_black", padding=(0, 1)))

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
    if args.transcriptome:
        steps.append("Salmon quantification")
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
    else:
        seqs = fp.sequences
    yield from _tick(t)

    # ── 2. ORF finding ────────────────────────────────────────────────────────
    t        = time.time()
    analyzer = OrfAnalyzer()
    for seq in seqs:
        analyzer.find_n_save_orfs(seq)
    yield from _tick(t)

    # ── 3. Diamond RefSeq filter ──────────────────────────────────────────────
    t    = time.time()
    dmnd = DiamondRunner(db_path=str(_DB["viral_dmnd"]), threads=args.cpu,
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
        for db_key, json_path, db_name in _HMM_FILTER:
            meta      = HmmMetadataLoader.load(str(json_path), db_name)
            hmm_hits  = searcher.search(seq_block, str(_DB[db_key]))
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
                str(_INDEX_DIR / "Pfam-index.json"), "Pfam"
            )
            pfam_hits = HmmSearcher(cpus=args.cpu).search(
                pfam_block, str(_DB["pfam"])
            )
            HmmResultAttacher.attach(pfam_hits, pfam_orf_map, pfam_meta, "Pfam")
            organizer.save_hmm_table(pfam_hits, "Pfam")
    yield from _tick(t)

    # ── 8. Taxonomy annotation ────────────────────────────────────────────────
    t = time.time()
    TaxonomyAnnotator(str(_DB["viral_tax"])).annotate(seqs)
    ViralFamilyAnnotator(
        str(_DB["fam_high"]), str(_DB["fam_low"])
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
        t = time.time()
        from .score_ai import SequenceScorer, LlmMode
        mode = LlmMode.HIGH if args.llm_tokens == "high" else LlmMode.LOW
        SequenceScorer(model_type=args.model_type,
                       model_name=args.model_name,
                       mode=mode,
                       api_key=args.api_key).score(seqs)
        yield from _tick(t)

    # ── 11. Salmon quantification ─────────────────────────────────────────────
    salmon_report = None
    if args.transcriptome and args.reads:
        t = time.time()
        from .salmon_quant import SalmonQuantPipeline
        salmon_report = SalmonQuantPipeline().run(
            user_transcriptome=Path(args.transcriptome),
            reads=args.reads,
            viral_seqs=seqs,
            outdir=outdir / "salmon",
        )
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
        title=f"[bold white]ViralQuest — {Path(args.input).name} complete[/bold white]",
        border_style="green", width=85,
    ))


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    from rich.console import Console
    console = Console(stderr=True)
    parser  = _build_parser()

    if len(sys.argv) == 1:
        console.print(
            "[bold red]ERROR:[/bold red] No arguments provided. "
            "Use [bold white]-h[/] or [bold white]--help[/] for usage."
        )
        sys.exit(1)

    args = parser.parse_args()

    if args.help:
        _show_rich_help()
        sys.exit(0)

    _validate_args(args, console)
    if not _check_databases(console):
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
