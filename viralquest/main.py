import argparse, sys, os
#import pyfiglet
import time

# rich dependencies
from rich.console import Console
from rich.panel import Panel
from rich.align import Align
from rich.columns import Columns
from rich.text import Text
from rich import box

from .info import __author__, __version__, __date__
from .cli import final_colored_ascii_banner


parser = argparse.ArgumentParser(
    description=final_colored_ascii_banner,
    add_help=False,  
    formatter_class=argparse.RawDescriptionHelpFormatter
)

# custom help
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit.')

# args
parser.add_argument("-in", "--input", type=str, dest="input")
parser.add_argument("-ref", "--viralRef", type=str, dest="viralRef")
parser.add_argument("-out", "--outdir", type=str)
parser.add_argument("--cap3", action="store_true")
parser.add_argument("-n", "--blastn_local", type=str)
parser.add_argument("--blastn_online", type=str, dest="blastn_online")
parser.add_argument("--blastn_onlineDB", type=str, default="nt", dest="blastn_onlineDB")
parser.add_argument("-x", "--diamond_blastx", type=str, dest="diamond_blastx")
parser.add_argument("-rvdb", "--rvdb_hmm", type=str, dest="rvdb_hmm")
parser.add_argument("-eggnog", "--eggnog_hmm", type=str, dest="eggnog_hmm")
parser.add_argument("-vfam", "--vfam_hmm", type=str, dest="vfam_hmm")
parser.add_argument("-pfam", "--pfam_hmm", type=str, dest="pfam_hmm")
parser.add_argument("-orf", "--maxORFs", type=int, dest="maxORFs", default=2)
parser.add_argument("-cpu", "--cpu", type=int, dest="cpu", default=2)
parser.add_argument("-dmnd_path", "--diamond_path", type=str, dest="diamond_path", default="None")
parser.add_argument("-v", "--version", action="version", version=f"ViralQuest v{__version__}")
parser.add_argument("--merge-json", type=str, dest="merge_json")
parser.add_argument('--model-type', required=False, choices=['ollama', 'openai', 'anthropic', 'google'])
parser.add_argument('--model-name', required=False, type=str)
parser.add_argument('--api-key', required=False, type=str)


# rich console
console = Console()

def show_rich_help():
    
    # ASCII banner
    #ascii_banner = pyfiglet.figlet_format("ViralQuest")
    
    # help content
    help_text = Text()
    # pyfiglet -  help_text.append(ascii_banner, style="bold blue")
    print('')
    print(final_colored_ascii_banner)
    help_text.append("\nA tool for viral diversity analysis and characterization.\n", style="italic")
    help_text.append("More info: https://github.com/gabrielvpina/viralquest\n\n", style="dim blue underline")
    #console.print(final_colored_ascii_banner)
    console.print(help_text)
    # eequired arguments section
    required_panel = Panel(
        "[bold white]--input/-in[/bold white]\n"
        "  .fasta file input. It's recomended a short name file (e.g.'CTL3.fasta') \n\n"
        "[bold white]--outdir/-out[/bold white]\n"
        "  Directory where the output files will be saved.\n\n"
        "[bold white]--viralRef/-ref[/bold white]\n"
        "  RefSeq Viral Protein Release file. Path to .dmnd (diamond db) file\n\n"
        "[bold white]--blastn_online[/bold white]\n"
        "  NCBI email to execute online BLASTn search using NCBI BLAST web service.\n\n"
        "[bold white]--diamond_blastx/-x[/bold white]\n"
        "  Path to the Diamond BLASTx database (.dmnd) for protein sequence comparison.\n\n"
        "[bold white]--pfam_hmm/-pfam[/bold white]\n"
        "  Path to the Pfam hmm for conserved domain analysis.",
        title="[bold red]REQUIRED ARGUMENTS[/bold red]",
        border_style="red",
        width=85, 
        box=box.ROUNDED
    )

    # Viral HMM 
    viral_hmm_panel = Panel(
        "[bold white]--rvdb_hmm/-rvdb[/bold white]\n"
        "  Path to the RVDB hmm for conserved domain analysis.\n\n"
        "[bold white]--eggnog_hmm/-eggnog[/bold white]\n"
        "  Path to the EggNOG hmm for conserved domain analysis.\n\n"
        "[bold white]--vfam_hmm/-vfam[/bold white]\n"
        "  Path to the Vfam hmm for conserved domain analysis.\n\n"
        "[bold yellow]Note:[/bold yellow] At least one of these is required.",
        title="[bold yellow]VIRAL HMM DATABASES[/bold yellow]",
        border_style="yellow",
        width=85, 
        box=box.ROUNDED
    )
    
    # optional arguments 
    optional_panel = Panel(
        "[bold white]--blastn_local/-n[/bold white]\n"
        "  Path to the BLASTn database for nucleotide sequence comparison.\n\n"
        "[bold white]--blastn_onlineDB[/bold white]\n"
        "  NCBI Nucleotide database for online BLASTn web service (DEFAULT='nt').\n\n"
        "[bold white]--maxORFs/-orf[/bold white]\n"
        "  Number of max largest non-overlapping ORFs from sequence (DEFAULT=2).\n\n"
        "[bold white]--cpu/-cpu[/bold white]\n"
        "  Number of CPU threads (DEFAULT=2).\n\n"
        "[bold white]--cap3[/bold white]\n"
        "  Activate CAP3 fasta assembly: Deactivated by default.",
        title="[bold green]OPTIONAL ARGUMENTS[/bold green]",
        border_style="green",
        width=85, 
        box=box.ROUNDED
    )
    
    # AI Summary section
    ai_panel = Panel(
        "[bold white]--model-type[/bold white]\n"
        "  Type of model to use for analysis (ollama, openai, anthropic, google).\n\n"
        "[bold white]--model-name[/bold white]\n"
        "  Name of the model (e.g., \"qwen3:4b\" for ollama, \"gpt-3.5-turbo\" for OpenAI).\n\n"
        "[bold white]--api-key[/bold white]\n"
        "  API key for cloud models (required for OpenAI, Anthropic, Google).",
        title="[bold magenta]AI SUMMARY (OPTIONAL)[/bold magenta]",
        border_style="magenta",
        width=85, 
        box=box.ROUNDED
    )
    
    # reports section
    merge_panel = Panel(
        "[bold white]--merge-json[/bold white]\n"
        "  Input Type: dir/ Merge JSON files in a directory to create a general ViralQuest HTML report. When used, other arguments are ignored.",
        title="[bold blue]MERGE REPORTS[/bold blue]",
        border_style="blue",
        width=85, 
        box=box.ROUNDED
    )
    
    # panels
    console.print(required_panel)
    console.print()
    console.print(viral_hmm_panel)
    console.print()
    console.print(optional_panel)
    console.print()
    console.print(ai_panel)
    console.print()
    console.print(merge_panel)
    console.print()
    
    # help
    help_footer = Panel(
        "[bold white]--help/-h[/bold white] Show this help message and exit\n"
        "[bold white]--version/-v[/bold white] Show program's version number and exit",
        title="[bold cyan]OTHER OPTIONS[/bold cyan]",
        border_style="cyan",
        width=85, 
        box=box.ROUNDED
    )
    console.print(help_footer)

# Show message if no arguments are provided
if len(sys.argv) == 1:
    console.print(
        "[bold red]ERROR:[/bold red] No arguments provided. Use '-h' or '--help' to see available options.",
        style="red"
    )
    sys.exit(1)

# Check for help before creating parser
if "-h" in sys.argv or "--help" in sys.argv:
    show_rich_help()
    sys.exit(0)

# Parse normal
args = parser.parse_args()



################################## MERGE JSON ########################################################

from .mergeJSON import merge_json_files

if args.merge_json:
    # warn the user
    if args.input or args.viralRef or args.outdir:
        print("Warning: --merge-json is exclusive. Other arguments will be ignored.")
    
    # merge_json function
    merge_json_files(args.merge_json)
    print("The merged report was saved as 'ViralQuestReport.html'")
    sys.exit(0) 
else:
    #  argument validation
    if not args.input or not args.viralRef or not args.outdir:
        parser.error("--input, --viralRef, and --outdir are required unless --merge-json is used")



############################## CHECK ARGUMENTS ###########################################################



####################################

from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.live import Live
from rich.spinner import Spinner
from rich.table import Table
from rich import box

def create_progress_table(steps, current_step, step_status):
    """Create a clean progress table with status indicators"""
    table = Table(
        box=box.ROUNDED,
        expand=True,
        show_header=True,
        border_style="bright_black",
        width=80,
        padding=(0, 1)
    )
    table.add_column("Status", style="cyan", width=15, justify="center")
    table.add_column("Current Process", style="white", width=50, justify="center")
    table.add_column("Time", style="dim", width=15, justify="center")

    # Create a spinner
    spinner = Spinner("dots", style="cyan")
    
    for i, step in enumerate(steps):
        if step is None:
            continue
            
        if i < current_step:
            # Completed step
            elapsed = step_status.get(i, {}).get('elapsed', '')
            table.add_row("✓", f"[blue]{step}[/]", f"[dim]{elapsed}[/]")
        elif i == current_step:
            # Current step
            table.add_row(spinner, f"[cyan bold]{step}[/]", "[dim]running...[/]")
        else:
            # Pending step
            table.add_row("⏳", f"[dim]{step}[/]", "")
    
    return table



def main():
    print(" ")
    print(final_colored_ascii_banner)

    from .cli import trackErrors
    trackErrors(args)

    from .inputFASTA import validateFasta, cap3, noCAP3
    from .processFASTA import copyfasta, filterfasta, renameFasta
    from .processORFs import findorf, ORFs
    from .initialBLAST import filterBLAST
    from .hmmsearch import process_VFAM, process_RVDB, process_EggNOG
    from .pfamDomains import pfamFasta, process_Pfam, generateFasta, mergeFASTA
    from .dmndBLASTx import diamond_blastx, generateFasta_blastn
    from .runBlastn import blastn, blastn_online
    from .ai_summary import analyze_viral_sequences
    from .finalAssets import finalTable, exec_cleaning
    from .finalReport import generate_html_report

    time_start = time.time()
    steps = [
        f"Check FASTA file - {args.input}",
        f"Finding ORFs - {args.maxORFs} maximum ORFs",
        f"Filter sequences - Viral RefSeq alignment",
        "HMMsearch RVDB - Detect viral elements" if args.rvdb_hmm else None,
        "HMMsearch Vfam - Detect viral elements" if args.vfam_hmm else None,
        "HMMsearch EggNOG - Detect viral elements" if args.eggnog_hmm else None,
        "HMMsearch Pfam - Functional Domains" if args.pfam_hmm else None,
        "Running BLASTx - Viral characterization" if args.diamond_blastx else None,
        "Running BLASTn - Viral characterization" if args.blastn_local else None,
        "Running Online BLASTn - Viral characterization" if args.blastn_online else None,
        "Generating final table",
        f"Generating AI Summary - {args.model_name}" if args.model_type and args.model_name else None,
        "Generating HTML report"
    ]
    steps = [step for step in steps if step is not None] 
    
    console = Console()
    console.print(Panel(Align.center(f"[bold white]Running ViralQuest[/bold white] 🔍"), border_style="blue", width=80))

    active_steps = [step for step in steps if step is not None]
    step_status = {}
    current_step = 0

    with Live(create_progress_table(active_steps, current_step, step_status), 
              console=console, refresh_per_second=8) as live:
        
        # Processing FASTA
        step_start = time.time()  # Start timing for this step
        is_valid = validateFasta(args.input)
        if not is_valid:
            print("FASTA validation failed. Check non-nuclotide characters in input file\n")
            sys.exit(1) 
        
        if args.cap3:
            cap3(args.input, args.outdir)
            copyfasta(args.outdir)
            filterfasta(args.outdir)
            renameFasta(args.outdir)
        else:
            noCAP3(args.input, args.outdir)
            copyfasta(args.outdir)
            filterfasta(args.outdir)
            renameFasta(args.outdir)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))
        
        # ORFs profile
        step_start = time.time()  # Reset timing for this step
        findorf(args.outdir)
        ORFs(args.outdir, args.maxORFs)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))
        
        # Filter w/ refseq sequences
        step_start = time.time()  # Reset timing for this step
        filterBLAST(args.outdir, args.viralRef, args.cpu, diamond_path=None, log_file=None)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))
        
        # HMM Search
        if args.pfam_hmm and not args.rvdb_hmm and not args.eggnog_hmm and not args.vfam_hmm:
            console.print("[bold red]Error:[/] --pfam_hmm requires --rvdb_hmm/--eggnog_hmm/--vfam_hmm to be specified.")
            sys.exit(1)

        if args.rvdb_hmm:
            step_start = time.time()  # Reset timing for this step
            process_RVDB(args.outdir, args.rvdb_hmm, args.cpu, score_threshold=10)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
            
        if args.vfam_hmm:
            step_start = time.time()  # Reset timing for this step
            process_VFAM(args.outdir, args.vfam_hmm, args.cpu, score_threshold=10)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
            
        if args.eggnog_hmm:
            step_start = time.time()  # Reset timing for this step
            process_EggNOG(args.outdir, args.eggnog_hmm, args.cpu, score_threshold=10)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
            
        step_start = time.time()  # Reset timing for this step
        pfamFasta(args.outdir)
        process_Pfam(args.outdir, args.pfam_hmm, args.cpu)
        generateFasta(args.outdir)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))
        
        mergeFASTA(args.outdir)
        
        # BLAST
        if args.diamond_blastx:
            step_start = time.time()  # Reset timing for this step
            diamond_blastx(args.outdir, args.diamond_blastx, args.cpu)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
            
        if args.blastn_local:
            step_start = time.time()  # Reset timing for this step
            generateFasta_blastn(args.outdir)
            blastn(args.outdir, args.blastn_local, args.cpu)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))

        if args.blastn_online:
            step_start = time.time()  # Reset timing for this step
            generateFasta_blastn(args.outdir)
            blastn_online(args.outdir, args.blastn_onlineDB, args.blastn_online)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
            
        # Final table
        step_start = time.time()  # Reset timing for this step
        finalTable(args.outdir, args)
        exec_cleaning(args.outdir)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))

        # AI summary
        if args.model_type and args.model_name:
            step_start = time.time()  # Reset timing for this step
            analyze_viral_sequences(args)
            step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
            current_step += 1
            live.update(create_progress_table(active_steps, current_step, step_status))
        
        # HTML report
        step_start = time.time()  # Reset timing for this step

        if args.cap3:
            cap3check = "true"
        else:
            cap3check = "false"

        input_repo = os.path.basename(args.input)
        outdir_repo = os.path.basename(args.outdir)

        if args.blastn_local:
            blastn_repo= os.path.basename(args.blastn_local)
        if args.blastn_online:
            blastn_repo= f"Online DB - {args.blastn_onlineDB}"

        diamond_blastx_repo = os.path.basename(args.diamond_blastx)
        pfam_hmm_repo = os.path.basename(args.pfam_hmm)
        cpu_repo = int(args.cpu)

        from .processFASTA import filterfasta, countOriginalFasta
        from .processORFs import countBiggestORFsFasta
        from .finalAssets import getViralSeqsNumber

        filteredSeqs = filterfasta(f"{args.outdir}/fasta-files")
        originalSeqs = countOriginalFasta(f"{args.outdir}/fasta-files")
        numberTotalORFs = countBiggestORFsFasta(f"{args.outdir}/fasta-files")
        number_viralSeqs = getViralSeqsNumber(args.outdir)

        generate_html_report(args.outdir, cap3check, input_repo, outdir_repo, blastn_repo, diamond_blastx_repo, pfam_hmm_repo, filteredSeqs, originalSeqs, numberTotalORFs, number_viralSeqs, cpu_repo)
        step_status[current_step] = {'elapsed': f"{time.time() - step_start:.1f}s"}
        current_step += 1
        live.update(create_progress_table(active_steps, current_step, step_status))
    
    time_end = time.time()
    console.print(f"[bold green]The pipeline took [/][bold yellow]{time_end - time_start:.2f}[/][bold green] seconds to run.[/]")
    
    # Log file generation (unchanged)
    name = os.path.basename(args.input)
    name = name.replace(".fasta", ".log")
    log_path = os.path.join(args.outdir, name)
    
    with open(log_path, "a") as log_file:
        log_file.write(f"""ViralQuest v{__version__} --- Start in: {time.strftime('%Y-%m-%d %H:%M:%S')}
#############################################
Input file: {args.input}
Output directory: {args.outdir}
Number of ORFs: {args.maxORFs}
CAP3 Assembly: {args.cap3}
CPU cores: {args.cpu}
BLASTn database: {args.blastn_local}
Diamond BLASTx database: {args.diamond_blastx}
RVDB hmm profile database: {args.rvdb_hmm}
Pfam hmm profile database: {args.pfam_hmm}
#############################################

The pipeline takes {time_end - time_start:.2f} seconds to run.
    """)


if __name__ == "__main__":
        main()
