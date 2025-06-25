import argparse, sys, os
import pyfiglet
import time

# rich dependencies
from rich.console import Console
from rich.panel import Panel
from rich.columns import Columns
from rich.text import Text
from rich import box

from .info import __author__, __version__, __date__


ascii_banner = pyfiglet.figlet_format("ViralQuest")

# deactivate help
parser = argparse.ArgumentParser(
    description=ascii_banner,
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
parser.add_argument("-N", "--blastn", type=str)
parser.add_argument("--blastn_online", type=str, dest="blastn_online")
parser.add_argument("--blastn_onlineDB", type=str, default="nt", dest="blastn_onlineDB")
parser.add_argument("-dX", "--diamond_blastx", type=str, dest="diamond_blastx")
parser.add_argument("-rvdb", "--rvdb_hmm", type=str, dest="rvdb_hmm")
parser.add_argument("-eggnog", "--eggnog_hmm", type=str, dest="eggnog_hmm")
parser.add_argument("-vfam", "--vfam_hmm", type=str, dest="vfam_hmm")
parser.add_argument("-pfam", "--pfam_hmm", type=str, dest="pfam_hmm")
parser.add_argument("-maxORFs", "--maxORFs", type=int, dest="maxORFs", default=2)
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
    """Display beautiful help using rich formatting"""
    
    # ASCII banner
    ascii_banner = pyfiglet.figlet_format("ViralQuest")
    
    # help content
    help_text = Text()
    help_text.append(ascii_banner, style="bold cyan")
    help_text.append("\nA tool for viral diversity analysis and characterization.\n", style="italic")
    help_text.append("More info: https://github.com/gabrielvpina/viralquest\n\n", style="dim blue underline")
    
    console.print(help_text)
    
    # eequired arguments section
    required_panel = Panel(
        "[bold white]-in/--input[/bold white]\n"
        "  Fasta file to be analyzed. It's recomended a short name file (e.g. 'CTL3.fasta') \n\n"
        "[bold white]-out/--outdir[/bold white]\n"
        "  Directory where the output files will be saved.\n\n"
        "[bold white]-ref/--viralRef[/bold white]\n"
        "  RefSeq Viral Protein Release file. Path to .dmnd file\n\n"
        "[bold white]--blastn_online[/bold white]\n"
        "  NCBI email to execute online BLASTn search using NCBI BLAST web service.\n\n"
        "[bold white]--diamond_blastx/-dX[/bold white]\n"
        "  Path to the Diamond BLASTx database (.dmnd) for protein sequence comparison.\n\n"
        "[bold white]-pfam/--pfam_hmm[/bold white]\n"
        "  Path to the Pfam hmm for conserved domain analysis.",
        title="[bold red]REQUIRED ARGUMENTS[/bold red]",
        border_style="red",
        box=box.ROUNDED
    )

    # Viral HMM 
    viral_hmm_panel = Panel(
        "[bold white]-rvdb/--rvdb_hmm[/bold white]\n"
        "  Path to the RVDB hmm for conserved domain analysis.\n\n"
        "[bold white]-eggnog/--eggnog_hmm[/bold white]\n"
        "  Path to the EggNOG hmm for conserved domain analysis.\n\n"
        "[bold white]-vfam/--vfam_hmm[/bold white]\n"
        "  Path to the Vfam hmm for conserved domain analysis.\n\n"
        "[bold yellow]Note:[/bold yellow] At least one of these is required.",
        title="[bold yellow]VIRAL HMM DATABASES[/bold yellow]",
        border_style="yellow",
        box=box.ROUNDED
    )
    
    # optional arguments 
    optional_panel = Panel(
        "[bold white]-N/--blastn[/bold white]\n"
        "  Path to the BLASTn database for nucleotide sequence comparison.\n\n"
        "[bold white]--blastn_onlineDB[/bold white]\n"
        "  NCBI Nucleotide database for online BLASTn web service (DEFAULT='nt').\n\n"
        "[bold white]-maxORFs/--maxORFs[/bold white]\n"
        "  Number of largest ORFs to select from the input sequences (DEFAULT=2).\n\n"
        "[bold white]-cpu/--cpu[/bold white]\n"
        "  Number of CPU threads (DEFAULT=2).\n\n"
        "[bold white]-dmnd_path/--diamond_path[/bold white]\n"
        "  Diamond executable application path for BLAST databases: path/to/diamond\n\n"
        "[bold white]--cap3[/bold white]\n"
        "  Activate CAP3 fasta assembly: Deactivated by default.",
        title="[bold green]OPTIONAL ARGUMENTS[/bold green]",
        border_style="green",
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
        box=box.ROUNDED
    )
    
    # reports section
    merge_panel = Panel(
        "[bold white]--merge-json[/bold white]\n"
        "  Merge JSON files in a directory to create a general ViralQuest HTML report.When used, other arguments are ignored.",
        title="[bold blue]MERGE REPORTS[/bold blue]",
        border_style="blue",
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
        "[bold white]-h/--help[/bold white]    Show this help message and exit\n"
        "[bold white]-v/--version[/bold white] Show program's version number and exit",
        title="[bold cyan]OTHER OPTIONS[/bold cyan]",
        border_style="cyan",
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


###################### Exclude directory ########################################################
if os.path.exists(args.outdir):
    error_footer = Panel(
        f"[bold white]ERROR:[/bold white] Output directory '{args.outdir}' already exists. Please choose a different directory or remove it.",
        title="[bold red]Output Error[/bold red]",
        border_style="red",
        box=box.ROUNDED
    )
    console.print(error_footer)
    sys.exit(1)
#################################################################################################





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

# Check BLASTx and BLASTn arguments
if not args.diamond_blastx:
    #parser.error("--diamond_blastx must be specified.")
    error_footer2 = Panel(
        f"[bold white]ERROR:[/bold white] --diamond_blastx must be specified.",
        title="[bold red]BLASTx Argument Required[/bold red]",
        border_style="red",
        box=box.ROUNDED
    )
    console.print(error_footer2)
    sys.exit(1)

if args.blastn and args.blastn_online:
    #parser.error("--blastn or --blastn_online can't be executed in the same time.")
    error_footer3 = Panel(
        f"[bold white]ERROR:[/bold white] --blastn or --blastn_online can't be executed in the same time.",
        title="[bold red]BLASTn Argument Error[/bold red]",
        border_style="red",
        box=box.ROUNDED
    )   
    console.print(error_footer3)
    sys.exit(1)

# HMM Files Check
if args.rvdb_hmm: 
    file_path = args.rvdb_hmm
    if not os.path.isfile(file_path):
        console.print(
            f"\n[bold red]ERROR:[/bold red] No RVDB HMM model file found. The file specified with {file_path} could not be located.",
            style="red"
        )
        sys.exit(1)

if args.eggnog_hmm: 
    file_path = args.eggnog_hmm
    if not os.path.isfile(file_path):
        console.print(
            f"\n[bold red]ERROR:[/bold red] No eggNOG HMM model file found. The file specified with {file_path} could not be located.",
            style="red"
        )
        sys.exit(1)

if args.vfam_hmm: 
    file_path = args.vfam_hmm
    if not os.path.isfile(file_path):
        console.print(
            f"\n[bold red]ERROR:[/bold red] No Vfam HMM model file found. The file specified with {file_path} could not be located.",
            style="red"
        )
        sys.exit(1)

if args.pfam_hmm: 
    file_path = args.pfam_hmm
    if not os.path.isfile(file_path):
        console.print(
            f"\n[bold red]ERROR:[/bold red] No Pfam HMM model file found. The file specified with {file_path} could not be located.",
            style="red"
        )
        sys.exit(1)

####################################

from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.live import Live
from rich.spinner import Spinner
from rich.table import Table
from rich import box

def main():

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
        "Processing FASTA",
        f"Finding ORFs - {args.maxORFs} maximum ORFs",
        f"Filter sequences - Viral RefSeq alignment",
        "HMMsearch RVDB - Detect viral elements" if args.rvdb_hmm else None,
        "HMMsearch Vfam - Detect viral elements" if args.vfam_hmm else None,
        "HMMsearch EggNOG - Detect viral elements" if args.eggnog_hmm else None,
        "HMMsearch Pfam - General characterization" if args.pfam_hmm else None,
        "Running BLASTx" if args.diamond_blastx else None,
        "Running BLASTn" if args.blastn else None,
        "Running Online BLASTn" if args.blastn_online else None,
        "Generating final table",
        "Generating AI Summary" if args.model_type and args.model_name else None,
        "Generating HTML report"
    ]
    steps = [step for step in steps if step is not None] 
    
    console = Console()
    console.print(Panel(Text(f"""Running ViralQuest 🔍""", style="bold green"), subtitle="Viral sequence analysis"))
    
    # Create a spinner
    spinner = Spinner("dots", style="cyan")
    
    # Initialize current_step before defining the function that uses it
    current_step = 0
    
    # Create a status table for all the steps
    def get_status_table():
        table = Table(box=box.ROUNDED, expand=True, show_header=False, border_style="bright_black")
        table.add_column("Status", style="cyan", width=3)
        table.add_column("Task", style="white")
        
        for i, step in enumerate(steps):
            if i < current_step:
                # Completed step
                table.add_row("✓", f"[green]{step}[/]")
            elif i == current_step:
                # Current step with spinner
                table.add_row(spinner, f"[cyan]{step}[/]")
            else:
                # Pending step
                table.add_row("○", f"[dim white]{step}[/]")
        
        return table
    
    # Start the live display
    with Live(get_status_table(), console=console, refresh_per_second=8) as live:
        
        # Processing FASTA
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
        current_step += 1
        live.update(get_status_table())
        
        # ORFs profile
        findorf(args.outdir)
        ORFs(args.outdir, args.maxORFs)
        current_step += 1
        live.update(get_status_table())
        
        # Filter w/ refseq sequences
        filterBLAST(args.outdir, args.viralRef, args.cpu, diamond_path=None, log_file=None)
        current_step += 1
        live.update(get_status_table())
        
        # HMM Search
        if args.pfam_hmm and not args.rvdb_hmm and not args.eggnog_hmm and not args.vfam_hmm:
            console.print("[bold red]Error:[/] --pfam_hmm requires --rvdb_hmm/--eggnog_hmm/--vfam_hmm to be specified.")
            sys.exit(1)

        if args.rvdb_hmm:
            process_RVDB(args.outdir, args.rvdb_hmm, args.cpu, score_threshold=10)
            current_step += 1
            live.update(get_status_table())
            
        if args.vfam_hmm:
            process_VFAM(args.outdir, args.vfam_hmm, args.cpu, score_threshold=10)
            current_step += 1
            live.update(get_status_table())
            
        if args.eggnog_hmm:
            process_EggNOG(args.outdir, args.eggnog_hmm, args.cpu, score_threshold=10)
            current_step += 1
            live.update(get_status_table())
            
        pfamFasta(args.outdir)
        process_Pfam(args.outdir, args.pfam_hmm, args.cpu)
        generateFasta(args.outdir)
        current_step += 1
        live.update(get_status_table())
        
        mergeFASTA(args.outdir)
        
        # BLAST
        if args.diamond_blastx:
            diamond_blastx(args.outdir, args.diamond_blastx, args.cpu)
            current_step += 1
            live.update(get_status_table())
            
        if args.blastn:
            generateFasta_blastn(args.outdir)
            blastn(args.outdir, args.blastn, args.cpu)
            current_step += 1
            live.update(get_status_table())

        if args.blastn_online:
            generateFasta_blastn(args.outdir)
            blastn_online(args.outdir, args.blastn_onlineDB, args.blastn_online)
            current_step += 1
            live.update(get_status_table())
            
        # Final table
        finalTable(args.outdir)
        exec_cleaning(args.outdir)
        current_step += 1
        live.update(get_status_table())

        # AI summary
        if args.model_type and args.model_name:
            current_step += 1
            analyze_viral_sequences(args)
        
        # HTML report #########################################################

        if args.cap3:
            cap3check = "true"
        else:
            cap3check = "false"

        input_repo = os.path.basename(args.input)
        outdir_repo = os.path.basename(args.outdir)

        if args.blastn:
            blastn_repo= os.path.basename(args.blastn)
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
        current_step += 1
        live.update(get_status_table())
    
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
BLASTn database: {args.blastn}
Diamond BLASTx database: {args.diamond_blastx}
RVDB hmm profile database: {args.rvdb_hmm}
Pfam hmm profile database: {args.pfam_hmm}
#############################################

The pipeline takes {time_end - time_start:.2f} seconds to run.
    """)


if __name__ == "__main__":
        main()