import sys, os
# rich dependencies
from rich.panel import Panel
from rich.console import Console
from rich import box
from rich.align import Align

console = Console()

def trackErrors(args):
    """Check if all args are correct.
    deal with empty blastx tables.
    """

    if os.path.exists(args.outdir):
        error_footer = Panel(
            f"[bold white]ERROR:[/bold white] Output directory '{args.outdir}' already exists. Please choose a different directory or remove it.",
            title="[bold red]Output Error[/bold red]",
            border_style="red",
            width=80,
            box=box.ROUNDED
        )
        console.print(error_footer)
        sys.exit(1)

    # Check BLASTx and BLASTn arguments
    if not args.diamond_blastx:
        #parser.error("--diamond_blastx must be specified.")
        error_footer2 = Panel(
            f"[bold white]ERROR:[/bold white] --diamond_blastx must be specified.",
            title="[bold red]BLASTx Argument Required[/bold red]",
            border_style="red",
            width=80,
            box=box.ROUNDED
        )
        console.print(error_footer2)
        sys.exit(1)

    if args.blastn_local and args.blastn_online:
        #parser.error("--blastn or --blastn_online can't be executed in the same time.")
        error_footer3 = Panel(
            f"[bold white]ERROR:[/bold white] --blastn or --blastn_online can't be executed in the same time.",
            title="[bold red]BLASTn Argument Error[/bold red]",
            border_style="red",
            width=80,
            box=box.ROUNDED
        )   
        console.print(error_footer3)
        sys.exit(1)

    # Check if BLASTx file exists
    if args.diamond_blastx: 
        file_path = args.diamond_blastx
        if not os.path.isfile(file_path):
            console.print(
                f"\n[bold red]ERROR:[/bold red] No file found. The file specified with {file_path} could not be located.",
                style="red"
            )
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


    # check NCBI connection
    import requests
    from rich.live import Live

    if args.blastn_online:
        try:
            # Create initial box
            request_box = Panel(
                Align.center(f"[bold white]Checking connection with NCBI...[/bold white]"),
                border_style="blue",
                width=80,
                box=box.ROUNDED
            )
            
            with Live(request_box, console=console, refresh_per_second=10) as live:
                response = requests.head("https://www.ncbi.nlm.nih.gov/", timeout=5)
                response.raise_for_status()
                
                # Update the same box with OK status
                success_box = Panel(
                    Align.center(f"[bold white]Checking connection with NCBI... [bold green]OK[/bold green][/bold white]"),
                    border_style="green",
                    width=80,
                    box=box.ROUNDED
                )
                live.update(success_box)    
            
        except requests.exceptions.RequestException as e:
            error_footer4 = Panel(
                f"[bold white]ERROR:[/bold white] No connection to NCBI. Please check your internet connection or try again later.\n{e}",
                title="[bold red]Connection FAILED[/bold red]",
                border_style="red",
                width=80,
                box=box.ROUNDED
            )
            console.print(error_footer4)
            sys.exit(1)


    import time

    # API Configuration and Key Validation
    # Check if model_type is provided
    if hasattr(args, 'model_type') and args.model_type:
        # Check if model_name is provided
        if not hasattr(args, 'model_name') or not args.model_name:
            error_footer5 = Panel(
                f"[bold white]ERROR:[/bold white] --model_name parameter is required when using --model_type.",
                title="[bold red]Configuration FAILED[/bold red]",
                border_style="red",
                width=80,
                box=box.ROUNDED
            )
            console.print(error_footer5)
            sys.exit(1)
        
        # Test OpenAI API
        if args.model_type == "openai":
            if not hasattr(args, 'api_key') or not args.api_key:
                error_footer6 = Panel(
                    f"[bold white]ERROR:[/bold white] --api_key is required for OpenAI models.",
                    title="[bold red]OpenAI Configuration FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer6)
                sys.exit(1)
            
            try:
                # Create initial box
                request_box = Panel(
                    Align.center(f"[bold white]Checking OpenAI API key:[/bold white] {args.model_name} 🔄"),
                    border_style="blue",
                    width=80,
                    box=box.ROUNDED
                )
                
                with Live(request_box, console=console, refresh_per_second=10) as live:
                    headers = {
                        'Authorization': f'Bearer {args.api_key}',
                        'Content-Type': 'application/json'
                    }
                    response = requests.get(
                        'https://api.openai.com/v1/models',
                        headers=headers,
                        timeout=10
                    )
                    response.raise_for_status()
                    
                    # Update the same box with OK status
                    success_box = Panel(
                        Align.center(f"[bold white]Checking OpenAI API key:[/bold white] {args.model_name} [bold green]OK[/bold green]"),
                        border_style="green",
                        width=80,
                        box=box.ROUNDED
                    )
                    live.update(success_box)
                    time.sleep(0.5)
                    
            except requests.exceptions.RequestException as e:
                error_footer7 = Panel(
                    f"[bold white]ERROR:[/bold white] Invalid OpenAI API key or connection failed.\nPlease check your API key and internet connection.\n{e}",
                    title="[bold red]OpenAI API FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer7)
                sys.exit(1)
        
        # Test Anthropic API
        elif args.model_type == "anthropic":
            if not hasattr(args, 'api_key') or not args.api_key:
                error_footer8 = Panel(
                    f"[bold white]ERROR:[/bold white] --api_key is required for Anthropic models.",
                    title="[bold red]Anthropic Configuration FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer8)
                sys.exit(1)
            
            try:
                # Create initial box
                request_box = Panel(
                    Align.center(f"[bold white]Checking Anthropic API key:[/bold white] {args.model_name} 🔄"),
                    border_style="blue",
                    width=80,
                    box=box.ROUNDED
                )
                
                with Live(request_box, console=console, refresh_per_second=10) as live:
                    headers = {
                        'x-api-key': args.api_key,
                        'Content-Type': 'application/json',
                        'anthropic-version': '2023-06-01'
                    }
                    data = {
                        'model': 'claude-3-haiku-20240307',
                        'max_tokens': 1,
                        'messages': [{'role': 'user', 'content': 'Hi'}]
                    }
                    response = requests.post(
                        'https://api.anthropic.com/v1/messages',
                        headers=headers,
                        json=data,
                        timeout=10
                    )
                    response.raise_for_status()
                    
                    # Update the same box with OK status
                    success_box = Panel(
                        Align.center(f"[bold white]Checking Anthropic API key:[/bold white] {args.model_name} [bold green]OK[/bold green]"),
                        border_style="green",
                        width=80,
                        box=box.ROUNDED
                    )
                    live.update(success_box)
                    time.sleep(0.5)
                    
            except requests.exceptions.RequestException as e:
                error_footer9 = Panel(
                    f"[bold white]ERROR:[/bold white] Invalid Anthropic API key or connection failed.\nPlease check your API key and internet connection.\n{e}",
                    title="[bold red]Anthropic API FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer9)
                sys.exit(1)
        
        # Test Google API
        elif args.model_type == "google":
            if not hasattr(args, 'api_key') or not args.api_key:
                error_footer10 = Panel(
                    f"[bold white]ERROR:[/bold white] --api_key is required for Google models.",
                    title="[bold red]Google Configuration FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer10)
                sys.exit(1)
            
            try:
                # Create initial box
                request_box = Panel(
                    Align.center(f"[bold white]Checking Google API key:[/bold white] {args.model_name} 🔄"),
                    border_style="blue",
                    width=80,
                    box=box.ROUNDED
                )
                
                with Live(request_box, console=console, refresh_per_second=10) as live:
                    response = requests.get(
                        f'https://generativelanguage.googleapis.com/v1/models?key={args.api_key}',
                        timeout=10
                    )
                    response.raise_for_status()
                    
                    # Update the same box with OK status
                    success_box = Panel(
                        Align.center(f"[bold white]Checking Google API key:[/bold white] {args.model_name} [bold green]OK[/bold green]"),
                        border_style="green",
                        width=80,
                        box=box.ROUNDED
                    )
                    live.update(success_box)
                    time.sleep(0.5)
                    
            except requests.exceptions.RequestException as e:
                error_footer11 = Panel(
                    f"[bold white]ERROR:[/bold white] Invalid Google API key or connection failed.\nPlease check your API key and internet connection.\n{e}",
                    title="[bold red]Google API FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer11)
                sys.exit(1)
        
        # Ollama doesn't need API key validation, but check if it's accessible
        elif args.model_type == "ollama":
            try:
                # Create initial box
                request_box = Panel(
                    Align.center(f"[bold white]Checking Ollama model:[/bold white] {args.model_name} 🔄"),
                    border_style="blue",
                    width=80,
                    box=box.ROUNDED
                )
                
                with Live(request_box, console=console, refresh_per_second=10) as live:
                    import ollama
                    ollama.show(args.model_name)
                    
                    # Update the same box with OK status
                    success_box = Panel(
                        Align.center(f"[bold white]Checking Ollama model:[/bold white] {args.model_name} [bold green]OK[/bold green]"),
                        border_style="green",
                        width=80,
                        box=box.ROUNDED
                    )
                    live.update(success_box)
                    time.sleep(0.5)
                    
            except ImportError:
                error_footer12 = Panel(
                    f"[bold white]ERROR:[/bold white] Ollama package not installed. Please install it with: pip install ollama",
                    title="[bold red]Ollama Package FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer12)
                sys.exit(1)
            except Exception as e:
                error_footer13 = Panel(
                    f"[bold white]ERROR:[/bold white] Ollama connection failed or model '{args.model_name}' not found.\nPlease ensure Ollama is running and model is downloaded: ollama pull {args.model_name}\n{e}",
                    title="[bold red]Ollama Connection FAILED[/bold red]",
                    border_style="red",
                    width=80,
                    box=box.ROUNDED
                )
                console.print(error_footer13)
                sys.exit(1)
        
        # Unsupported model type
        else:
            error_footer14 = Panel(
                f"[bold white]ERROR:[/bold white] Unsupported model type: {args.model_type}\nSupported types: openai, anthropic, google, ollama",
                title="[bold red]Model Type FAILED[/bold red]",
                border_style="red",
                width=80,
                box=box.ROUNDED
            )
            console.print(error_footer14)
            sys.exit(1)






# this is pure LLM code 
def get_gradient_color(progress, colors):
    """
    Interpolates a color along a gradient defined by a list of RGB colors.
    """
    if progress <= 0:
        return colors[0]
    if progress >= 1:
        return colors[-1]

    num_segments = len(colors) - 1
    segment_index = int(progress * num_segments)
    
    if segment_index == num_segments:
        segment_index -= 1 

    start_color = colors[segment_index]
    end_color = colors[segment_index + 1]

    # Calculate progress within the current segment
    segment_progress = (progress * num_segments) - segment_index

    # Interpolate each RGB component
    r = int(start_color[0] + segment_progress * (end_color[0] - start_color[0]))
    g = int(start_color[1] + segment_progress * (end_color[1] - start_color[1]))
    b = int(start_color[2] + segment_progress * (end_color[2] - start_color[2]))
    
    return (r, g, b)


# The ASCII art banner provided by the user
ascii_banner = f"""
██╗   ██╗██╗██████╗  █████╗ ██╗      ██████╗ ██╗   ██╗███████╗███████╗████████╗
██║   ██║██║██╔══██╗██╔══██╗██║     ██╔═══██╗██║   ██║██╔════╝██╔════╝╚══██╔══╝
██║   ██║██║██████╔╝███████║██║     ██║   ██║██║   ██║█████╗  ███████╗   ██║   
╚██╗ ██╔╝██║██╔══██╗██╔══██║██║     ██║▄▄ ██║██║   ██║██╔══╝  ╚════██║   ██║   
 ╚████╔╝ ██║██║  ██║██║  ██║███████╗╚██████╔╝╚██████╔╝███████╗███████║   ██║   
  ╚═══╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚══▀▀═╝  ╚═════╝ ╚══════╝╚══════╝   ╚═╝  
"""

# Split the banner into individual lines
lines = ascii_banner.strip().split('\n')

# Determine the maximum width of the banner to properly apply the gradient
max_width = max(len(line) for line in lines)

# Define the gradient colors (RGB tuples)
# These colors are chosen to approximate the blue -> purple -> pink/red gradient from the image.
gradient_colors = [
    (0, 115, 230),  # Bright Blue (start)
    (100, 50, 200), # Medium Purple (middle)
    (230, 0, 450)   # Pinkish Red (end)
]

# List to store the fully colored lines
colored_banner_lines = []

# Iterate through each line and each character to apply the gradient
for line in lines:
    colored_line_chars = []
    for i, char in enumerate(line):
        # Do not color spaces; they should remain transparent to show the terminal background
        if char == ' ': 
            colored_line_chars.append(' ')
            continue
        
        # Calculate the progress for the current character across the entire banner width
        # This determines where in the gradient the character's color should be
        progress = i / (max_width - 1) if max_width > 1 else 0
        
        # get the interpolated RGB color for the current progress
        r, g, b = get_gradient_color(progress, gradient_colors)
        
        # Apply ANSI 24-bit true color escape code to the character
        # \x1b[38;2;R;G;Bm sets the foreground color
        colored_line_chars.append(f"\x1b[38;2;{r};{g};{b}m{char}")
    
    # reset color at the end of each line to avoid bleeding into subsequent terminal output
    # \x1b[0m resets all SGR (Select Graphic Rendition) attributes
    colored_banner_lines.append("".join(colored_line_chars) + "\x1b[0m") 

final_colored_ascii_banner = "\n".join(colored_banner_lines)
#final_colored_ascii_banner = print("\n".join(colored_banner_lines))

