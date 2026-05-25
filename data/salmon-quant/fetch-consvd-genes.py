#!/usr/bin/env python3
"""
fetch_conserved_genes.py
========================
Downloads mRNA/CDS FASTA sequences for kingdom-specific conserved housekeeping
genes from NCBI Entrez (Nucleotide + Gene databases).

Usage
-----
    pip install biopython tqdm
    python fetch_conserved_genes.py --email your@email.com --outdir ./gene_fastas

    # Fetch only specific kingdoms:
    python fetch_conserved_genes.py --email your@email.com --kingdoms mammals plants fungi

    # Limit hits per gene (faster testing):
    python fetch_conserved_genes.py --email your@email.com --max_per_gene 5

Output
------
    gene_fastas/
        mammals/
            ACTB.fasta
            GAPDH.fasta
            ...
        plants/
            ACT2.fasta
            ...
        combined/
            mammals_combined.fasta      <- ready for salmon index
            plants_combined.fasta
            all_kingdoms_combined.fasta <- full multi-kingdom decoy-aware index
"""

import os
import sys
import time
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional

# ── third-party ──────────────────────────────────────────────────────────────
try:
    from Bio import Entrez, SeqIO
    from Bio.Entrez import efetch, esearch, read as entrez_read
except ImportError:
    sys.exit("Biopython not found. Install with:  pip install biopython")

try:
    from tqdm import tqdm
except ImportError:
    # Graceful fallback if tqdm not installed
    def tqdm(iterable, **kwargs):
        return iterable


# ─────────────────────────────────────────────────────────────────────────────
#  GENE CATALOGUE
#  Each entry: gene_symbol -> (ncbi_gene_name_query, representative_taxon_hint)
#  The taxon_hint is appended to the Entrez query to bias towards well-annotated
#  reference sequences while still allowing cross-species hits.
# ─────────────────────────────────────────────────────────────────────────────

KINGDOM_GENES: Dict[str, Dict[str, dict]] = {

    "mammals": {
        "ACTB":   {"query": "beta-actin",           "taxon": "Homo sapiens",       "db": "nuccore"},
        "GAPDH":  {"query": "glyceraldehyde-3-phosphate dehydrogenase", "taxon": "Homo sapiens", "db": "nuccore"},
        "B2M":    {"query": "beta-2-microglobulin",  "taxon": "Homo sapiens",       "db": "nuccore"},
        "RPLP0":  {"query": "60S acidic ribosomal protein P0", "taxon": "Homo sapiens", "db": "nuccore"},
        "RPS18":  {"query": "40S ribosomal protein S18", "taxon": "Homo sapiens",   "db": "nuccore"},
        "EEF1A1": {"query": "eukaryotic translation elongation factor 1 alpha 1", "taxon": "Homo sapiens", "db": "nuccore"},
        "HPRT1":  {"query": "hypoxanthine phosphoribosyltransferase 1", "taxon": "Homo sapiens", "db": "nuccore"},
        "TBP":    {"query": "TATA-box binding protein", "taxon": "Homo sapiens",   "db": "nuccore"},
    },

    "plants": {
        "ACT2":   {"query": "actin 2",               "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "TUB4":   {"query": "tubulin beta-4",        "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "EF1A":   {"query": "elongation factor 1-alpha", "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "UBQ10":  {"query": "polyubiquitin 10",      "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "PP2A":   {"query": "serine/threonine protein phosphatase 2A", "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "GAPDH":  {"query": "glyceraldehyde-3-phosphate dehydrogenase", "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "CYP":    {"query": "cyclophilin",           "taxon": "Arabidopsis thaliana", "db": "nuccore"},
        "CLATHRIN": {"query": "clathrin adaptor complex",  "taxon": "Arabidopsis thaliana", "db": "nuccore"},
    },

    "fungi": {
        "ACT1":   {"query": "actin",                 "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "TEF1":   {"query": "translation elongation factor EF-1 alpha", "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "TUB2":   {"query": "beta-tubulin",          "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "GPD1":   {"query": "glycerol-3-phosphate dehydrogenase", "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "PGK1":   {"query": "phosphoglycerate kinase 1", "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "UBC":    {"query": "ubiquitin conjugating enzyme", "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
        "RDN18":  {"query": "18S ribosomal RNA",     "taxon": "Saccharomyces cerevisiae", "db": "nuccore"},
    },

    "arthropods": {
        "RpL32":      {"query": "ribosomal protein L32",   "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "EF1alpha":   {"query": "elongation factor 1-alpha", "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "GAPDH":      {"query": "glyceraldehyde 3-phosphate dehydrogenase", "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "Act5C":      {"query": "actin 5C",               "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "RpS20":      {"query": "ribosomal protein S20",  "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "alphaTub84B":{"query": "alpha-Tubulin 84B",      "taxon": "Drosophila melanogaster", "db": "nuccore"},
        "RpL13A":     {"query": "ribosomal protein L13A", "taxon": "Drosophila melanogaster", "db": "nuccore"},
    },

    "bacteria": {
        "rpoB":   {"query": "DNA-directed RNA polymerase subunit beta", "taxon": "Escherichia coli", "db": "nuccore"},
        "gyrB":   {"query": "DNA gyrase subunit B",  "taxon": "Escherichia coli",   "db": "nuccore"},
        "recA":   {"query": "recombinase A",          "taxon": "Escherichia coli",   "db": "nuccore"},
        "dnaK":   {"query": "molecular chaperone DnaK", "taxon": "Escherichia coli", "db": "nuccore"},
        "groEL":  {"query": "chaperonin GroEL",       "taxon": "Escherichia coli",   "db": "nuccore"},
        "16S":    {"query": "16S ribosomal RNA",      "taxon": "Escherichia coli",   "db": "nuccore"},
    },

    "archaea": {
        "rpoB":   {"query": "DNA-directed RNA polymerase subunit B", "taxon": "Methanocaldococcus jannaschii", "db": "nuccore"},
        "EF2":    {"query": "elongation factor 2",   "taxon": "Methanocaldococcus jannaschii", "db": "nuccore"},
        "aIF2":   {"query": "translation initiation factor 2 alpha", "taxon": "Sulfolobus solfataricus", "db": "nuccore"},
        "16S":    {"query": "16S ribosomal RNA",     "taxon": "Methanocaldococcus jannaschii", "db": "nuccore"},
        "radA":   {"query": "RadA recombinase",      "taxon": "Methanocaldococcus jannaschii", "db": "nuccore"},
    },
}


# ─────────────────────────────────────────────────────────────────────────────
#  HELPERS
# ─────────────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

RATE_LIMIT_DELAY = 0.4   # seconds between NCBI requests (stay under 3 req/s)
RETRY_DELAY      = 10    # seconds to wait after a failed request
MAX_RETRIES      = 3


def entrez_search(db: str, query: str, max_results: int = 20) -> List[str]:
    """Run esearch and return a list of UIDs."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            handle = esearch(db=db, term=query, retmax=max_results, usehistory="y")
            record = entrez_read(handle)
            handle.close()
            time.sleep(RATE_LIMIT_DELAY)
            return record["IdList"]
        except Exception as exc:
            log.warning("esearch attempt %d failed: %s", attempt, exc)
            time.sleep(RETRY_DELAY)
    log.error("esearch gave up after %d attempts for query: %s", MAX_RETRIES, query)
    return []


def entrez_fetch_fasta(db: str, ids: List[str]) -> Optional[str]:
    """Fetch sequences in FASTA format for a list of UIDs."""
    if not ids:
        return None
    id_str = ",".join(ids)
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            handle = efetch(db=db, id=id_str, rettype="fasta", retmode="text")
            fasta  = handle.read()
            handle.close()
            time.sleep(RATE_LIMIT_DELAY)
            return fasta
        except Exception as exc:
            log.warning("efetch attempt %d failed: %s", attempt, exc)
            time.sleep(RETRY_DELAY)
    log.error("efetch gave up after %d attempts for ids: %s...", MAX_RETRIES, id_str[:60])
    return None


def build_query(gene_info: dict) -> str:
    """Build a targeted Entrez query for mRNA/CDS sequences."""
    name   = gene_info["query"]
    taxon  = gene_info["taxon"]
    # Restrict to mRNA or CDS records; exclude whole-genome shotgun contigs
    return (
        f'"{name}"[Gene Name] AND "{taxon}"[Organism] '
        f'AND (mRNA[Filter] OR CDS) '
        f'NOT wgs[Filter] NOT patent[Filter]'
    )


def fetch_gene(kingdom: str, symbol: str, gene_info: dict,
               outdir: Path, max_per_gene: int) -> int:
    """Fetch sequences for a single gene and save to outdir/kingdom/symbol.fasta.
    Returns the number of sequences saved."""

    dest_dir = outdir / kingdom
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_file = dest_dir / f"{symbol}.fasta"

    if out_file.exists() and out_file.stat().st_size > 0:
        log.info("  ✓ %s/%s already exists — skipping", kingdom, symbol)
        # Count existing records
        return sum(1 for _ in SeqIO.parse(str(out_file), "fasta"))

    query = build_query(gene_info)
    log.info("  Searching: %s", query)

    ids = entrez_search(db=gene_info["db"], query=query, max_results=max_per_gene)
    if not ids:
        log.warning("  No results for %s/%s", kingdom, symbol)
        return 0

    log.info("  Found %d UIDs — fetching FASTA …", len(ids))
    fasta = entrez_fetch_fasta(db=gene_info["db"], ids=ids)
    if not fasta or not fasta.strip():
        log.warning("  Empty FASTA for %s/%s", kingdom, symbol)
        return 0

    out_file.write_text(fasta)
    n = fasta.count(">")
    log.info("  Saved %d sequences → %s", n, out_file)
    return n


def combine_fastas(kingdom: str, outdir: Path) -> Path:
    """Concatenate all per-gene FASTAs for a kingdom into one combined file."""
    combined_dir = outdir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_path = combined_dir / f"{kingdom}_combined.fasta"

    king_dir = outdir / kingdom
    if not king_dir.exists():
        return combined_path

    total = 0
    with open(combined_path, "w") as out_fh:
        for fasta_file in sorted(king_dir.glob("*.fasta")):
            content = fasta_file.read_text()
            if content.strip():
                out_fh.write(content)
                out_fh.write("\n")
                total += content.count(">")

    log.info("Combined %s → %d sequences in %s", kingdom, total, combined_path)
    return combined_path


def build_all_kingdoms_fasta(outdir: Path, kingdoms: List[str]) -> Path:
    """Merge all kingdom-combined FASTAs into a single multi-kingdom FASTA."""
    combined_dir = outdir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    all_path = combined_dir / "all_kingdoms_combined.fasta"

    total = 0
    with open(all_path, "w") as out_fh:
        for kingdom in kingdoms:
            king_combined = combined_dir / f"{kingdom}_combined.fasta"
            if king_combined.exists():
                content = king_combined.read_text()
                if content.strip():
                    out_fh.write(content)
                    out_fh.write("\n")
                    total += content.count(">")

    log.info("All-kingdoms FASTA → %d sequences in %s", total, all_path)
    return all_path


def write_salmon_index_script(outdir: Path, kingdoms: List[str]) -> Path:
    """Generate a ready-to-run bash script for salmon index + quant."""
    script_path = outdir / "run_salmon.sh"
    combined_fasta = outdir / "combined" / "all_kingdoms_combined.fasta"

    lines = [
        "#!/usr/bin/env bash",
        "# Auto-generated by fetch_conserved_genes.py",
        "# Runs salmon index on conserved housekeeping genes, then quant on all samples.",
        "",
        "set -euo pipefail",
        "",
        "FASTA=\"" + str(combined_fasta) + "\"",
        "INDEX_DIR=\"./salmon_index_housekeeping\"",
        "THREADS=8",
        "",
        "# ── 1. Build the salmon index ──────────────────────────────────────────────",
        "echo '>>> Building salmon index …'",
        "salmon index \\",
        "    -t \"${FASTA}\" \\",
        "    -i \"${INDEX_DIR}\" \\",
        "    --threads \"${THREADS}\" \\",
        "    --keepDuplicates",
        "",
        "echo 'Index built at '\"${INDEX_DIR}\"",
        "",
        "# ── 2. Quantify each sample ────────────────────────────────────────────────",
        "# Edit SAMPLES array: each entry is a sample name; reads assumed at",
        "#   ./reads/<sample>_R1.fastq.gz  ./reads/<sample>_R2.fastq.gz",
        "SAMPLES=(sample1 sample2 sample3)",
        "READS_DIR=\"./reads\"",
        "QUANT_DIR=\"./salmon_quant\"",
        "",
        "mkdir -p \"${QUANT_DIR}\"",
        "",
        "for SAMPLE in \"${SAMPLES[@]}\"; do",
        "    echo \">>> Quantifying ${SAMPLE} …\"",
        "    salmon quant \\",
        "        -i  \"${INDEX_DIR}\" \\",
        "        -l  A \\",
        "        -1  \"${READS_DIR}/${SAMPLE}_R1.fastq.gz\" \\",
        "        -2  \"${READS_DIR}/${SAMPLE}_R2.fastq.gz\" \\",
        "        -p  \"${THREADS}\" \\",
        "        --validateMappings \\",
        "        --gcBias \\",
        "        --seqBias \\",
        "        -o  \"${QUANT_DIR}/${SAMPLE}\"",
        "done",
        "",
        "echo '>>> All quantifications complete.'",
        "echo 'Results are in '\"${QUANT_DIR}\"",
    ]
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)
    log.info("Salmon run script written → %s", script_path)
    return script_path


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Fetch conserved housekeeping genes per kingdom from NCBI Entrez."
    )
    p.add_argument(
        "--email", required=True,
        help="Your e-mail address (required by NCBI Entrez policy)."
    )
    p.add_argument(
        "--outdir", default="./gene_fastas",
        help="Root output directory (default: ./gene_fastas)."
    )
    p.add_argument(
        "--kingdoms", nargs="+",
        choices=list(KINGDOM_GENES.keys()),
        default=list(KINGDOM_GENES.keys()),
        help="Kingdoms to fetch (default: all)."
    )
    p.add_argument(
        "--max_per_gene", type=int, default=20,
        help="Max sequences to download per gene symbol (default: 20)."
    )
    p.add_argument(
        "--api_key", default=None,
        help="NCBI API key (optional but raises rate limit from 3 to 10 req/s)."
    )
    return p.parse_args()


def main():
    args = parse_args()

    # Configure Entrez
    Entrez.email   = args.email
    Entrez.tool    = "fetch_conserved_genes"
    if args.api_key:
        Entrez.api_key = args.api_key
        global RATE_LIMIT_DELAY
        RATE_LIMIT_DELAY = 0.12   # 10 req/s with API key

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, Dict[str, int]] = {}

    for kingdom in args.kingdoms:
        log.info("━━━ Kingdom: %s ━━━", kingdom.upper())
        genes = KINGDOM_GENES[kingdom]
        summary[kingdom] = {}

        for symbol, gene_info in tqdm(genes.items(), desc=kingdom):
            n = fetch_gene(
                kingdom=kingdom,
                symbol=symbol,
                gene_info=gene_info,
                outdir=outdir,
                max_per_gene=args.max_per_gene,
            )
            summary[kingdom][symbol] = n

        combine_fastas(kingdom, outdir)

    build_all_kingdoms_fasta(outdir, args.kingdoms)
    write_salmon_index_script(outdir, args.kingdoms)

    # ── Print summary table ───────────────────────────────────────────────────
    print("\n" + "═" * 55)
    print(f"{'KINGDOM':<15} {'GENE':<15} {'SEQS':>6}")
    print("═" * 55)
    total_seqs = 0
    for kingdom, genes in summary.items():
        for symbol, n in genes.items():
            print(f"{kingdom:<15} {symbol:<15} {n:>6}")
            total_seqs += n
        print("─" * 55)
    print(f"{'TOTAL':<31} {total_seqs:>6}")
    print("═" * 55)
    print(f"\nAll FASTAs saved to:  {outdir}/")
    print(f"Combined FASTAs:      {outdir}/combined/")
    print(f"Salmon run script:    {outdir}/run_salmon.sh")


if __name__ == "__main__":
    main()
