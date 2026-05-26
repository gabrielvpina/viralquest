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

    # ── Mammals ──────────────────────────────────────────────────────────────
    # Representative taxa: human, mouse, rat, cow, pig, dog, macaque, zebrafish
    "mammals": {
        "ACTB": {
            "query": "beta-actin",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa", "Canis lupus familiaris",
                     "Macaca mulatta", "Oryctolagus cuniculus"],
            "db": "nuccore",
        },
        "GAPDH": {
            "query": "glyceraldehyde-3-phosphate dehydrogenase",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa", "Canis lupus familiaris",
                     "Ovis aries", "Equus caballus"],
            "db": "nuccore",
        },
        "B2M": {
            "query": "beta-2-microglobulin",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa"],
            "db": "nuccore",
        },
        "RPLP0": {
            "query": "60S acidic ribosomal protein P0",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa", "Macaca mulatta"],
            "db": "nuccore",
        },
        "RPS18": {
            "query": "40S ribosomal protein S18",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Canis lupus familiaris"],
            "db": "nuccore",
        },
        "EEF1A1": {
            "query": "eukaryotic translation elongation factor 1 alpha 1",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa", "Macaca mulatta",
                     "Canis lupus familiaris"],
            "db": "nuccore",
        },
        "HPRT1": {
            "query": "hypoxanthine phosphoribosyltransferase 1",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Sus scrofa", "Oryctolagus cuniculus"],
            "db": "nuccore",
        },
        "TBP": {
            "query": "TATA-box binding protein",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Macaca mulatta"],
            "db": "nuccore",
        },
        "RPL13A": {
            "query": "ribosomal protein L13a",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa"],
            "db": "nuccore",
        },
        "PPIA": {
            "query": "peptidylprolyl isomerase A",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa", "Canis lupus familiaris"],
            "db": "nuccore",
        },
        "SDHA": {
            "query": "succinate dehydrogenase complex subunit A",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus"],
            "db": "nuccore",
        },
        "UBC": {
            "query": "ubiquitin C",
            "taxa": ["Homo sapiens", "Mus musculus", "Rattus norvegicus",
                     "Bos taurus", "Sus scrofa"],
            "db": "nuccore",
        },
    },

    # ── Plants ────────────────────────────────────────────────────────────────
    # Representative taxa: Arabidopsis, rice, maize, tomato, soybean, tobacco,
    #                       wheat, potato, barley, poplar
    "plants": {
        "ACT2": {
            "query": "actin 2",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max", "Nicotiana tabacum",
                     "Triticum aestivum", "Hordeum vulgare"],
            "db": "nuccore",
        },
        "TUB4": {
            "query": "tubulin beta-4",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max"],
            "db": "nuccore",
        },
        "EF1A": {
            "query": "elongation factor 1-alpha",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max", "Nicotiana tabacum",
                     "Triticum aestivum", "Hordeum vulgare", "Phaseolus vulgaris"],
            "db": "nuccore",
        },
        "UBQ10": {
            "query": "polyubiquitin 10",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max", "Nicotiana tabacum"],
            "db": "nuccore",
        },
        "PP2A": {
            "query": "serine/threonine protein phosphatase 2A",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Nicotiana tabacum"],
            "db": "nuccore",
        },
        "GAPDH": {
            "query": "glyceraldehyde-3-phosphate dehydrogenase",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max", "Triticum aestivum",
                     "Hordeum vulgare", "Nicotiana tabacum"],
            "db": "nuccore",
        },
        "CYP": {
            "query": "cyclophilin",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max"],
            "db": "nuccore",
        },
        "CLATHRIN": {
            "query": "clathrin adaptor complex",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Solanum lycopersicum",
                     "Glycine max", "Nicotiana tabacum"],
            "db": "nuccore",
        },
        "eIF4A": {
            "query": "eukaryotic initiation factor 4A",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Zea mays",
                     "Solanum lycopersicum", "Glycine max", "Nicotiana tabacum"],
            "db": "nuccore",
        },
        "SAND": {
            "query": "SAND domain protein",
            "taxa": ["Arabidopsis thaliana", "Oryza sativa", "Solanum lycopersicum"],
            "db": "nuccore",
        },
    },

    # ── Fungi ─────────────────────────────────────────────────────────────────
    # Representative taxa: S. cerevisiae, S. pombe, C. albicans, A. niger,
    #                       A. fumigatus, N. crassa, Trichoderma reesei
    "fungi": {
        "ACT1": {
            "query": "actin",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger",
                     "Aspergillus fumigatus", "Neurospora crassa",
                     "Trichoderma reesei", "Cryptococcus neoformans"],
            "db": "nuccore",
        },
        "TEF1": {
            "query": "translation elongation factor EF-1 alpha",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger",
                     "Aspergillus fumigatus", "Neurospora crassa",
                     "Cryptococcus neoformans"],
            "db": "nuccore",
        },
        "TUB2": {
            "query": "beta-tubulin",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger",
                     "Aspergillus fumigatus", "Neurospora crassa"],
            "db": "nuccore",
        },
        "GPD1": {
            "query": "glycerol-3-phosphate dehydrogenase",
            "taxa": ["Saccharomyces cerevisiae", "Candida albicans",
                     "Aspergillus niger", "Neurospora crassa"],
            "db": "nuccore",
        },
        "PGK1": {
            "query": "phosphoglycerate kinase 1",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger",
                     "Aspergillus fumigatus", "Neurospora crassa"],
            "db": "nuccore",
        },
        "UBC": {
            "query": "ubiquitin conjugating enzyme",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger", "Neurospora crassa"],
            "db": "nuccore",
        },
        "RDN18": {
            "query": "18S ribosomal RNA",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus niger",
                     "Aspergillus fumigatus", "Neurospora crassa",
                     "Trichoderma reesei", "Cryptococcus neoformans"],
            "db": "nuccore",
        },
        "CYP2": {
            "query": "cyclophilin",
            "taxa": ["Saccharomyces cerevisiae", "Candida albicans",
                     "Aspergillus fumigatus", "Neurospora crassa"],
            "db": "nuccore",
        },
        "FKS1": {
            "query": "glucan synthase",
            "taxa": ["Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
                     "Candida albicans", "Aspergillus fumigatus"],
            "db": "nuccore",
        },
    },

    # ── Arthropods ────────────────────────────────────────────────────────────
    # Representative taxa: Drosophila, mosquitoes, honey bee, silkworm,
    #                       red flour beetle, shrimp, spider mite
    "arthropods": {
        "RpL32": {
            "query": "ribosomal protein L32",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera",
                     "Tribolium castaneum", "Bombyx mori",
                     "Culex quinquefasciatus"],
            "db": "nuccore",
        },
        "EF1alpha": {
            "query": "elongation factor 1-alpha",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera",
                     "Tribolium castaneum", "Bombyx mori",
                     "Penaeus vannamei", "Tetranychus urticae"],
            "db": "nuccore",
        },
        "GAPDH": {
            "query": "glyceraldehyde 3-phosphate dehydrogenase",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera",
                     "Tribolium castaneum", "Bombyx mori",
                     "Penaeus vannamei"],
            "db": "nuccore",
        },
        "ACT": {
            "query": "actin",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera",
                     "Tribolium castaneum", "Bombyx mori",
                     "Penaeus vannamei", "Tetranychus urticae"],
            "db": "nuccore",
        },
        "RpS20": {
            "query": "ribosomal protein S20",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Apis mellifera", "Tribolium castaneum", "Bombyx mori"],
            "db": "nuccore",
        },
        "alphaTub": {
            "query": "alpha-tubulin",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera",
                     "Tribolium castaneum", "Bombyx mori", "Penaeus vannamei"],
            "db": "nuccore",
        },
        "RpL13A": {
            "query": "ribosomal protein L13A",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Apis mellifera", "Tribolium castaneum",
                     "Bombyx mori", "Culex quinquefasciatus"],
            "db": "nuccore",
        },
        "RPS3": {
            "query": "ribosomal protein S3",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Apis mellifera", "Tribolium castaneum", "Bombyx mori"],
            "db": "nuccore",
        },
        "HSP70": {
            "query": "heat shock protein 70",
            "taxa": ["Drosophila melanogaster", "Aedes aegypti",
                     "Anopheles gambiae", "Apis mellifera", "Bombyx mori"],
            "db": "nuccore",
        },
    },

    # ── Bacteria ──────────────────────────────────────────────────────────────
    # Representative taxa: E. coli, B. subtilis, S. aureus, M. tuberculosis,
    #                       P. aeruginosa, Salmonella, Streptococcus, Listeria
    "bacteria": {
        "rpoB": {
            "query": "DNA-directed RNA polymerase subunit beta",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica",
                     "Streptococcus pneumoniae", "Listeria monocytogenes"],
            "db": "nuccore",
        },
        "gyrB": {
            "query": "DNA gyrase subunit B",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica",
                     "Streptococcus pneumoniae"],
            "db": "nuccore",
        },
        "recA": {
            "query": "recombinase A",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica"],
            "db": "nuccore",
        },
        "dnaK": {
            "query": "molecular chaperone DnaK",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica",
                     "Listeria monocytogenes"],
            "db": "nuccore",
        },
        "groEL": {
            "query": "chaperonin GroEL",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica"],
            "db": "nuccore",
        },
        "16S": {
            "query": "16S ribosomal RNA",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Staphylococcus aureus", "Mycobacterium tuberculosis",
                     "Pseudomonas aeruginosa", "Salmonella enterica",
                     "Streptococcus pneumoniae", "Listeria monocytogenes",
                     "Helicobacter pylori", "Clostridioides difficile"],
            "db": "nuccore",
        },
        "atpD": {
            "query": "ATP synthase subunit beta",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Pseudomonas aeruginosa", "Mycobacterium tuberculosis",
                     "Salmonella enterica"],
            "db": "nuccore",
        },
        "infB": {
            "query": "translation initiation factor IF-2",
            "taxa": ["Escherichia coli", "Bacillus subtilis",
                     "Pseudomonas aeruginosa", "Staphylococcus aureus"],
            "db": "nuccore",
        },
    },

    # ── Archaea ───────────────────────────────────────────────────────────────
    # Representative taxa: M. jannaschii, S. solfataricus, T. kodakarensis,
    #                       P. furiosus, H. volcanii, A. fulgidus
    "archaea": {
        "rpoB": {
            "query": "DNA-directed RNA polymerase subunit B",
            "taxa": ["Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus",
                     "Haloferax volcanii", "Archaeoglobus fulgidus",
                     "Methanobacterium thermoautotrophicum"],
            "db": "nuccore",
        },
        "EF2": {
            "query": "elongation factor 2",
            "taxa": ["Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus",
                     "Haloferax volcanii", "Archaeoglobus fulgidus"],
            "db": "nuccore",
        },
        "aIF2": {
            "query": "translation initiation factor 2 alpha",
            "taxa": ["Sulfolobus solfataricus", "Methanocaldococcus jannaschii",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus",
                     "Haloferax volcanii"],
            "db": "nuccore",
        },
        "16S": {
            "query": "16S ribosomal RNA",
            "taxa": ["Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus",
                     "Haloferax volcanii", "Archaeoglobus fulgidus",
                     "Methanobacterium thermoautotrophicum",
                     "Halobacterium salinarum"],
            "db": "nuccore",
        },
        "radA": {
            "query": "RadA recombinase",
            "taxa": ["Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus",
                     "Haloferax volcanii", "Archaeoglobus fulgidus"],
            "db": "nuccore",
        },
        "MCM": {
            "query": "minichromosome maintenance protein",
            "taxa": ["Methanocaldococcus jannaschii", "Sulfolobus solfataricus",
                     "Thermococcus kodakarensis", "Pyrococcus furiosus"],
            "db": "nuccore",
        },
        "atp-synthase": {
            "query": "archaeal ATP synthase subunit B",
            "taxa": ["Methanocaldococcus jannaschii", "Thermococcus kodakarensis",
                     "Pyrococcus furiosus", "Archaeoglobus fulgidus",
                     "Haloferax volcanii"],
            "db": "nuccore",
        },
    },

    # ── Fish (Teleosts) ───────────────────────────────────────────────────────
    # Representative taxa: zebrafish, medaka, tilapia, Atlantic salmon,
    #                       rainbow trout, common carp, fugu
    "fish": {
        "actb1": {
            "query": "beta-actin 1",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss", "Cyprinus carpio",
                     "Takifugu rubripes", "Sparus aurata"],
            "db": "nuccore",
        },
        "eef1a1": {
            "query": "eukaryotic translation elongation factor 1 alpha 1",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss", "Cyprinus carpio"],
            "db": "nuccore",
        },
        "gapdh": {
            "query": "glyceraldehyde-3-phosphate dehydrogenase",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss", "Cyprinus carpio",
                     "Sparus aurata"],
            "db": "nuccore",
        },
        "rpl13a": {
            "query": "ribosomal protein L13a",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss"],
            "db": "nuccore",
        },
        "rps18": {
            "query": "40S ribosomal protein S18",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss"],
            "db": "nuccore",
        },
        "b2m": {
            "query": "beta-2-microglobulin",
            "taxa": ["Danio rerio", "Oryzias latipes", "Oreochromis niloticus",
                     "Salmo salar", "Oncorhynchus mykiss"],
            "db": "nuccore",
        },
    },

    # ── Nematodes ─────────────────────────────────────────────────────────────
    # Representative taxa: C. elegans, C. briggsae, C. brenneri, Brugia malayi,
    #                       Ascaris suum, Haemonchus contortus
    "nematodes": {
        "act-1": {
            "query": "actin",
            "taxa": ["Caenorhabditis elegans", "Caenorhabditis briggsae",
                     "Brugia malayi", "Ascaris suum",
                     "Haemonchus contortus", "Pristionchus pacificus"],
            "db": "nuccore",
        },
        "eef-1A": {
            "query": "elongation factor 1-alpha",
            "taxa": ["Caenorhabditis elegans", "Caenorhabditis briggsae",
                     "Brugia malayi", "Ascaris suum", "Haemonchus contortus"],
            "db": "nuccore",
        },
        "gpd-2": {
            "query": "glyceraldehyde-3-phosphate dehydrogenase",
            "taxa": ["Caenorhabditis elegans", "Caenorhabditis briggsae",
                     "Brugia malayi", "Haemonchus contortus"],
            "db": "nuccore",
        },
        "rpl-1": {
            "query": "ribosomal protein L1",
            "taxa": ["Caenorhabditis elegans", "Caenorhabditis briggsae",
                     "Brugia malayi", "Ascaris suum"],
            "db": "nuccore",
        },
        "pmp-5": {
            "query": "peptidylprolyl isomerase",
            "taxa": ["Caenorhabditis elegans", "Caenorhabditis briggsae",
                     "Brugia malayi", "Haemonchus contortus"],
            "db": "nuccore",
        },
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


def build_query(symbol: str, gene_info: dict) -> str:
    name = gene_info["query"]
    taxa = gene_info["taxa"]
    taxa_clause = " OR ".join(f'"{t}"[Organism]' for t in taxa)
    
    # 1. Mantemos a busca ampla por Símbolo ou Nome
    # 2. Trocamos o (mRNA[Filter] OR CDS) por biomol_mrna[prop]
    # 3. Excluímos explicitamente DNA genômico por garantia
    return (
        f'({symbol}[Gene] OR "{name}"[Protein Name]) AND ({taxa_clause}) '
        f'AND biomol_mrna[prop] '
        f'NOT "genomic dna"[Filter] NOT wgs[Filter] NOT patent[Filter]'
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

    query = build_query(symbol, gene_info)
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
