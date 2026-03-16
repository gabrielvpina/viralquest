<br>

<div align="center">

<img src="https://github.com/gabrielvpina/viralquest/blob/main/misc/headerLogo.png?raw=true" width="430" height="140">
  
  <p align="center">
    <strong>A pipeline for viral diversity analysis</strong>
    <br>
    <br>
      <a href="https://pypi.org/project/viralquest/">
        <img alt="Static Badge" src="https://img.shields.io/badge/ViralQuest-v2.6.22-COLOR2%3Fcolor%3DCOLOR1">
        <img alt="Static Badge" src="https://img.shields.io/badge/Linux64-8A2BE2">
      </a>
  </p>
</div>


<p align="center">
  <a href="#setup">
    <img src="https://img.shields.io/badge/Setup-informational" alt="Setup">
  </a>
  <a href="#install-databases">
    <img src="https://img.shields.io/badge/Install_Databases-informational" alt="Install Databases">
  </a>
  <a href="#viral-hmm-models">
    <img src="https://img.shields.io/badge/Viral_HMM_Models-informational" alt="Viral HMM Models">
  </a>
  <a href="#install-pfam-model">
    <img src="https://img.shields.io/badge/Install_Pfam_Model-informational" alt="Install Pfam Model">
  </a>
  <a href="#ai-summary">
    <img src="https://img.shields.io/badge/AI_Summary-informational" alt="AI Summary">
  </a>
  <a href="#usage">
    <img src="https://img.shields.io/badge/Usage-informational" alt="Usage">
  </a>
  <a href="#output-files">
    <img src="https://img.shields.io/badge/Output_Files-informational" alt="Output Files">
  </a>
</p>



## Introduction
ViralQuest is a Python-based bioinformatics pipeline designed to detect, identify, and characterize viral sequences from assembled contig datasets. It streamlines the analysis of metagenomic or transcriptomic data by integrating multiple steps—such as sequence alignment, taxonomic classification, and annotation—into a cohesive and automated workflow. ViralQuest is particularly useful for virome studies, enabling researchers to uncover viral diversity, assess potential host-virus interactions, and explore the ecological or clinical significance of detected viruses.

Link to the paper: [https://link.springer.com/article/10.1186/s12859-026-06391-6](https://link.springer.com/article/10.1186/s12859-026-06391-6)


<img src="https://github.com/gabrielvpina/viralquest/blob/main/misc/figure1.png?raw=true" width="850" height="550">

### HTML Output
[Example of HTML Viral Report Output (Click Here)](https://aqua-cristi-28.tiiny.site)
> ⚠️ **Warning:** The HTML file may have some bugs in resolutions below 1920x1080p.
<img src="https://github.com/gabrielvpina/viralquest/blob/main/misc/screenshot_vq_COV.png?raw=true" width="850" height="550">

## Setup

### Install via PyPI (Recommended)

Use pip to install the latest stable version of ViralQuest for `linux64` machines.
```
pip install viralquest
```

### Install via Docker

```
# Clone the repository from GitHub:
git clone https://github.com/gabrielvpina/viralquest.git
cd viralquest

# Build the Dockerfile:
docker build -t viralquest .
```

## Install Databases

### RefSeq Viral release
The RefSeq viral release is a curated collection of viral genome and protein sequences provided by the NCBI Reference Sequence (RefSeq) database. It includes high-quality, non-redundant, and well-annotated reference sequences for viruses, maintained and updated regularly by NCBI. The required file is `viral.1.protein.faa.gz`, download via [this link](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz).
- Convert the fasta file to a Diamond Database (.dmnd):
```
diamond makedb --in viral.1.protein.faa --db viralDB.dmnd
```

The final `viralDB.dmnd` file has a size of **~219 MB**.

### BLAST nr/nt Databases
The BLAST nr (non-redundant protein) and nt (nucleotide) databases are essential resources for viral identification. The nt database is useful for identifying viral genomes or transcripts using nucleotide similarity, while nr is especially powerful for detecting and annotating viral proteins, even in divergent or novel viruses, through translated searches like blastx.
Download the nr/nt databases in fasta format via [this link](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
### nr database
1) The file `nr.gz` is the nr database in FASTA
```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
```
2) Decompress the file with `gunzip nr.gz` command.

3) Convert the fasta file to a Diamond Database (.dmnd):
```
diamond makedb --in nr --db nr.dmnd
```
> ⚠️ **Warning:** Check the version of diamond, make sure that is the same version or higher then the used to build the RefSeq Viral Release `.dmnd` file.

The final `nr.dmnd` file has a size of **346 GB**.

### nt database (optional)
1) The `nt.gz` file correspond to nt.fasta
```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz 
```
2) Decompress the file with `gunzip nt.gz` command.

## Viral HMM Models
### Important note
Hidden Markov Model (HMM) models are essential for identifying divergent viral sequences and refining sequence selection.

For this task, three models are available:

- RVDB (Reference Viral DataBase) Protein
- Vfam
- eggNOG

At least one of these models is necessary to run the pipeline. However, it's recommended to use all three concurrently.

### Vfam HMM
The VFam HMM models are profile Hidden Markov Models (HMMs) specifically designed for the identification of viral proteins. 

### eggNOG Viral HMM
The eggNOG viral OGs HMM models are part of the eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) resource and are designed to identify and annotate viral genes and proteins based on orthologous groups (OGs).

### RVDB Viral HMM
The Reference Viral Database (RVDB) is a curated collection of viral sequences, and its protein HMM models—RVDB-prot and RVDB-prot-HMM—are designed to enhance the detection and annotation of viral proteins.

### Pfam Model
Pfam is a widely used database of protein families, each represented by a profile Hidden Markov Model (HMM). These models are built from curated multiple sequence alignments and represent conserved domains or full-length protein families. Download the **version 37.2**.

## Download all models

It's possible to download all models via this script:

```
#!/bin/bash

# create new dir
mkdir vq-hmms
cd vq-hmms

# download models
wget -O EggNOG-4.5.hmm.xz https://zenodo.org/records/18715455/files/EggNOG-4.5.hmm.xz?download=1
wget -O U-RVDBv29.0-prot.hmm.xz https://zenodo.org/records/18715455/files/U-RVDBv29.0-prot.hmm.xz?download=1
wget -O Vfam-228.hmm.xz https://zenodo.org/records/18715455/files/Vfam-228.hmm.xz?download=1
wget -O Pfam-A.hmm.xz https://zenodo.org/records/18715455/files/Pfam-A.hmm.xz?download=1

# decompress models
unxz -v *xz
```

Save the script as `getModels.sh` and execute the command `chmod +x getModels.sh`. This script will create a new directory called `vq-hmms` with all hmms required. 

## AI Summary
You can use either a local LLM (via Ollama) or an API key to process and integrate viral data — such as BLAST results and HMM characterizations — with the internal ViralQuest database, which includes viral family information from ICTV (International Committee on Taxonomy of Viruses) and ViralZone. This database contains information on over 200 viral families, including details such as host range, geographic distribution, viral vectors, and more. The LLM can summarize this information to provide a broader and more insightful perspective on the viral data.

### Local LLM (via Ollama)
You can run a local LLM on your machine using Ollama. However, it is important to select a model that is well-suited for processing the data. In our tests, the smallest model that provided acceptable performance was `qwen3:4b`. Therefore, we recommend using this model as a minimum requirement for running this type of analysis.

### LLM Assistance via API
ViralQuest supports API-based LLMs from `Google`, `OpenAI`, and `Anthropic`, corresponding to the Gemini, ChatGPT, and Claude models, respectively. Please review the usage terms of each service, as a high number of requests in a short period (e.g., 3 to 15 requests per minute, depending on the number of viral sequences) may be subject to rate limits or usage restrictions.

### LLM in ViralQuest
The arguments available to use local or API LLMs are:
```
--model-type 
    Type of model to use for analysis (ollama, openai, anthropic, google).
--model-name
    Name of the model (e.g., "qwen3:4b" for ollama, "gpt-3.5-turbo" for OpenAI).
--api-key
    API key for cloud models (required for OpenAI, Anthropic, Google).
```
This is a use of the arguments with a **Local LLM (Ollama)**:
```
--model-type ollama --model-name "qwen3:8b"
```
Now using an **API key**:
```
--model-type google --model-name "gemini-2.0-flash" --api-key "12345-My-API-Key_HERE67890"
```

A tutorial to install a local LLM via ollama or Google Gemini free API is available in the [wiki](https://github.com/gabrielvpina/viralquest/wiki/Setup-AI-Summary-resource) page.

## Usage
### Query example
This is a structure of viralquest query (without AI summary resource):
```
viralquest -in SAMPLE.fasta \
-ref viral/release/viralDB.dmnd \
--blastn_online yourNCBI@email.com \
--diamond_blastx path/to/nr/diamond/database/nr.dmnd \
-rvdb /path/to/RVDB/hmm/U-RVDBv29.0-prot.hmm \
-eggnog /path/to/eggNOG/hmm/eggNOG.hmm \
-vfam /path/to/Vfam/hmm/Vfam228.hmm \
-pfam /path/to/Pfam/hmm/Pfam-A.hmm \
-cpu 4 -maxORFs 4 \
-out SAMPLE
```
> ⚠️ **Warning:** Check the version of Diamond aligner with `diamond --version` to ensure that the databases use the same version of the diamond blastx executable. The argument `dmnd_path` can be used to select a specific version of a diamond binary to be used in the pipeline.


## Output Files
This is the output directory structure:
```
INPUT: SAMPLE.fasta

OUTPUT_sample/
├── fasta-files
│   ├── SAMPLE_all_ORFs.fasta
│   ├── SAMPLE_biggest_ORFs.fasta
│   ├── SAMPLE_filtered.fasta
│   ├── SAMPLE_orig.fasta
│   ├── SAMPLE_pfam_ORFs.fasta
│   ├── SAMPLE_viralHMM.fasta
│   ├── SAMPLE_viralSeq.fasta
│   └── SAMPLE_vq.fasta
├── hit_tables
│   ├── SAMPLE_all-BLAST.csv
│   ├── SAMPLE_blastn.tsv
│   ├── SAMPLE_blastx.tsv
│   ├── SAMPLE_EggNOG.csv
│   ├── SAMPLE_hmm.csv
│   └── SAMPLE_ref.csv
├── SAMPLE_bestSeqs.json      # JSON with BLAST, HMM and ORFs information
├── SAMPLE.log                # Some parameters used in the execution of the pipeline
├── SAMPLE_viral-BLAST.csv    # BLAST result of viral sequences found
├── SAMPLE_viral.fa           # FASTA of viral sequences found
└── SAMPLE_visualization.html # HTML report
```
