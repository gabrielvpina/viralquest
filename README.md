<div align="center">
  <img src="viralquest/components/assets/logo-bg-text.png" width="380">
  <br><br>
  <strong>A pipeline for viral diversity analysis from assembled contigs</strong>
  <br><br>
  <img src="https://img.shields.io/badge/version-3.0.0--beta-orange">
  <img src="https://img.shields.io/badge/status-BETA-red">
  <img src="https://img.shields.io/badge/platform-Linux64-8A2BE2">
  <img src="https://img.shields.io/badge/python-3.12-blue">
</div>

---

> **BETA / TEST VERSION** — This is a development release. APIs, arguments, and output formats may change without notice. Not recommended for production use.

## Overview

ViralQuest v3 detects and characterizes viral sequences from assembled metagenomics or transcriptomics contigs. The pipeline integrates:

- **Diamond BLASTx** against a curated viral RefSeq database (filter) and NCBI NR (confirmation)
- **HMM profiling** with RVDB, Vfam, and EggNOG (viral detection) + Pfam (functional annotation)
- **BLASTn** — local or online via NCBI qblast
- **Taxonomy annotation** from NCBI viral taxonomy + ICTV
- **Sequence clustering** by species
- **Salmon quantification** with host transcriptome integration (optional)
- **LLM scoring** via Ollama (local) or OpenAI / Anthropic / Google APIs (optional)
- **Self-contained HTML report** with interactive genome viewer, cluster analysis, taxonomy, and Salmon plots

Paper: [https://link.springer.com/article/10.1186/s12859-026-06391-6](https://link.springer.com/article/10.1186/s12859-026-06391-6)

---

## Installation

### 1. Install pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Clone and set up the environment

```bash
git clone https://github.com/gabrielvpina/viralquest.git
cd viralquest
pixi install
pip install -e .
```

### 3. Install bioinformatics tools and databases

```bash
viralquest-setup
```

This command handles:
- Bioinformatics tools (`diamond`, `blast`, `salmon`, `cap3`) via pixi/bioconda
- Reference databases from Zenodo (HMM profiles + viralDB.dmnd)

To re-download or update databases only:

```bash
viralquest-download           # skips files already present
viralquest-download --force   # overwrites all
```

---

## External databases (user-provided)

These are optional but required for full pipeline runs.

| Database | Purpose | Approx. size |
|---|---|---|
| `nr.dmnd` | NR BLASTx confirmation | ~346 GB |
| BLAST nt | BLASTn local search | ~1 TB |

Build NR Diamond database:

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
diamond makedb --in nr --db nr.dmnd
```

BLASTn can also run online (`--blastn-online`), and the NR step can be omitted entirely.

---

## Usage

Minimal run:

```bash
viralquest \
  -in   SAMPLE.fasta \
  -out  SAMPLE_output \
  -cpu  8
```

Full run with NR, Salmon, and local LLM:

```bash
viralquest \
  -in            SAMPLE.fasta \
  -out           SAMPLE_output \
  -nr            /path/to/nr.dmnd \
  --blastn-online your@email.com \
  --transcriptome host_transcriptome.fasta \
  --reads         R1.fastq R2.fastq \
  --model-type    ollama \
  --model-name    qwen3:4b \
  --llm-tokens    low \
  -cpu            8
```

### Arguments

| Argument | Description |
|---|---|
| `-in` | Input FASTA (assembled contigs) |
| `-out` | Output directory |
| `-nr` | NR Diamond database (optional) |
| `-n` / `--blastn-local` | Local BLAST nt database path |
| `--blastn-online EMAIL` | NCBI qblast instead of local BLASTn |
| `-cpu` | Threads (default: 2) |
| `--cap3` | Run CAP3 assembly before the pipeline |
| `--min-identity` | Clustering identity threshold (default: 70%) |
| `--force` | Export all sequences ignoring confirmation filters |
| `--transcriptome` | Host transcriptome FASTA for Salmon |
| `--reads` | RNA-seq reads for Salmon quantification |
| `--model-type` | LLM provider: `ollama`, `openai`, `anthropic`, `google` |
| `--model-name` | Model name (e.g. `qwen3:4b`, `gemini-2.0-flash`) |
| `--llm-tokens` | Prompt mode: `high` (detailed) or `low` (compact) |
| `--api-key` | API key for cloud LLM providers |

---

## Output

```
SAMPLE_output/
├── diamond/
│   ├── refseq.tsv               # Diamond BLASTx vs viral RefSeq
│   └── nr.tsv                   # Diamond BLASTx vs NR
├── hmm/
│   ├── RVDB.tsv
│   ├── Vfam.tsv
│   ├── EggNOG.tsv
│   └── Pfam.tsv
├── blastn/
│   └── blastn.tsv
├── SAMPLE_viral_contigs.fasta   # confirmed viral sequences (FASTA)
├── SAMPLE_viralquest.json       # full structured report (JSON)
└── SAMPLE_viralquest.html       # interactive HTML report
```

The HTML report is fully self-contained (no external dependencies) and includes:
- General statistics dashboard
- Cluster analysis with alignment view
- Sequence viewer with genome map, ORFs, and HMM domains
- Taxonomy tree
- Salmon quantification plots (if run)

---

## LLM scoring

Runs only on NR-confirmed sequences. Minimum recommended model for Ollama: `qwen3:4b`.

```bash
# Local
--model-type ollama --model-name qwen3:4b --llm-tokens low

# API
--model-type google --model-name gemini-2.0-flash --api-key YOUR_KEY --llm-tokens low
```

Setup guide: [wiki/Setup-AI-Summary-resource](https://github.com/gabrielvpina/viralquest/wiki/Setup-AI-Summary-resource)
