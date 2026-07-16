# Salmon Normalization — Data-Driven HK Gene Selection

## Context

The current `SalmonQuantPipeline` uses pre-built kingdom FASTA files
(`data/salmon-quant/*.fasta`) and a fixed set of 22 Pfam HK domain profiles to
identify housekeeping gene anchors for normalization. The weakness: if those
pre-built genes are not assembled or expressed in the sample, the normalization
anchors are phantom references.

## The Idea

Replace the static HK gene lists with a **data-driven, two-pass approach** that
selects normalization candidates directly from the Salmon `quant.sf` output,
then characterizes them with BLASTn (online or local) before accepting them as
anchors. This makes normalization adaptive to whatever is actually expressed in
the sample.

## Implementation Plan

### Pass 1 — Stratified candidate selection from `quant.sf`

After `SalmonQuantRunner` produces `quant.sf`, parse all entries by **TPM**
(not NumReads — TPM already corrects for transcript length and sequencing depth).

Select candidates in two groups:

1. **Median cluster** — 10 sequences closest to the median TPM value among
   non-viral entries. These are the strongest normalization anchors: constitutively
   expressed, low noise, unlikely to be condition-specific.

2. **Upper-moderate cluster** — 10 sequences in the 75th–90th TPM percentile.
   Likely structural or metabolic genes (ribosomal proteins, translation factors).
   Useful secondary anchors.

Exclusion rules before selection:
- Drop any sequence already tagged `VQ_VIRAL_` (obvious — never use viral
  transcripts as normalization references).
- Drop the bottom 10% of TPM (noise floor — too few reads, high Poisson variance).
- Drop the top 1% of TPM (extreme outliers — may be viral or condition-specific
  highly induced genes).

Do **not** select the bottom 10 expressed sequences. Low expression = low read
counts = unreliable normalization. These are only useful as contamination or
quality markers, not as anchors.

### Pass 2 — BLASTn characterization and filtering

BLAST the selected ~20 candidates against NCBI `nt` (online via `-remote`) or a
local `nt` database. Use the existing `blastn` binary already required by the
pipeline.

Parse the top hit per candidate and classify:
- **Keep**: ribosomal protein, translation factor, metabolic enzyme, cytoskeletal
  gene, chaperone — i.e., core cellular functions present in all conditions.
- **Discard**: viral gene, transposable element, stress-response gene, unknown
  with no hit, or any sequence where the top hit e-value is above threshold.

Survivors become the **sample-derived normalization anchors** for that run,
replacing the static kingdom FASTA approach.

### LLM integration (optional, if user has an LLM configured)

After BLASTn, pass the top 3 hits per candidate to the LLM with a prompt like:

> "Given these BLASTn hits for a transcript in a [kingdom] metatranscriptomic
> sample, is this transcript a reliable housekeeping gene for normalization?
> Explain why or why not."

The LLM synthesizes the BLAST annotation and flags edge cases (e.g., a gene that
is technically cellular but known to be stress-regulated). Results are included
in the report.

## Key Design Decisions

| Decision | Choice | Reason |
|---|---|---|
| Expression metric | TPM | Length- and depth-normalized; comparable across transcripts |
| Anchor tier | Median cluster (primary) | Most stable, least condition-specific |
| Anchor tier | Upper-moderate cluster (secondary) | Likely structural/metabolic genes |
| Excluded tier | Bottom 10 expressed | Too noisy; unreliable as references |
| Excluded tier | Top 1% expressed | Outliers; risk of viral or induced genes |
| Characterization | BLASTn + annotation parsing | Confirms biological role of each candidate |
| LLM role | Interpret BLAST hits, flag edge cases | Adds biological reasoning layer |

## Where to Implement

- New class `AdaptiveHkSelector` in `viralquest/salmon_quant.py`
- Called after `QuantsfParser.parse()` in both `_run_reference` and `_run_de_novo`
- Needs access to: the `quant.sf` path, the assembled FASTA (to extract sequences
  for BLASTn), and optionally the LLM client
- Results added to `SalmonQuantReport` as a new field (e.g. `adaptive_hk_anchors`)

## Caveats

- **De-novo pathway**: assembled contigs are often partial/fragmented. Require
  minimum query coverage ≥ 60% and stricter e-value (e.g. 1e-10) before trusting
  BLASTn annotation on contig sequences.
- **Online BLASTn**: slow for 20 sequences, but acceptable as a post-processing
  step. Local `nt` database is preferred if available — add `--db-dir` support.
- **Small samples**: if fewer than 30 non-viral entries exist in `quant.sf`,
  fall back to the current static HK approach rather than selecting from noise.
- **Single-sample runs**: this approach cannot measure stability across conditions.
  The median/upper-moderate selection is a proxy for stability, not a guarantee.
  True cross-condition stability analysis requires multiple samples.
