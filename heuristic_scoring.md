# Heuristic Scoring (`vq_score` without an LLM)

ViralQuest produces a 0–100 confidence score (`vq_score`) and a classification
(`viral-known` / `viral-unknown` / `non-viral`) for every confirmed viral
sequence. Two independent engines can produce this metric:

| Engine | Module | Field on `NucSequence` | Requirements |
|--------|--------|------------------------|--------------|
| **Heuristic** (this document) | `score_heuristic.py` | `heuristic_output` | none — runs by default on every run |
| LLM | `score_ai.py` | `llm_output` | `--model-type`, an NR database, and (cloud) an API key |

The heuristic engine is **deterministic and rule-based**: identical input always
yields the identical score. It needs no API key, no NR database, and runs
instantly, so it scores every exported sequence unconditionally. The LLM engine,
when enabled, is an *additional* layer — the two coexist and are reported
side-by-side; neither overwrites the other.

This document describes exactly how the heuristic `vq_score` is computed. Every
number quoted here is a named constant in the `CALIBRATION` block at the top of
`viralquest/score_heuristic.py` and is meant to be tuned.

---

## 1. Overview

The score is built from **three components**, each independently mapped to a
0–100 scale, then combined by a **weighted average**. Two rules shape that
average, and they are deliberately kept separate because they answer different
questions:

```
                 ┌─────────────────────────────────────────────┐
  evidence  ───► │  per-component scores  (blastn / blastx / hmm) │
                 └─────────────────────────────────────────────┘
                                   │
                 (1) renormalise   │   absent components are dropped,
                     over present  ▼   their weight redistributed
                 ┌─────────────────────────────────────────────┐
                 │            weighted average  → raw_score       │
                 └─────────────────────────────────────────────┘
                                   │
                 (2) penalty       │   strong NON-viral dominant BLASTn hit
                     after sum     ▼   damps the score multiplicatively
                 ┌─────────────────────────────────────────────┐
                 │        vq_score = round(raw_score × penalty)  │
                 └─────────────────────────────────────────────┘
```

- **Renormalisation** handles *missing* evidence (e.g. BLASTn was never run).
  An absent component is excluded and its weight is spread over the components
  that are present, so the score is never artificially capped. **Absent is not
  the same as zero.**
- **The false-positive penalty** handles *present but contradictory* evidence: a
  contig flagged viral that actually aligns almost perfectly to a non-viral
  subject (typically the host genome) in BLASTn. This is applied *after* the
  weighted sum and can override otherwise strong BLASTx/HMM signal.

---

## 2. The ramp function

Every raw metric (percent identity, coverage, bit score, hit count) is mapped to
0–100 with a clamped linear ramp:

```
ramp(x, lo, hi) =  0                       if x ≤ lo
                   100                      if x ≥ hi
                   (x − lo)/(hi − lo)·100   otherwise
```

`lo` is the noise floor (below it, no credit) and `hi` is the saturation point
(at or above it, full credit). Choosing `lo`/`hi` per metric is the main
calibration lever.

---

## 3. The three components

### 3.1 Alignment components (BLASTx and BLASTn)

Both BLAST components use the same shape: the *best hit*'s percent identity and
query coverage are each ramped, then mixed:

```
alignment = ID_FRAC · ramp(pident, id_lo, id_hi)
          + COV_FRAC · ramp(qcov,  cov_lo, cov_hi)
```

| Constant | Value | Meaning |
|----------|-------|---------|
| `ID_FRAC` | 0.6 | weight of identity within the component |
| `COV_FRAC` | 0.4 | weight of coverage within the component |
| `BX_ID_LO`, `BX_ID_HI` | 20, 90 | BLASTx identity ramp (protein-level — meaningful identity is lower) |
| `BN_ID_LO`, `BN_ID_HI` | 70, 98 | BLASTn identity ramp (nucleotide-level) |
| `COV_LO`, `COV_HI` | 10, 80 | shared query-coverage ramp |

**BLASTx** uses `best_blastx` (the NR hit if present, else the RefSeq hit). Both
pools are viral by construction — RefSeq is a viral-only database, and NR hits
are pre-filtered to viral subjects upstream (see §6) — so BLASTx always counts
as positive viral evidence when present.

**BLASTn** is different: the database is general, so a hit can be to *anything*.
The positive BLASTn component is computed **only from the best hit whose subject
title looks viral** (§5). If no BLASTn hit looks viral, the BLASTn component is
treated as **absent** (`None`) and renormalised out — it contributes no positive
score. The non-viral hits are not forgotten, though: they feed the penalty (§4).

### 3.2 HMM component

The HMM component blends a **qualitative** signal (how strong the best domain
is) with a **quantitative** signal (how many domains were found), restricted to
the viral-confirmation databases:

```
FILTER_DBS = { RVDB, Vfam, EggNOG }          # Pfam is functional, not counted here

qual  = ramp(best_filter_domain_bitscore, HMM_SCORE_LO, HMM_SCORE_HI)
quant = ramp(total_filter_hit_count,       HMM_COUNT_LO, HMM_COUNT_HI)
hmm   = HMM_QUAL_W · qual + HMM_QUANT_W · quant
```

| Constant | Value | Meaning |
|----------|-------|---------|
| `HMM_SCORE_LO`, `HMM_SCORE_HI` | 50, 150 | best FILTER-domain bit-score ramp |
| `HMM_COUNT_LO`, `HMM_COUNT_HI` | 1, 6 | FILTER-hit multiplicity ramp |
| `HMM_QUAL_W` | 0.7 | weight of the best-score signal |
| `HMM_QUANT_W` | 0.3 | weight of the hit-count signal |

The hit count comes from `NucSequence.hmm_raw_counts`, which aggregates
`Orf.raw_hmm_counts` — the number of HMM hits **above the search threshold
(score ≥ 50)** recorded *before* positional de-duplication in
`HmmResultAttacher.attach`. This makes multiplicity (e.g. RVDB matching an ORF
many times) a usable complementary signal that the de-duplicated `hmm_domains`
list would otherwise hide.

> **Calibration note.** The count is over *domain hits*, not distinct models, so
> a multi-domain ORF inflates it. Adjust `HMM_COUNT_HI` if this over-rewards.

If a sequence has no FILTER-database evidence at all (no domain and zero raw
count), the HMM component is **absent** and renormalised out.

---

## 4. Combining components

### 4.1 Base weights and renormalisation

The base weights follow the project rubric:

| Component | Base weight |
|-----------|-------------|
| `blastn` | 0.30 |
| `blastx` | 0.30 |
| `hmm` | 0.40 |

Only the **present** components participate. Their base weights are renormalised
to sum to 1, then applied:

```
present  = { components that are not None }
total_w  = Σ base_weight[k]  for k in present
applied  = base_weight[k] / total_w
raw_score = Σ applied[k] · component[k]
```

This is why a run with no BLASTn database is **not** capped at 70: the missing
30% is redistributed across BLASTx and HMM rather than counted as zero. The
weights actually applied are stored in `components.weights` for transparency.

If *no* component is present, `raw_score = 0`.

### 4.2 False-positive penalty

After the weighted sum, the **dominant** BLASTn hit (highest bit score,
regardless of subject) is examined. If it is **non-viral** and **strong**, the
whole score is multiplied by a heavy damping factor:

```
if dominant_blastn is non-viral
   and dominant_blastn.pident ≥ FP_ID
   and dominant_blastn.qcovhsp ≥ FP_COV:
       penalty = NON_VIRAL_PENALTY      # 0.2
else:
       penalty = 1.0                    # no penalty

vq_score = clamp(round(raw_score × penalty), 0, 100)
```

| Constant | Value | Meaning |
|----------|-------|---------|
| `FP_ID` | 90 | identity at/above which a non-viral hit is "strong" |
| `FP_COV` | 80 | coverage at/above which a non-viral hit is "strong" |
| `NON_VIRAL_PENALTY` | 0.2 | multiplicative damping when the penalty fires |

A `×0.2` factor turns a score of 56 into 11 — strong enough to override
otherwise convincing BLASTx/HMM evidence and flag a likely host or contaminant
contig. The applied factor is stored in `components.penalty` (1.0 = none).

---

## 5. The viral-subject regex

Whether a BLASTn subject "looks viral" is decided by searching its title
*anywhere* (not only inside NCBI `[brackets]`) for viral / phage terms:

```
vir(us|al|idae|ales|inae|oid|ion|aceae)  |  phage  |  bacteriophage
```

Searching the whole title — rather than only bracketed organism names — is what
lets a **non-standard or general BLASTn database** still be classified
correctly. This matters in two directions:

- a **viral** match drives the positive BLASTn component and can reach
  `viral-known`;
- a **non-viral** dominant hit at high identity/coverage triggers the penalty.

---

## 6. Classification

The classification is **presence-of-evidence logic**, computed independently of
the numeric bands:

```
if penalty < 1.0:                          → "non-viral"   # false positive flagged
elif no component is present:              → "non-viral"   # no viral evidence at all
elif (viral BLASTn  with pident ≥ 90 and qcov ≥ 70)
  or (best BLASTx   with pident ≥ 90 and qcov ≥ 70):
                                           → "viral-known"
else:                                      → "viral-unknown"
```

| Constant | Value | Meaning |
|----------|-------|---------|
| `KNOWN_ID` | 90 | identity required for species-level ("known") evidence |
| `KNOWN_COV` | 70 | coverage required for species-level ("known") evidence |

> **Design note / calibration choice.** Classification reflects *evidence
> presence*, while `vq_score` reflects *evidence strength*. They can therefore
> disagree — e.g. a sequence with weak-but-present viral signal can be
> `viral-unknown` with a low score (see Example D). If you prefer classification
> tied to the numeric bands (as in the LLM rubric), that change lives entirely
> in `_classify`.

`blastn_species` is a best-effort, deterministic extraction of the species name
from the best viral BLASTn subject title (leading accession dropped, cut at the
first descriptor keyword / comma / parenthesis). It is intentionally minimal.

---

## 7. Worked examples

These four cases are taken directly from the scorer's smoke test; the arithmetic
matches the implementation exactly.

### A — strong viral evidence on all fronts → **91, viral-known**

| | identity | coverage | component |
|--|--|--|--|
| BLASTx (viral) | 95% | 90% | `0.6·100 + 0.4·100` = **100** |
| BLASTn (viral) | 96% | 90% | `0.6·92.9 + 0.4·100` = **95.7** |
| HMM | RVDB score 140, 4 raw hits | | `0.7·90 + 0.3·60` = **81** |

All present → weights stay 0.3 / 0.3 / 0.4:
`0.3·95.7 + 0.3·100 + 0.4·81 = 91.1` → **91**. No penalty. Both BLAST channels
clear 90/70 → `viral-known`.

### B — HMM-only, BLASTn not run → **53, viral-unknown** (not capped at 70)

| | value | component |
|--|--|--|
| BLASTx | 45% id / 60% cov | **50** |
| BLASTn | *absent* | — |
| HMM | Vfam score 120, 2 raw hits | **55** |

BLASTn absent → renormalise over {blastx 0.3, hmm 0.4}, `total_w = 0.7`:
applied = 0.43 / 0.57. `0.43·50 + 0.57·55 = 52.9` → **53**. This demonstrates
renormalisation: the missing BLASTn weight does not cap the score at 70.

### C — false positive: dominant BLASTn hit is the host genome → **11, non-viral**

| | value | component |
|--|--|--|
| BLASTx | 40% id / 50% cov | **40** |
| BLASTn (viral) | *none — only hit is non-viral* | absent |
| HMM | RVDB score 130, 3 raw hits | **68** |

Weighted sum over {blastx, hmm}: `0.43·40 + 0.57·68 = 56.0`. The dominant BLASTn
hit is *Homo sapiens chromosome 7* at 99% id / 95% cov → non-viral and strong →
`penalty = 0.2`. `56.0 × 0.2 = 11.2` → **11**, and the penalty forces
`non-viral`. Note how the penalty overrode the strong HMM signal.

### D — weak evidence everywhere → **6, viral-unknown**

BLASTx 25% id / 20% cov → component 10; HMM EggNOG score 55, 1 raw hit →
component 3.5. Renormalised sum ≈ 6, no penalty. Viral evidence is present but
far below the known threshold → `viral-unknown` with a low score (see the
classification design note in §6).

---

## 8. Calibration constants — quick reference

All in the `CALIBRATION` block of `viralquest/score_heuristic.py`:

```python
WEIGHTS            = {"blastn": 0.30, "blastx": 0.30, "hmm": 0.40}
ID_FRAC, COV_FRAC  = 0.6, 0.4
BX_ID_LO, BX_ID_HI = 20.0, 90.0     # BLASTx identity ramp (protein)
BN_ID_LO, BN_ID_HI = 70.0, 98.0     # BLASTn identity ramp (nucleotide)
COV_LO,   COV_HI   = 10.0, 80.0     # shared coverage ramp
HMM_SCORE_LO, HMM_SCORE_HI = 50.0, 150.0
HMM_COUNT_LO, HMM_COUNT_HI = 1.0,  6.0
HMM_QUAL_W,   HMM_QUANT_W  = 0.7,  0.3
KNOWN_ID, KNOWN_COV        = 90.0, 70.0
FP_ID,    FP_COV           = 90.0, 80.0
NON_VIRAL_PENALTY          = 0.2
```

**Suggested calibration workflow.** The LLM engine writes a `vq_score` to the
same exported sequences whenever it runs. After a pipeline run that used both
engines, correlate `heuristic_output.vq_score` against `llm_output.vq_score`
across sequences and adjust the ramps/weights so the two broadly agree on the
clear cases, while keeping the heuristic engine's deterministic edge on the
false-positive cases the LLM may miss.

---

## 9. Output schema

Serialised under each sequence as `heuristic_output`:

```json
{
  "seq_id": "contig_42",
  "vq_score": 91,
  "classification": "viral-known",
  "blastn_species": "Tomato bushy stunt virus",
  "analysis": "BLASTx best hit 95.0% identity / 90.0% coverage ...",
  "components": {
    "blastn": 95.71,
    "blastx": 100.0,
    "hmm": 81.0,
    "weights": { "blastn": 0.3, "blastx": 0.3, "hmm": 0.4 },
    "penalty": 1.0
  }
}
```

A `component` value of `null` means that component was absent and excluded from
the average. `weights` are the weights *actually applied* after renormalisation.
`penalty` is the multiplicative factor (`1.0` = none).

Run-level aggregates appear under `pipeline_stats.heuristic` (count, per-class
totals, average score), mirroring `pipeline_stats.llm`.
```
