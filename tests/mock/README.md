# tests/mock — small mock dataset

Tiny, deterministic fixture for unit/integration tests of the `--reads` modules
(read QC, coverage, sequence quality). Not biologically realistic.

## Files

| File | Content |
|------|---------|
| `contigs.fasta` | 4 contigs (500–600 bp) with distinct structural features |
| `reads_R1.fastq` / `reads_R2.fastq` | 31 paired-end reads (150 bp) simulated from the contigs |
| `generate.py` | deterministic generator (seed `20260714`) — regenerates all of the above |

## Contigs (features to exercise the analysis)

| ID | Feature | Exercises |
|----|---------|-----------|
| `contig_random` | plain random sequence, no repeats | baseline / negative control |
| `contig_lowcomplexity` | `AT`-dinucleotide stretch + poly-A run | dustmasker low-complexity |
| `contig_direct_repeat` | 80 bp block duplicated downstream | self-BLASTn direct repeat + dot plot (slope +1) |
| `contig_inverted_repeat` | 80 bp block + its reverse complement | self-BLASTn inverted repeat + dot plot (slope −1) |

The reads tile ~230 bp fragments across every contig (so all four get coverage);
every 7th R2 read carries a quality dip in its 3′ half to exercise the
base-quality colouring of the coverage track.

## Regenerate

```bash
python tests/mock/generate.py
```
