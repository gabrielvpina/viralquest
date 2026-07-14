#!/usr/bin/env python3
"""
generate.py — build the small mock test dataset in tests/mock/.

Produces four contigs with distinct, testable structural features and a matching
paired-end read set simulated from them (so coverage / read-QC / seq-quality
modules all have real-ish input). Deterministic (fixed seed) and intentionally
tiny — for unit/integration tests only, not biological realism.

Run from the repo root:  python tests/mock/generate.py
"""

from __future__ import annotations

import random
from pathlib import Path

SEED    = 20260714
OUT_DIR = Path(__file__).parent
RNG     = random.Random(SEED)

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def rand_seq(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def build_contigs() -> list[tuple[str, str, str]]:
    """Return [(id, description, sequence), ...]."""
    # 1) plain random contig — no repeats, normal complexity.
    c1 = rand_seq(600)

    # 2) low-complexity contig — AT dinucleotide stretch + poly-A run embedded.
    low = "AT" * 45 + "A" * 40                       # ~130 bp dustmasker target
    c2  = rand_seq(180) + low + rand_seq(190)

    # 3) direct-repeat contig — an 80 bp block duplicated downstream.
    block = rand_seq(80)
    c3    = rand_seq(120) + block + rand_seq(140) + block + rand_seq(120)

    # 4) inverted-repeat contig — a block followed later by its reverse complement.
    inv = rand_seq(80)
    c4  = rand_seq(130) + inv + rand_seq(150) + revcomp(inv) + rand_seq(110)

    return [
        ("contig_random",         "plain random sequence, no repeats",        c1),
        ("contig_lowcomplexity",  "AT-rich + poly-A low-complexity region",   c2),
        ("contig_direct_repeat",  "80 bp block duplicated (direct repeat)",   c3),
        ("contig_inverted_repeat","80 bp block + its reverse complement",     c4),
    ]


def write_fasta(contigs, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for cid, desc, seq in contigs:
            fh.write(f">{cid} {desc}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def qual_string(n: int, dip: bool) -> str:
    """Phred+33 quality string; `dip` lowers quality in the read's 3' half."""
    out = []
    for i in range(n):
        q = 38 if i < n * 0.6 else 30
        if dip and i > n * 0.5:
            q = 14
        q += RNG.randint(-2, 2)
        q = max(2, min(40, q))
        out.append(chr(q + 33))
    return "".join(out)


def simulate_reads(contigs, read_len=150, frag=230, step=45):
    """Tile paired-end fragments across each contig. Returns (r1_records, r2_records)."""
    r1, r2 = [], []
    pair_id = 0
    for cid, _desc, seq in contigs:
        for start in range(0, max(1, len(seq) - frag + 1), step):
            fragment = seq[start:start + frag]
            if len(fragment) < read_len:
                continue
            read1 = fragment[:read_len]
            read2 = revcomp(fragment[-read_len:])
            pair_id += 1
            name = f"{cid}_frag{pair_id:04d}"
            # every 7th pair gets a quality dip in R2 to exercise the QC colouring
            dip = (pair_id % 7 == 0)
            r1.append((name, read1, qual_string(len(read1), False)))
            r2.append((name, read2, qual_string(len(read2), dip)))
    return r1, r2


def write_fastq(records, path: Path, mate: int) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}/{mate}\n{seq}\n+\n{qual}\n")


def main() -> None:
    contigs = build_contigs()
    write_fasta(contigs, OUT_DIR / "contigs.fasta")
    r1, r2 = simulate_reads(contigs)
    write_fastq(r1, OUT_DIR / "reads_R1.fastq", 1)
    write_fastq(r2, OUT_DIR / "reads_R2.fastq", 2)
    print(f"contigs: {len(contigs)}  |  read pairs: {len(r1)}")
    for cid, _d, seq in contigs:
        print(f"  {cid}: {len(seq)} bp")


if __name__ == "__main__":
    main()
