"""
seq_quality.py — Structural quality signals for viral sequences.

Complementary to the read-coverage / base-quality track, this module flags
signals that suggest mis-assembly or sequencing artefacts, using the confirmed
viral sequences as input:

  * **dustmasker** — low-complexity regions (per-base masked runs);
  * **self-BLASTn** — internal repeats (direct and inverted), with the trivial
    full-length self-diagonal filtered out;
  * **jellyfish** — a k-mer repetitiveness score (fraction of k-mer instances
    that recur), computed per sequence;
  * a **self-similarity dot plot**, precomputed in numpy and emitted as line
    segments for rendering as a native ``<svg>`` (off-diagonal = direct repeat,
    anti-diagonal = inverted repeat).

Every external tool degrades gracefully: if a binary is missing the
corresponding signal is simply omitted, mirroring ``CoveragePipeline``.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from loguru import logger

from viralquest.biodata import DotPlot, NucSequence, SelfRepeat, SequenceQuality


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class SequenceQualityPipeline:
    """Build a ``SequenceQuality`` for each confirmed viral sequence."""

    _DOTPLOT_MAX_DIM      = 500     # sample the sequence to ~this many anchors per axis
    _DOTPLOT_MAX_SEGMENTS = 4000    # cap emitted SVG segments per sequence
    _SELF_BLAST_EVALUE    = "1e-5"

    def __init__(
        self,
        threads:        int = 4,
        kmer_size:      int = 15,   # jellyfish k
        dotplot_word:   int = 12,   # dot plot seed word
        dustmasker_bin: str = "dustmasker",
        blastn_bin:     str = "blastn",
        jellyfish_bin:  str = "jellyfish",
        cleanup:        bool = True,
    ):
        self.threads        = threads
        self.kmer_size      = kmer_size
        self.dotplot_word   = dotplot_word
        self.dustmasker_bin = dustmasker_bin
        self.blastn_bin     = blastn_bin
        self.jellyfish_bin  = jellyfish_bin
        self.cleanup        = cleanup

    # ── public ────────────────────────────────────────────────────────────────

    def run(
        self,
        viral_seqs: list[NucSequence],
        outdir:     Path,
    ) -> dict[str, SequenceQuality]:
        """
        Returns ``{seq_id: SequenceQuality}`` for every sequence.  Missing tools
        only drop their own signal; a sequence always gets an entry (possibly
        with empty repeat / low-complexity lists).
        """
        seqs = [s for s in viral_seqs if s.sequence]
        if not seqs:
            logger.warning("SequenceQualityPipeline: no viral sequences — skipping.")
            return {}

        outdir.mkdir(parents=True, exist_ok=True)
        ref_fasta = outdir / "viral_ref.fasta"
        self._write_fasta(seqs, ref_fasta)

        try:
            low_cx   = self._dustmask(ref_fasta)             # {seq_id: [[start,end], ...]}
            repeats  = self._self_blast(ref_fasta)           # {seq_id: [SelfRepeat, ...]}
            profiles = self._build(seqs, low_cx, repeats, outdir)
        finally:
            if self.cleanup:
                ref_fasta.unlink(missing_ok=True)

        self._export(profiles, outdir)
        logger.success(f"Sequence quality signals built for {len(profiles)} sequence(s).")
        return profiles

    # ── fasta ─────────────────────────────────────────────────────────────────

    @staticmethod
    def _write_fasta(seqs: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in seqs:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    # ── persistent export ───────────────────────────────────────────────────────

    @staticmethod
    def _export(profiles: "dict[str, SequenceQuality]", outdir: Path) -> None:
        """
        Persist the sequence-quality signals to disk (in-memory results also go to
        the JSON/HTML report). Writes:

          * ``<seq_id>.tsv`` — one row per structural feature (low-complexity
            region or self-repeat), for every sequence that has at least one.
            Coordinates are 1-based inclusive.
          * ``summary.tsv``  — one row per sequence with the scalar signals.

        Mirrors the per-sequence export style of ``coverage.py``.
        """
        _safe = re.compile(r"[^\w\-.]")
        _STR  = {"plus": "+", "minus": "-"}

        for prof in profiles.values():
            rows: list[tuple] = []
            # feature  q_start  q_end  s_start  s_end  strand  pct_id  length
            for a, b in prof.low_complexity_regions:
                rows.append(("low_complexity", a + 1, b + 1, "", "", "", "", b - a + 1))
            for r in prof.self_repeats:
                kind = "repeat_direct" if r.strand == "plus" else "repeat_inverted"
                rows.append((kind, r.q_start + 1, r.q_end + 1,
                             r.s_start + 1, r.s_end + 1,
                             _STR.get(r.strand, r.strand), r.pct_id, r.length))
            if not rows:
                continue
            safe = _safe.sub("_", prof.seq_id)
            with open(outdir / f"{safe}.tsv", "w", encoding="utf-8") as fh:
                fh.write("feature\tq_start\tq_end\ts_start\ts_end\tstrand\tpct_id\tlength\n")
                for row in rows:
                    fh.write("\t".join(str(x) for x in row) + "\n")

        with open(outdir / "summary.tsv", "w", encoding="utf-8") as fh:
            fh.write("seq_id\tlength\tlow_complexity_frac\tn_low_complexity\t"
                     "n_repeats\tn_direct\tn_inverted\tkmer_size\tkmer_distinct\t"
                     "kmer_total\tkmer_repeat_score\n")
            for prof in profiles.values():
                n_direct   = sum(1 for r in prof.self_repeats if r.strand == "plus")
                n_inverted = len(prof.self_repeats) - n_direct
                fh.write("\t".join(str(x) for x in (
                    prof.seq_id, prof.length, prof.low_complexity_frac,
                    len(prof.low_complexity_regions), len(prof.self_repeats),
                    n_direct, n_inverted, prof.kmer_size, prof.kmer_distinct,
                    prof.kmer_total, prof.kmer_repeat_score,
                )) + "\n")

    # ── dustmasker ──────────────────────────────────────────────────────────────

    def _dustmask(self, ref_fasta: Path) -> dict[str, list[list[int]]]:
        """
        Run ``dustmasker -outfmt interval`` and return per-sequence low-complexity
        runs as ``[[start, end], ...]`` (0-based, inclusive).
        """
        if not shutil.which(self.dustmasker_bin):
            logger.warning("SequenceQualityPipeline: dustmasker not found — low-complexity signal skipped.")
            return {}

        proc = subprocess.run(
            [self.dustmasker_bin, "-in", str(ref_fasta), "-outfmt", "interval"],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            logger.error(f"dustmasker failed: {proc.stderr.strip()}")
            return {}

        out: dict[str, list[list[int]]] = {}
        current: str | None = None
        for line in proc.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                out.setdefault(current, [])
            elif current is not None and " - " in line:
                a, b = line.split(" - ")
                out[current].append([int(a), int(b)])
        return out

    # ── self-BLASTn ─────────────────────────────────────────────────────────────

    def _self_blast(self, ref_fasta: Path) -> dict[str, list[SelfRepeat]]:
        """
        Align the sequences against themselves and keep only same-sequence HSPs
        that are not the trivial full-length self-diagonal — i.e. internal direct
        (plus) or inverted (minus) repeats.
        """
        if not shutil.which(self.blastn_bin):
            logger.warning("SequenceQualityPipeline: blastn not found — self-repeat signal skipped.")
            return {}

        fmt = "6 qseqid sseqid qstart qend sstart send pident length sstrand"
        proc = subprocess.run(
            [self.blastn_bin,
             "-query", str(ref_fasta), "-subject", str(ref_fasta),
             "-dust", "yes", "-evalue", self._SELF_BLAST_EVALUE,
             "-word_size", "11", "-outfmt", fmt],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            logger.error(f"self-blastn failed: {proc.stderr.strip()}")
            return {}

        out: dict[str, list[SelfRepeat]] = {}
        seen: set[tuple] = set()
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            qseqid, sseqid = parts[0], parts[1]
            if qseqid != sseqid:
                continue                      # cross-sequence hit — not a self-repeat
            qs, qe = int(parts[2]) - 1, int(parts[3]) - 1   # → 0-based inclusive
            ss, se = int(parts[4]) - 1, int(parts[5]) - 1
            pid    = float(parts[6])
            length = int(parts[7])
            strand = parts[8]                 # "plus" | "minus"

            # blastn reports sstart>send on the minus strand — normalise for storage.
            s_lo, s_hi = (ss, se) if ss <= se else (se, ss)

            if strand == "plus" and qs == s_lo and qe == s_hi:
                continue                      # trivial full self-diagonal

            # Symmetric duplicate (A→B and B→A) — keep one, upper triangle by q start.
            key = (min(qs, s_lo), max(qs, s_lo), strand)
            if key in seen:
                continue
            seen.add(key)

            out.setdefault(qseqid, []).append(SelfRepeat(
                q_start=qs, q_end=qe, s_start=s_lo, s_end=s_hi,
                strand=strand, pct_id=round(pid, 2), length=length,
            ))
        return out

    # ── jellyfish ───────────────────────────────────────────────────────────────

    def _jellyfish(self, seq: NucSequence, outdir: Path) -> tuple[int, int, float]:
        """
        Return ``(distinct, total, repeat_score)`` where repeat_score is the
        fraction of k-mer *instances* that recur (multiplicity > 1).  Returns
        ``(0, 0, 0.0)`` if jellyfish is unavailable or the sequence is too short.
        """
        if not shutil.which(self.jellyfish_bin):
            return (0, 0, 0.0)
        if seq.length < self.kmer_size:
            return (0, 0, 0.0)

        safe   = re.sub(r"[^\w\-.]", "_", seq.id)
        fa     = outdir / f"{safe}.jf.fa"
        jf     = outdir / f"{safe}.jf"
        try:
            fa.write_text(f">{seq.id}\n{seq.sequence}\n", encoding="utf-8")
            count = subprocess.run(
                [self.jellyfish_bin, "count", "-m", str(self.kmer_size),
                 "-s", "10M", "-C", "-t", str(self.threads), "-o", str(jf), str(fa)],
                capture_output=True, text=True,
            )
            if count.returncode != 0:
                logger.debug(f"jellyfish count failed for {seq.id}: {count.stderr.strip()}")
                return (0, 0, 0.0)
            histo = subprocess.run(
                [self.jellyfish_bin, "histo", str(jf)],
                capture_output=True, text=True,
            )
            if histo.returncode != 0:
                logger.debug(f"jellyfish histo failed for {seq.id}: {histo.stderr.strip()}")
                return (0, 0, 0.0)
        finally:
            if self.cleanup:
                fa.unlink(missing_ok=True)
                jf.unlink(missing_ok=True)

        # histo lines: "<multiplicity> <num_distinct_kmers_with_that_multiplicity>"
        distinct = total = singletons = 0
        for line in histo.stdout.splitlines():
            bits = line.split()
            if len(bits) != 2:
                continue
            mult, d = int(bits[0]), int(bits[1])
            distinct += d
            total    += mult * d
            if mult == 1:
                singletons += d
        score = round((total - singletons) / total, 4) if total else 0.0
        return (distinct, total, score)

    # ── dot plot (native SVG segments) ──────────────────────────────────────────

    @classmethod
    def _revcomp(cls, s: str) -> str:
        return s.translate(_COMPLEMENT)[::-1]

    _MAX_KMER_HITS = 200   # skip ultra-repetitive k-mers (low-complexity noise)

    def _dotplot(self, seq: NucSequence) -> DotPlot | None:
        """
        Self-similarity dot plot as extended line segments in bp coordinates.

        A **full** k-mer index over every position is built, then seeded from
        sampled anchors (every ``stride`` bases, so each axis has at most
        ``_DOTPLOT_MAX_DIM`` seed points).  Because the index holds all positions,
        a repeat is detected whatever its offset — sampling only the *seeds*, not
        the index, avoids missing repeats whose two copies don't share the
        sampling lattice.  Forward matches on offset ``d = y - x`` become segments
        along a diagonal (``d = 0`` is the identity line); reverse-complement
        matches on anti-diagonal ``s = x + y`` mark inverted repeats.  Total
        segments are capped for a lightweight SVG.
        """
        w = self.dotplot_word
        s = seq.sequence.upper()
        L = len(s)
        if L < w:
            return None

        stride  = max(1, L // self._DOTPLOT_MAX_DIM)
        anchors = range(0, L - w + 1, stride)

        # Full forward k-mer index — every start position, not just the anchors.
        fwd: dict[str, list[int]] = {}
        for i in range(0, L - w + 1):
            fwd.setdefault(s[i:i + w], []).append(i)

        # Direct repeats: from each anchor, all forward copies of its k-mer
        # (offset d = j - i ≥ 0; d = 0 is the identity diagonal).
        by_offset: dict[int, list[int]] = {}
        for i in anchors:
            hits = fwd.get(s[i:i + w])
            if not hits or len(hits) > self._MAX_KMER_HITS:
                continue
            for j in hits:
                if j >= i:
                    by_offset.setdefault(j - i, []).append(i)

        segments: list[list[int]] = []
        self._emit_segments(by_offset, w, stride, anti=False, out=segments)

        # Inverted repeats: seq[i:] matches revcomp(seq[j:]) ⇒ anti-diagonal.
        by_sum: dict[int, list[int]] = {}
        for i in anchors:
            hits = fwd.get(self._revcomp(s[i:i + w]))
            if not hits or len(hits) > self._MAX_KMER_HITS:
                continue
            for j in hits:
                if j > i:                    # dedupe symmetric anti-diagonal matches
                    by_sum.setdefault(i + j, []).append(i)
        self._emit_segments(by_sum, w, stride, anti=True, out=segments)

        if not segments:
            return None
        return DotPlot(word=w, length=L, segments=segments[: self._DOTPLOT_MAX_SEGMENTS])

    @classmethod
    def _emit_segments(
        cls, groups: dict[int, list[int]], w: int, stride: int,
        anti: bool, out: list[list[int]],
    ) -> None:
        """
        Merge collinear anchor matches into segments.  ``groups`` maps a constant
        (offset for diagonals, sum for anti-diagonals) to the list of x anchors;
        contiguous x runs (gap <= 2·stride) become one [x1, y1, x2, y2] segment.
        """
        for const, xs in groups.items():
            xs = sorted(set(xs))
            if not xs:
                continue
            run_start = prev = xs[0]
            for x in xs[1:] + [None]:
                if x is not None and x - prev <= stride * 2:
                    prev = x
                    continue
                x1, x2 = run_start, prev + w
                if anti:                     # y = const - x
                    y1, y2 = const - run_start, const - (prev + w)
                else:                        # y = x + const
                    y1, y2 = run_start + const, (prev + w) + const
                out.append([x1, y1, x2, y2])
                if len(out) >= cls._DOTPLOT_MAX_SEGMENTS:
                    return
                if x is not None:
                    run_start = prev = x

    # ── assembly ────────────────────────────────────────────────────────────────

    def _build(
        self,
        seqs:    list[NucSequence],
        low_cx:  dict[str, list[list[int]]],
        repeats: dict[str, list[SelfRepeat]],
        outdir:  Path,
    ) -> dict[str, SequenceQuality]:
        out: dict[str, SequenceQuality] = {}
        for seq in seqs:
            regions = low_cx.get(seq.id, [])
            masked  = sum(b - a + 1 for a, b in regions)
            distinct, total, score = self._jellyfish(seq, outdir)
            out[seq.id] = SequenceQuality(
                seq_id                 = seq.id,
                length                 = seq.length,
                low_complexity_regions = regions,
                self_repeats           = repeats.get(seq.id, []),
                kmer_size              = self.kmer_size,
                kmer_distinct          = distinct,
                kmer_total             = total,
                kmer_repeat_score      = score,
                low_complexity_frac    = round(masked / seq.length, 4) if seq.length else 0.0,
                dotplot                = self._dotplot(seq),
            )
        return out
