"""
coverage.py — Per-base read-coverage profiling of viral sequences.

When the user supplies sequencing reads (``--reads``), the confirmed viral
sequences are used as a small reference; the reads are aligned with minimap2 and
per-base depth **and mean base quality** are read from the *same* BAM in a single
``samtools mpileup`` pass (no second alignment).  The resulting
``CoverageProfile`` (one per sequence) is attached to each ``NucSequence`` and
rendered as a coverage track in the genome viewer, coloured by base quality.

The biological motive: a correctly assembled viral contig has roughly uniform
coverage, while a chimeric / mis-assembled contig tends to show a coverage
discontinuity (a step change, or a near-zero internal drop) at the artificial
junction.  The track makes that signal visible under the ORF lanes.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
from loguru import logger

from viralquest.biodata import CoverageProfile, NucSequence


# ---------------------------------------------------------------------------
# CoveragePipeline
# ---------------------------------------------------------------------------

class CoveragePipeline:
    """
    Align reads to the viral sequences and build a per-sequence coverage profile.

    Single-end (one reads file) and paired-end (two files) are both supported,
    selected by the number of entries in ``reads`` — same convention as
    ``SalmonQuantRunner``.
    """

    _N_BINS         = 600     # max plotted points per sequence
    _DROP_FRACTION  = 0.20    # depth < mean * this → "low coverage"
    _MIN_RUN_FRAC   = 0.01    # ignore low-cov runs shorter than length * this …
    _MIN_RUN_ABS    = 20      # … or this many bases, whichever is larger

    # minimap2 preset per read technology (--read-type CLI flag)
    _MM2_PRESET: dict[str, str] = {
        "sr":   "sr",        # Illumina short reads
        "ont":  "map-ont",   # Oxford Nanopore
        "pb":   "map-pb",    # PacBio CLR
        "hifi": "map-hifi",  # PacBio HiFi / CCS
    }

    def __init__(
        self,
        threads:       int  = 4,
        low_memory:    bool = False,
        read_type:     str  = "sr",
        minimap2_bin:  str  = "minimap2",
        samtools_bin:  str  = "samtools",
        cleanup:       bool = True,
    ):
        self.threads      = threads
        self.low_memory   = low_memory
        self.read_type    = read_type
        self.minimap2_bin = minimap2_bin
        self.samtools_bin = samtools_bin
        self.cleanup      = cleanup

    # ── public ────────────────────────────────────────────────────────────────

    def run(
        self,
        reads:      list[str],
        viral_seqs: list[NucSequence],
        outdir:     Path,
    ) -> dict[str, CoverageProfile]:
        """
        Returns ``{seq_id: CoverageProfile}`` for every sequence with mapped reads.
        Sequences with no coverage (or on any tooling/IO failure) are simply
        absent from the returned dict.
        """
        if len(reads) not in (1, 2):
            logger.warning(f"CoveragePipeline: expected 1 or 2 reads files, got {len(reads)} — skipping.")
            return {}

        seqs = [s for s in viral_seqs if s.sequence]
        if not seqs:
            logger.warning("CoveragePipeline: no viral sequences — skipping coverage.")
            return {}

        if not shutil.which(self.minimap2_bin) or not shutil.which(self.samtools_bin):
            logger.warning(
                "CoveragePipeline: minimap2 and/or samtools not found on PATH — "
                "coverage track will be unavailable."
            )
            return {}

        outdir.mkdir(parents=True, exist_ok=True)
        ref_fasta = outdir / "viral_ref.fasta"
        bam       = outdir / "coverage.bam"

        try:
            self._write_fasta(seqs, ref_fasta)
            if not self._align(ref_fasta, reads, bam):
                return {}
            depth, quality, fwd, rev = self._compute_pileup(seqs, bam)
            profiles = self._build_profiles(seqs, depth, quality, fwd, rev)
        except Exception as exc:
            logger.error(f"CoveragePipeline failed — skipping coverage: {exc}")
            return {}
        finally:
            if self.cleanup:
                for p in (ref_fasta, bam, bam.with_suffix(".bam.bai")):
                    p.unlink(missing_ok=True)

        self._export_tsv(profiles, outdir)
        logger.success(f"Coverage profiles built for {len(profiles)}/{len(seqs)} viral sequence(s).")
        return profiles

    # ── private ───────────────────────────────────────────────────────────────

    @staticmethod
    def _write_fasta(seqs: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in seqs:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    def _align(self, ref_fasta: Path, reads: list[str], bam: Path) -> bool:
        """minimap2 -ax sr | samtools sort → indexed BAM.  Returns success."""
        preset  = self._MM2_PRESET.get(self.read_type, "sr")
        mm2_cmd = [
            self.minimap2_bin, "-ax", preset,
            "-t", str(self.threads),
            str(ref_fasta), *reads,
        ]
        sort_cmd = [self.samtools_bin, "sort", "-@", str(self.threads), "-o", str(bam)]

        mode = "paired" if len(reads) == 2 else "single"
        logger.info(f"minimap2 [{mode}-end, preset={preset}, {self.threads} threads] → '{bam.name}' ...")

        mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sort = subprocess.run(sort_cmd, stdin=mm2.stdout, capture_output=True, text=True)
        if mm2.stdout:
            mm2.stdout.close()
        mm2_err = mm2.stderr.read().decode("utf-8", "replace") if mm2.stderr else ""
        mm2.wait()

        if mm2.returncode != 0:
            logger.error(f"minimap2 failed: {mm2_err.strip().splitlines()[-1] if mm2_err.strip() else 'unknown error'}")
            return False
        if sort.returncode != 0:
            logger.error(f"samtools sort failed: {sort.stderr.strip()}")
            return False

        idx = subprocess.run(
            [self.samtools_bin, "index", str(bam)], capture_output=True, text=True
        )
        if idx.returncode != 0:
            logger.error(f"samtools index failed: {idx.stderr.strip()}")
            return False
        return True

    # Forward-strand base symbols in the mpileup read-bases column are uppercase
    # (``.`` match, ``ACGTN`` mismatch); reverse-strand are lowercase (``,`` /
    # ``acgtn``). Deletions (``*``) and reference skips (``><``) carry no strand.
    _FWD_BASES = frozenset(".ACGTN")
    _REV_BASES = frozenset(",acgtn")

    @classmethod
    def _count_strands(cls, read_bases: str) -> tuple[int, int]:
        """
        Count forward vs reverse aligned bases in one mpileup read-bases string,
        correctly skipping read-start markers (``^x``), read ends (``$``) and
        insertion/deletion runs (``[+-]<len><bases>``).
        """
        fwd = rev = 0
        i, n = 0, len(read_bases)
        while i < n:
            c = read_bases[i]
            if c == "^":            # read start: next char is mapping quality
                i += 2
                continue
            if c == "$":            # read end marker
                i += 1
                continue
            if c in "+-":           # indel: sign, length, then <length> bases
                i += 1
                num = ""
                while i < n and read_bases[i].isdigit():
                    num += read_bases[i]
                    i += 1
                i += int(num) if num else 0
                continue
            if c in cls._FWD_BASES:
                fwd += 1
            elif c in cls._REV_BASES:
                rev += 1
            i += 1
        return fwd, rev

    def _compute_pileup(
        self, seqs: list[NucSequence], bam: Path
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray],
               dict[str, np.ndarray], dict[str, np.ndarray]]:
        """
        Read per-base depth, mean base quality **and strand-split depth** from the
        same BAM in one ``samtools mpileup`` pass — no second alignment.

        Returns ``(depth, quality, depth_fwd, depth_rev)`` keyed by seq_id. Strand
        is taken from the case of the read-bases column (uppercase = forward /
        sense, lowercase = reverse / antisense). The pileup filters are all
        disabled (``-A -B -Q 0 -q 0``) so depth matches the previous
        ``samtools depth -a`` count, and every position is emitted
        (``-a``, ``-d 0`` removes the depth cap).
        """
        depth   = {s.id: np.zeros(s.length, dtype=np.float64) for s in seqs}
        quality = {s.id: np.zeros(s.length, dtype=np.float64) for s in seqs}
        fwd     = {s.id: np.zeros(s.length, dtype=np.float64) for s in seqs}
        rev     = {s.id: np.zeros(s.length, dtype=np.float64) for s in seqs}

        proc = subprocess.Popen(
            [self.samtools_bin, "mpileup",
             "-a", "-A", "-B", "-Q", "0", "-q", "0", "-d", "0",
             str(bam)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            # mpileup: 0=name 1=pos 2=ref 3=depth 4=read_bases 5=base_quals
            if len(parts) < 4:
                continue
            arr_d = depth.get(parts[0])
            if arr_d is None:
                continue
            pos = int(parts[1]) - 1          # mpileup is 1-based
            if not (0 <= pos < arr_d.size):
                continue
            arr_d[pos] = float(parts[3])
            # read-bases column → forward / reverse strand split
            if len(parts) >= 5 and parts[4] and parts[4] != "*":
                f, r = self._count_strands(parts[4])
                fwd[parts[0]][pos] = f
                rev[parts[0]][pos] = r
            # base-quality column: one ASCII char (Phred+33) per stacked base
            if len(parts) >= 6 and parts[5] and parts[5] != "*":
                quals = [ord(c) - 33 for c in parts[5]]
                if quals:
                    quality[parts[0]][pos] = sum(quals) / len(quals)
        err = proc.stderr.read() if proc.stderr else ""
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"samtools mpileup failed: {err.strip()}")
        return depth, quality, fwd, rev

    def _build_profiles(
        self,
        seqs:    list[NucSequence],
        depth:   dict[str, np.ndarray],
        quality: dict[str, np.ndarray],
        fwd:     dict[str, np.ndarray],
        rev:     dict[str, np.ndarray],
    ) -> dict[str, CoverageProfile]:
        profiles: dict[str, CoverageProfile] = {}
        for seq in seqs:
            arr = depth.get(seq.id)
            if arr is None or arr.size == 0:
                continue
            mean = float(arr.mean())
            if mean <= 0:
                continue   # no reads mapped → no track

            std        = float(arr.std())
            cv         = std / mean if mean > 0 else 0.0
            breadth    = float((arr >= 1).mean())
            bins       = self._downsample(arr)
            low_runs   = self._low_cov_regions(arr, mean)

            # Base quality is only meaningful where reads map; average over
            # covered positions so uncovered gaps don't dilute the colour.
            qarr        = quality.get(seq.id)
            covered     = arr > 0
            qual_bins   = self._downsample_quality(qarr, covered) if qarr is not None else []
            mean_qual   = float(qarr[covered].mean()) if (qarr is not None and covered.any()) else 0.0

            # Strand-split depth (sense = forward, antisense = reverse).
            farr = fwd.get(seq.id)
            rarr = rev.get(seq.id)
            bins_fwd = self._downsample(farr) if farr is not None else []
            bins_rev = self._downsample(rarr) if rarr is not None else []
            mean_fwd = float(farr.mean()) if farr is not None else 0.0
            mean_rev = float(rarr.mean()) if rarr is not None else 0.0

            profiles[seq.id] = CoverageProfile(
                seq_id          = seq.id,
                length          = int(arr.size),
                mean_depth      = round(mean, 2),
                max_depth       = round(float(arr.max()), 2),
                cv              = round(cv, 3),
                breadth_1x      = round(breadth, 4),
                bins            = bins,
                low_cov_regions = low_runs,
                quality_bins    = qual_bins,
                mean_quality    = round(mean_qual, 2),
                bins_fwd        = bins_fwd,
                bins_rev        = bins_rev,
                mean_depth_fwd  = round(mean_fwd, 2),
                mean_depth_rev  = round(mean_rev, 2),
            )
        return profiles

    @staticmethod
    def _export_tsv(profiles: "dict[str, CoverageProfile]", outdir: Path) -> None:
        """
        Write one TSV per sequence: bin index, genomic position, total mean depth
        and the strand-split (sense / antisense) depth per bin.
        """
        import re
        _safe = re.compile(r'[^\w\-.]')
        for profile in profiles.values():
            safe_name = _safe.sub('_', profile.seq_id)
            tsv_path  = outdir / f"{safe_name}.tsv"
            nb        = len(profile.bins)
            fwd = profile.bins_fwd or []
            rev = profile.bins_rev or []
            with open(tsv_path, "w", encoding="utf-8") as fh:
                fh.write("bin\tposition_nt\tmean_depth\tsense_depth\tantisense_depth\n")
                for i, depth in enumerate(profile.bins):
                    pos = int(((i + 0.5) / nb) * profile.length)
                    f = fwd[i] if i < len(fwd) else ""
                    r = rev[i] if i < len(rev) else ""
                    fh.write(f"{i + 1}\t{pos}\t{depth}\t{f}\t{r}\n")

    @classmethod
    def _downsample(cls, arr: np.ndarray) -> list[float]:
        """Mean depth per bin, capped at _N_BINS points."""
        if arr.size <= cls._N_BINS:
            return [round(float(v), 2) for v in arr]
        chunks = np.array_split(arr, cls._N_BINS)
        return [round(float(c.mean()), 2) for c in chunks]

    @classmethod
    def _downsample_quality(cls, qarr: np.ndarray, covered: np.ndarray) -> list[float]:
        """
        Mean base quality per bin, aligned 1:1 with ``_downsample`` bins, but
        averaged only over covered positions (``covered``) so gaps read as 0
        rather than dragging the colour down.
        """
        if qarr.size <= cls._N_BINS:
            return [round(float(q), 2) if c else 0.0 for q, c in zip(qarr, covered)]
        q_chunks = np.array_split(qarr,    cls._N_BINS)
        c_chunks = np.array_split(covered, cls._N_BINS)
        out: list[float] = []
        for q, c in zip(q_chunks, c_chunks):
            out.append(round(float(q[c].mean()), 2) if c.any() else 0.0)
        return out

    @classmethod
    def _low_cov_regions(cls, arr: np.ndarray, mean: float) -> list[list[int]]:
        """
        Internal [start, end] runs (0-based, inclusive) where depth falls below
        ``mean * _DROP_FRACTION``.  Runs touching either end are dropped (natural
        coverage taper), as are runs shorter than the min-run threshold.
        """
        threshold = mean * cls._DROP_FRACTION
        below     = arr < threshold
        if not below.any():
            return []

        min_run = max(cls._MIN_RUN_ABS, int(arr.size * cls._MIN_RUN_FRAC))
        regions: list[list[int]] = []

        # Find contiguous True runs.
        idx   = np.flatnonzero(np.diff(np.concatenate(([0], below.view(np.int8), [0]))))
        edges = idx.reshape(-1, 2)   # [start, end_exclusive] pairs
        last  = arr.size - 1
        for start, end_excl in edges:
            end = int(end_excl) - 1
            if start == 0 or end == last:
                continue                       # touches a sequence end
            if (end - start + 1) < min_run:
                continue                       # too short to flag
            regions.append([int(start), end])
        return regions
