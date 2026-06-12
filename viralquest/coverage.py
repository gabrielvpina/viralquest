"""
coverage.py — Per-base read-coverage profiling of viral sequences.

When the user supplies sequencing reads (``--reads``), the confirmed viral
sequences are used as a small reference; the reads are aligned with minimap2 and
per-base depth is computed with ``samtools depth``.  The resulting
``CoverageProfile`` (one per sequence) is attached to each ``NucSequence`` and
rendered as a coverage track in the genome viewer.

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
            depth = self._compute_depth(seqs, bam)
            profiles = self._build_profiles(seqs, depth)
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

    def _compute_depth(
        self, seqs: list[NucSequence], bam: Path
    ) -> dict[str, np.ndarray]:
        """Run ``samtools depth -a`` and return {seq_id: per-base depth array}."""
        depth = {s.id: np.zeros(s.length, dtype=np.float64) for s in seqs}

        proc = subprocess.Popen(
            [self.samtools_bin, "depth", "-a", "-@", str(self.threads), str(bam)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            arr = depth.get(parts[0])
            if arr is None:
                continue
            pos = int(parts[1]) - 1          # samtools depth is 1-based
            if 0 <= pos < arr.size:
                arr[pos] = float(parts[2])
        err = proc.stderr.read() if proc.stderr else ""
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"samtools depth failed: {err.strip()}")
        return depth

    def _build_profiles(
        self, seqs: list[NucSequence], depth: dict[str, np.ndarray]
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

            profiles[seq.id] = CoverageProfile(
                seq_id          = seq.id,
                length          = int(arr.size),
                mean_depth      = round(mean, 2),
                max_depth       = round(float(arr.max()), 2),
                cv              = round(cv, 3),
                breadth_1x      = round(breadth, 4),
                bins            = bins,
                low_cov_regions = low_runs,
            )
        return profiles

    @staticmethod
    def _export_tsv(profiles: "dict[str, CoverageProfile]", outdir: Path) -> None:
        """Write one TSV per sequence with bin index, genomic position and mean depth."""
        import re
        _safe = re.compile(r'[^\w\-.]')
        for profile in profiles.values():
            safe_name = _safe.sub('_', profile.seq_id)
            tsv_path  = outdir / f"{safe_name}.tsv"
            nb        = len(profile.bins)
            with open(tsv_path, "w", encoding="utf-8") as fh:
                fh.write("bin\tposition_nt\tmean_depth\n")
                for i, depth in enumerate(profile.bins):
                    pos = int(((i + 0.5) / nb) * profile.length)
                    fh.write(f"{i + 1}\t{pos}\t{depth}\n")

    @classmethod
    def _downsample(cls, arr: np.ndarray) -> list[float]:
        """Mean depth per bin, capped at _N_BINS points."""
        if arr.size <= cls._N_BINS:
            return [round(float(v), 2) for v in arr]
        chunks = np.array_split(arr, cls._N_BINS)
        return [round(float(c.mean()), 2) for c in chunks]

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
