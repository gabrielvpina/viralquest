"""
read_qc.py — Sequencing-read quality control with fastp.

When the user supplies sequencing reads (``--reads``), fastp is run once over the
input read file(s) to produce a global QC snapshot: read/base counts before and
after filtering, Q20/Q30 rates, GC content, duplication and adapter rates, plus
per-cycle quality and GC curves for plotting.  The result is a single
``ReadQcReport`` attached to the report as an optional section.

fastp is designed for short reads; for long-read technologies (``ont``/``pb``/
``hifi``) adapter trimming is disabled and ``adapter_rate`` is reported as 0.

fastp runs in **QC-only** mode: the filtered reads are discarded (written to
``os.devnull``) — the goal is to inform the user about read quality, not to trim.
Only the ``fastp.json`` / ``fastp.html`` QC reports are written to the output
directory; no cleaned FASTQ is produced.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

from loguru import logger

from viralquest.biodata import ReadQcReport


# Long-read technologies — fastp's adapter model does not apply.
_LONG_READ_TYPES = ("ont", "pb", "hifi")


class ReadQcPipeline:
    """
    Run fastp over the input reads and parse its JSON report into a ``ReadQcReport``.

    Single-end (one reads file) and paired-end (two files) are both supported,
    selected by the number of entries in ``reads`` — same convention as
    ``CoveragePipeline`` / ``SalmonQuantRunner``.
    """

    _MAX_CURVE_POINTS = 300     # downsample per-cycle curves to at most this many points

    def __init__(
        self,
        threads:    int  = 4,
        read_type:  str  = "sr",
        fastp_bin:  str  = "fastp",
    ):
        self.threads   = threads
        self.read_type = read_type
        self.fastp_bin = fastp_bin

    # ── public ────────────────────────────────────────────────────────────────

    def run(
        self,
        reads:  list[str],
        outdir: Path,
    ) -> ReadQcReport | None:
        """
        Returns a ``ReadQcReport`` for the supplied reads, or ``None`` on any
        tooling / IO failure (fastp missing, bad input, unparseable JSON).
        """
        if len(reads) not in (1, 2):
            logger.warning(f"ReadQcPipeline: expected 1 or 2 reads files, got {len(reads)} — skipping.")
            return None

        if not shutil.which(self.fastp_bin):
            logger.warning(
                "ReadQcPipeline: fastp not found on PATH — read QC will be unavailable."
            )
            return None

        outdir.mkdir(parents=True, exist_ok=True)
        json_path = outdir / "fastp.json"
        html_path = outdir / "fastp.html"

        try:
            self._run_fastp(reads, json_path, html_path)
            data = json.loads(json_path.read_text(encoding="utf-8"))
            report = self._parse(reads, data)
        except Exception as exc:
            logger.error(f"ReadQcPipeline failed — skipping read QC: {exc}")
            return None

        logger.success(
            f"Read QC: {report.reads_after:,}/{report.reads_before:,} reads passed "
            f"(Q30 {report.q30_rate * 100:.1f}%, GC {report.gc_pct:.1f}%)."
        )
        return report

    # ── private ───────────────────────────────────────────────────────────────

    def _run_fastp(
        self, reads: list[str], json_path: Path, html_path: Path
    ) -> None:
        mode = "paired" if len(reads) == 2 else "single"
        cmd = [
            self.fastp_bin,
            "--thread", str(self.threads),
            "-j", str(json_path),
            "-h", str(html_path),
        ]
        if self.read_type in _LONG_READ_TYPES:
            cmd.append("--disable_adapter_trimming")

        # QC only — discard the filtered reads (we report quality, we don't trim).
        if len(reads) == 2:
            cmd += ["-i", reads[0], "-I", reads[1], "-o", os.devnull, "-O", os.devnull]
        else:
            cmd += ["-i", reads[0], "-o", os.devnull]

        logger.info(f"fastp [{mode}-end, {self.threads} threads] → '{json_path.name}' ...")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            last = proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else "unknown error"
            raise RuntimeError(f"fastp failed: {last}")

    def _parse(self, reads: list[str], data: dict) -> ReadQcReport:
        before = data.get("summary", {}).get("before_filtering", {}) or {}
        after  = data.get("summary", {}).get("after_filtering", {})  or {}
        dup    = data.get("duplication", {}) or {}
        adapt  = data.get("adapter_cutting", {}) or {}

        reads_before = int(before.get("total_reads", 0))
        reads_after  = int(after.get("total_reads", 0))

        # Adapter rate: fraction of input reads that had an adapter trimmed.
        if self.read_type in _LONG_READ_TYPES:
            adapter_rate = 0.0
        else:
            trimmed      = float(adapt.get("adapter_trimmed_reads", 0))
            adapter_rate = round(trimmed / reads_before, 4) if reads_before else 0.0

        # Per-cycle curves come from read1 after filtering (read2 similar for PE).
        r1     = data.get("read1_after_filtering", {}) or {}
        qcurve = (r1.get("quality_curves", {}) or {}).get("mean", []) or []
        gcurve = (r1.get("content_curves", {}) or {}).get("GC", [])   or []

        return ReadQcReport(
            reads            = list(reads),
            read_type        = self.read_type,
            reads_before     = reads_before,
            reads_after      = reads_after,
            bases_before     = int(before.get("total_bases", 0)),
            bases_after      = int(after.get("total_bases", 0)),
            q20_rate         = round(float(after.get("q20_rate", 0.0)), 4),
            q30_rate         = round(float(after.get("q30_rate", 0.0)), 4),
            gc_pct           = round(float(after.get("gc_content", 0.0)) * 100, 2),
            mean_length      = float(after.get("read1_mean_length", 0.0)),
            dup_rate         = round(float(dup.get("rate", 0.0)), 4),
            adapter_rate     = adapter_rate,
            per_base_quality = self._downsample(qcurve),
            per_base_gc      = self._downsample(gcurve),
        )

    @classmethod
    def _downsample(cls, curve: list) -> list[float]:
        """Downsample a per-cycle curve to at most _MAX_CURVE_POINTS block means."""
        vals = [float(v) for v in curve]
        n    = len(vals)
        if n <= cls._MAX_CURVE_POINTS:
            return [round(v, 3) for v in vals]
        step = n / cls._MAX_CURVE_POINTS
        out: list[float] = []
        for i in range(cls._MAX_CURVE_POINTS):
            lo = int(i * step)
            hi = int((i + 1) * step) or lo + 1
            block = vals[lo:hi]
            if block:
                out.append(round(sum(block) / len(block), 3))
        return out
