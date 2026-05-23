# viralquest/diamond.py
import csv
import math
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from loguru import logger
from .biodata import NucSequence, BlastxResult


class DiamondRunner:
    
    _OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    def __init__(
        self,
        diamond_bin: str = "diamond",
        db_path: str = "",
        threads: int = 4,
        e_value: float = 1e-5,
        max_target_seqs: int = 5,
        outdir: str | None = None,
        batch_size: int = 1000,       # <-- new
    ):
        self.diamond_bin = diamond_bin
        self.db_path = db_path
        self.threads = threads
        self.e_value = e_value
        self.max_target_seqs = max_target_seqs
        self.outdir = Path(outdir) if outdir else Path(tempfile.gettempdir())
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size

    @staticmethod
    def _chunk(seqs: list[NucSequence], size: int) -> list[list[NucSequence]]:
        """Split a flat list into sublists of at most `size` items."""
        return [seqs[i : i + size] for i in range(0, len(seqs), size)]

    def _write_batch_fasta(self, batch: list[NucSequence], path: Path) -> None:
        """Write multiple sequences into one FASTA file."""
        with open(path, "w", encoding="utf-8") as fh:
            for seq in batch:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    def run_blastx(self, batch: list[NucSequence], batch_index: int) -> Path:
        """
        Runs diamond blastx for a batch of sequences.
        Returns the path to the TSV output file.
        """
        fasta_path = self.outdir / f"batch_{batch_index}.fa"
        out_path   = self.outdir / f"batch_{batch_index}.tsv"
        self._write_batch_fasta(batch, fasta_path)

        cmd = [
            self.diamond_bin, "blastx",
            "--db",              self.db_path,
            "--query",           str(fasta_path),
            "--out",             str(out_path),
            "--outfmt",          *self._OUTFMT.split(),
            "--threads",         str(self.threads),
            "--evalue",          str(self.e_value),
            "--max-target-seqs", str(self.max_target_seqs),
            "--quiet",
        ]
        logger.debug(f"Running diamond — batch {batch_index} ({len(batch)} sequences)")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"Diamond failed on batch {batch_index}: {result.stderr.strip()}")
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    def run_all(
        self, nuc_seqs: list[NucSequence], max_workers: int = 4
    ) -> list[Path]:
        """
        Chunks sequences into batches of `self.batch_size`, runs each batch
        in parallel, and returns the ordered list of TSV output paths.
        """
        batches = self._chunk(nuc_seqs, self.batch_size)
        n = len(batches)
        logger.info(
            f"Dispatching {len(nuc_seqs)} sequences across {n} batch(es) "
            f"(batch_size={self.batch_size}, workers={max_workers})"
        )

        tsv_paths: list[Path | None] = [None] * n

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(self.run_blastx, batch, i): i
                for i, batch in enumerate(batches)
            }
            for future in as_completed(futures):
                i = futures[future]
                try:
                    tsv_paths[i] = future.result()
                    logger.success(f"Batch {i} done → {tsv_paths[i].name}")
                except Exception as exc:
                    logger.error(f"Batch {i} failed: {exc}")

        return [p for p in tsv_paths if p is not None]


class DiamondOutputParser:
    """Parses Diamond TSV output into BlastxResult objects."""

    @staticmethod
    def parse(tsv_path: Path) -> list[BlastxResult]:
        hits = []
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return hits
        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 12:
                    continue
                try:
                    hits.append(BlastxResult(
                        query_id=row[0],
                        subject_id=row[1],
                        pct_identity=float(row[2]),
                        aln_length=int(row[3]),
                        mismatches=int(row[4]),
                        gap_opens=int(row[5]),
                        query_start=int(row[6]),
                        query_end=int(row[7]),
                        subject_start=int(row[8]),
                        subject_end=int(row[9]),
                        e_value=float(row[10]),
                        bit_score=float(row[11]),
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(f"Skipping malformed row in {tsv_path}: {exc}")
        return hits


class DiamondResultAttacher:
    """Attaches parsed BlastxResult hits to the right NucSequence."""

    @staticmethod
    def attach(
        hits: list[BlastxResult],
        seq: NucSequence,
    ) -> None:
        for hit in hits:
            hit.compute_coverage(seq.length)
            seq.blast_hits.append(hit)
        logger.debug(f"{len(hits)} hits attached to {seq.id}")


class DiamondFilterer:
    """Post-hoc filtering of BlastxResult lists."""

    @staticmethod
    def filter(
        hits: list[BlastxResult],
        min_identity: float = 0.0,
        max_evalue: float = 1e-5,
        min_coverage: float = 0.0,
    ) -> list[BlastxResult]:
        return [
            h for h in hits
            if h.pct_identity >= min_identity
            and h.e_value <= max_evalue
            and h.query_coverage >= min_coverage
        ]