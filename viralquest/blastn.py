import csv
import glob
import lzma
import os
import subprocess
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO, BiopythonWarning
from loguru import logger

from .biodata import NucSequence


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class BlastnResult:
    """One BLASTn hit — compatible with both local and online search."""
    query_id: str
    query_length: int
    subject_length: int
    query_coverage: float
    pct_identity: float
    e_value: float
    subject_title: str


# ---------------------------------------------------------------------------
# Binary resolver (extracted from the original runBlastn.py logic)
# ---------------------------------------------------------------------------

class BlastnBinaryResolver:
    """
    Finds or extracts the blastn binary.

    Search order:
      1. An explicit path passed by the caller.
      2. The bundled binary at viralquest/bin/blastn.
      3. The compressed fallback at viralquest/bin/blastn.xz.
    """

    def __init__(self, explicit_path: str | None = None):
        self.explicit_path = explicit_path
        self._resolved: Path | None = None

    def resolve(self) -> Path:
        if self._resolved is not None:
            return self._resolved

        if self.explicit_path:
            path = Path(self.explicit_path)
            if not path.exists():
                raise FileNotFoundError(f"BLASTn binary not found at: {path}")
            self._make_executable(path)
            self._resolved = path
            return self._resolved

        # locate relative to this file → project_root/viralquest/bin/
        bin_dir = Path(__file__).resolve().parent / "bin"
        binary  = bin_dir / "blastn"
        compressed = bin_dir / "blastn.xz"

        if not binary.exists() and compressed.exists():
            logger.info(f"Extracting {compressed} …")
            try:
                with lzma.open(compressed, "rb") as src, open(binary, "wb") as dst:
                    dst.write(src.read())
            except Exception as exc:
                raise RuntimeError(f"Failed to extract {compressed}: {exc}") from exc

        if not binary.exists():
            raise FileNotFoundError(
                f"BLASTn binary not found. Expected: {binary} (or {compressed})"
            )

        self._make_executable(binary)
        self._resolved = binary
        return self._resolved

    @staticmethod
    def _make_executable(path: Path) -> None:
        if not os.access(path, os.X_OK):
            path.chmod(0o755)


# ---------------------------------------------------------------------------
# Output parser (shared between local and online results)
# ---------------------------------------------------------------------------

class BlastnOutputParser:
    """Parses a BLASTn TSV (outfmt 6 variant) into BlastnResult objects."""

    @staticmethod
    def parse(tsv_path: Path) -> list[BlastnResult]:
        hits: list[BlastnResult] = []
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return hits

        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 7:
                    continue
                try:
                    hits.append(BlastnResult(
                        query_id=row[0],
                        query_length=int(row[1]),
                        subject_length=int(row[2]),
                        query_coverage=float(row[3]),
                        pct_identity=float(row[4]),
                        e_value=float(row[5]),
                        subject_title=row[6],
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(f"Skipping malformed row in {tsv_path.name}: {exc}")
        return hits


# ---------------------------------------------------------------------------
# Result attacher
# ---------------------------------------------------------------------------

class BlastnResultAttacher:
    """Attaches parsed BlastnResult hits to the matching NucSequence objects."""

    @staticmethod
    def attach(hits: list[BlastnResult], nuc_seqs: list[NucSequence]) -> None:
        seq_map = {seq.id: seq for seq in nuc_seqs}
        for hit in hits:
            seq = seq_map.get(hit.query_id)
            if seq is None:
                logger.warning(f"No NucSequence found for query '{hit.query_id}' — skipping.")
                continue
            seq.blast_hits.append(hit)
        logger.debug(f"{len(hits)} BLASTn hits attached.")


# ---------------------------------------------------------------------------
# Local BLASTn runner
# ---------------------------------------------------------------------------

class BlastnRunner:
    """
    Runs the local blastn binary against batches of NucSequence objects.

    Batching strategy mirrors DiamondRunner: up to `batch_size` sequences
    are written into one temporary FASTA, then a single subprocess is
    launched per batch. Batches run in parallel via ThreadPoolExecutor.

    Columns written (outfmt 6 custom):
        qseqid  qlen  slen  qcovs  pident  evalue  stitle
    """

    _OUTFMT = "6 qseqid qlen slen qcovs pident evalue stitle"

    def __init__(
        self,
        database: str,
        blastn_path: str | None = None,
        threads: int = 4,
        e_value: float = 1e-5,
        max_target_seqs: int = 1,
        outdir: str | None = None,
        batch_size: int = 1000,
    ):
        self.database = database
        self.resolver = BlastnBinaryResolver(blastn_path)
        self.threads = threads
        self.e_value = e_value
        self.max_target_seqs = max_target_seqs
        self.outdir = Path(outdir) if outdir else Path(tempfile.gettempdir())
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _chunk(seqs: list[NucSequence], size: int) -> list[list[NucSequence]]:
        return [seqs[i : i + size] for i in range(0, len(seqs), size)]

    def _write_batch_fasta(self, batch: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in batch:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    def _run_batch(self, batch: list[NucSequence], batch_index: int) -> Path:
        """Runs blastn for one batch. Returns path to the TSV output."""
        binary    = self.resolver.resolve()
        fasta_path = self.outdir / f"blastn_batch_{batch_index}.fa"
        out_path   = self.outdir / f"blastn_batch_{batch_index}.tsv"

        self._write_batch_fasta(batch, fasta_path)

        cmd = [
            str(binary),
            "-query",           str(fasta_path),
            "-db",              self.database,
            "-out",             str(out_path),
            "-outfmt",          self._OUTFMT,
            "-num_threads",     str(self.threads),
            "-evalue",          str(self.e_value),
            "-max_target_seqs", str(self.max_target_seqs),
        ]

        logger.debug(f"blastn — batch {batch_index} ({len(batch)} sequences)")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"blastn failed on batch {batch_index}: {result.stderr.strip()}")
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_all(
        self, nuc_seqs: list[NucSequence], max_workers: int = 4
    ) -> list[Path]:
        """
        Chunks sequences, dispatches batches in parallel, returns ordered
        list of TSV paths. Pass results to BlastnOutputParser.parse().
        """
        batches = self._chunk(nuc_seqs, self.batch_size)
        logger.info(
            f"blastn: {len(nuc_seqs)} sequences → {len(batches)} batch(es) "
            f"(batch_size={self.batch_size}, workers={max_workers})"
        )

        tsv_paths: list[Path | None] = [None] * len(batches)

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(self._run_batch, batch, i): i
                for i, batch in enumerate(batches)
            }
            for future in as_completed(futures):
                i = futures[future]
                try:
                    tsv_paths[i] = future.result()
                    logger.success(f"blastn batch {i} done → {tsv_paths[i].name}")
                except Exception as exc:
                    logger.error(f"blastn batch {i} failed: {exc}")

        return [p for p in tsv_paths if p is not None]


# ---------------------------------------------------------------------------
# Online BLASTn runner (NCBI qblast)
# ---------------------------------------------------------------------------

class BlastnOnlineRunner:
    """
    Submits sequences to NCBI qblast one at a time (NCBI rate-limit policy).

    Results are written in the same TSV format as BlastnRunner so that
    BlastnOutputParser.parse() works identically for both backends.
    """

    def __init__(
        self,
        database: str = "nt",
        email: str = "",
        hitlist_size: int = 1,
        sleep_interval: float = 1.0,
        outdir: str | None = None,
    ):
        if not email:
            raise ValueError("NCBI requires an e-mail address for online BLAST.")
        self.database = database
        self.email = email
        self.hitlist_size = hitlist_size
        self.sleep_interval = sleep_interval
        self.outdir = Path(outdir) if outdir else Path(tempfile.gettempdir())
        self.outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _query_one(self, seq: NucSequence) -> BlastnResult | None:
        """Runs qblast for a single sequence. Returns a result or None."""
        from Bio.Blast import NCBIXML, NCBIWWW

        NCBIWWW.email = self.email

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonWarning)
                handle = NCBIWWW.qblast(
                    program="blastn",
                    database=self.database,
                    sequence=seq.sequence,
                    format_type="XML",
                    hitlist_size=self.hitlist_size,
                )

            record = NCBIXML.read(handle)
            handle.close()

            if not record.alignments:
                logger.debug(f"No hits for {seq.id}")
                return BlastnResult(
                    query_id=seq.id,
                    query_length=seq.length,
                    subject_length=0,
                    query_coverage=0.0,
                    pct_identity=0.0,
                    e_value=float("nan"),
                    subject_title="no_hit",
                )

            alignment = record.alignments[0]
            hsp = alignment.hsps[0]
            coverage = (hsp.align_length / record.query_length) * 100
            identity = (hsp.identities / hsp.align_length) * 100

            return BlastnResult(
                query_id=seq.id,
                query_length=record.query_length,
                subject_length=alignment.length,
                query_coverage=round(coverage, 2),
                pct_identity=round(identity, 2),
                e_value=hsp.expect,
                subject_title=alignment.title,
            )

        except Exception as exc:
            logger.error(f"qblast failed for {seq.id}: {exc}")
            return None

    def _write_tsv(self, results: list[BlastnResult], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for r in results:
                fh.write(
                    f"{r.query_id}\t{r.query_length}\t{r.subject_length}\t"
                    f"{r.query_coverage}\t{r.pct_identity}\t{r.e_value}\t"
                    f"{r.subject_title}\n"
                )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_all(self, nuc_seqs: list[NucSequence]) -> Path:
        """
        Queries NCBI sequentially (rate-limit safe). Returns path to TSV.
        Note: no parallelism here — NCBI enforces per-IP request limits.
        """
        out_path = self.outdir / "blastn_online.tsv"
        results: list[BlastnResult] = []

        logger.info(f"Online BLASTn: {len(nuc_seqs)} sequences → NCBI '{self.database}'")

        for i, seq in enumerate(nuc_seqs):
            logger.debug(f"  [{i+1}/{len(nuc_seqs)}] querying {seq.id} …")
            result = self._query_one(seq)
            if result is not None:
                results.append(result)
            time.sleep(self.sleep_interval)

        self._write_tsv(results, out_path)
        logger.success(f"Online BLASTn done — {len(results)} results → {out_path.name}")
        return out_path