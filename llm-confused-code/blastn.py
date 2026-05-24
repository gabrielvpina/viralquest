import csv
import subprocess
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import BiopythonWarning
from loguru import logger

from .biodata import NucSequence, BlastnResult


# ---------------------------------------------------------------------------
# Output parser  (shared between local and online — same TSV schema)
# ---------------------------------------------------------------------------

class BlastnOutputParser:
    """
    Parses BLASTn TSV output into BlastnResult objects.

    Expected columns (outfmt 6 custom):
        qseqid  qlen  slen  qcovs  pident  evalue  stitle
    """

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
                        query_coverage=int(float(row[3])),  # qcovs is integer; float() guards "99.0"
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
    """Attaches BlastnResult hits to matching NucSequence objects."""

    @staticmethod
    def attach(hits: list[BlastnResult], nuc_seqs: list[NucSequence]) -> None:
        seq_map = {seq.id: seq for seq in nuc_seqs}
        attached = 0
        for hit in hits:
            seq = seq_map.get(hit.query_id)
            if seq is None:
                logger.warning(f"No NucSequence for query '{hit.query_id}' — skipping.")
                continue
            seq.blastn_hits.append(hit)
            attached += 1
        logger.debug(f"blastn: {attached} hits attached.")


# ---------------------------------------------------------------------------
# Local BLASTn runner
# ---------------------------------------------------------------------------

class BlastnRunner:
    """
    Runs the system blastn binary (expected on PATH via pixi) in batches.

    outfmt columns:
        qseqid  qlen  slen  qcovs  pident  evalue  stitle
    """

    _OUTFMT = "6 qseqid qlen slen qcovs pident evalue stitle"

    def __init__(
        self,
        database: str,
        threads: int = 4,
        e_value: float = 1e-5,
        max_target_seqs: int = 1,
        outdir: str | None = None,
        batch_size: int = 1000,
    ):
        self.database = database
        self.threads = threads
        self.e_value = e_value
        self.max_target_seqs = max_target_seqs
        self.outdir = Path(outdir) if outdir else Path(tempfile.gettempdir())
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size

    @staticmethod
    def _chunk(seqs: list[NucSequence], size: int) -> list[list[NucSequence]]:
        return [seqs[i : i + size] for i in range(0, len(seqs), size)]

    def _write_batch_fasta(self, batch: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in batch:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    def _run_batch(self, batch: list[NucSequence], batch_index: int) -> Path:
        fasta_path = self.outdir / f"blastn_batch{batch_index}.fa"
        out_path   = self.outdir / f"blastn_batch{batch_index}.tsv"
        self._write_batch_fasta(batch, fasta_path)

        cmd = [
            "blastn",
            "-query",           str(fasta_path),
            "-db",              self.database,
            "-out",             str(out_path),
            "-outfmt",          self._OUTFMT,
            "-num_threads",     str(self.threads),
            "-evalue",          str(self.e_value),
            "-max_target_seqs", str(self.max_target_seqs),
        ]

        logger.debug(f"blastn batch {batch_index} — {len(batch)} sequences")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"blastn failed batch {batch_index}: {result.stderr.strip()}")
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    def run_all(
        self, nuc_seqs: list[NucSequence], max_workers: int = 1
    ) -> list[Path]:
        if not nuc_seqs:
            logger.warning("blastn: no sequences to search.")
            return []

        batches = self._chunk(nuc_seqs, self.batch_size)
        logger.info(
            f"blastn: {len(nuc_seqs)} sequences → "
            f"{len(batches)} batch(es), workers={max_workers}"
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
# Online BLASTn runner  (NCBI qblast — always sequential)
# ---------------------------------------------------------------------------

class BlastnOnlineRunner:
    """
    Submits sequences to NCBI qblast one at a time (NCBI rate-limit policy).
    Writes the same 7-column TSV as BlastnRunner so BlastnOutputParser works
    for both backends.
    """

    def __init__(
        self,
        database: str = "nt",
        email: str = "",
        hitlist_size: int = 1,
        sleep_interval: float = 1.5,
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

    def _query_one(self, seq: NucSequence) -> BlastnResult | None:
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
                logger.debug(f"blastn online: no hits for {seq.id}")
                return BlastnResult(
                    query_id=seq.id,
                    query_length=seq.length,
                    subject_length=0,
                    query_coverage=0,
                    pct_identity=0.0,
                    e_value=float("nan"),
                    subject_title="no_hit",
                )

            alignment = record.alignments[0]
            hsp = alignment.hsps[0]
            coverage = int((hsp.align_length / record.query_length) * 100)
            identity = round((hsp.identities / hsp.align_length) * 100, 2)

            return BlastnResult(
                query_id=seq.id,
                query_length=record.query_length,
                subject_length=alignment.length,
                query_coverage=coverage,
                pct_identity=identity,
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

    def run_all(self, nuc_seqs: list[NucSequence]) -> Path:
        """Queries NCBI one sequence at a time. Returns path to TSV."""
        if not nuc_seqs:
            logger.warning("blastn online: no sequences to search.")
            out_path = self.outdir / "blastn_online.tsv"
            out_path.write_text("")
            return out_path

        out_path = self.outdir / "blastn_online.tsv"
        results: list[BlastnResult] = []

        logger.info(
            f"blastn online: {len(nuc_seqs)} sequences → NCBI '{self.database}'"
        )

        for i, seq in enumerate(nuc_seqs):
            logger.debug(f"  [{i+1}/{len(nuc_seqs)}] {seq.id}")
            result = self._query_one(seq)
            if result is not None:
                results.append(result)
            time.sleep(self.sleep_interval)   # respect NCBI rate limits

        self._write_tsv(results, out_path)
        logger.success(
            f"blastn online done — {len(results)} results → {out_path.name}"
        )
        return out_path