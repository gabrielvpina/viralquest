import csv
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from pathlib import Path

from loguru import logger
from .biodata import NucSequence, BlastxResult


class DiamondPhase(Enum):
    REFSEQ_FILTER  = "refseq"   # phase-1: small viral RefSeq → select candidates
    NR_CHARACTERIZE = "nr"      # phase-2: NR global → full characterisation


# Compiled once — matches viral subject titles from NR
_VIRAL_PATTERN = re.compile(
    r"vir(us|al|idae|ales|iform|oid)|phage|bacteriophage",
    re.IGNORECASE,
)


def is_viral_subject(title: str) -> bool:
    """
    Returns True if the subject title looks like a virus.
    Catches: 'virus', 'viral', 'Viridae', 'Virales', 'viroid',
             'phage', 'bacteriophage' — case-insensitive.
    """
    return bool(_VIRAL_PATTERN.search(title))


class DiamondRunner:
    """
    Runs diamond blastx in batches of up to `batch_size` sequences.

    outfmt columns (13):
        qseqid sseqid stitle pident length mismatch gapopen
        qstart qend sstart send evalue bitscore
    """

    _OUTFMT = (
        "6 qseqid sseqid stitle pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore"
    )

    def __init__(
        self,
        db_path: str,
        diamond_bin: str = "diamond",
        threads: int = 2,
        e_value: float = 1e-5,
        max_target_seqs: int = 5,
        outdir: str | None = None,
        batch_size: int = 500,
    ):
        self.db_path = db_path
        self.diamond_bin = diamond_bin
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

    def _run_batch(self, batch: list[NucSequence], batch_index: int, tag: str) -> Path:
        fasta_path = self.outdir / f"dmnd_{tag}_batch{batch_index}.fa"
        out_path   = self.outdir / f"dmnd_{tag}_batch{batch_index}.tsv"
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

        logger.debug(f"diamond [{tag}] batch {batch_index} — {len(batch)} seqs")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"diamond [{tag}] batch {batch_index} failed: {result.stderr.strip()}")
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    def run_all(
        self,
        nuc_seqs: list[NucSequence],
        phase: DiamondPhase,
        max_workers: int = 1,
    ) -> list[Path]:
        if not nuc_seqs:
            logger.warning(f"diamond [{phase.value}]: no sequences.")
            return []

        tag = phase.value
        batches = self._chunk(nuc_seqs, self.batch_size)
        logger.info(
            f"diamond [{tag}]: {len(nuc_seqs)} seqs → "
            f"{len(batches)} batch(es), {self.threads} threads"
        )

        tsv_paths: list[Path | None] = [None] * len(batches)

        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(self._run_batch, batch, i, tag): i
                for i, batch in enumerate(batches)
            }
            for future in as_completed(futures):
                i = futures[future]
                try:
                    tsv_paths[i] = future.result()
                    logger.success(f"diamond [{tag}] batch {i} done")
                except Exception as exc:
                    logger.error(f"diamond [{tag}] batch {i} failed: {exc}")

        return [p for p in tsv_paths if p is not None]


class DiamondOutputParser:
    """Parses Diamond TSV (13-column) into BlastxResult objects."""

    @staticmethod
    def parse(tsv_path: Path) -> list[BlastxResult]:
        hits: list[BlastxResult] = []
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return hits

        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 13:
                    continue
                try:
                    hits.append(BlastxResult(
                        query_id=row[0],
                        subject_id=row[1],
                        subject_title=row[2],
                        pct_identity=float(row[3]),
                        aln_length=int(row[4]),
                        mismatches=int(row[5]),
                        gap_opens=int(row[6]),
                        query_start=int(row[7]),
                        query_end=int(row[8]),
                        subject_start=int(row[9]),
                        subject_end=int(row[10]),
                        e_value=float(row[11]),
                        bit_score=float(row[12]),
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(f"Skipping malformed row in {tsv_path.name}: {exc}")
        return hits


class DiamondResultAttacher:
    """Attaches BlastxResult hits to NucSequence objects."""

    @staticmethod
    def attach(
        hits: list[BlastxResult],
        nuc_seqs: list[NucSequence],
        phase: DiamondPhase,
    ) -> int:
        seq_map = {seq.id: seq for seq in nuc_seqs}
        attached = 0
        for hit in hits:
            seq = seq_map.get(hit.query_id)
            if seq is None:
                continue
            hit.compute_coverage(seq.length)
            if phase == DiamondPhase.REFSEQ_FILTER:
                seq.blastx_hits.append(hit)
                seq.is_viral = True          # refseq db is already viral-only
            else:
                # phase-2 NR: only attach if the subject title looks viral
                if is_viral_subject(hit.subject_title):
                    seq.blastx_nr_hits.append(hit)
            attached += 1
        logger.debug(f"diamond [{phase.value}]: {attached} hits processed.")
        return attached


class DiamondFilterer:
    """Post-hoc quality filtering of BlastxResult lists."""

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