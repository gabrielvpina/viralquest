import csv
import re
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Generator
from enum import Enum
from pathlib import Path

from loguru import logger
from .biodata import NucSequence, BlastxResult, merge_query_intervals


class DiamondPhase(Enum):
    REFSEQ_FILTER   = "refseq"  # phase-1: small viral RefSeq → select candidates
    NR_CHARACTERIZE = "nr"      # phase-2: NR global → full characterisation


# Matches viral organism names inside brackets in NCBI subject titles,
# e.g. "[Tomato bushy stunt virus]" or "[Escherichia phage T4]"
_VIRAL_PATTERN = re.compile(
    r"\[[^\]]*(?:vir(?:us|al|idae|ales|iform|oid)|phage|bacteriophage)[^\]]*\]",
    re.IGNORECASE,
)


def is_viral_subject(title: str) -> bool:
    """True if the subject title contains a viral organism name in brackets."""
    return bool(_VIRAL_PATTERN.search(title))


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

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
        batch_size: int = 5000,
        block_size: float | None = None,
        index_chunks: int | None = None,
        tmpdir: str | None = None,
    ):
        self.db_path         = db_path
        self.diamond_bin     = diamond_bin
        self.threads         = threads
        self.e_value         = e_value
        self.max_target_seqs = max_target_seqs
        self.outdir          = Path(outdir) if outdir else Path(tempfile.gettempdir())
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.batch_size      = batch_size
        self.block_size      = block_size
        self.index_chunks    = index_chunks
        self.tmpdir          = tmpdir

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
        if self.block_size is not None:
            cmd += ["--block-size", str(self.block_size)]
        if self.index_chunks is not None:
            cmd += ["--index-chunks", str(self.index_chunks)]
        if self.tmpdir is not None:
            cmd += ["--tmpdir", str(self.tmpdir)]

        logger.debug(f"diamond [{tag}] batch {batch_index} — {len(batch)} seqs")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(
                f"diamond [{tag}] batch {batch_index} failed: {result.stderr.strip()}"
            )
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    def run_batched(
        self,
        nuc_seqs: list[NucSequence],
        phase: DiamondPhase,
    ) -> Generator[tuple[int, int, Path], None, None]:
        """
        Sequential batch runner. Yields (batch_num_1based, total_batches, tsv_path)
        after each batch completes so callers can show progress.
        """
        if not nuc_seqs:
            logger.warning(f"diamond [{phase.value}]: no sequences.")
            return
        tag     = phase.value
        batches = self._chunk(nuc_seqs, self.batch_size)
        n       = len(batches)
        logger.info(
            f"diamond [{tag}]: {len(nuc_seqs)} seqs → {n} batch(es), {self.threads} threads"
        )
        for i, batch in enumerate(batches):
            path = self._run_batch(batch, i, tag)
            logger.success(f"diamond [{tag}] batch {i + 1}/{n} done")
            yield i + 1, n, path

    def run_single(
        self,
        nuc_seqs: list[NucSequence],
        phase: DiamondPhase,
    ) -> Path | None:
        """Run all sequences in a single Diamond call — no batching."""
        if not nuc_seqs:
            logger.warning(f"diamond [{phase.value}]: no sequences.")
            return None
        tag = phase.value
        logger.info(
            f"diamond [{tag}]: {len(nuc_seqs)} seqs → single run, {self.threads} threads"
        )
        path = self._run_batch(nuc_seqs, 0, tag)
        logger.success(f"diamond [{tag}] done")
        return path


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Attacher  (merged-interval coverage)
# ---------------------------------------------------------------------------

class DiamondResultAttacher:
    """
    Attaches BlastxResult hits to NucSequence objects and computes
    merged query coverage across all HSPs for each query sequence.

    Coverage logic:
      - Each hit gets its own single-interval coverage via compute_coverage().
      - The best hit (highest bit_score) additionally receives the *merged*
        coverage across all hits for that query, which is what the reporter
        uses for the JSON output.  This avoids double-counting overlapping
        aligned regions when the same query matches the subject in multiple
        non-contiguous segments.
    """

    @staticmethod
    def attach(
        hits: list[BlastxResult],
        nuc_seqs: list[NucSequence],
        phase: DiamondPhase,
    ) -> int:
        seq_map = {seq.id: seq for seq in nuc_seqs}

        # group hits by query_id so we can compute merged coverage once
        hits_by_query: dict[str, list[BlastxResult]] = defaultdict(list)
        for hit in hits:
            hits_by_query[hit.query_id].append(hit)

        attached = 0

        for query_id, query_hits in hits_by_query.items():
            seq = seq_map.get(query_id)
            if seq is None:
                continue

            # per-hit single-interval coverage
            for hit in query_hits:
                hit.compute_coverage(seq.length)

            # merged coverage across all HSPs → assigned to the best hit
            merged_bases    = merge_query_intervals(query_hits)
            merged_coverage = (
                round((merged_bases / seq.length) * 100, 2) if seq.length > 0 else 0.0
            )
            best_hit = max(query_hits, key=lambda h: h.bit_score)
            best_hit.query_coverage = merged_coverage

            if phase == DiamondPhase.REFSEQ_FILTER:
                seq.blastx_hits.extend(query_hits)
                seq.is_viral = True          # RefSeq db is already viral-only

            else:
                # NR phase: only keep hits whose subject title looks viral
                viral_hits = [h for h in query_hits if is_viral_subject(h.subject_title)]
                if viral_hits:
                    seq.blastx_nr_hits.extend(viral_hits)

            attached += len(query_hits)

        logger.debug(f"diamond [{phase.value}]: {attached} hits processed.")
        return attached


# ---------------------------------------------------------------------------
# Filterer
# ---------------------------------------------------------------------------

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
            if h.pct_identity  >= min_identity
            and h.e_value      <= max_evalue
            and h.query_coverage >= min_coverage
        ]