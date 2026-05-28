import csv
import subprocess
import tempfile
import time
from collections import defaultdict
from enum import Enum
from pathlib import Path

from loguru import logger

from viralquest.biodata import BlastnResult, NucSequence


class BlastnMode(Enum):
    LOCAL  = "local"   # subprocess against a local BLAST+ nucleotide DB
    ONLINE = "online"  # NCBI qblast via Biopython NCBIWWW (no local DB needed)


# outfmt columns for local mode (order matters for the parser)
_OUTFMT = (
    "6 qseqid sseqid stitle pident length qlen slen qcovhsp evalue bitscore sacc qstart qend"
)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

class BlastnRunner:
    """
    Runs BLASTn on viral-flagged NucSequence objects only (is_viral=True).

    Modes
    -----
    LOCAL  — subprocess blastn against a local nucleotide database.
             Sequences are submitted in batches for efficiency.
    ONLINE — NCBI qblast via Biopython (no local installation needed).
             Submits one sequence at a time with a configurable delay to
             respect NCBI rate limits (3 req/s without key, 10 with key).

    Results are returned as a flat list of BlastnResult objects.
    Use BlastnResultAttacher to attach them to NucSequence objects.

    Parameters
    ----------
    mode            : BlastnMode.LOCAL or BlastnMode.ONLINE
    db_path         : path to local BLAST DB prefix (LOCAL only, required)
    blastn_bin      : blastn executable name or full path (LOCAL only)
    threads         : -num_threads value passed to blastn (LOCAL only)
    batch_size      : sequences per blastn subprocess call (LOCAL only)
    ncbi_api_key    : NCBI API key — raises rate limit to 10 req/s (ONLINE only)
    request_delay   : seconds between NCBI requests; default 0.4 s (ONLINE only)
    e_value         : e-value cutoff (both modes)
    max_target_seqs : maximum hits returned per query (both modes)
    outdir          : directory for temporary files (LOCAL only)
    """

    def __init__(
        self,
        mode: BlastnMode = BlastnMode.LOCAL,
        db_path: str | None = None,
        blastn_bin: str = "blastn",
        threads: int = 2,
        e_value: float = 1e-5,
        max_target_seqs: int = 5,
        outdir: str | None = None,
        batch_size: int = 200,
        ncbi_api_key: str | None = None,
        request_delay: float = 0.4,
    ):
        if mode == BlastnMode.LOCAL and not db_path:
            raise ValueError("db_path is required for BlastnMode.LOCAL.")

        self.mode            = mode
        self.db_path         = db_path
        self.blastn_bin      = blastn_bin
        self.threads         = threads
        self.e_value         = e_value
        self.max_target_seqs = max_target_seqs
        self.outdir          = Path(outdir) if outdir else Path(tempfile.gettempdir()) / "vq_blastn"
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.batch_size      = batch_size
        self.ncbi_api_key    = ncbi_api_key
        self.request_delay   = request_delay

    # --- public ---

    def run(self, nuc_seqs: list[NucSequence]) -> list[BlastnResult]:
        """
        Filter to viral sequences, run BLASTn, and return all hits.
        Attaching results to NucSequence objects is done separately by
        BlastnResultAttacher.
        """
        if not nuc_seqs:
            logger.warning("BLASTn: no sequences to search — skipping.")
            return []

        logger.info(
            f"BLASTn [{self.mode.value}]: searching {len(nuc_seqs)} sequence(s)."
        )

        if self.mode == BlastnMode.LOCAL:
            return self._run_local(nuc_seqs)
        return self._run_online(nuc_seqs)

    # --- local mode ----------------------------------------------------------

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
            self.blastn_bin,
            "-db",              self.db_path,
            "-query",           str(fasta_path),
            "-out",             str(out_path),
            "-outfmt",          _OUTFMT,
            "-num_threads",     str(self.threads),
            "-evalue",          str(self.e_value),
            "-max_target_seqs", str(self.max_target_seqs),
        ]

        logger.debug(f"BLASTn [local] batch {batch_index} — {len(batch)} seqs")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(
                f"BLASTn [local] batch {batch_index} failed "
                f"(exit {result.returncode}): {result.stderr.strip()}"
            )
            raise RuntimeError(result.stderr.strip())

        fasta_path.unlink(missing_ok=True)
        return out_path

    def _run_local(self, seqs: list[NucSequence]) -> list[BlastnResult]:
        batches  = self._chunk(seqs, self.batch_size)
        all_hits: list[BlastnResult] = []

        logger.info(
            f"BLASTn [local]: {len(batches)} batch(es), "
            f"db='{self.db_path}', threads={self.threads}"
        )

        for i, batch in enumerate(batches):
            try:
                tsv_path = self._run_batch(batch, i)
                hits     = BlastnOutputParser.parse_tsv(tsv_path)
                all_hits.extend(hits)
                logger.success(f"BLASTn [local] batch {i}: {len(hits)} hit(s)")
            except Exception as exc:
                logger.error(f"BLASTn [local] batch {i} failed: {exc}")

        return all_hits

    # --- online mode ---------------------------------------------------------

    def _run_online(self, seqs: list[NucSequence]) -> list[BlastnResult]:
        try:
            from Bio.Blast import NCBIWWW, NCBIXML
        except ImportError:
            logger.error(
                "Biopython is required for online BLASTn. "
                "Install it with: pip install biopython"
            )
            return []

        all_hits: list[BlastnResult] = []

        for i, seq in enumerate(seqs):
            logger.info(
                f"BLASTn [online] [{i + 1}/{len(seqs)}] querying '{seq.id}' ..."
            )
            try:
                kwargs: dict = dict(
                    program="blastn",
                    database="nt",
                    sequence=seq.sequence,
                    hitlist_size=self.max_target_seqs,
                    expect=self.e_value,
                    format_type="XML",
                )
                if self.ncbi_api_key:
                    kwargs["api_key"] = self.ncbi_api_key

                result_handle = NCBIWWW.qblast(**kwargs)
                blast_record  = next(NCBIXML.parse(result_handle))
                hits = BlastnOutputParser.parse_xml_record(blast_record, seq.id)
                all_hits.extend(hits)
                logger.success(f"BLASTn [online] '{seq.id}': {len(hits)} hit(s)")

            except Exception as exc:
                logger.error(f"BLASTn [online] '{seq.id}' failed: {exc}")

            # rate-limit: skip delay after the last sequence
            if i < len(seqs) - 1:
                time.sleep(self.request_delay)

        return all_hits


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

class BlastnOutputParser:
    """Parses BLASTn output (local TSV or online XML) into BlastnResult objects."""

    # Column indices for the custom outfmt defined in _OUTFMT
    # 0:qseqid  1:sseqid  2:stitle  3:pident  4:length
    # 5:qlen    6:slen    7:qcovhsp 8:evalue  9:bitscore  10:sacc  11:qstart  12:qend
    _QSEQID   = 0
    _STITLE   = 2
    _PIDENT   = 3
    _QLEN     = 5
    _SLEN     = 6
    _QCOVHSP  = 7
    _EVALUE   = 8
    _BITSCORE = 9
    _SACC     = 10
    _QSTART   = 11
    _QEND     = 12

    @classmethod
    def parse_tsv(cls, tsv_path: Path) -> list[BlastnResult]:
        """Parse a local blastn outfmt-6 TSV file."""
        hits: list[BlastnResult] = []
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return hits

        with open(tsv_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 10:
                    continue
                try:
                    hits.append(BlastnResult(
                        qseqid      = row[cls._QSEQID],
                        qlen        = int(row[cls._QLEN]),
                        slen        = int(row[cls._SLEN]),
                        qcovhsp     = int(float(row[cls._QCOVHSP])),
                        pident      = float(row[cls._PIDENT]),
                        evalue      = float(row[cls._EVALUE]),
                        bit_score   = float(row[cls._BITSCORE]),
                        stitle      = row[cls._STITLE],
                        accession   = row[cls._SACC]   if len(row) > cls._SACC  else None,
                        query_start = int(row[cls._QSTART]) if len(row) > cls._QSTART else None,
                        query_end   = int(row[cls._QEND])   if len(row) > cls._QEND   else None,
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(
                        f"Skipping malformed BLASTn row in {tsv_path.name}: {exc}"
                    )
        return hits

    @staticmethod
    def parse_xml_record(blast_record, query_id: str) -> list[BlastnResult]:
        """
        Parse one Biopython NCBIXML blast_record into BlastnResult objects.

        For each alignment, the best HSP (highest bit score) is selected.
        query_coverage is computed from the HSP query coordinates.
        """
        hits: list[BlastnResult] = []
        qlen = blast_record.query_length or 0

        for alignment in blast_record.alignments:
            if not alignment.hsps:
                continue

            best_hsp = max(alignment.hsps, key=lambda h: h.bits)

            aln_span = best_hsp.query_end - best_hsp.query_start + 1
            qcovhsp  = int(aln_span / qlen * 100) if qlen > 0 else 0
            pident   = (
                round(best_hsp.identities / best_hsp.align_length * 100, 2)
                if best_hsp.align_length > 0 else 0.0
            )
            # hit_def is the description without the accession prefix
            stitle = alignment.hit_def or alignment.title.lstrip(">").strip()

            try:
                hits.append(BlastnResult(
                    qseqid      = query_id,
                    qlen        = qlen,
                    slen        = alignment.length,
                    qcovhsp     = qcovhsp,
                    pident      = pident,
                    evalue      = best_hsp.expect,
                    bit_score   = float(best_hsp.bits),
                    stitle      = stitle,
                    accession   = getattr(alignment, "accession", None) or None,
                    query_start = best_hsp.query_start,
                    query_end   = best_hsp.query_end,
                ))
            except Exception as exc:
                logger.warning(
                    f"Skipping alignment '{alignment.title[:60]}': {exc}"
                )

        return hits


# ---------------------------------------------------------------------------
# Attacher
# ---------------------------------------------------------------------------

class BlastnResultAttacher:
    """
    Attaches BlastnResult hits to NucSequence objects.

    Hits are sorted by bit_score descending before attaching, so
    NucSequence.best_blastn always returns the top hit.

    Parameters
    ----------
    max_hits_per_query : maximum hits to keep per sequence (default 5)
    """

    @staticmethod
    def attach(
        hits: list[BlastnResult],
        nuc_seqs: list[NucSequence],
        max_hits_per_query: int = 5,
    ) -> int:
        """
        Attach hits to the matching NucSequence objects.
        Returns the total number of hits attached.
        """
        seq_map: dict[str, NucSequence]       = {s.id: s for s in nuc_seqs}
        hits_by_query: dict[str, list[BlastnResult]] = defaultdict(list)

        for hit in hits:
            hits_by_query[hit.qseqid].append(hit)

        attached = 0
        for query_id, query_hits in hits_by_query.items():
            seq = seq_map.get(query_id)
            if seq is None:
                continue
            ranked = sorted(query_hits, key=lambda h: h.bit_score, reverse=True)
            kept   = ranked[:max_hits_per_query]
            seq.blastn_hits.extend(kept)
            attached += len(kept)

        logger.debug(f"BLASTn: {attached} hit(s) attached across {len(hits_by_query)} sequence(s).")
        return attached
