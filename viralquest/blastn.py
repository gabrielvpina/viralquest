import csv
import io
import re
import subprocess
import tempfile
import time
from collections import defaultdict
from enum import Enum
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from loguru import logger

from viralquest.biodata import BlastnResult, NucSequence


class BlastnMode(Enum):
    LOCAL  = "local"   # subprocess against a local BLAST+ nucleotide DB
    ONLINE = "online"  # NCBI BLAST URL API (no local DB needed)


# outfmt columns for local mode (order matters for the parser)
_OUTFMT = (
    "6 qseqid sseqid stitle pident length qlen slen qcovhsp evalue bitscore sacc qstart qend"
)


# ---------------------------------------------------------------------------
# NCBI BLAST URL API client (online mode)
# ---------------------------------------------------------------------------

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


class OnlineSearchTimeout(RuntimeError):
    """An NCBI search did not finish within its time limit."""


class NcbiBlastClient:
    """
    Minimal client for NCBI's BLAST URL API (CMD=Put / CMD=Get) with a time
    limit per search. Biopython's qblast polls until the search finishes, with
    no limit, and loops forever on a FAILED status; this client gives up after
    `timeout` seconds (raising OnlineSearchTimeout) and raises on FAILED /
    UNKNOWN at once.

    NCBI usage rules are kept, as in qblast:
      * >= MIN_CONTACT_GAP s between any two requests to the server;
      * first status poll after FIRST_POLL s, then every POLL_EVERY s
        (never more often than once a minute per search);
      * `email` and `tool` sent with every submission.
    https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html

    `http`, `clock` and `sleep` are injectable for tests.
    """

    MIN_CONTACT_GAP = 10.0
    FIRST_POLL      = 20.0
    POLL_EVERY      = 60.0

    _RID_RE    = re.compile(r"^\s*RID = (\S+)", re.MULTILINE)
    _STATUS_RE = re.compile(r"Status=(\w+)")

    def __init__(self, email: str | None = None, tool: str = "viralquest",
                 url: str = NCBI_BLAST_URL, http=None, clock=None, sleep=None):
        self.email = email
        self.tool  = tool
        self.url   = url
        # Resolved at call time (not bound as defaults), so tests can patch
        # time.sleep / time.monotonic module-wide.
        self._http  = http or self._urlopen
        self._clock = clock or (lambda: time.monotonic())
        self._sleep = sleep or (lambda s: time.sleep(s))
        self._last  = float("-inf")       # time of the last request (any search)

    # -- transport --------------------------------------------------------
    def _urlopen(self, params: dict) -> str:
        data = urlencode(params).encode()
        req = Request(self.url, data, {"User-Agent": f"{self.tool} (BLAST URL API client)"})
        with urlopen(req, timeout=120) as resp:
            return resp.read().decode()

    def _request(self, params: dict) -> str:
        wait = self._last + self.MIN_CONTACT_GAP - self._clock()
        if wait > 0:
            self._sleep(wait)
        self._last = self._clock()
        return self._http(params)

    # -- API --------------------------------------------------------------
    def submit(self, **put) -> str:
        """CMD=Put — send the search, return its RID."""
        params = {k.upper(): v for k, v in put.items() if v is not None}
        params["CMD"] = "Put"
        if self.email:
            params["email"] = self.email
        params["tool"] = self.tool
        page = self._request(params)
        m = self._RID_RE.search(page)
        if not m:
            snippet = " ".join(re.sub(r"<[^>]+>", " ", page).split())[:200]
            raise RuntimeError(f"NCBI did not accept the search: {snippet or 'empty reply'}")
        return m.group(1)

    def search(self, timeout: float, **put) -> str:
        """
        Submit and wait for the XML result. Raises OnlineSearchTimeout when the
        next poll would fall past `timeout` seconds after submission.
        """
        rid = self.submit(**put)
        start = self._clock()
        get = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}
        delay = self.FIRST_POLL
        while True:
            if self._clock() - start + delay > timeout:
                raise OnlineSearchTimeout(
                    f"search {rid} not finished after {self._clock() - start:.0f} s "
                    f"(limit {timeout:.0f} s)")
            self._sleep(delay)
            page = self._request(get)
            delay = self.POLL_EVERY
            if not page.strip():
                continue                                  # in progress, empty page
            m = self._STATUS_RE.search(page)
            if not m:
                return page                               # finished: the XML itself
            status = m.group(1).upper()
            if status == "WAITING":
                continue
            if status == "READY":
                delay = 0.0                               # fetch the XML next round
                continue
            raise RuntimeError(f"NCBI search {rid} ended with status {status}")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

class BlastnRunner:
    """
    Runs BLASTn on the NucSequence objects it is given. It does not filter:
    the caller (cli.py) selects them — NR-confirmed sequences with --nr-db,
    is_viral sequences otherwise.

    Modes
    -----
    LOCAL  — subprocess blastn against a local nucleotide database.
             Sequences are submitted in batches for efficiency.
    ONLINE — NCBI BLAST URL API via NcbiBlastClient (no local installation needed).
             Sequences are sent in batches (one multi-FASTA search per
             `online_batch_size` sequences), as NCBI recommends: fewer searches
             run faster and stay under the 100-searches-per-day threshold
             above which NCBI moves a user to a slower queue. NCBI's pacing
             rules are kept (see NcbiBlastClient). A search that exceeds
             `online_timeout` is resubmitted up to `online_retries` times; a
             batch that still times out, or fails, is split in half and
             retried, so one bad query cannot sink the others.

    Results are returned as a flat list of BlastnResult objects.
    Use BlastnResultAttacher to attach them to NucSequence objects.

    Parameters
    ----------
    mode            : BlastnMode.LOCAL or BlastnMode.ONLINE
    db_path         : path to local BLAST DB prefix (LOCAL only, required)
    blastn_bin      : blastn executable name or full path (LOCAL only)
    threads         : -num_threads value passed to blastn (LOCAL only)
    batch_size      : sequences per blastn subprocess call (LOCAL only)
    online_batch_size : sequences per NCBI search; default 20 (ONLINE only)
    online_timeout  : seconds one NCBI search may take before it is resubmitted;
                      default 900 (15 min) (ONLINE only)
    online_retries  : resubmissions of a timed-out search before the batch is
                      split; default 2 (ONLINE only)
    online_db       : NCBI nucleotide database to query; default "nt" (ONLINE only)
    email           : contact e-mail sent to NCBI with each request (ONLINE only)
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
        online_batch_size: int = 20,
        online_timeout: float = 900.0,
        online_retries: int = 2,
        online_db: str = "nt",
        email: str | None = None,
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
        self.online_db       = online_db
        self.email           = email
        self.failed_ids: list[str] = []   # queries whose search errored (reset per run)
        self.online_batch_size = max(1, int(online_batch_size))
        self.online_timeout    = float(online_timeout)
        self.online_retries    = max(0, int(online_retries))
        self.ncbi              = NcbiBlastClient(email=email)

    # --- public ---

    def run(self, nuc_seqs: list[NucSequence]) -> list[BlastnResult]:
        """
        Run BLASTn on every given sequence and return all hits.
        Queries that error are listed in ``self.failed_ids`` (reset per call).
        Attaching results to NucSequence objects is done separately by
        BlastnResultAttacher.
        """
        self.failed_ids = []
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
                self.failed_ids.extend(s.id for s in batch)

        return all_hits

    # --- online mode ---------------------------------------------------------

    def _run_online(self, seqs: list[NucSequence]) -> list[BlastnResult]:
        try:
            from Bio.Blast import NCBIXML  # noqa: F401  (imported here: optional dependency)
        except ImportError:
            logger.error(
                "Biopython is required for online BLASTn. "
                "Install it with: pip install biopython"
            )
            self.failed_ids = [s.id for s in seqs]
            return []

        batches = self._chunk(seqs, self.online_batch_size)
        logger.info(
            f"BLASTn [online]: {len(seqs)} sequence(s) in {len(batches)} search(es) "
            f"of up to {self.online_batch_size} against '{self.online_db}' "
            f"(each search: >=1 min; limit {self.online_timeout / 60:.0f} min, "
            f"{self.online_retries} resubmission(s))."
        )
        all_hits: list[BlastnResult] = []
        for i, batch in enumerate(batches, 1):
            label = f"batch {i}/{len(batches)}"
            all_hits.extend(self._online_batch(batch, label))
        return all_hits

    def _online_batch(self, batch: list[NucSequence], label: str) -> list[BlastnResult]:
        """
        One NCBI search for a multi-FASTA batch. On failure the batch is split in
        half and each half retried, down to single sequences, which are then
        recorded in failed_ids.
        """
        from Bio.Blast import NCBIXML

        total_bp = sum(len(s.sequence) for s in batch)
        logger.info(f"BLASTn [online] {label}: {len(batch)} seq(s), {total_bp:,} bp — submitting ...")
        try:
            fasta = "".join(f">{s.id}\n{s.sequence}\n" for s in batch)
            xml = self._search_with_resubmission(fasta, label)
            records = list(NCBIXML.parse(io.StringIO(xml)))
        except Exception as exc:
            if len(batch) == 1:
                logger.error(f"BLASTn [online] '{batch[0].id}' failed: {exc}")
                self.failed_ids.append(batch[0].id)
                return []
            half = len(batch) // 2
            logger.warning(
                f"BLASTn [online] {label} failed ({exc}) — splitting into "
                f"{half} + {len(batch) - half} and retrying."
            )
            return (self._online_batch(batch[:half], f"{label}a")
                    + self._online_batch(batch[half:], f"{label}b"))

        hits = self._attach_records(records, batch)
        logger.success(f"BLASTn [online] {label}: {len(hits)} hit(s) across {len(batch)} seq(s)")
        return hits

    def _search_with_resubmission(self, fasta: str, label: str) -> str:
        """
        Run one NCBI search; when it exceeds online_timeout, abandon it and
        resubmit (up to online_retries times). The last timeout is re-raised so
        the caller splits the batch.
        """
        attempts = self.online_retries + 1
        for attempt in range(1, attempts + 1):
            try:
                return self.ncbi.search(
                    self.online_timeout,
                    program="blastn",
                    database=self.online_db,
                    query=fasta,
                    hitlist_size=self.max_target_seqs,
                    expect=self.e_value,
                )
            except OnlineSearchTimeout as exc:
                if attempt == attempts:
                    raise
                logger.warning(
                    f"BLASTn [online] {label}: {exc} — resubmitting "
                    f"(attempt {attempt + 1}/{attempts})."
                )
        raise AssertionError("unreachable")

    def _attach_records(self, records: list, batch: list[NucSequence]) -> list[BlastnResult]:
        """
        Map each XML record back to its query. NCBI echoes the FASTA defline
        (the sequence id) in <Iteration_query-def>; records come back in input
        order, which is the fallback when the id cannot be read. A query with
        no record at all is counted as failed.
        """
        ids = {s.id for s in batch}
        by_id: dict[str, object] = {}
        for pos, rec in enumerate(records):
            qid = (getattr(rec, "query", "") or "").split()[0] if getattr(rec, "query", "") else ""
            if qid not in ids and pos < len(batch):
                qid = batch[pos].id
            by_id.setdefault(qid, rec)

        hits: list[BlastnResult] = []
        for s in batch:
            rec = by_id.get(s.id)
            if rec is None:
                logger.error(f"BLASTn [online] '{s.id}': no result returned for this query.")
                self.failed_ids.append(s.id)
                continue
            hits.extend(BlastnOutputParser.parse_xml_record(rec, s.id))
        return hits


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

    # NCBI sequence-id database tags — the accession follows one of these in the
    # pipe-delimited Hit_id (e.g. "gi|123|ref|NC_045512.2|").
    _DB_TAGS = frozenset({
        "ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd",
        "sp", "pir", "prf", "pdb", "gnl", "lcl",
    })

    @classmethod
    def _accession_from_alignment(cls, alignment) -> str | None:
        """
        Resolve the subject accession for one alignment.

        Prefer the explicit ``<Hit_accession>`` (present in local BLAST+/XML v2
        output). Online (NCBI BLAST URL API) XML omits it, exposing the accession only
        inside ``Hit_id`` (e.g. ``gi|1798174254|ref|NC_045512.2|``), so fall back
        to parsing it from there.
        """
        acc = getattr(alignment, "accession", None)
        if acc:
            return acc

        hit_id = getattr(alignment, "hit_id", "") or ""
        parts  = [p for p in hit_id.split("|") if p]
        if not parts:
            return hit_id or None
        for i, part in enumerate(parts):
            if part in cls._DB_TAGS and i + 1 < len(parts):
                return parts[i + 1]
        return parts[-1]

    @classmethod
    def parse_xml_record(cls, blast_record, query_id: str) -> list[BlastnResult]:
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
                    accession   = cls._accession_from_alignment(alignment),
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
