import csv
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from viralquest.biodata import BlastnResult, BlastxResult, NucSequence
from viralquest.blastn import (
    NcbiBlastClient,
    OnlineSearchTimeout,
    BlastnMode,
    BlastnOutputParser,
    BlastnResultAttacher,
    BlastnRunner,
    BlastnSubjectFetcher,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_seq(
    seq_id: str = "seq1",
    sequence: str = "ATGC" * 300,
    is_viral: bool = False,
    with_nr_hit: bool = False,
) -> NucSequence:
    seq = NucSequence(id=seq_id, sequence=sequence)
    seq.is_viral = is_viral
    if with_nr_hit:
        seq.blastx_nr_hits.append(
            BlastxResult(
                query_id=seq_id, subject_id="NR_001",
                subject_title=f"polymerase [Influenza A virus]",
                pct_identity=92.0, aln_length=300, mismatches=8,
                gap_opens=0, query_start=1, query_end=300,
                subject_start=1, subject_end=300,
                e_value=1e-90, bit_score=450.0,
            )
        )
    return seq


def _tsv_row(
    qseqid: str  = "seq1",
    sseqid: str  = "NC_002017",
    stitle: str  = "Influenza A virus segment 1 RNA polymerase PB2",
    pident: str  = "92.5",
    length: str  = "1200",
    qlen: str    = "1200",
    slen: str    = "13600",
    qcovhsp: str = "80",
    evalue: str  = "1e-60",
    bitscore: str = "200.0",
) -> str:
    return "\t".join([qseqid, sseqid, stitle, pident, length, qlen, slen,
                      qcovhsp, evalue, bitscore])


def _write_tsv(tmp_path: Path, rows: list[str]) -> Path:
    p = tmp_path / "out.tsv"
    p.write_text("\n".join(rows), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# BlastnMode
# ---------------------------------------------------------------------------

class TestBlastnMode:
    def test_values(self):
        assert BlastnMode.LOCAL.value  == "local"
        assert BlastnMode.ONLINE.value == "online"


# ---------------------------------------------------------------------------
# BlastnRunner  — init
# ---------------------------------------------------------------------------

class TestBlastnRunnerInit:
    def test_local_without_db_path_raises(self):
        with pytest.raises(ValueError, match="db_path"):
            BlastnRunner(mode=BlastnMode.LOCAL)

    def test_local_with_db_path_succeeds(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path=str(tmp_path / "db"))
        assert runner.db_path == str(tmp_path / "db")

    def test_online_without_db_path_succeeds(self):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        assert runner.mode == BlastnMode.ONLINE

    def test_default_params(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/db")
        assert runner.blastn_bin      == "blastn"
        assert runner.threads         == 2
        assert runner.e_value         == 1e-5
        assert runner.max_target_seqs == 5
        assert runner.batch_size      == 200
        assert runner.online_batch_size == 20

    def test_outdir_created(self, tmp_path):
        outdir = tmp_path / "blast_out"
        BlastnRunner(mode=BlastnMode.LOCAL, db_path="/db", outdir=str(outdir))
        assert outdir.exists()


# ---------------------------------------------------------------------------
# BlastnRunner.run  — dispatch (sequence selection happens in cli.py)
# ---------------------------------------------------------------------------

class TestBlastnRunnerRun:
    def _runner(self, tmp_path) -> BlastnRunner:
        return BlastnRunner(
            mode=BlastnMode.LOCAL,
            db_path="/fake/db",
            outdir=str(tmp_path),
        )

    def test_passes_every_given_sequence(self, tmp_path):
        # The runner no longer filters: cli.py picks the sequences
        # (NR-confirmed with --nr-db, is_viral otherwise).
        seqs = [
            make_seq("v1", is_viral=True,  with_nr_hit=True),
            make_seq("v2", is_viral=True,  with_nr_hit=False),
            make_seq("v3", is_viral=False, with_nr_hit=False),
        ]
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local", return_value=[]) as mock:
            runner.run(seqs)
        mock.assert_called_once()
        assert [s.id for s in mock.call_args[0][0]] == ["v1", "v2", "v3"]

    def test_empty_input_returns_empty(self, tmp_path):
        runner = self._runner(tmp_path)
        assert runner.run([]) == []

    def test_empty_input_does_not_search(self, tmp_path):
        runner = self._runner(tmp_path)
        with patch.object(runner, "_run_local") as mock:
            runner.run([])
        mock.assert_not_called()

    def test_routes_online_mode(self, tmp_path):
        seq = make_seq(is_viral=True, with_nr_hit=True)
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path))
        with patch.object(runner, "_run_online", return_value=[]) as mock:
            runner.run([seq])
        mock.assert_called_once()


# ---------------------------------------------------------------------------
# BlastnRunner — failed query tracking
# ---------------------------------------------------------------------------

class TestBlastnFailedIds:
    def test_local_failed_batch_records_its_ids(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/fake/db",
                              outdir=str(tmp_path), batch_size=1)
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        with patch.object(runner, "_run_batch",
                          side_effect=[RuntimeError("boom"), empty]):
            runner.run([make_seq("a"), make_seq("b")])
        assert runner.failed_ids == ["a"]

    def test_local_all_ok_leaves_failed_ids_empty(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/fake/db",
                              outdir=str(tmp_path))
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        with patch.object(runner, "_run_batch", return_value=empty):
            runner.run([make_seq("a")])
        assert runner.failed_ids == []

    def test_failed_ids_reset_between_runs(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/fake/db",
                              outdir=str(tmp_path))
        with patch.object(runner, "_run_batch", side_effect=RuntimeError("x")):
            runner.run([make_seq("a")])
            runner.run([make_seq("b")])
        assert runner.failed_ids == ["b"]

    def test_online_failed_query_records_its_id(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path))

        def boom(timeout, **put):
            raise RuntimeError("server error")

        runner.ncbi.search = boom
        runner.run([make_seq("a"), make_seq("b")])
        assert runner.failed_ids == ["a", "b"]


# ---------------------------------------------------------------------------
# NCBI BLAST URL API — fake server + client tests
# ---------------------------------------------------------------------------

class _FakeRecord:
    """Stand-in for a Biopython NCBIXML record: just the query defline."""
    def __init__(self, query):
        self.query = query


class FakeNcbi:
    """
    Scripted NCBI server for NcbiBlastClient. `polls` is the list of pages the
    Get requests return, in order (e.g. ["Status=WAITING", "<xml>"]). A virtual
    clock advances on every sleep, so timeouts are tested without waiting.
    """
    def __init__(self, polls=("<BlastOutput/>",), put_page="    RID = RID123\n    RTOE = 30\n"):
        self.polls = list(polls)
        self.put_page = put_page
        self.requests = []          # (time, params)
        self.now = 0.0

    def http(self, params):
        self.requests.append((self.now, dict(params)))
        if params["CMD"] == "Put":
            return self.put_page
        return self.polls.pop(0) if self.polls else "<BlastOutput/>"

    def clock(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds

    def client(self, **kw):
        return NcbiBlastClient(http=self.http, clock=self.clock, sleep=self.sleep, **kw)


class TestNcbiBlastClient:
    def test_submit_sends_put_with_email_and_tool(self):
        ncbi = FakeNcbi()
        rid = ncbi.client(email="me@x.org").submit(program="blastn", database="nt", query=">a\nACGT\n")
        put = ncbi.requests[0][1]
        assert rid == "RID123"
        assert put["CMD"] == "Put" and put["DATABASE"] == "nt" and put["PROGRAM"] == "blastn"
        assert put["QUERY"] == ">a\nACGT\n" and put["email"] == "me@x.org" and put["tool"] == "viralquest"

    def test_rejected_submission_raises_with_ncbi_message(self):
        ncbi = FakeNcbi(put_page="<p>Message ID#24 Error: Failed to read the Blast query</p>")
        with pytest.raises(RuntimeError, match="Failed to read the Blast query"):
            ncbi.client().submit(program="blastn", database="nt", query="x")

    def test_waits_until_results_ready(self):
        ncbi = FakeNcbi(polls=["Status=WAITING", "Status=WAITING", "<BlastOutput>done</BlastOutput>"])
        xml = ncbi.client().search(900, program="blastn", database="nt", query="x")
        assert xml == "<BlastOutput>done</BlastOutput>"
        gets = [r for r in ncbi.requests if r[1]["CMD"] == "Get"]
        assert len(gets) == 3 and gets[0][1]["RID"] == "RID123" and gets[0][1]["FORMAT_TYPE"] == "XML"

    def test_respects_ncbi_polling_rules(self):
        ncbi = FakeNcbi(polls=["Status=WAITING", "Status=WAITING", "<BlastOutput/>"])
        ncbi.client().search(900, program="blastn", database="nt", query="x")
        times = [t for t, _ in ncbi.requests]
        gaps = [b - a for a, b in zip(times, times[1:])]
        assert gaps[0] >= 20                      # first poll after 20 s
        assert all(g >= 60 for g in gaps[1:])     # then once a minute per search
        assert all(g >= 10 for g in gaps)         # never < 10 s between contacts

    def test_empty_page_means_still_running(self):
        ncbi = FakeNcbi(polls=["\n\n", "<BlastOutput/>"])
        assert ncbi.client().search(900, program="blastn", database="nt", query="x") == "<BlastOutput/>"

    def test_ready_status_fetches_the_xml(self):
        ncbi = FakeNcbi(polls=["Status=READY", "<BlastOutput/>"])
        assert ncbi.client().search(900, program="blastn", database="nt", query="x") == "<BlastOutput/>"

    @pytest.mark.parametrize("status", ["FAILED", "UNKNOWN"])
    def test_failed_status_raises_immediately(self, status):
        ncbi = FakeNcbi(polls=[f"Status={status}"])
        with pytest.raises(RuntimeError, match=status):
            ncbi.client().search(900, program="blastn", database="nt", query="x")

    def test_timeout_raised_when_next_poll_would_pass_the_limit(self):
        ncbi = FakeNcbi(polls=["Status=WAITING"] * 100)
        with pytest.raises(OnlineSearchTimeout, match="RID123"):
            ncbi.client().search(300, program="blastn", database="nt", query="x")
        # 20 s + 4 x 60 s = 260 s of polling; the next poll (320 s) would pass 300 s
        assert ncbi.now - ncbi.requests[0][0] <= 300


# ---------------------------------------------------------------------------
# BlastnRunner — online database, e-mail and resubmission on timeout
# ---------------------------------------------------------------------------

class TestBlastnOnlineOptions:
    def _put(self, tmp_path, **runner_kw) -> dict:
        ncbi = FakeNcbi()
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path), **runner_kw)
        runner.ncbi = ncbi.client(email=runner_kw.get("email"))
        runner.run([make_seq("a")])
        return ncbi.requests[0][1]

    def test_default_database_is_nt(self, tmp_path):
        assert self._put(tmp_path)["DATABASE"] == "nt"

    def test_chosen_database_is_used(self, tmp_path):
        assert self._put(tmp_path, online_db="refseq_viruses_rep_genomes")["DATABASE"] == "refseq_viruses_rep_genomes"

    def test_email_sent_to_ncbi(self, tmp_path):
        assert self._put(tmp_path, email="me@example.org")["email"] == "me@example.org"

    def test_runner_builds_client_with_email(self, tmp_path):
        r = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path), email="me@x.org")
        assert isinstance(r.ncbi, NcbiBlastClient) and r.ncbi.email == "me@x.org"

    def test_defaults(self, tmp_path):
        r = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path))
        assert (r.online_batch_size, r.online_timeout, r.online_retries) == (20, 900.0, 2)


class TestBlastnOnlineResubmission:
    def _runner(self, tmp_path, outcomes, **kw):
        """outcomes: per search() call, an exception to raise or an XML string."""
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path), **kw)
        calls = []

        def search(timeout, **put):
            calls.append((timeout, put["query"].count(">")))
            out = outcomes.pop(0)
            if isinstance(out, Exception):
                raise out
            return out

        runner.ncbi.search = search
        return runner, calls

    def test_timed_out_search_is_resubmitted(self, tmp_path):
        runner, calls = self._runner(tmp_path, [OnlineSearchTimeout("slow"), "<x/>"], online_timeout=600)
        with patch("Bio.Blast.NCBIXML.parse", return_value=[_FakeRecord("a")]), \
             patch.object(BlastnOutputParser, "parse_xml_record", return_value=["hit"]):
            hits = runner.run([make_seq("a")])
        assert hits == ["hit"] and runner.failed_ids == []
        assert calls == [(600.0, 1), (600.0, 1)]           # same batch, sent twice

    def test_retries_exhausted_then_batch_is_split(self, tmp_path):
        timeouts = [OnlineSearchTimeout("slow")] * 3        # 1 try + 2 retries for the full batch
        runner, calls = self._runner(tmp_path, timeouts + ["<x/>", "<x/>"], online_retries=2)
        with patch("Bio.Blast.NCBIXML.parse", side_effect=[[_FakeRecord("a")], [_FakeRecord("b")]]), \
             patch.object(BlastnOutputParser, "parse_xml_record", side_effect=lambda r, q: [q]):
            hits = runner.run([make_seq("a"), make_seq("b")])
        assert [n for _, n in calls] == [2, 2, 2, 1, 1]
        assert sorted(hits) == ["a", "b"] and runner.failed_ids == []

    def test_zero_retries_splits_at_once(self, tmp_path):
        runner, calls = self._runner(tmp_path, [OnlineSearchTimeout("slow"), "<x/>", "<x/>"], online_retries=0)
        with patch("Bio.Blast.NCBIXML.parse", side_effect=[[_FakeRecord("a")], [_FakeRecord("b")]]), \
             patch.object(BlastnOutputParser, "parse_xml_record", side_effect=lambda r, q: [q]):
            runner.run([make_seq("a"), make_seq("b")])
        assert [n for _, n in calls] == [2, 1, 1]

    def test_single_sequence_timing_out_every_time_is_failed(self, tmp_path):
        runner, calls = self._runner(tmp_path, [OnlineSearchTimeout("slow")] * 3, online_retries=2)
        assert runner.run([make_seq("a")]) == []
        assert runner.failed_ids == ["a"] and len(calls) == 3


# ---------------------------------------------------------------------------
# BlastnRunner — online batching (one multi-FASTA NCBI search per batch)
# ---------------------------------------------------------------------------

class TestBlastnOnlineBatching:
    """NCBI is simulated: the client's search() records the FASTA it receives, NCBIXML.parse
    returns one fake record per query, parse_xml_record tags the hit with the
    query id it was given — so the tests see exactly how records map back."""

    def _run(self, tmp_path, seqs, qblast_effect=None, records_for=None, **kw):
        from Bio.Blast import NCBIXML
        sent = []
        runner = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path), **kw)

        def search(timeout, **put):
            fasta = put["query"]
            sent.append(fasta)
            if qblast_effect:
                qblast_effect(fasta)
            return fasta                                   # the "XML"

        def parse(handle):
            ids = [l[1:] for l in handle.getvalue().splitlines() if l.startswith(">")]
            return [_FakeRecord(q) for q in (records_for(ids) if records_for else ids)]

        def parse_record(rec, query_id):
            return [f"hit:{query_id}<-{rec.query}"]

        runner.ncbi.search = search
        with patch.object(NCBIXML, "parse", side_effect=parse), \
             patch.object(BlastnOutputParser, "parse_xml_record", side_effect=parse_record):
            hits = runner.run(seqs)
        return runner, hits, sent

    def _seqs(self, n):
        return [make_seq(f"s{i}", "ACGT" * 25) for i in range(n)]

    def test_sequences_sent_in_batches_of_20(self, tmp_path):
        _, hits, sent = self._run(tmp_path, self._seqs(45))
        assert [f.count(">") for f in sent] == [20, 20, 5]
        assert len(hits) == 45

    def test_batch_size_is_configurable(self, tmp_path):
        _, _, sent = self._run(tmp_path, self._seqs(10), online_batch_size=4)
        assert [f.count(">") for f in sent] == [4, 4, 2]

    def test_one_search_is_a_multi_fasta(self, tmp_path):
        _, _, sent = self._run(tmp_path, self._seqs(2))
        assert sent == [">s0\n" + "ACGT" * 25 + "\n>s1\n" + "ACGT" * 25 + "\n"]

    def test_records_mapped_by_query_id_even_out_of_order(self, tmp_path):
        _, hits, _ = self._run(tmp_path, self._seqs(3), records_for=lambda ids: ids[::-1])
        assert sorted(hits) == ["hit:s0<-s0", "hit:s1<-s1", "hit:s2<-s2"]

    def test_defline_with_description_still_maps(self, tmp_path):
        _, hits, _ = self._run(tmp_path, self._seqs(2),
                               records_for=lambda ids: [f"{i} len=100 some text" for i in ids])
        assert sorted(h.split("<-")[0] for h in hits) == ["hit:s0", "hit:s1"]

    def test_query_without_record_is_failed(self, tmp_path):
        runner, hits, _ = self._run(tmp_path, self._seqs(3), records_for=lambda ids: ids[:2])
        assert runner.failed_ids == ["s2"]
        assert len(hits) == 2

    def test_failed_batch_is_split_and_retried(self, tmp_path):
        calls = {"n": 0}

        def flaky(fasta):
            calls["n"] += 1
            if calls["n"] == 1:
                raise RuntimeError("CPU usage limit was exceeded")

        runner, hits, sent = self._run(tmp_path, self._seqs(6), qblast_effect=flaky)
        assert [f.count(">") for f in sent] == [6, 3, 3]
        assert runner.failed_ids == [] and len(hits) == 6

    def test_one_bad_query_is_isolated(self, tmp_path):
        def reject_bad(fasta):
            if ">s5\n" in fasta:
                raise RuntimeError("bad query")

        runner, hits, _ = self._run(tmp_path, self._seqs(8), qblast_effect=reject_bad)
        assert runner.failed_ids == ["s5"]
        assert len(hits) == 7


# ---------------------------------------------------------------------------
# BlastnRunner._chunk  (static)
# ---------------------------------------------------------------------------

class TestChunk:
    def _runner(self) -> BlastnRunner:
        return BlastnRunner(mode=BlastnMode.ONLINE)

    def test_exact_multiple(self):
        seqs = [object()] * 6
        chunks = BlastnRunner._chunk(seqs, 3)
        assert len(chunks) == 2
        assert all(len(c) == 3 for c in chunks)

    def test_remainder(self):
        seqs = [object()] * 7
        chunks = BlastnRunner._chunk(seqs, 3)
        assert len(chunks) == 3
        assert len(chunks[-1]) == 1

    def test_single_chunk(self):
        seqs = [object()] * 3
        chunks = BlastnRunner._chunk(seqs, 10)
        assert len(chunks) == 1

    def test_empty(self):
        assert BlastnRunner._chunk([], 5) == []


# ---------------------------------------------------------------------------
# BlastnRunner._write_batch_fasta
# ---------------------------------------------------------------------------

class TestWriteBatchFasta:
    def test_fasta_format(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        seq1 = make_seq("seq1", "ATGCATGC")
        seq2 = make_seq("seq2", "TTTTCCCC")
        out = tmp_path / "batch.fa"
        runner._write_batch_fasta([seq1, seq2], out)
        text = out.read_text()
        assert ">seq1\nATGCATGC\n" in text
        assert ">seq2\nTTTTCCCC\n" in text

    def test_empty_list_creates_empty_file(self, tmp_path):
        runner = BlastnRunner(mode=BlastnMode.ONLINE)
        out = tmp_path / "empty.fa"
        runner._write_batch_fasta([], out)
        assert out.read_text() == ""


# ---------------------------------------------------------------------------
# BlastnOutputParser.parse_tsv
# ---------------------------------------------------------------------------

class TestParseTsv:
    def test_parses_valid_row(self, tmp_path):
        tsv = _write_tsv(tmp_path, [_tsv_row()])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 1
        h = hits[0]
        assert h.qseqid    == "seq1"
        assert h.pident    == 92.5
        assert h.qlen      == 1200
        assert h.slen      == 13600
        assert h.qcovhsp   == 80
        assert h.evalue    == pytest.approx(1e-60)
        assert h.bit_score == 200.0
        assert "Influenza" in h.stitle

    def test_parses_multiple_rows(self, tmp_path):
        rows = [_tsv_row("seq1"), _tsv_row("seq2")]
        tsv  = _write_tsv(tmp_path, rows)
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 2
        assert {h.qseqid for h in hits} == {"seq1", "seq2"}

    def test_nonexistent_file_returns_empty(self, tmp_path):
        hits = BlastnOutputParser.parse_tsv(tmp_path / "missing.tsv")
        assert hits == []

    def test_empty_file_returns_empty(self, tmp_path):
        p = tmp_path / "empty.tsv"
        p.write_text("")
        assert BlastnOutputParser.parse_tsv(p) == []

    def test_malformed_row_skipped(self, tmp_path):
        rows = ["only\ttwo\tcolumns", _tsv_row()]
        tsv  = _write_tsv(tmp_path, rows)
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert len(hits) == 1

    def test_qcovhsp_float_truncated_to_int(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row(qcovhsp="79.9")])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert hits[0].qcovhsp == 79

    def test_scientific_evalue(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row(evalue="1.2e-100")])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert hits[0].evalue == pytest.approx(1.2e-100)

    def test_returns_blastn_result_instances(self, tmp_path):
        tsv  = _write_tsv(tmp_path, [_tsv_row()])
        hits = BlastnOutputParser.parse_tsv(tsv)
        assert all(isinstance(h, BlastnResult) for h in hits)


# ---------------------------------------------------------------------------
# BlastnOutputParser.parse_xml_record
# ---------------------------------------------------------------------------

def _make_hsp(bits=200.0, identities=90, align_length=100,
              query_start=1, query_end=100, expect=1e-50):
    hsp = MagicMock()
    hsp.bits         = bits
    hsp.identities   = identities
    hsp.align_length = align_length
    hsp.query_start  = query_start
    hsp.query_end    = query_end
    hsp.expect       = expect
    return hsp


def _make_alignment(title="NC_002017 Influenza PB2", hit_def="Influenza PB2",
                    length=13600, hsps=None):
    aln = MagicMock()
    aln.title   = title
    aln.hit_def = hit_def
    aln.length  = length
    aln.hsps    = hsps if hsps is not None else [_make_hsp()]
    return aln


def _make_blast_record(query_length=1200, alignments=None):
    rec = MagicMock()
    rec.query_length = query_length
    rec.alignments   = alignments if alignments is not None else [_make_alignment()]
    return rec


class TestParseXmlRecord:
    def test_basic_parse(self):
        rec  = _make_blast_record()
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert len(hits) == 1
        h = hits[0]
        assert h.qseqid    == "seq1"
        assert h.qlen      == 1200
        assert h.slen      == 13600
        assert h.bit_score == 200.0
        assert h.evalue    == pytest.approx(1e-50)
        assert h.stitle    == "Influenza PB2"

    def test_pident_computed_from_hsp(self):
        hsp = _make_hsp(identities=95, align_length=100)
        rec = _make_blast_record(alignments=[_make_alignment(hsps=[hsp])])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].pident == 95.0

    def test_qcovhsp_computed(self):
        # query_start=1, query_end=600, qlen=1200 → span=600, cov=50%
        hsp = _make_hsp(query_start=1, query_end=600)
        rec = _make_blast_record(query_length=1200, alignments=[_make_alignment(hsps=[hsp])])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].qcovhsp == 50

    def test_best_hsp_selected_by_bitscore(self):
        hsp_low  = _make_hsp(bits=100.0, identities=60)
        hsp_high = _make_hsp(bits=400.0, identities=95)
        aln = _make_alignment(hsps=[hsp_low, hsp_high])
        hits = BlastnOutputParser.parse_xml_record(_make_blast_record(alignments=[aln]), "seq1")
        assert hits[0].bit_score == 400.0
        assert hits[0].pident    == 95.0

    def test_multiple_alignments(self):
        rec = _make_blast_record(alignments=[_make_alignment(), _make_alignment(hit_def="Other virus")])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert len(hits) == 2

    def test_empty_alignments_returns_empty(self):
        rec = _make_blast_record(alignments=[])
        assert BlastnOutputParser.parse_xml_record(rec, "seq1") == []

    def test_alignment_without_hsps_skipped(self):
        aln = _make_alignment(hsps=[])
        rec = _make_blast_record(alignments=[aln])
        assert BlastnOutputParser.parse_xml_record(rec, "seq1") == []

    def test_hit_def_none_falls_back_to_title(self):
        hsp = _make_hsp()
        aln = _make_alignment(title=">NC_002017 Fallback title", hit_def=None, hsps=[hsp])
        aln.hit_def = None
        rec = _make_blast_record(alignments=[aln])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert "Fallback title" in hits[0].stitle

    def test_zero_query_length_gives_zero_coverage(self):
        rec = _make_blast_record(query_length=0)
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].qcovhsp == 0

    def test_returns_blastn_result_instances(self):
        rec  = _make_blast_record()
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert all(isinstance(h, BlastnResult) for h in hits)


# ---------------------------------------------------------------------------
# BlastnOutputParser._accession_from_alignment  (online Hit_id fallback)
# ---------------------------------------------------------------------------

def _aln(hit_id=None, accession="__unset__", hsps=None):
    """Alignment stub. accession is omitted from the object when '__unset__',
    mirroring Biopython alignments parsed from online qblast XML (no Hit_accession)."""
    from types import SimpleNamespace
    ns = SimpleNamespace(title="t", hit_def="Some virus", length=1000,
                         hsps=hsps if hsps is not None else [_make_hsp()])
    if hit_id is not None:
        ns.hit_id = hit_id
    if accession != "__unset__":
        ns.accession = accession
    return ns


class TestAccessionResolution:
    _resolve = staticmethod(BlastnOutputParser._accession_from_alignment)

    def test_explicit_hit_accession_preferred(self):
        # Local BLAST+/XML v2 sets alignment.accession — it wins over Hit_id.
        aln = _aln(hit_id="gi|123|ref|X99999.9|", accession="NC_045512")
        assert self._resolve(aln) == "NC_045512"

    def test_ref_tag_in_pipe_hit_id(self):
        aln = _aln(hit_id="gi|1798174254|ref|NC_045512.2|")
        assert self._resolve(aln) == "NC_045512.2"

    def test_gb_tag_in_pipe_hit_id(self):
        aln = _aln(hit_id="gi|555|gb|MN908947.3|")
        assert self._resolve(aln) == "MN908947.3"

    def test_plain_accession_hit_id_no_pipes(self):
        aln = _aln(hit_id="NC_045512.2")
        assert self._resolve(aln) == "NC_045512.2"

    def test_no_tag_falls_back_to_last_token(self):
        # Pipe-delimited but no recognised db tag → last token.
        aln = _aln(hit_id="12345|NC_045512.2")
        assert self._resolve(aln) == "NC_045512.2"

    def test_no_accession_and_no_hit_id_returns_none(self):
        aln = _aln(hit_id=None)
        assert self._resolve(aln) is None

    def test_parse_xml_record_uses_hit_id_fallback(self):
        # Online path: no Hit_accession → accession pulled from Hit_id.
        aln  = _aln(hit_id="gi|1798174254|ref|NC_045512.2|")
        rec  = _make_blast_record(alignments=[aln])
        hits = BlastnOutputParser.parse_xml_record(rec, "seq1")
        assert hits[0].accession == "NC_045512.2"


# ---------------------------------------------------------------------------
# BlastnResultAttacher.attach
# ---------------------------------------------------------------------------

class TestBlastnResultAttacher:
    def _hit(self, qseqid: str, bit_score: float) -> BlastnResult:
        return BlastnResult(
            qseqid=qseqid, qlen=1200, slen=13600,
            qcovhsp=80, pident=92.5, evalue=1e-50,
            bit_score=bit_score,
            stitle="Influenza A",
        )

    def test_hits_attached_to_matching_seq(self):
        seq  = make_seq("seq1")
        hit  = self._hit("seq1", 200.0)
        BlastnResultAttacher.attach([hit], [seq])
        assert len(seq.blastn_hits) == 1

    def test_hits_sorted_by_bitscore_descending(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", b) for b in [100.0, 300.0, 200.0]]
        BlastnResultAttacher.attach(hits, [seq])
        scores = [h.bit_score for h in seq.blastn_hits]
        assert scores == sorted(scores, reverse=True)

    def test_max_hits_capped(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i)) for i in range(10)]
        BlastnResultAttacher.attach(hits, [seq], max_hits_per_query=3)
        assert len(seq.blastn_hits) == 3

    def test_top_hits_kept_after_capping(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i * 10)) for i in range(5)]
        BlastnResultAttacher.attach(hits, [seq], max_hits_per_query=2)
        scores = [h.bit_score for h in seq.blastn_hits]
        assert scores[0] == 40.0
        assert scores[1] == 30.0

    def test_unknown_query_id_ignored(self):
        seq = make_seq("seq1")
        hit = self._hit("seq_unknown", 200.0)
        attached = BlastnResultAttacher.attach([hit], [seq])
        assert len(seq.blastn_hits) == 0
        assert attached == 0

    def test_returns_total_attached_count(self):
        seq1 = make_seq("seq1")
        seq2 = make_seq("seq2")
        hits = [self._hit("seq1", 100.0), self._hit("seq2", 200.0), self._hit("seq2", 300.0)]
        count = BlastnResultAttacher.attach(hits, [seq1, seq2])
        assert count == 3

    def test_empty_hits_returns_zero(self):
        seq = make_seq("seq1")
        assert BlastnResultAttacher.attach([], [seq]) == 0
        assert seq.blastn_hits == []

    def test_empty_seqs_returns_zero(self):
        hit = self._hit("seq1", 200.0)
        assert BlastnResultAttacher.attach([hit], []) == 0

    def test_multiple_seqs_independent(self):
        seq1 = make_seq("seq1")
        seq2 = make_seq("seq2")
        hits = [self._hit("seq1", 100.0), self._hit("seq2", 200.0)]
        BlastnResultAttacher.attach(hits, [seq1, seq2])
        assert len(seq1.blastn_hits) == 1
        assert len(seq2.blastn_hits) == 1
        assert seq1.blastn_hits[0].bit_score == 100.0
        assert seq2.blastn_hits[0].bit_score == 200.0

    def test_default_max_is_five(self):
        seq  = make_seq("seq1")
        hits = [self._hit("seq1", float(i)) for i in range(10)]
        BlastnResultAttacher.attach(hits, [seq])
        assert len(seq.blastn_hits) == 5


# ---------------------------------------------------------------------------
# BlastnSubjectFetcher — full subject sequences (blastdbcmd → efetch)
# ---------------------------------------------------------------------------

def _bn_hit(acc, title="Dengue virus 2, complete genome", slen=10723):
    return BlastnResult(qseqid="q", qlen=1000, slen=slen, qcovhsp=90, pident=95.0,
                        evalue=0.0, bit_score=900.0, stitle=title, accession=acc)


class _Proc:
    def __init__(self, stdout="", stderr="", returncode=0):
        self.stdout, self.stderr, self.returncode = stdout, stderr, returncode


class TestBlastnSubjectFetcher:

    def _fetcher(self, run=None, http=None, db="/db/nt", **kw):
        calls = {"run": [], "http": []}

        def _run(cmd, **_):
            calls["run"].append(cmd)
            return run(cmd) if run else _Proc()

        def _http(params):
            calls["http"].append(params)
            return http(params) if http else ""

        f = BlastnSubjectFetcher(db_path=db, http=_http, run=_run, sleep=lambda s: None, **kw)
        return f, calls

    def test_non_viral_titles_are_never_fetched(self):
        f, calls = self._fetcher()
        hit = _bn_hit("NC_000001", title="Homo sapiens chromosome 1")
        counts = f.fetch([hit])
        assert counts["requested"] == 0
        assert calls == {"run": [], "http": []}
        assert hit.subject_seq is None

    def test_oversized_subjects_are_skipped(self):
        f, calls = self._fetcher(max_len=1000)
        assert f.fetch([_bn_hit("MN1", slen=5000)])["requested"] == 0
        assert calls["http"] == []

    def test_blastdbcmd_result_is_used(self):
        with patch("viralquest.blastn.shutil.which", return_value="/bin/blastdbcmd"):
            f, calls = self._fetcher(run=lambda cmd: _Proc(">MN908947\nACGT\n"))
            hit = _bn_hit("MN908947")
            counts = f.fetch([hit])
        assert hit.subject_seq == "ACGT"
        assert counts["blastdbcmd"] == 1 and counts["efetch"] == 0
        assert calls["http"] == []
        assert "-entry_batch" in calls["run"][0]

    def test_missing_entries_fall_back_to_efetch(self):
        # blastdbcmd exits non-zero for a missing entry but prints the found ones.
        with patch("viralquest.blastn.shutil.which", return_value="/bin/blastdbcmd"):
            f, calls = self._fetcher(
                run=lambda cmd: _Proc(">A1\nAAAA\n", "Entry not found: B2", 1),
                http=lambda p: ">B2.1 Some virus, complete genome\nCC\nGG\n",
            )
            a, b = _bn_hit("A1"), _bn_hit("B2")
            counts = f.fetch([a, b])
        assert (a.subject_seq, b.subject_seq) == ("AAAA", "CCGG")
        assert calls["http"][0]["id"] == "B2"
        assert calls["http"][0]["db"] == "nuccore"
        assert counts == {"requested": 2, "blastdbcmd": 1, "efetch": 1, "missing": 0}

    def test_missing_blastdbcmd_binary_uses_efetch(self):
        with patch("viralquest.blastn.shutil.which", return_value=None):
            f, calls = self._fetcher(http=lambda p: ">A1.2\nTT\n")
            hit = _bn_hit("A1")
            f.fetch([hit])
        assert calls["run"] == []
        assert hit.subject_seq == "TT"

    def test_online_mode_uses_efetch_only(self):
        f, calls = self._fetcher(db=None, http=lambda p: ">A1.1\nGG\n", email="me@x.org")
        hit = _bn_hit("A1")
        f.fetch([hit])
        assert calls["run"] == []
        assert calls["http"][0]["email"] == "me@x.org"
        assert hit.subject_seq == "GG"

    def test_duplicate_accessions_are_fetched_once(self):
        f, calls = self._fetcher(db=None, http=lambda p: ">A1.1\nGG\n")
        hits = [_bn_hit("A1"), _bn_hit("A1")]
        f.fetch(hits)
        assert calls["http"][0]["id"] == "A1"
        assert [h.subject_seq for h in hits] == ["GG", "GG"]

    def test_efetch_is_batched(self):
        f, calls = self._fetcher(db=None)
        f.EFETCH_BATCH = 2
        f.fetch([_bn_hit(f"A{i}") for i in range(5)])
        assert [p["id"] for p in calls["http"]] == ["A0,A1", "A2,A3", "A4"]

    def test_efetch_failure_does_not_raise(self):
        def boom(p):
            raise OSError("network down")
        f, _ = self._fetcher(db=None, http=boom)
        hit = _bn_hit("A1")
        counts = f.fetch([hit])
        assert hit.subject_seq is None
        assert counts["missing"] == 1

    def test_runner_builds_fetcher_for_its_mode(self, tmp_path):
        local = BlastnRunner(mode=BlastnMode.LOCAL, db_path="/db/nt",
                             blastn_bin="/opt/blast/bin/blastn", outdir=str(tmp_path))
        f = local.subject_fetcher()
        assert f.db_path == "/db/nt"
        assert f.blastdbcmd_bin == "/opt/blast/bin/blastdbcmd"
        online = BlastnRunner(mode=BlastnMode.ONLINE, outdir=str(tmp_path), email="me@x.org")
        g = online.subject_fetcher()
        assert g.db_path is None and g.email == "me@x.org"
