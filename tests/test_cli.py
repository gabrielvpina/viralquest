"""Tests for cli.py helpers: step labels, workflow options and the --live screen."""
import json
from collections import deque

import pytest
from rich.console import Console

from viralquest import cli


def _args(*extra: str):
    return cli._build_parser().parse_args(["-in", "sample.fasta", "-out", "out", *extra])


def _render(renderable, width: int = 140, height: int = 30) -> str:
    console = Console(width=width, height=height, record=True,
                      force_terminal=False, color_system=None)
    console.print(renderable, height=height)
    return console.export_text()


# ---------------------------------------------------------------------------
# _build_steps — labels must follow the real execution order in _run_pipeline
# ---------------------------------------------------------------------------

class TestBuildSteps:
    def _index(self, steps, prefix):
        return next(i for i, s in enumerate(steps) if s.startswith(prefix))

    def test_nr_runs_before_blastn(self):
        steps = cli._build_steps(_args("--nr-db", "nr.dmnd", "--blastn-online", "me@x.org"))
        assert self._index(steps, "Diamond BLASTx NR") < self._index(steps, "BLASTn online")

    def test_blastn_online_label_uses_chosen_db(self):
        steps = cli._build_steps(_args("--blastn-online", "me@x.org",
                                       "--blastn-online-db", "refseq_viruses_rep_genomes"))
        assert any(s == "BLASTn online  —  refseq_viruses_rep_genomes" for s in steps)

    def test_llm_step_needs_model_type_and_name(self):
        only_type = cli._build_steps(_args("--model-type", "ollama"))
        both      = cli._build_steps(_args("--model-type", "ollama", "--model-name", "qwen3:4b"))
        assert not any(s.startswith("LLM scoring") for s in only_type)
        assert any(s.startswith("LLM scoring") for s in both)

    def test_minimal_run_ends_with_export_and_html(self):
        steps = cli._build_steps(_args())
        assert steps[-2:] == ["Export JSON report", "Build HTML report"]
        assert not any(s.startswith(("Diamond BLASTx NR", "BLASTn", "Salmon", "LLM"))
                       for s in steps)


# ---------------------------------------------------------------------------
# _workflow_options — shown in the report's workflow card
# ---------------------------------------------------------------------------

class TestWorkflowOptions:
    def test_never_includes_the_api_key(self):
        args = _args("--model-type", "openai", "--model-name", "gpt-4o",
                     "--api-key", "sk-secret-123")
        opts = cli._workflow_options(args)
        assert "sk-secret-123" not in json.dumps(opts)
        assert opts["llm"] == "openai / gpt-4o"

    def test_paths_reduced_to_file_names(self):
        args = _args("--nr-db", "/data/dbs/nr.dmnd", "--blastn-local", "/data/dbs/nt")
        opts = cli._workflow_options(args)
        assert opts["nr_db"] == "nr.dmnd"
        assert opts["blastn"] == "local · nt"
        assert opts["input"] == "sample.fasta"

    def test_disabled_modules_are_falsy(self):
        opts = cli._workflow_options(_args())
        assert not opts["nr_db"] and not opts["blastn"] and not opts["llm"]
        assert opts["reads"] == [] and opts["read_type"] is None

    def test_is_json_serialisable(self):
        json.dumps(cli._workflow_options(_args("--cap3", "--reads", "r1.fq", "r2.fq")))


# ---------------------------------------------------------------------------
# Time formatting on the live screen
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("sec, text", [(5, "0:05"), (75, "1:15"), (3725, "1:02:05"), (-3, "0:00")])
def test_fmt_clock(sec, text):
    assert cli._fmt_clock(sec) == text


@pytest.mark.parametrize("sec, text", [(0.1, "0.1s"), (39.44, "39.4s"), (552.7, "9:12")])
def test_fmt_step_time(sec, text):
    assert cli._fmt_step_time(sec) == text


# ---------------------------------------------------------------------------
# _LiveState
# ---------------------------------------------------------------------------

class TestLiveState:
    def test_finish_step_advances_and_records(self):
        st = cli._LiveState(["a", "b"])
        st.finish_step(0, 1.5)
        assert st.current == 1 and st.timings == {0: 1.5} and not st.done

    def test_last_step_marks_done_and_stops_clock(self):
        st = cli._LiveState(["a"])
        st.finish_step(0, 1.0)
        assert st.done and st.end is not None
        assert st.elapsed() == pytest.approx(st.end - st.run_start)

    def test_fail_marks_current_step(self):
        st = cli._LiveState(["a", "b", "c"])
        st.finish_step(0, 1.0)
        st.fail()
        assert st.failed == 1 and st.end is not None

    def test_batch_cleared_when_step_finishes(self):
        st = cli._LiveState(["a", "b"])
        st.batch = "2/5"
        st.finish_step(0, 1.0)
        assert st.batch is None


# ---------------------------------------------------------------------------
# _LogTail / _LiveScreen rendering
# ---------------------------------------------------------------------------

class TestLiveScreen:
    STEPS = ["Parse FASTA", "Find ORFs", "Export JSON report"]

    def _screen(self, records=(), **state_kw):
        st = cli._LiveState(self.STEPS)
        for k, v in state_kw.items():
            setattr(st, k, v)
        info = cli._workflow_options(_args()) | {"outdir": "out"}
        return cli._LiveScreen(st, deque(records), info), st

    def test_log_text_with_brackets_is_not_parsed_as_markup(self):
        # Regression: "[refseq]" was swallowed and "[/x]" raised MarkupError.
        records = [("15:00:01", "INFO", "diamond [refseq] batch [/weird] done")]
        text = _render(cli._LogTail(deque(records)))
        assert "diamond [refseq] batch [/weird] done" in text

    def test_log_tail_shows_latest_lines_that_fit(self):
        records = deque(("15:00:00", "INFO", f"line {i}") for i in range(50))
        text = _render(cli._LogTail(records), height=5)
        assert "line 49" in text and "line 0 " not in text

    def test_empty_log_shows_placeholder(self):
        assert "Waiting for output" in _render(cli._LogTail(deque()))

    @pytest.mark.parametrize("width", [70, 110, 150])
    def test_renders_at_common_widths(self, width):
        screen, st = self._screen([("15:00:01", "WARNING", "careful")])
        st.finish_step(0, 0.4)
        text = _render(screen, width=width)
        for label in self.STEPS:
            assert label in text
        assert "1/3" in text

    def test_failed_step_is_reported(self):
        screen, st = self._screen()
        st.fail()
        assert "failed" in _render(screen)

    def test_completed_run_is_reported(self):
        screen, st = self._screen()
        for i in range(len(self.STEPS)):
            st.finish_step(i, 1.0)
        text = _render(screen)
        assert "complete" in text and "3/3" in text

    def test_warning_and_error_counters(self):
        screen, _ = self._screen(warnings=2, errors=1)
        text = _render(screen)
        assert "2 warnings" in text and "1 error" in text


# ---------------------------------------------------------------------------
# _log_summary — final results as COMPLETE log lines (no Rich panel)
# ---------------------------------------------------------------------------

class TestLogSummary:
    def _capture(self, args):
        from loguru import logger
        lines = []
        sink = logger.add(lambda m: lines.append((m.record["level"].name, m.record["message"])),
                          level="TRACE")
        try:
            cli._log_summary(args, {0: 1.5, 1: 2.0})
        finally:
            logger.remove(sink)
        return lines

    def _seqs(self, n_viral, n_nr):
        from viralquest.biodata import BlastxResult, NucSequence
        seqs = []
        for i in range(n_viral):
            s = NucSequence(id=f"s{i}", sequence="ATGC" * 20)
            s.is_viral = True
            if i < n_nr:
                s.blastx_nr_hits.append(BlastxResult(
                    query_id=s.id, subject_id="x", subject_title="p [Virus]",
                    pct_identity=90.0, aln_length=10, mismatches=0, gap_opens=0,
                    query_start=1, query_end=10, subject_start=1, subject_end=10,
                    e_value=1e-9, bit_score=50.0))
            seqs.append(s)
        return seqs

    def test_every_line_uses_the_complete_level(self):
        args = _args()
        args._result_seqs, args._result_clusters = self._seqs(3, 0), []
        lines = self._capture(args)
        assert lines and {lvl for lvl, _ in lines} == {"COMPLETE"}
        assert lines[0][1] == "ViralQuest — sample.fasta complete"

    def test_confirmed_count_follows_nr_rule(self):
        # 5 flagged viral, 2 confirmed by NR → the summary reports 2, like the report.
        args = _args("--nr-db", "nr.dmnd")
        args._result_seqs, args._result_clusters = self._seqs(5, 2), []
        msgs = [m for _, m in self._capture(args)]
        assert any(m.startswith("Confirmed viral sequences:") and m.endswith(" 2") for m in msgs)

    def test_unset_outputs_are_omitted(self):
        args = _args()
        args._result_seqs, args._result_clusters = [], []
        msgs = [m for _, m in self._capture(args)]
        assert not any(m.startswith("BLASTn results:") for m in msgs)
        assert any(m.startswith("Total time:") and m.endswith("3.5s") for m in msgs)


# ---------------------------------------------------------------------------
# Hidden dev flag: --reload rebuilds the HTML from an existing report JSON
# ---------------------------------------------------------------------------

_MINI_REPORT = {"meta": {}, "sequences": [], "clusters": [], "pipeline_stats": {}}


class TestReload:
    def _json(self, path, mtime=None):
        import os
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(_MINI_REPORT), encoding="utf-8")
        if mtime is not None:
            os.utime(path, (mtime, mtime))
        return path

    def test_flag_is_hidden_from_help(self):
        assert "reload" not in cli._build_parser().format_help()

    def test_flag_parses_with_and_without_target(self):
        p = cli._build_parser()
        assert p.parse_args(["--reload"]).reload == ""
        assert p.parse_args(["--reload", "x.json"]).reload == "x.json"
        assert p.parse_args(["-in", "a", "-out", "b"]).reload is None

    def test_finds_json_in_current_directory(self, tmp_path, monkeypatch):
        target = self._json(tmp_path / "AT39_viralquest.json")
        monkeypatch.chdir(tmp_path)
        assert cli._find_report_json("") == target

    def test_newest_json_wins(self, tmp_path):
        self._json(tmp_path / "old_viralquest.json", mtime=1_000)
        new = self._json(tmp_path / "new_viralquest.json", mtime=2_000)
        assert cli._find_report_json(str(tmp_path)) == new

    def test_explicit_file_is_used_as_is(self, tmp_path):
        f = self._json(tmp_path / "any_name.json")
        assert cli._find_report_json(str(f)) == f

    def test_other_json_files_are_ignored(self, tmp_path):
        self._json(tmp_path / "settings.json")
        with pytest.raises(FileNotFoundError, match="_viralquest.json"):
            cli._find_report_json(str(tmp_path))

    def test_missing_path_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            cli._find_report_json(str(tmp_path / "nope"))

    def test_rebuilds_html_next_to_json(self, tmp_path):
        from unittest.mock import patch
        src = self._json(tmp_path / "AT39_viralquest.json")
        with patch("urllib.request.urlopen", side_effect=OSError("offline")):
            out = cli._reload_report(str(src))
        assert out == tmp_path / "AT39_viralquest.html"
        assert "d3js.org" in out.read_text(encoding="utf-8")

    def test_main_reload_skips_tool_checks(self, tmp_path, monkeypatch):
        from unittest.mock import patch
        self._json(tmp_path / "AT39_viralquest.json")
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr("sys.argv", ["viralquest", "--reload"])
        with patch("viralquest.setup_env.missing_tools") as tools, \
             pytest.raises(SystemExit) as exit_:
            cli.main()
        assert exit_.value.code == 0
        tools.assert_not_called()
        assert (tmp_path / "AT39_viralquest.html").is_file()

    def test_main_reload_without_json_exits_1(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr("sys.argv", ["viralquest", "--reload"])
        with pytest.raises(SystemExit) as exit_:
            cli.main()
        assert exit_.value.code == 1


class TestBlastnOnlineBatchOption:
    def test_default_and_custom_batch(self):
        assert _args("--blastn-online", "me@x.org").blastn_online_batch == 20
        assert _args("--blastn-online", "me@x.org", "--blastn-online-batch", "5").blastn_online_batch == 5


class TestBlastnOnlineTimeoutOptions:
    def test_defaults(self):
        a = _args("--blastn-online", "me@x.org")
        assert (a.blastn_online_timeout, a.blastn_online_retries) == (15.0, 2)

    def test_custom_values(self):
        a = _args("--blastn-online", "me@x.org",
                  "--blastn-online-timeout", "5", "--blastn-online-retries", "0")
        assert (a.blastn_online_timeout, a.blastn_online_retries) == (5.0, 0)
