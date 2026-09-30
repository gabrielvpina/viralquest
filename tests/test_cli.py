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
