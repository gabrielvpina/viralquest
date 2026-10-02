import json
import subprocess
import sys
import warnings
from pathlib import Path
from unittest.mock import patch

import pytest

from viralquest.html_report import (
    _build_taxonomy_tree,
    _enrich_report,
    _fetch_d3,
    _render_template,
    write_report,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

STUB_D3 = "/* d3-stub */"

MINIMAL_REPORT: dict = {
    "meta": {
        "viralquest_version": "1.0.0",
        "timestamp": "2026-01-01T00:00:00",
        "input_file": {"name": "sample.fasta", "size_bytes": 1024},
    },
    "sequences": [],
    "clusters": [],
}


def _seq(seq_id: str, phylum: str, order: str, family: str, genus: str) -> dict:
    return {
        "id": seq_id,   # exported sequences carry their id under "id"
        "taxonomy": {
            "phylum": phylum,
            "order": order,
            "family": family,
            "genus": genus,
        },
    }


# ---------------------------------------------------------------------------
# _build_taxonomy_tree
# ---------------------------------------------------------------------------


class TestBuildTaxonomyTree:
    def test_empty_sequences_returns_root(self):
        tree = _build_taxonomy_tree([])
        assert tree["name"] == "Viruses"
        assert tree["children"] == []

    def test_single_sequence_full_path(self):
        seqs = [_seq("seq1", "Negarnaviricota", "Amarillovirales", "Flaviviridae", "Flavivirus")]
        tree = _build_taxonomy_tree(seqs)

        phylum_node = tree["children"][0]
        assert phylum_node["name"] == "Negarnaviricota"

        order_node = phylum_node["children"][0]
        assert order_node["name"] == "Amarillovirales"

        family_node = order_node["children"][0]
        assert family_node["name"] == "Flaviviridae"

        genus_node = family_node["children"][0]
        assert genus_node["name"] == "Flavivirus"

        leaf = genus_node["children"][0]
        assert leaf["seq_id"] == "seq1"
        assert leaf["name"] == "seq1"
        assert leaf["family"] == "Flaviviridae"
        assert leaf["is_leaf"] is True

    def test_two_sequences_same_family_share_nodes(self):
        seqs = [
            _seq("seq1", "Negarnaviricota", "Amarillovirales", "Flaviviridae", "Flavivirus"),
            _seq("seq2", "Negarnaviricota", "Amarillovirales", "Flaviviridae", "Flavivirus"),
        ]
        tree = _build_taxonomy_tree(seqs)
        # Root → phylum → order → family → genus → 2 leaves
        genus_node = tree["children"][0]["children"][0]["children"][0]["children"][0]
        assert len(genus_node["children"]) == 2
        leaf_ids = {l["seq_id"] for l in genus_node["children"]}
        assert leaf_ids == {"seq1", "seq2"}

    def test_two_sequences_different_families(self):
        seqs = [
            _seq("seq1", "Negarnaviricota", "Amarillovirales", "Flaviviridae", "Flavivirus"),
            _seq("seq2", "Negarnaviricota", "Amarillovirales", "Parvoviridae",  "Amdoparvovirus"),
        ]
        tree = _build_taxonomy_tree(seqs)
        order_node = tree["children"][0]["children"][0]
        family_names = {n["name"] for n in order_node["children"]}
        assert family_names == {"Flaviviridae", "Parvoviridae"}

    def test_two_sequences_different_phyla(self):
        seqs = [
            _seq("seq1", "Negarnaviricota", "A", "B", "C"),
            _seq("seq2", "Pisuviricota",    "X", "Y", "Z"),
        ]
        tree = _build_taxonomy_tree(seqs)
        phylum_names = {n["name"] for n in tree["children"]}
        assert phylum_names == {"Negarnaviricota", "Pisuviricota"}

    def test_missing_taxonomy_fields_use_unclassified(self):
        seqs = [{"seq_id": "seq1", "taxonomy": {}}]
        tree = _build_taxonomy_tree(seqs)
        phylum_node = tree["children"][0]
        assert phylum_node["name"] == "Unclassified"

    def test_null_taxonomy_uses_unclassified(self):
        seqs = [{"seq_id": "seq1", "taxonomy": None}]
        tree = _build_taxonomy_tree(seqs)
        phylum_node = tree["children"][0]
        assert phylum_node["name"] == "Unclassified"

    def test_missing_taxonomy_key_uses_unclassified(self):
        seqs = [{"seq_id": "seq1"}]
        tree = _build_taxonomy_tree(seqs)
        assert tree["children"][0]["name"] == "Unclassified"

    def test_leaf_family_falls_back_to_unclassified(self):
        seqs = [{"seq_id": "s1", "taxonomy": {"phylum": "P", "order": "O", "genus": "G"}}]
        tree = _build_taxonomy_tree(seqs)
        # Walk to leaf
        leaf = (
            tree["children"][0]
            ["children"][0]
            ["children"][0]   # family = Unclassified
            ["children"][0]   # genus
            ["children"][0]   # leaf
        )
        assert leaf["family"] == "Unclassified"

    def test_no_duplicate_intermediate_nodes(self):
        seqs = [
            _seq("s1", "P", "O", "F", "G"),
            _seq("s2", "P", "O", "F", "G"),
            _seq("s3", "P", "O", "F", "G"),
        ]
        tree = _build_taxonomy_tree(seqs)
        # Only one phylum node
        assert len(tree["children"]) == 1

    def test_seq_id_empty_string_handled(self):
        seqs = [{"seq_id": "", "taxonomy": {"phylum": "P", "order": "O", "family": "F", "genus": "G"}}]
        tree = _build_taxonomy_tree(seqs)
        leaf = tree["children"][0]["children"][0]["children"][0]["children"][0]["children"][0]
        assert leaf["seq_id"] == ""
        assert leaf["is_leaf"] is True


# ---------------------------------------------------------------------------
# _enrich_report
# ---------------------------------------------------------------------------


class TestEnrichReport:
    def test_adds_taxonomy_tree_key(self):
        report = {"sequences": []}
        enriched = _enrich_report(report)
        assert "_taxonomy_tree" in enriched

    def test_modifies_report_in_place(self):
        report = {"sequences": []}
        result = _enrich_report(report)
        assert result is report

    def test_tree_built_from_sequences(self):
        report = {
            "sequences": [_seq("s1", "P", "O", "F", "G")]
        }
        _enrich_report(report)
        assert report["_taxonomy_tree"]["children"][0]["name"] == "P"

    def test_no_sequences_key(self):
        report = {}
        _enrich_report(report)
        assert report["_taxonomy_tree"] == {"name": "Viruses", "children": []}


# ---------------------------------------------------------------------------
# _render_template
# ---------------------------------------------------------------------------


class TestRenderTemplate:
    def test_all_placeholders_replaced(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        remaining = [m for m in ["{{VQ_STYLES}}", "{{VQ_DATA}}", "{{VQ_D3}}",
                                  "{{VQ_EXPORT}}", "{{VQ_STATS}}", "{{VQ_VIRUSES}}", "{{VQ_CLUSTERS}}",
                                  "{{VQ_VIEWER}}", "{{VQ_TAXONOMY}}", "{{VQ_SALMON}}",
                                  "{{SAMPLE_NAME}}"]
                     if m in html]
        assert remaining == [], f"Unresolved placeholders: {remaining}"

    def test_run_and_virome_tabs_each_have_their_own_script(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        assert 'id="section-run"' in html and 'id="section-virome"' in html
        assert "window.vqInitStats = vqInitStats" in html          # Run tab
        assert "window.vqInitViruses = vqInitViruses" in html      # Virome tab
        assert "vqInitStats(VQ_REPORT);" in html and "vqInitViruses(VQ_REPORT);" in html

    def test_overview_helpers_shared_through_vq(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        for name in ("statChip", "fmtNum", "pct", "cardGroup", "redrawOnResize"):
            assert f"{name}:" in html, name

    def test_report_data_is_valid_json_in_output(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        # Extract the JSON assigned to VQ_REPORT
        import re
        match = re.search(r"const VQ_REPORT = ({.*?});", html, re.DOTALL)
        assert match, "VQ_REPORT not found in rendered HTML"
        parsed = json.loads(match.group(1))
        assert parsed["meta"]["viralquest_version"] == "1.0.0"

    def test_sample_name_from_meta(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        assert "sample.fasta" in html

    def test_sample_name_fallback_when_missing(self):
        report = {**MINIMAL_REPORT, "meta": {}}
        html = _render_template(report, d3_js=STUB_D3)
        assert "ViralQuest Report" in html

    def test_d3_stub_inlined(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        assert STUB_D3 in html

    def test_taxonomy_tree_included_in_json(self):
        report = _enrich_report({**MINIMAL_REPORT, "sequences": [_seq("s1", "P", "O", "F", "G")]})
        html = _render_template(report, d3_js=STUB_D3)
        assert "_taxonomy_tree" in html

    def test_is_valid_html_document(self):
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        assert html.strip().startswith("<!DOCTYPE html")
        assert "</html>" in html

    def test_no_vq_placeholders_remain(self):
        import re
        html = _render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        leftover = re.findall(r"\{\{VQ_\w+\}\}", html)
        assert leftover == []

    def test_warns_on_unresolved_placeholder(self, tmp_path, monkeypatch):
        # Inject a broken template with an unknown placeholder
        broken_template = (
            Path(__file__).parent.parent
            / "viralquest" / "components" / "shell.html"
        ).read_text()
        broken_template += "\n<!-- {{VQ_UNKNOWN}} -->"

        fake_shell = tmp_path / "shell.html"
        fake_shell.write_text(broken_template)

        import viralquest.html_report as hr
        monkeypatch.setattr(hr, "_COMPONENTS", tmp_path)
        # Copy real component files so other placeholders resolve
        import shutil
        real_comp = Path(__file__).parent.parent / "viralquest" / "components"
        for f in real_comp.iterdir():
            if f.name == "shell.html":
                continue
            if f.is_dir():                      # e.g. components/assets/
                shutil.copytree(f, tmp_path / f.name)
            else:
                shutil.copy(f, tmp_path / f.name)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            hr._render_template(MINIMAL_REPORT, d3_js=STUB_D3)
        assert any("VQ_UNKNOWN" in str(warning.message) for warning in w)


# ---------------------------------------------------------------------------
# D3 loading (d3_asset.load_d3, used by both report builders)
# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, data: bytes):
        self._data = data
    def read(self):
        return self._data
    def __enter__(self):
        return self
    def __exit__(self, *_):
        pass


class TestFetchD3:
    def test_bundled_d3_ships_with_the_package(self):
        from viralquest.d3_asset import D3_VENDORED, D3_VERSION
        assert D3_VENDORED.is_file()
        head = D3_VENDORED.read_text(encoding="utf-8")[:80]
        assert f"d3js.org v{D3_VERSION}" in head
        assert (D3_VENDORED.parent / "LICENSE-d3.txt").is_file()

    def test_bundled_copy_used_without_network(self, tmp_path):
        from viralquest.d3_asset import load_d3
        vendored = tmp_path / "d3.min.js"
        vendored.write_text("/* bundled-d3 */")
        with patch("urllib.request.urlopen") as net:
            assert load_d3(vendored, tmp_path / "cache.js") == "/* bundled-d3 */"
        net.assert_not_called()

    def test_user_cache_used_when_bundle_missing(self, tmp_path):
        from viralquest.d3_asset import load_d3
        cache = tmp_path / "cache.js"
        cache.write_text("/* cached-d3 */")
        assert load_d3(tmp_path / "missing.js", cache) == "/* cached-d3 */"

    def test_download_written_to_user_cache(self, tmp_path):
        from viralquest.d3_asset import load_d3
        cache = tmp_path / "sub" / "cache.js"
        with patch("urllib.request.urlopen", return_value=_FakeResp(b"/* fetched-d3 */")):
            assert load_d3(tmp_path / "missing.js", cache) == "/* fetched-d3 */"
        assert cache.read_text() == "/* fetched-d3 */"

    def test_unwritable_cache_does_not_break(self, tmp_path):
        from viralquest.d3_asset import load_d3
        blocker = tmp_path / "file"
        blocker.write_text("")                  # a file where a directory is needed
        with patch("urllib.request.urlopen", return_value=_FakeResp(b"/* d3 */")):
            assert load_d3(tmp_path / "missing.js", blocker / "cache.js") == "/* d3 */"

    def test_offline_and_not_bundled_raises(self, tmp_path):
        from viralquest.d3_asset import load_d3
        with patch("urllib.request.urlopen", side_effect=OSError("no network")):
            with pytest.raises(RuntimeError, match="d3_js="):
                load_d3(tmp_path / "missing.js", tmp_path / "cache.js")

    def test_user_cache_honours_xdg_cache_home(self, tmp_path, monkeypatch):
        from viralquest.d3_asset import user_cache_path
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
        assert user_cache_path().parent == tmp_path / "viralquest"

    def test_both_report_builders_use_the_bundled_d3(self):
        import viralquest.html_report as hr
        import viralquest.report_builder as rb
        from viralquest.d3_asset import D3_VENDORED
        bundled = D3_VENDORED.read_text(encoding="utf-8")
        with patch("urllib.request.urlopen") as net:
            assert hr._fetch_d3() == bundled
            assert rb._fetch_d3() == bundled
        net.assert_not_called()


# ---------------------------------------------------------------------------
# write_report
# ---------------------------------------------------------------------------


class TestWriteReport:
    def test_creates_output_file(self, tmp_path):
        out = tmp_path / "report.html"
        write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert out.exists()

    def test_returns_path_object(self, tmp_path):
        out = tmp_path / "report.html"
        result = write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert isinstance(result, Path)
        assert result == out

    def test_accepts_str_path(self, tmp_path):
        out = str(tmp_path / "report.html")
        result = write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert Path(result).exists()

    def test_creates_parent_dirs(self, tmp_path):
        out = tmp_path / "nested" / "deep" / "report.html"
        write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert out.exists()

    def test_output_is_utf8(self, tmp_path):
        report = {**MINIMAL_REPORT, "meta": {
            **MINIMAL_REPORT["meta"],
            "input_file": {"name": "mosquitões.fasta", "size_bytes": 0},
        }}
        out = tmp_path / "report.html"
        write_report(report, out, d3_js=STUB_D3)
        content = out.read_text(encoding="utf-8")
        assert "mosquitões" in content

    def test_enriches_report_with_taxonomy_tree(self, tmp_path):
        report = {
            **MINIMAL_REPORT,
            "sequences": [_seq("s1", "Negarnaviricota", "O", "Flaviviridae", "G")],
        }
        out = tmp_path / "report.html"
        write_report(report, out, d3_js=STUB_D3)
        html = out.read_text()
        assert "_taxonomy_tree" in html
        assert "Negarnaviricota" in html

    def test_salmon_quant_included_when_present(self, tmp_path):
        report = {
            **MINIMAL_REPORT,
            "salmon_quant": {"mapping_rate": 55.2, "total_reads": 1000000,
                             "viral_quant": [], "conserved_quant": [], "host_viral_hits": []},
        }
        out = tmp_path / "report.html"
        write_report(report, out, d3_js=STUB_D3)
        html = out.read_text()
        assert "salmon_quant" in html

    def test_no_placeholders_in_output(self, tmp_path):
        import re
        out = tmp_path / "report.html"
        write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        html = out.read_text()
        leftover = re.findall(r"\{\{VQ_\w+\}\}", html)
        assert leftover == []

    def test_uses_d3_js_kwarg_skips_fetch(self, tmp_path):
        out = tmp_path / "report.html"
        with patch("urllib.request.urlopen", side_effect=OSError("should not be called")):
            write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert out.exists()

    def test_overwrites_existing_file(self, tmp_path):
        out = tmp_path / "report.html"
        out.write_text("old content")
        write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        assert "old content" not in out.read_text()

    def test_json_in_output_is_parseable(self, tmp_path):
        import re
        out = tmp_path / "report.html"
        write_report(MINIMAL_REPORT, out, d3_js=STUB_D3)
        html = out.read_text()
        match = re.search(r"const VQ_REPORT = ({.*?});", html, re.DOTALL)
        assert match
        parsed = json.loads(match.group(1))
        assert parsed["meta"]["viralquest_version"] == "1.0.0"


# ---------------------------------------------------------------------------
# CLI (_cli)
# ---------------------------------------------------------------------------


class TestCli:
    def test_cli_writes_output(self, tmp_path):
        in_json = tmp_path / "report.json"
        in_json.write_text(json.dumps(MINIMAL_REPORT), encoding="utf-8")
        out_html = tmp_path / "out.html"
        d3_file  = tmp_path / "d3.min.js"
        d3_file.write_text(STUB_D3)

        result = subprocess.run(
            [sys.executable, "-m", "viralquest.html_report",
             str(in_json), str(out_html), "--d3", str(d3_file)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert out_html.exists()

    def test_cli_prints_path_to_stderr(self, tmp_path):
        in_json = tmp_path / "report.json"
        in_json.write_text(json.dumps(MINIMAL_REPORT), encoding="utf-8")
        out_html = tmp_path / "out.html"
        d3_file  = tmp_path / "d3.min.js"
        d3_file.write_text(STUB_D3)

        result = subprocess.run(
            [sys.executable, "-m", "viralquest.html_report",
             str(in_json), str(out_html), "--d3", str(d3_file)],
            capture_output=True, text=True,
        )
        assert "out.html" in result.stderr

    def test_without_d3_uses_bundled_copy_offline(self, tmp_path):
        import viralquest.html_report as hr
        from viralquest.d3_asset import D3_VENDORED, D3_VERSION
        out_html = tmp_path / "out.html"
        with patch("urllib.request.urlopen", side_effect=OSError("offline")) as net:
            hr.write_report(json.loads(json.dumps(MINIMAL_REPORT)), out_html)
        net.assert_not_called()
        assert f"d3js.org v{D3_VERSION}" in out_html.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# _enrich_report — workflow card data
# ---------------------------------------------------------------------------

class TestEnrichWorkflow:
    def test_workflow_exposed_from_pipeline_stats(self):
        from viralquest.html_report import _enrich_report
        wf = {"total_seconds": 12.5, "steps": [{"key": "parse", "status": "done"}]}
        report = _enrich_report({"sequences": [], "clusters": [],
                                 "pipeline_stats": {"workflow": wf}})
        assert report["workflow"] == wf

    def test_older_reports_without_workflow_get_none(self):
        from viralquest.html_report import _enrich_report
        report = _enrich_report({"sequences": [], "clusters": [], "pipeline_stats": {}})
        assert report["workflow"] is None
