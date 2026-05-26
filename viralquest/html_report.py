"""
html_report.py — Assemble the self-contained ViralQuest HTML report.

Reads a ViralQuest JSON report dict (from exporter.ReportExporter.export()),
builds the _taxonomy_tree structure, inlines D3, all CSS, and all JS
component files, then writes a single .html file that has no external deps.
"""

from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path
from typing import Any

_HERE        = Path(__file__).parent
_COMPONENTS  = _HERE / "components"
_D3_VERSION  = "7.9.0"
_D3_CDN_URL  = f"https://cdn.jsdelivr.net/npm/d3@{_D3_VERSION}/dist/d3.min.js"
_D3_CACHE    = _COMPONENTS / ".d3.min.js.cache"


# ── Public entry point ────────────────────────────────────────────────────────

def write_report(report: dict, output_path: str | Path, *, d3_js: str | None = None) -> Path:
    """
    Build and write the self-contained HTML report.

    Parameters
    ----------
    report      : dict produced by ReportExporter.export()
    output_path : destination .html file
    d3_js       : optional pre-fetched D3 source (skips network fetch)

    Returns
    -------
    Path to the written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = _enrich_report(report)
    html   = _render_template(report, d3_js=d3_js or _fetch_d3())
    output_path.write_text(html, encoding="utf-8")
    return output_path


# ── Report enrichment ─────────────────────────────────────────────────────────

def _enrich_report(report: dict) -> dict:
    """Add _taxonomy_tree built from sequences[].taxonomy."""
    seqs = report.get("sequences", [])
    report["_taxonomy_tree"] = _build_taxonomy_tree(seqs)
    return report


def _build_taxonomy_tree(sequences: list[dict]) -> dict:
    """
    Convert flat sequence taxonomy dicts into a D3-hierarchy-compatible tree.

    Shape:
        {name: "Viruses", children: [
            {name: <phylum>, children: [
                {name: <order>, children: [
                    {name: <family>, children: [
                        {name: <genus>, children: [
                            {name: <seq_id>, seq_id: <seq_id>, family: <family>, is_leaf: true}
                        ]}
                    ]}
                ]}
            ]}
        ]}
    """
    root: dict[str, Any] = {"name": "Viruses", "children": []}
    idx: dict[str, dict] = {}

    for seq in sequences:
        tax   = seq.get("taxonomy") or {}
        parts = [
            tax.get("phylum")  or "Unclassified",
            tax.get("order")   or "Unclassified",
            tax.get("family")  or "Unclassified",
            tax.get("genus")   or "Unclassified",
        ]
        parent: dict = root
        key    = ""
        for name in parts:
            key += "/" + name
            if key not in idx:
                node: dict[str, Any] = {"name": name, "children": []}
                parent["children"].append(node)
                idx[key] = node
            parent = idx[key]

        parent["children"].append({
            "name":    seq.get("seq_id", ""),
            "seq_id":  seq.get("seq_id", ""),
            "family":  tax.get("family") or "Unclassified",
            "is_leaf": True,
        })

    return root


# ── Template rendering ────────────────────────────────────────────────────────

def _render_template(report: dict, d3_js: str) -> str:
    template_path = _COMPONENTS / "shell.html"
    template      = template_path.read_text(encoding="utf-8")

    sample_name = (
        (report.get("meta") or {}).get("input_file", {}).get("name", "")
        or "ViralQuest Report"
    )

    replacements = {
        "{{SAMPLE_NAME}}":  sample_name,
        "{{VQ_STYLES}}":    (_COMPONENTS / "base.css").read_text(encoding="utf-8"),
        "{{VQ_DATA}}":      json.dumps(report, ensure_ascii=False),
        "{{VQ_D3}}":        d3_js,
        "{{VQ_EXPORT}}":    (_COMPONENTS / "export.js").read_text(encoding="utf-8"),
        "{{VQ_STATS}}":     (_COMPONENTS / "section_stats.js").read_text(encoding="utf-8"),
        "{{VQ_CLUSTERS}}":  (_COMPONENTS / "section_clusters.js").read_text(encoding="utf-8"),
        "{{VQ_VIEWER}}":    (_COMPONENTS / "section_viewer.js").read_text(encoding="utf-8"),
        "{{VQ_TAXONOMY}}":  (_COMPONENTS / "section_taxonomy.js").read_text(encoding="utf-8"),
        "{{VQ_SALMON}}":    (_COMPONENTS / "section_salmon.js").read_text(encoding="utf-8"),
    }

    html = template
    for placeholder, content in replacements.items():
        html = html.replace(placeholder, content, 1)

    # Warn if any placeholder was missed
    missed = re.findall(r"\{\{VQ_\w+\}\}", html)
    if missed:
        import warnings
        warnings.warn(f"Unresolved template placeholders: {missed}", stacklevel=2)

    return html


# ── D3 fetching / caching ─────────────────────────────────────────────────────

def _fetch_d3() -> str:
    """
    Return D3 source.  Uses a local cache file to avoid repeated network calls.
    Falls back to a minimal stub that raises a clear error if offline and uncached.
    """
    if _D3_CACHE.exists():
        return _D3_CACHE.read_text(encoding="utf-8")

    try:
        with urllib.request.urlopen(_D3_CDN_URL, timeout=15) as resp:
            js = resp.read().decode("utf-8")
        _D3_CACHE.write_text(js, encoding="utf-8")
        return js
    except Exception as exc:
        raise RuntimeError(
            f"Could not fetch D3 v{_D3_VERSION} from CDN and no cache found.\n"
            f"Supply the D3 source manually via the d3_js= parameter.\n"
            f"Original error: {exc}"
        ) from exc


# ── CLI convenience ───────────────────────────────────────────────────────────

def _cli() -> None:
    import argparse, sys

    parser = argparse.ArgumentParser(
        description="Build a self-contained ViralQuest HTML report from a JSON file."
    )
    parser.add_argument("input_json",  help="Path to ViralQuest JSON report")
    parser.add_argument("output_html", help="Destination HTML file")
    parser.add_argument(
        "--d3", metavar="d3.min.js",
        help="Local D3 source file (skips network fetch)"
    )
    args = parser.parse_args()

    report = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    d3_js  = Path(args.d3).read_text(encoding="utf-8") if args.d3 else None

    out = write_report(report, args.output_html, d3_js=d3_js)
    print(f"Report written to {out}", file=sys.stderr)


if __name__ == "__main__":
    _cli()
