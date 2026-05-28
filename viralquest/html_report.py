"""
html_report.py — Assemble the self-contained ViralQuest HTML report.

Reads a ViralQuest JSON report dict (from exporter.ReportExporter.export()),
builds the _taxonomy_tree structure, inlines D3, all CSS, and all JS
component files, then writes a single .html file that has no external deps.
"""

from __future__ import annotations

import base64
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
    """Add _taxonomy_tree and pre-aggregated stats to the report."""
    seqs     = report.get("sequences", [])
    clusters = report.get("clusters", [])
    ps       = report.get("pipeline_stats") or {}   # pre-computed by exporter when available

    report["_taxonomy_tree"] = _build_taxonomy_tree(seqs)
    report["summary"]        = _build_summary(seqs, clusters, ps)
    report["blast_stats"]    = _build_blast_stats(seqs, ps)
    report["hmm_stats"]      = _build_hmm_stats(seqs, ps)
    report["llm_stats"]      = _build_llm_stats(seqs, ps)
    report["salmon_stats"]   = _build_salmon_stats(report.get("salmon_quant"))
    return report


def _build_summary(seqs: list[dict], clusters: list[dict], ps: dict) -> dict:
    blast = ps.get("blast") or {}
    hmm   = ps.get("hmm")   or {}
    return {
        "total_sequences": blast.get("total_input",    len(seqs)),
        "confirmed_viral": blast.get("total_confirmed", sum(1 for s in seqs if s.get("is_viral"))),
        "total_orfs":      hmm.get("total_orfs",       sum(len(s.get("orfs") or []) for s in seqs)),
        "total_clusters":  ps.get("clusters",          len(clusters)),
        "cap3_used":       False,
    }


def _build_blast_stats(seqs: list[dict], ps: dict) -> dict:
    blast = ps.get("blast") or {}
    return {
        "refseq_unique_seqs": blast.get("refseq_unique_seqs", sum(1 for s in seqs if s.get("blastx_hits"))),
        "nr_unique_seqs":     blast.get("nr_unique_seqs",     sum(1 for s in seqs if s.get("blastx_nr_hits"))),
        "blastn_unique_seqs": blast.get("blastn_unique_seqs", sum(1 for s in seqs if s.get("blastn_hits"))),
    }


def _build_hmm_stats(seqs: list[dict], ps: dict) -> dict:
    hmm = ps.get("hmm") or {}
    if hmm:
        return {
            "rvdb_hits":   hmm.get("rvdb_hits",   0),
            "vfam_hits":   hmm.get("vfam_hits",   0),
            "eggnog_hits": hmm.get("eggnog_hits", 0),
            "pfam_hits":   hmm.get("pfam_hits",   0),
        }
    counts: dict[str, int] = {"RVDB": 0, "Vfam": 0, "EggNOG": 0, "Pfam": 0}
    for s in seqs:
        for orf in (s.get("orfs") or []):
            for dom in (orf.get("domains") or []):
                db = dom.get("database", "")
                if db in counts:
                    counts[db] += 1
    return {
        "rvdb_hits":   counts["RVDB"],
        "vfam_hits":   counts["Vfam"],
        "eggnog_hits": counts["EggNOG"],
        "pfam_hits":   counts["Pfam"],
    }


def _build_llm_stats(seqs: list[dict], ps: dict) -> dict:
    llm = ps.get("llm") or {}
    if llm and llm.get("present") is not None:
        return llm
    scored = [s for s in seqs if s.get("llm_output")]
    if not scored:
        return {"present": False}
    scores = [
        s["llm_output"]["vq_score"]
        for s in scored
        if s["llm_output"].get("vq_score") is not None
    ]
    first = scored[0]["llm_output"]
    return {
        "present":       True,
        "model":         first.get("model"),
        "mode":          first.get("mode"),
        "scored":        len(scored),
        "viral_known":   sum(1 for s in scored if s["llm_output"].get("classification") == "viral-known"),
        "viral_unknown": sum(1 for s in scored if s["llm_output"].get("classification") == "viral-unknown"),
        "avg_score":     round(sum(scores) / len(scores), 1) if scores else None,
    }


def _build_salmon_stats(salmon_quant: dict | None) -> dict:
    if not salmon_quant:
        return {"present": False}
    return {
        "present":         True,
        "pathway":         salmon_quant.get("pathway", "reference"),
        "mapping_rate":    salmon_quant.get("mapping_rate"),
        "total_reads":     salmon_quant.get("total_reads"),
        "host_viral_hits": len(salmon_quant.get("host_viral_hits") or []),
        "ref_hk_count":    len(salmon_quant.get("ref_hk_quant") or []),
    }


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
            "name":    seq.get("id", ""),
            "seq_id":  seq.get("id", ""),
            "family":  tax.get("family") or "Unclassified",
            "is_leaf": True,
        })

    return root


# ── Asset helpers ─────────────────────────────────────────────────────────────

def _asset_data_uri(filename: str) -> str:
    """Return a base64 data URI for an asset in components/assets/."""
    path = _COMPONENTS / "assets" / filename
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode()
    ext  = path.suffix.lower().lstrip(".")
    mime = {"png": "image/png", "svg": "image/svg+xml", "ico": "image/x-icon"}.get(ext, "image/png")
    return f"data:{mime};base64,{data}"


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
        "{{VQ_FAVICON}}":   _asset_data_uri("favicon.png"),
        "{{VQ_LOGO}}":      _asset_data_uri("logo-bg-text.png"),
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
