"""
report_builder.py — Assemble the self-contained consolidated HTML for
``viralquest-report`` from a ``VQ_REPORT`` data dict (see report_data.py).

Mirrors html_report.py but draws its template, CSS and JS from
``viralquest/components-report/`` and injects the multi-sample sections.
"""

from __future__ import annotations

import base64
import json
import re
import urllib.request
from pathlib import Path

_HERE       = Path(__file__).parent
_COMPONENTS = _HERE / "components-report"
_D3_VERSION = "7.9.0"
_D3_CDN_URL = f"https://cdn.jsdelivr.net/npm/d3@{_D3_VERSION}/dist/d3.min.js"
_D3_CACHE   = _COMPONENTS / ".d3.min.js.cache"


# ── Public entry point ─────────────────────────────────────────────────────

def write_report(report_data: dict, output_html: str | Path, *, d3_js: str | None = None) -> Path:
    """Build and write the consolidated self-contained HTML report."""
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    html = _render_template(report_data, d3_js=d3_js or _fetch_d3())
    output_html.write_text(html, encoding="utf-8")
    return output_html


# ── Asset helpers ──────────────────────────────────────────────────────────

def _asset_data_uri(filename: str) -> str:
    path = _COMPONENTS / "assets" / filename
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode()
    ext  = path.suffix.lower().lstrip(".")
    mime = {"png": "image/png", "svg": "image/svg+xml", "ico": "image/x-icon"}.get(ext, "image/png")
    return f"data:{mime};base64,{data}"


def _read(name: str) -> str:
    return (_COMPONENTS / name).read_text(encoding="utf-8")


# ── Template rendering ─────────────────────────────────────────────────────

def _render_template(report_data: dict, d3_js: str) -> str:
    template = _read("shell.html")

    n = report_data.get("meta", {}).get("n_samples", 0)
    title = f"ViralQuest Report — {n} sample{'s' if n != 1 else ''}"

    replacements = {
        "{{TITLE}}":        title,
        "{{VQ_FAVICON}}":   _asset_data_uri("favicon.png"),
        "{{VQ_LOGO}}":      _asset_data_uri("logo-bg-text.png"),
        "{{VQ_STYLES}}":    _read("base.css"),
        "{{VQ_DATA}}":      json.dumps(report_data, ensure_ascii=False),
        "{{VQ_D3}}":        d3_js,
        "{{VQ_EXPORT}}":    _read("export.js"),
        "{{VQ_ABOUT}}":     _read("section_about.js"),
        "{{VQ_OVERVIEW}}":  _read("section_overview.js"),
        "{{VQ_CLUSTERS}}":  _read("section_clusters.js"),
        "{{VQ_VIEWER}}":         _read("section_viewer.js"),
        "{{VQ_VIEWER_REPORT}}":  _read("section_viewer_report.js"),
        "{{VQ_QUANT}}":          _read("section_quant.js"),
    }

    html = template
    for placeholder, content in replacements.items():
        html = html.replace(placeholder, content, 1)

    missed = re.findall(r"\{\{VQ_\w+\}\}", html)
    if missed:
        import warnings
        warnings.warn(f"Unresolved template placeholders: {missed}", stacklevel=2)

    return html


# ── D3 fetching / caching (shared cache file with the per-sample report) ───

def _fetch_d3() -> str:
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
