"""
d3_asset.py — Locate the D3 source inlined into every HTML report.

D3 v7.9.0 (ISC licence, see vendor/LICENSE-d3.txt) ships inside the package,
so report generation works offline and on read-only installs. The CDN is only
a fallback for broken installs, and its download is cached in the user's
cache directory — never inside the (possibly read-only) package directory.
"""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

D3_VERSION  = "7.9.0"
D3_CDN_URL  = f"https://cdn.jsdelivr.net/npm/d3@{D3_VERSION}/dist/d3.min.js"
D3_VENDORED = Path(__file__).parent / "vendor" / f"d3-{D3_VERSION}.min.js"


def user_cache_path() -> Path:
    """Per-user download cache ($XDG_CACHE_HOME or ~/.cache)."""
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "viralquest" / f"d3-{D3_VERSION}.min.js"


def load_d3(vendored: Path | None = None, cache: Path | None = None) -> str:
    """
    Return the D3 source: bundled copy → user cache → CDN download.
    Raises RuntimeError when none is available (offline and not bundled).
    """
    vendored = D3_VENDORED if vendored is None else vendored
    cache    = user_cache_path() if cache is None else cache

    for path in (vendored, cache):
        if path.is_file():
            return path.read_text(encoding="utf-8")

    try:
        with urllib.request.urlopen(D3_CDN_URL, timeout=15) as resp:
            js = resp.read().decode("utf-8")
    except Exception as exc:
        raise RuntimeError(
            f"D3 v{D3_VERSION} is missing from the package and could not be "
            f"downloaded from the CDN.\n"
            f"Reinstall viralquest, or supply the D3 source via the d3_js= parameter.\n"
            f"Original error: {exc}"
        ) from exc

    try:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(js, encoding="utf-8")
    except OSError:
        pass   # an unwritable cache must not break report generation
    return js
