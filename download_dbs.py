#!/usr/bin/env python3
"""
download_dbs.py — Download and decompress ViralQuest reference databases.

Fetches all HMM profiles and the Diamond viral reference database from Zenodo
and places them in the expected locations:

  data/hmm-dbs/EggNOG-4.5.hmm
  data/hmm-dbs/Pfam-A.hmm
  data/hmm-dbs/U-RVDBv29.0-prot.hmm
  data/hmm-dbs/Vfam-228.hmm
  data/viralDB.dmnd

Already-downloaded targets are skipped automatically.
Run with  --force  to re-download everything.
"""

from __future__ import annotations

import argparse
import gzip
import lzma
import shutil
import sys
import tempfile
from pathlib import Path

import requests
from tqdm import tqdm

# ── Database manifest ────────────────────────────────────────────────────────

_ROOT    = Path(__file__).parent
_HMM_DIR = _ROOT / "data" / "hmm-dbs"
_DATA    = _ROOT / "data"

_ZENODO  = "https://zenodo.org/records/18715455/files"

DATABASES = [
    {
        "name":   "EggNOG HMM profiles",
        "url":    f"{_ZENODO}/EggNOG-4.5.hmm.xz?download=1",
        "target": _HMM_DIR / "EggNOG-4.5.hmm",
        "fmt":    "xz",
    },
    {
        "name":   "Pfam-A HMM profiles",
        "url":    f"{_ZENODO}/Pfam-A.hmm.xz?download=1",
        "target": _HMM_DIR / "Pfam-A.hmm",
        "fmt":    "xz",
    },
    {
        "name":   "U-RVDB viral HMM profiles",
        "url":    f"{_ZENODO}/U-RVDBv29.0-prot.hmm.xz?download=1",
        "target": _HMM_DIR / "U-RVDBv29.0-prot.hmm",
        "fmt":    "xz",
    },
    {
        "name":   "Vfam viral HMM profiles",
        "url":    f"{_ZENODO}/Vfam-228.hmm.xz?download=1",
        "target": _HMM_DIR / "Vfam-228.hmm",
        "fmt":    "xz",
    },
    {
        "name":   "ViralDB Diamond database",
        "url":    f"{_ZENODO}/viralDB.dmnd.gz?download=1",
        "target": _DATA / "viralDB.dmnd",
        "fmt":    "gz",
    },
]

# ── Download helpers ─────────────────────────────────────────────────────────

_CHUNK = 1 << 20   # 1 MiB


def _download_to_temp(url: str, label: str) -> Path:
    """Stream URL to a temp file; show a tqdm progress bar. Returns temp path."""
    resp = requests.get(url, stream=True, timeout=30)
    resp.raise_for_status()

    total = int(resp.headers.get("content-length", 0)) or None

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".download")
    try:
        with tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=f"  downloading {label}",
            leave=False,
        ) as bar:
            for chunk in resp.iter_content(chunk_size=_CHUNK):
                tmp.write(chunk)
                bar.update(len(chunk))
        tmp.flush()
        return Path(tmp.name)
    except Exception:
        tmp.close()
        Path(tmp.name).unlink(missing_ok=True)
        raise
    finally:
        tmp.close()


def _decompress(src: Path, dst: Path, fmt: str) -> None:
    """Decompress src → dst using lzma (.xz) or gzip (.gz)."""
    opener = lzma.open if fmt == "xz" else gzip.open
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp_dst = dst.with_suffix(".tmp")
    try:
        with opener(src, "rb") as fin, tmp_dst.open("wb") as fout:
            with tqdm(
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=f"  decompressing",
                leave=False,
            ) as bar:
                while True:
                    block = fin.read(_CHUNK)
                    if not block:
                        break
                    fout.write(block)
                    bar.update(len(block))
        tmp_dst.replace(dst)
    except Exception:
        tmp_dst.unlink(missing_ok=True)
        raise


# ── Main logic ───────────────────────────────────────────────────────────────

def download_all(*, force: bool = False) -> None:
    _HMM_DIR.mkdir(parents=True, exist_ok=True)

    ok = skipped = failed = 0

    for db in DATABASES:
        target: Path = db["target"]
        name: str    = db["name"]

        print(f"\n[{name}]")

        if target.exists() and not force:
            size_mb = target.stat().st_size / 1e6
            print(f"  already present ({size_mb:.1f} MB) — skipping. Use --force to re-download.")
            skipped += 1
            continue

        tmp: Path | None = None
        try:
            tmp = _download_to_temp(db["url"], target.name)
            _decompress(tmp, target, db["fmt"])
            size_mb = target.stat().st_size / 1e6
            print(f"  done — {size_mb:.1f} MB written to {target.relative_to(_ROOT)}")
            ok += 1
        except requests.HTTPError as exc:
            print(f"  ERROR: HTTP {exc.response.status_code} — {exc}", file=sys.stderr)
            failed += 1
        except Exception as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            failed += 1
        finally:
            if tmp is not None:
                tmp.unlink(missing_ok=True)

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    total = len(DATABASES)
    print(f"Done: {ok} downloaded, {skipped} skipped, {failed} failed  (of {total} total)")

    if failed:
        sys.exit(1)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-download and overwrite files that already exist.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    download_all(force=args.force)
