#!/usr/bin/env python3
"""
download_dbs.py — Download and decompress ViralQuest reference databases.

Fetches all HMM profiles and the Diamond viral reference database from Zenodo
and places them in the specified database directory:

  <db-dir>/hmm-dbs/EggNOG-4.5.hmm
  <db-dir>/hmm-dbs/Pfam-A.hmm
  <db-dir>/hmm-dbs/U-RVDBv29.0-prot.hmm
  <db-dir>/hmm-dbs/Vfam-228.hmm
  <db-dir>/viral-db/viralDB.dmnd

Already-downloaded targets are skipped automatically.
Run with  --force  to re-download everything.
"""

from __future__ import annotations

import argparse
import gzip
import lzma
import subprocess
import sys
import tempfile
from pathlib import Path

from tqdm import tqdm

# ── Database manifest ─────────────────────────────────────────────────────────

_ZENODO = "https://zenodo.org/records/18715455/files"
_CHUNK  = 1 << 20   # 1 MiB


def _default_db_dir() -> Path:
    return Path(__file__).parent.parent / "data"


def _get_databases(db_dir: Path, flat: bool = False) -> list[dict]:
    """Build database manifest.

    flat=False (default): nested layout — HMM files go into db_dir/hmm-dbs/.
    flat=True: all five files land directly in db_dir with no subdirectory.
    This matches the layout expected by viralquest --db-dir.
    """
    hmm_dir  = db_dir if flat else (db_dir / "hmm-dbs")
    dmnd_dir = db_dir if flat else (db_dir / "viral-db")
    return [
        {
            "name":   "EggNOG HMM profiles",
            "url":    f"{_ZENODO}/EggNOG-4.5.hmm.xz?download=1",
            "target": hmm_dir / "EggNOG-4.5.hmm",
            "fmt":    "xz",
        },
        {
            "name":   "Pfam-A HMM profiles",
            "url":    f"{_ZENODO}/Pfam-A.hmm.xz?download=1",
            "target": hmm_dir / "Pfam-A.hmm",
            "fmt":    "xz",
        },
        {
            "name":   "U-RVDB viral HMM profiles",
            "url":    f"{_ZENODO}/U-RVDBv29.0-prot.hmm.xz?download=1",
            "target": hmm_dir / "U-RVDBv29.0-prot.hmm",
            "fmt":    "xz",
        },
        {
            "name":   "Vfam viral HMM profiles",
            "url":    f"{_ZENODO}/Vfam-228.hmm.xz?download=1",
            "target": hmm_dir / "Vfam-228.hmm",
            "fmt":    "xz",
        },
        {
            "name":   "ViralDB Diamond database",
            "url":    f"{_ZENODO}/viralDB.dmnd.gz?download=1",
            "target": dmnd_dir / "viralDB.dmnd",
            "fmt":    "gz",
        },
    ]


# ── Download helpers ──────────────────────────────────────────────────────────

def _download_to_temp(url: str, label: str) -> Path:
    """Download URL to a temp file using curl. Returns temp path."""
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".download")
    tmp.close()
    tmp_path = Path(tmp.name)
    print(f"  downloading {label}")
    try:
        result = subprocess.run(
            ["curl", "-fSL", "--progress-bar", "-o", str(tmp_path), url],
            check=False,
        )
        if result.returncode != 0:
            tmp_path.unlink(missing_ok=True)
            raise RuntimeError(f"curl exited with code {result.returncode} for {label}")
        return tmp_path
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


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
                desc="  decompressing",
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


# ── Pre-flight check ──────────────────────────────────────────────────────────

def check_databases(db_dir: Path | None = None, flat: bool = False) -> dict[str, bool]:
    """Return a mapping of database name → True (present) / False (missing)."""
    if db_dir is None:
        db_dir = _default_db_dir()
        flat   = False
    return {db["name"]: db["target"].exists() for db in _get_databases(db_dir, flat=flat)}


# ── Main logic ────────────────────────────────────────────────────────────────

def download_all(*, force: bool = False, db_dir: Path | None = None, flat: bool = False) -> None:
    if db_dir is None:
        db_dir = _default_db_dir()
        flat   = False

    databases = _get_databases(db_dir, flat=flat)
    if flat:
        db_dir.mkdir(parents=True, exist_ok=True)
    else:
        (db_dir / "hmm-dbs").mkdir(parents=True, exist_ok=True)
        (db_dir / "viral-db").mkdir(parents=True, exist_ok=True)

    ok = skipped = failed = 0

    for db in databases:
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
            print(f"  done — {size_mb:.1f} MB written to {target}")
            ok += 1
        except Exception as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            failed += 1
        finally:
            if tmp is not None:
                tmp.unlink(missing_ok=True)

    print()
    total = len(databases)
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
    p.add_argument(
        "--db-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="Directory where databases will be stored "
             "(default: data/ inside the package).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args   = _parse_args()
    db_dir = Path(args.db_dir) if args.db_dir else None
    flat   = db_dir is not None
    download_all(force=args.force, db_dir=db_dir, flat=flat)
