#!/usr/bin/env python3
"""
download_dbs.py — Download and decompress ViralQuest reference databases.

Fetches all HMM profiles and the Diamond viral reference database from either
Zenodo (default) or Google Drive, and places them in the database directory:

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
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from tqdm import tqdm

# ── Database manifest ─────────────────────────────────────────────────────────

_ZENODO = "https://zenodo.org/records/18715455/files"
_CHUNK  = 1 << 20   # 1 MiB

_GDRIVE_IDS: dict[str, str] = {
    "EggNOG-4.5.hmm.xz":      "1NnbDE_z8BsziIlXilMfs69hJ1mep27GW",
    "Pfam-A.hmm.xz":           "1757Itpc4t-gLVrKg2nmN-EBrTI0nfsR9",
    "U-RVDBv29.0-prot.hmm.xz": "154ASdC9zIi9CnRnF0Yp4jjDWmS78aPaS",
    "Vfam-228.hmm.xz":         "1GpaN_BpRd2dXtoAsqrwGwh4Mmf5-4jtv",
    "viralDB.dmnd.gz":         "1TH78rAg-31eOY4bbtJakXjgNS_Zwm25y",
}


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
            "name":      "EggNOG HMM profiles",
            "filename":  "EggNOG-4.5.hmm.xz",
            "url":       f"{_ZENODO}/EggNOG-4.5.hmm.xz?download=1",
            "target":    hmm_dir / "EggNOG-4.5.hmm",
            "fmt":       "xz",
        },
        {
            "name":      "Pfam-A HMM profiles",
            "filename":  "Pfam-A.hmm.xz",
            "url":       f"{_ZENODO}/Pfam-A.hmm.xz?download=1",
            "target":    hmm_dir / "Pfam-A.hmm",
            "fmt":       "xz",
        },
        {
            "name":      "U-RVDB viral HMM profiles",
            "filename":  "U-RVDBv29.0-prot.hmm.xz",
            "url":       f"{_ZENODO}/U-RVDBv29.0-prot.hmm.xz?download=1",
            "target":    hmm_dir / "U-RVDBv29.0-prot.hmm",
            "fmt":       "xz",
        },
        {
            "name":      "Vfam viral HMM profiles",
            "filename":  "Vfam-228.hmm.xz",
            "url":       f"{_ZENODO}/Vfam-228.hmm.xz?download=1",
            "target":    hmm_dir / "Vfam-228.hmm",
            "fmt":       "xz",
        },
        {
            "name":      "ViralDB Diamond database",
            "filename":  "viralDB.dmnd.gz",
            "url":       f"{_ZENODO}/viralDB.dmnd.gz?download=1",
            "target":    dmnd_dir / "viralDB.dmnd",
            "fmt":       "gz",
        },
    ]


# ── Download helpers ──────────────────────────────────────────────────────────

def _download_zenodo(url: str, label: str) -> Path:
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


def _download_gdrive(file_id: str, label: str) -> Path:
    """Download a file from Google Drive using gdown. Returns temp path."""
    import gdown

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".download")
    tmp.close()
    tmp_path = Path(tmp.name)
    print(f"  downloading {label} (Google Drive)")
    try:
        result = gdown.download(id=file_id, output=str(tmp_path), quiet=False)
        if result is None:
            tmp_path.unlink(missing_ok=True)
            raise RuntimeError(f"gdown failed to download {label}")
        return tmp_path
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


_MAGIC = {
    "xz": b"\xfd7zXZ\x00",
    "gz": b"\x1f\x8b",
}


def _is_compressed(path: Path, fmt: str) -> bool:
    """Return True if the file header matches the expected compression format."""
    magic = _MAGIC[fmt]
    with path.open("rb") as f:
        header = f.read(len(magic))
    return header == magic


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


# ── Source selection prompt ───────────────────────────────────────────────────

def prompt_source() -> str:
    """Ask the user to pick a download source. Returns 'zenodo' or 'gdrive'."""
    print()
    print("  Download source:")
    print("    [1] Zenodo       (official, always available)")
    print("    [2] Google Drive (may be faster, daily quota applies)")
    print()
    while True:
        try:
            answer = input("  Choice [1/2] (default: 1): ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return "zenodo"
        if answer in ("", "1"):
            return "zenodo"
        if answer == "2":
            return "gdrive"
        print("  Please enter 1 or 2.")


# ── Pre-flight check ──────────────────────────────────────────────────────────

def check_databases(db_dir: Path | None = None, flat: bool = False) -> dict[str, bool]:
    """Return a mapping of database name → True (present) / False (missing)."""
    if db_dir is None:
        db_dir = _default_db_dir()
        flat   = False
    return {db["name"]: db["target"].exists() for db in _get_databases(db_dir, flat=flat)}


# ── Main logic ────────────────────────────────────────────────────────────────

def download_all(
    *,
    force: bool = False,
    db_dir: Path | None = None,
    flat: bool = False,
    source: str = "zenodo",
) -> None:
    if db_dir is None:
        db_dir = _default_db_dir()
        flat   = False

    if source not in ("zenodo", "gdrive"):
        raise ValueError(f"Unknown source '{source}'. Choose 'zenodo' or 'gdrive'.")

    databases = _get_databases(db_dir, flat=flat)
    if flat:
        db_dir.mkdir(parents=True, exist_ok=True)
    else:
        (db_dir / "hmm-dbs").mkdir(parents=True, exist_ok=True)
        (db_dir / "viral-db").mkdir(parents=True, exist_ok=True)

    source_label = "Zenodo" if source == "zenodo" else "Google Drive"
    print(f"  Source: {source_label}\n")

    ok = skipped = failed = 0

    for db in databases:
        target:   Path = db["target"]
        name:     str  = db["name"]
        filename: str  = db["filename"]

        print(f"\n[{name}]")

        if target.exists() and not force:
            size_mb = target.stat().st_size / 1e6
            print(f"  already present ({size_mb:.1f} MB) — skipping. Use --force to re-download.")
            skipped += 1
            continue

        tmp: Path | None = None
        try:
            if source == "gdrive":
                gdrive_id = _GDRIVE_IDS[filename]
                tmp = _download_gdrive(gdrive_id, target.name)
            else:
                tmp = _download_zenodo(db["url"], target.name)

            if _is_compressed(tmp, db["fmt"]):
                _decompress(tmp, target, db["fmt"])
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(tmp), str(target))
                tmp = None  # already moved, skip unlink in finally
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
    p.add_argument(
        "--source",
        choices=["zenodo", "gdrive"],
        default=None,
        metavar="SOURCE",
        help="Download source: 'zenodo' (default) or 'gdrive' (Google Drive). "
             "If omitted, you will be prompted interactively.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args   = _parse_args()
    db_dir = Path(args.db_dir) if args.db_dir else None
    flat   = db_dir is not None
    source = args.source if args.source else prompt_source()
    download_all(force=args.force, db_dir=db_dir, flat=flat, source=source)
