"""
setup_env.py — viralquest-setup / viralquest-download entry points.

viralquest-setup:
    1. Checks for required bioinformatics binaries.
    2. Installs pixi if absent, then runs `pixi install`.
    3. Downloads all reference databases from Zenodo.

viralquest-download:
    Re-runs only Step 3 (database download). Use --force to overwrite.

Run viralquest-setup once after `pip install viralquest`.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Binaries provided by conda/bioconda that pip cannot install.
_REQUIRED_TOOLS: list[tuple[str, str]] = [
    ("diamond",  "diamond  >=2.1.24  (bioconda)"),
    ("blastn",   "blast    >=2.16    (bioconda)"),
    ("salmon",   "salmon   >=1.11.4  (bioconda)"),
    ("cap3",     "cap3     >=10.2011 (bioconda)"),
    ("hmmbuild", "hmmer    >=3.4     (bioconda)"),
]

_PIXI_INSTALL_CMD  = "curl -fsSL https://pixi.sh/install.sh | bash"
_PIXI_TOML_NAME    = "pixi.toml"


# ---------------------------------------------------------------------------
# Public helpers (also imported by cli.py for the startup check)
# ---------------------------------------------------------------------------

def missing_tools() -> list[str]:
    """Return names of required binaries not found on PATH."""
    return [name for name, _ in _REQUIRED_TOOLS if not shutil.which(name)]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_pixi_toml() -> Path | None:
    """Walk up from this file looking for pixi.toml (works in editable and wheel installs)."""
    search = Path(__file__).resolve().parent
    for _ in range(4):
        candidate = search / _PIXI_TOML_NAME
        if candidate.exists():
            return candidate
        search = search.parent
    return None


def _pixi_available() -> bool:
    return shutil.which("pixi") is not None


def _install_pixi() -> bool:
    print("  pixi not found — installing via official script...")
    result = subprocess.run(_PIXI_INSTALL_CMD, shell=True)
    if result.returncode != 0:
        return False
    pixi_bin = Path.home() / ".pixi" / "bin"
    os.environ["PATH"] = str(pixi_bin) + ":" + os.environ.get("PATH", "")
    return _pixi_available()


def _run_pixi_install(pixi_toml: Path) -> bool:
    print(f"  Running: pixi install  (cwd: {pixi_toml.parent})")
    return subprocess.run(["pixi", "install"], cwd=pixi_toml.parent).returncode == 0


def _print_header(text: str) -> None:
    bar = "─" * len(text)
    print(f"\n{bar}\n{text}\n{bar}")


# ---------------------------------------------------------------------------
# Entry point 1 — viralquest-setup
# ---------------------------------------------------------------------------

def setup() -> None:
    _print_header("ViralQuest environment setup")

    # ── Step 1: check binaries ───────────────────────────────────────────────
    absent = missing_tools()
    if absent:
        print("Missing bioinformatics tools:")
        absent_set = set(absent)
        for name, desc in _REQUIRED_TOOLS:
            if name in absent_set:
                print(f"  ✗  {desc}")

        pixi_toml = _find_pixi_toml()
        if pixi_toml is None:
            _print_header("Manual installation required")
            print(
                "Could not find pixi.toml. Install the missing tools manually:\n\n"
                "  conda install -c conda-forge -c bioconda \\\n"
                "    'diamond>=2.1.24' 'blast>=2.16' 'salmon>=1.11.4' \\\n"
                "    'cap3>=10.2011'   'hmmer>=3.4'\n"
            )
            sys.exit(1)

        print(f"\nFound pixi.toml → {pixi_toml}")

        # ── Step 2: install pixi ─────────────────────────────────────────────
        _print_header("Step 1 / 3 — pixi")
        if _pixi_available():
            print("  ✓ pixi is already installed.")
        else:
            if not _install_pixi():
                print(
                    "\n✗ pixi installation failed.\n"
                    "  Install manually: https://pixi.sh\n"
                    "  Then re-run: viralquest-setup\n"
                )
                sys.exit(1)
            print("  ✓ pixi installed successfully.")

        # ── Step 3: pixi install ─────────────────────────────────────────────
        _print_header("Step 2 / 3 — conda / bioconda packages")
        if not _run_pixi_install(pixi_toml):
            print(
                "\n✗ pixi install failed.\n"
                "  Check the output above, then re-run: viralquest-setup\n"
            )
            sys.exit(1)

        still_missing = missing_tools()
        if still_missing:
            print(
                f"\n✗ Still missing after setup: {', '.join(still_missing)}\n"
                "  Open a new terminal (PATH may not be updated yet) and\n"
                "  re-run: viralquest-setup\n"
            )
            sys.exit(1)

        print("  ✓ All bioinformatics tools installed.")
    else:
        print("✓ All bioinformatics tools already available on PATH.")

    # ── Step 4: check and download databases ────────────────────────────────
    _print_header("Step 3 / 3 — Reference databases")
    from viralquest.download_dbs import check_databases, download_all

    db_status = check_databases()
    present   = [name for name, ok in db_status.items() if ok]
    missing   = [name for name, ok in db_status.items() if not ok]

    for name in present:
        print(f"  ✓  {name}")
    for name in missing:
        print(f"  ✗  {name}  ← will be downloaded")

    if not missing:
        print("\n  All databases already present — skipping download.")
    else:
        print(f"\n  Downloading {len(missing)} missing database(s)...")
        download_all(force=False)

    _print_header("Setup complete")
    print("✓ ViralQuest is ready to use.\n")
    print("  Run:  viralquest --help\n")


# ---------------------------------------------------------------------------
# Entry point 2 — viralquest-download
# ---------------------------------------------------------------------------

def download() -> None:
    """Standalone database downloader (re-runs only the download step)."""
    parser = argparse.ArgumentParser(
        prog="viralquest-download",
        description="Download ViralQuest reference databases from Zenodo.",
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-download and overwrite files that already exist.",
    )
    args = parser.parse_args()

    _print_header("ViralQuest — database download")
    from viralquest.download_dbs import check_databases, download_all

    db_status = check_databases()
    for name, ok in db_status.items():
        mark = "✓" if ok else "✗"
        note = "" if ok else "  ← will be downloaded"
        print(f"  {mark}  {name}{note}")

    if not args.force and all(db_status.values()):
        print("\n  All databases already present. Use --force to re-download.\n")
        return

    print()
    download_all(force=args.force)
