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
import platform
import shutil
import subprocess
import sys
from pathlib import Path

# Binaries provided by conda/bioconda that pip cannot install.
# Note: HMM searches use pyhmmer (Python binding) — no hmmer binary needed.
# CAP3 ships no macOS build on bioconda (neither osx-arm64 nor osx-64), so the
# optional --cap3 assembly step is Linux-only; it's excluded from the macOS
# requirement list.
_REQUIRED_TOOLS: list[tuple[str, str]] = [
    ("diamond",   "diamond        >=2.1.24  (bioconda)"),
    ("blastn",    "blast          >=2.16    (bioconda)"),  # also provides dustmasker (self-BLAST / low-complexity)
    ("salmon",    "salmon         >=1.11.4  (bioconda)"),
    ("minimap2",  "minimap2       >=2.28    (bioconda)"),
    ("samtools",  "samtools       >=1.21    (bioconda)"),
    ("jellyfish", "kmer-jellyfish >=2.3     (bioconda)"),  # k-mer repetitiveness
] + ([] if platform.system() == "Darwin" else [
    ("cap3", "cap3      >=10.2011 (bioconda)"),
])

_PIXI_INSTALL_CMD  = "curl -fsSL https://pixi.sh/install.sh | bash"
_PIXI_TOML_NAME    = "pixi.toml"


# ---------------------------------------------------------------------------
# Public helpers (also imported by cli.py for the startup check)
# ---------------------------------------------------------------------------

def _pixi_env_bin() -> Path | None:
    """Return the pixi default-env bin dir if it exists."""
    pixi_toml = _find_pixi_toml()
    if pixi_toml is None:
        return None
    candidate = pixi_toml.parent / ".pixi" / "envs" / "default" / "bin"
    return candidate if candidate.exists() else None


def activate_pixi_env() -> Path | None:
    """
    Prepend the pixi default-env bin dir to PATH for this process (idempotent).

    Every entry point that shells out to a bioinformatics binary should call
    this first, so `blastn`, `diamond`, … resolve to the pixi installation
    instead of whatever happens to sit on the user's PATH. Returns the bin dir,
    or None when there is no pixi env to activate.
    """
    env_bin = _pixi_env_bin()
    if env_bin is None:
        return None
    entries = os.environ.get("PATH", "").split(os.pathsep)
    if str(env_bin) not in entries:
        os.environ["PATH"] = os.pathsep.join([str(env_bin), *entries]).strip(os.pathsep)
    return env_bin


def resolve_tool(name: str, override: str | None = None) -> str | None:
    """
    Full path to a bioinformatics binary, preferring the pixi installation.

    *override* is an explicit user-supplied binary: a path (used as-is when it
    exists) or a bare name looked up on PATH. Without an override the pixi env
    is activated and *name* is resolved from PATH. Returns None if not found.
    """
    activate_pixi_env()
    target = override or name
    if os.sep in target or (os.altsep and os.altsep in target):
        candidate = Path(target).expanduser()
        return str(candidate) if candidate.is_file() else None
    return shutil.which(target)


def missing_tools() -> list[str]:
    """Return names of required binaries not found on PATH or the pixi env."""
    activate_pixi_env()
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
    ok = subprocess.run(["pixi", "install"], cwd=pixi_toml.parent).returncode == 0
    if ok:
        activate_pixi_env()
    return ok


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
            cap3_pin = "" if platform.system() == "Darwin" else "    'cap3>=10.2011'   "
            print(
                "Could not find pixi.toml. Install the missing tools manually:\n\n"
                "  conda install -c conda-forge -c bioconda \\\n"
                "    'diamond>=2.1.24' 'blast>=2.16' 'salmon>=1.11.4' \\\n"
                f"{cap3_pin}'minimap2>=2.28' 'samtools>=1.21'\n"
                + ("" if platform.system() != "Darwin" else
                   "\n  Note: CAP3 has no native macOS build on bioconda,\n"
                   "  so the optional --cap3 assembly step is unavailable on macOS.\n")
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

    # ── Step 4: offer database download ─────────────────────────────────────
    _print_header("Step 3 / 3 — Reference databases")
    from viralquest.download_dbs import check_databases, download_all

    db_status = check_databases()
    present   = [name for name, ok in db_status.items() if ok]
    missing   = [name for name, ok in db_status.items() if not ok]

    for name in present:
        print(f"  ✓  {name}")
    for name in missing:
        print(f"  ✗  {name}")

    if not missing:
        print("\n  All databases already present.")
    else:
        print(f"\n  {len(missing)} database(s) are missing (~several GB total).")
        print("  They will be downloaded via curl and stored inside the package.")
        print()
        try:
            answer = input("  Download now? [Y/n]: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            answer = "n"
            print()

        if answer in ("", "y", "yes"):
            from viralquest.download_dbs import prompt_source
            source = prompt_source()
            print()
            download_all(force=False, source=source)
        else:
            print()
            print("  Skipped. You can download later with:")
            print("    viralquest-download")
            print()
            print("  Or place the database files (.hmm + viralDB.dmnd) in any flat")
            print("  directory and point viralquest at it at runtime:")
            print("    viralquest --db-dir /path/to/dbs  -in input.fasta  -out results/")
            print()

    _print_header("Setup complete")
    print("✓ ViralQuest is ready to use.\n")
    print("  Run:  viralquest --help\n")


# ---------------------------------------------------------------------------
# Entry point 2 — viralquest-download
# ---------------------------------------------------------------------------

def download() -> None:
    """Standalone database downloader (re-runs only the download step)."""
    from pathlib import Path

    parser = argparse.ArgumentParser(
        prog="viralquest-download",
        description="Download ViralQuest reference databases from Zenodo or Google Drive.",
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Re-download and overwrite files that already exist.",
    )
    parser.add_argument(
        "--db-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="Directory where databases will be stored "
             "(default: data/ inside the package).",
    )
    parser.add_argument(
        "--source",
        choices=["zenodo", "gdrive"],
        default=None,
        metavar="SOURCE",
        help="Download source: 'zenodo' (default) or 'gdrive' (Google Drive). "
             "If omitted, you will be prompted interactively.",
    )
    args = parser.parse_args()

    db_dir = Path(args.db_dir) if args.db_dir else None
    flat   = db_dir is not None   # custom dir → flat layout (no hmm-dbs/ subdir)

    _print_header("ViralQuest — database download")
    from viralquest.download_dbs import check_databases, download_all, prompt_source

    db_status = check_databases(db_dir=db_dir, flat=flat)
    for name, ok in db_status.items():
        mark = "✓" if ok else "✗"
        note = "" if ok else "  ← will be downloaded"
        print(f"  {mark}  {name}{note}")

    if not args.force and all(db_status.values()):
        print("\n  All databases already present. Use --force to re-download.\n")
        return

    source = args.source if args.source else prompt_source()
    print()
    download_all(force=args.force, db_dir=db_dir, flat=flat, source=source)
