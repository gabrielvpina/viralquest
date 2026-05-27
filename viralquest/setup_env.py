"""
setup_env.py — viralquest-setup entry point.

Checks for required bioinformatics binaries, installs pixi if absent,
then runs `pixi install` to pull conda/bioconda packages.
Run this once after `pip install viralquest`.
"""

from __future__ import annotations

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

_PIXI_INSTALL_CMD = (
    'curl -fsSL https://pixi.sh/install.sh | bash'
)

_PIXI_TOML_RELATIVE = "pixi.toml"   # relative to the repo / install root


# ---------------------------------------------------------------------------
# Public helpers (imported by cli.py for the startup check)
# ---------------------------------------------------------------------------

def missing_tools() -> list[str]:
    """Return names of required binaries not found on PATH."""
    return [name for name, _ in _REQUIRED_TOOLS if not shutil.which(name)]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _find_pixi_toml() -> Path | None:
    """
    Walk up from this file's location looking for pixi.toml.
    Works both in editable installs (repo root) and installed wheels
    (the pixi.toml is bundled as package data one level above the package).
    """
    search = Path(__file__).resolve().parent
    for _ in range(4):
        candidate = search / _PIXI_TOML_RELATIVE
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
    # After install the binary lands in ~/.pixi/bin — refresh PATH in-process.
    pixi_bin = Path.home() / ".pixi" / "bin"
    import os
    os.environ["PATH"] = str(pixi_bin) + ":" + os.environ.get("PATH", "")
    return _pixi_available()


def _run_pixi_install(pixi_toml: Path) -> bool:
    print(f"  Running: pixi install  (cwd: {pixi_toml.parent})")
    result = subprocess.run(
        ["pixi", "install"],
        cwd=pixi_toml.parent,
    )
    return result.returncode == 0


def _print_header(text: str) -> None:
    bar = "─" * len(text)
    print(f"\n{bar}\n{text}\n{bar}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def setup() -> None:
    _print_header("ViralQuest environment setup")

    # ── 1. Check which tools are missing ─────────────────────────────────────
    absent = missing_tools()

    if not absent:
        print("✓ All required tools are already available on PATH.")
        print("  No action needed.\n")
        return

    print("Missing bioinformatics tools:")
    absent_set = set(absent)
    for name, desc in _REQUIRED_TOOLS:
        if name in absent_set:
            print(f"  ✗  {desc}")

    # ── 2. Locate pixi.toml ──────────────────────────────────────────────────
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

    print(f"\nFound pixi.toml → {pixi_toml}\n")

    # ── 3. Ensure pixi is installed ──────────────────────────────────────────
    _print_header("Step 1 / 2 — pixi")
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

    # ── 4. Run pixi install ──────────────────────────────────────────────────
    _print_header("Step 2 / 2 — conda / bioconda packages")
    if not _run_pixi_install(pixi_toml):
        print(
            "\n✗ pixi install failed.\n"
            "  Check the output above for details, then re-run: viralquest-setup\n"
        )
        sys.exit(1)

    # ── 5. Verify ────────────────────────────────────────────────────────────
    still_missing = missing_tools()
    if still_missing:
        print(
            f"\n✗ Still missing after setup: {', '.join(still_missing)}\n"
            "  Try opening a new terminal (PATH may not be updated yet) and\n"
            "  re-running: viralquest-setup\n"
        )
        sys.exit(1)

    _print_header("Setup complete")
    print("✓ All tools are installed and available.\n")
    print("  Run:  viralquest --help\n")
