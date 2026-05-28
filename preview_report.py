#!/usr/bin/env python3
"""
preview_report.py — Generate a ViralQuest HTML report from a pipeline JSON.

Usage:
    python preview_report.py                        # reads example.json, writes preview.html
    python preview_report.py out.html               # custom output path
    python preview_report.py data.json out.html     # custom input + output
    python preview_report.py out.html --open        # open in browser automatically
"""

import json
import sys
import webbrowser
from pathlib import Path

_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT))

from viralquest.html_report import write_report

_DEFAULT_JSON = _ROOT / "example.json"


def main() -> None:
    args = [a for a in sys.argv[1:] if a != "--open"]
    open_browser = "--open" in sys.argv

    if len(args) == 0:
        json_path = _DEFAULT_JSON
        out = Path("preview.html")
    elif len(args) == 1:
        p = Path(args[0])
        if p.suffix == ".json":
            json_path = p
            out = Path("preview.html")
        else:
            json_path = _DEFAULT_JSON
            out = p
    else:
        json_path = Path(args[0])
        out = Path(args[1])

    if not json_path.exists():
        sys.exit(f"Error: JSON file not found: {json_path}")

    print(f"Loading {json_path}")
    report = json.loads(json_path.read_text(encoding="utf-8"))

    print(f"Building report → {out}")
    write_report(report, out)
    print(f"Done. Open in browser: file://{out.resolve()}")

    if open_browser:
        webbrowser.open(f"file://{out.resolve()}")


if __name__ == "__main__":
    main()
