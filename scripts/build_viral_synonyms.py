#!/usr/bin/env python3
"""
build_viral_synonyms.py — Build data/viralSynonyms.json.xz from the NCBI taxdump.

Maintainer tool, run once whenever the bundled taxonomy is refreshed; the
pipeline itself never touches the network.  It keeps, for every taxid of the
bundled viral taxonomy (data/viralTax.json.xz), the alternative names NCBI
knows it by — old names, equivalent names and common names — so a BLAST title
such as "Orchid fleck virus genomic RNA ..." still lands on its current taxon
(*Dichorhavirus orchidaceae*, taxid 152177).

Usage
-----
    python scripts/build_viral_synonyms.py                      # downloads taxdump
    python scripts/build_viral_synonyms.py --taxdump taxdump.tar.gz

Rules
-----
* name classes kept: synonym, equivalent name, genbank common name, common name
  (acronyms are left out: short tokens like "OFV" collide with unrelated text);
* only taxids present in the bundled viral taxonomy (merged taxids followed
  through merged.dmp), so every synonym resolves to a lineage;
* a name already used as a scientific / species name in the bundle is skipped
  (the primary name wins);
* a name pointing to several taxids is kept only when they all share the same
  species (the lowest taxid is used), otherwise it is dropped as ambiguous;
* names shorter than 4 characters, without a letter, or generic ("virus") are skipped.

Output: {"meta": {...}, "names": {lower-case name: taxid}} as xz-compressed JSON.
"""

from __future__ import annotations

import argparse
import io
import json
import lzma
import subprocess
import sys
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
KEEP_CLASSES = ("synonym", "equivalent name", "genbank common name", "common name")
GENERIC = {"virus", "viruses", "phage", "phages", "viroid", "viroids", "satellite", "satellites"}

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_VIRAL_TAX = ROOT / "data" / "viralTax.json.xz"
DEFAULT_OUT = ROOT / "data" / "viralSynonyms.json.xz"


def _dmp_rows(fh):
    """Rows of an NCBI .dmp file ("\\t|\\t"-separated, "\\t|" terminated)."""
    for line in io.TextIOWrapper(fh, encoding="utf-8"):
        yield [f.strip() for f in line.rstrip("\n").rstrip("|").split("\t|\t")]


def build(taxdump: Path, viral_tax: Path) -> dict:
    with lzma.open(viral_tax, "rt", encoding="utf-8") as fh:
        records = json.load(fh)
    by_taxid = {int(r["TaxId"]): r for r in records if r.get("TaxId")}
    primary = set()
    for r in records:
        for key in ("ScientificName", "Species"):
            if r.get(key):
                primary.add(r[key].lower())

    with tarfile.open(taxdump, "r:gz") as tar:
        # merged.dmp: old_taxid → new_taxid.  The bundle may carry either id, so
        # map every current id back to the bundle's id(s).
        current_to_bundle: dict[int, set[int]] = defaultdict(set)
        for t in by_taxid:
            current_to_bundle[t].add(t)
        for old, new, *_ in _dmp_rows(tar.extractfile("merged.dmp")):
            old, new = int(old), int(new)
            if old in by_taxid:
                current_to_bundle[new].add(old)

        names: dict[str, set[int]] = defaultdict(set)
        kept_by_class: dict[str, int] = defaultdict(int)
        for tid, name, _unique, cls, *_ in _dmp_rows(tar.extractfile("names.dmp")):
            if cls not in KEEP_CLASSES:
                continue
            bundle_ids = current_to_bundle.get(int(tid))
            if not bundle_ids:
                continue
            key = " ".join(name.split()).lower()
            if (len(key) < 4 or not any(c.isalpha() for c in key)
                    or key in GENERIC or key in primary):
                continue
            names[key].update(bundle_ids)
            kept_by_class[cls] += 1

    out: dict[str, int] = {}
    ambiguous = 0
    for key, tids in names.items():
        if len(tids) == 1:
            out[key] = next(iter(tids))
            continue
        species = {by_taxid[t].get("Species") or by_taxid[t].get("ScientificName") for t in tids}
        if len(species) == 1:
            out[key] = min(tids)
        else:
            ambiguous += 1

    return {
        "meta": {
            "source": TAXDUMP_URL,
            "built": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
            "taxdump_mtime": datetime.fromtimestamp(taxdump.stat().st_mtime, timezone.utc).strftime("%Y-%m-%d"),
            "name_classes": list(KEEP_CLASSES),
            "names": len(out),
            "taxids": len(set(out.values())),
            "dropped_ambiguous": ambiguous,
            "kept_by_class": dict(kept_by_class),
        },
        "names": dict(sorted(out.items())),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--taxdump", type=Path, help="local taxdump.tar.gz (downloaded when omitted)")
    ap.add_argument("--viral-tax", type=Path, default=DEFAULT_VIRAL_TAX)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    taxdump = args.taxdump
    if taxdump is None:
        taxdump = Path("taxdump.tar.gz")
        print(f"Downloading {TAXDUMP_URL} …", file=sys.stderr)
        subprocess.run(["curl", "-fSL", "--progress-bar", "-o", str(taxdump), TAXDUMP_URL], check=True)

    data = build(taxdump, args.viral_tax)
    with lzma.open(args.out, "wt", encoding="utf-8", preset=9) as fh:
        json.dump(data, fh, separators=(",", ":"))
    meta = data["meta"]
    print(f"{args.out}: {meta['names']} names → {meta['taxids']} taxids "
          f"({meta['dropped_ambiguous']} ambiguous dropped), "
          f"{args.out.stat().st_size / 1e6:.2f} MB", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
