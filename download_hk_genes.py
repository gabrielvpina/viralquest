#!/usr/bin/env python3
"""
Download housekeeping gene mRNA sequences from NCBI Nucleotide and merge
into the ViralQuest kingdom FASTA files, maximising species diversity.

Strategy
--------
For each query the script fetches a large pool of candidate UIDs (--pool,
default 500), retrieves their document summaries to extract TaxId (species
identifier) and accession prefix, then keeps at most --per-species records
per species — preferring curated RefSeq (NM_) over predicted (XM_) over
GenBank. This avoids the default NCBI relevance ranking that floods results
with a handful of model organisms.

Usage
-----
    python download_hk_genes.py --email your@email.com
    python download_hk_genes.py --email your@email.com --api-key YOUR_KEY
    python download_hk_genes.py --email your@email.com --kingdoms fish human mouse
    python download_hk_genes.py --email your@email.com --pool 1000 --per-species 2
    python download_hk_genes.py --email your@email.com --dry-run

NCBI API key (free): https://www.ncbi.nlm.nih.gov/account/
With key: 10 req/s.  Without: 3 req/s (slower but works fine).
"""

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR = Path(__file__).parent / "data" / "salmon-quant"

# ── Query catalogue ───────────────────────────────────────────────────────────
# biomol_mrna[Properties]  — mRNA molecules only (most reliable filter).
# [SLEN]                   — exclude fragments < 300 bp and genomic > N bp.
# Bacteria: no mRNA filter (prokaryotic records are mostly annotated as
#           genomic DNA/CDS). NOT "complete genome" guards against chromosomes.

KINGDOMS: dict[str, dict] = {

    "arthropods": {
        "file": DATA_DIR / "arthropods.fasta",
        "queries": [
            (
                "beta-tubulin[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:6000[SLEN]",
                "beta-tubulin",
            ),
            (
                "ubiquitin[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "ubiquitin",
            ),
            (
                "cyclophilin[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "cyclophilin",
            ),
            (
                "TATA-box binding protein[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "TBP",
            ),
            (
                "ribosomal protein S18[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPS18",
            ),
            (
                "ribosomal protein S3[Title] AND Arthropoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPS3",
            ),
        ],
    },

    "bacteria": {
        "file": DATA_DIR / "bacteria.fasta",
        "queries": [
            (
                "DNA gyrase subunit B[Title] AND Bacteria[Organism]"
                ' AND 300:4000[SLEN] NOT "complete genome"[Title]',
                "gyrB",
            ),
            (
                "recA[Title] AND Bacteria[Organism]"
                ' AND 300:3000[SLEN] NOT "complete genome"[Title]',
                "recA",
            ),
            (
                "dnaK[Title] AND Bacteria[Organism]"
                ' AND 300:4000[SLEN] NOT "complete genome"[Title]',
                "dnaK",
            ),
            (
                "superoxide dismutase[Title] AND Bacteria[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "sodA",
            ),
            (
                "translation initiation factor IF-2[Title] AND Bacteria[Organism]"
                ' AND 300:4000[SLEN] NOT "complete genome"[Title]',
                "infB",
            ),
        ],
    },

    "fish": {
        "file": DATA_DIR / "fish.fasta",
        "queries": [
            (
                "hypoxanthine phosphoribosyltransferase[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:5000[SLEN]",
                "HPRT1",
            ),
            (
                "TATA-box binding protein[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "TBP",
            ),
            (
                "succinate dehydrogenase complex subunit A[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:5000[SLEN]",
                "SDHA",
            ),
            (
                "peptidylprolyl isomerase A[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "PPIA",
            ),
            (
                "ubiquitin[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "ubiquitin",
            ),
            (
                "14-3-3 protein zeta[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "YWHAZ",
            ),
            (
                "ribosomal protein L32[Title] AND Actinopterygii[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPL32",
            ),
        ],
    },

    "fungi": {
        "file": DATA_DIR / "fungi.fasta",
        "queries": [
            (
                "glyceraldehyde-3-phosphate dehydrogenase[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "GAPDH",
            ),
            (
                "calmodulin[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "calmodulin",
            ),
            (
                "histone H3[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:2000[SLEN]",
                "histone H3",
            ),
            (
                "ubiquitin[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "ubiquitin",
            ),
            (
                "TATA-box binding protein[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "TBP",
            ),
            (
                "phosphoglycerate kinase[Title] AND Fungi[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "PGK",
            ),
        ],
    },

    "mammals": {
        "file": DATA_DIR / "mammals.fasta",
        "queries": [
            (
                "hypoxanthine phosphoribosyltransferase 1[Title] AND Mammalia[Organism]"
                " AND biomol_mrna[Properties] AND 300:5000[SLEN]"
                ' NOT "Homo sapiens"[Organism] NOT "Mus musculus"[Organism]',
                "HPRT1",
            ),
            (
                "14-3-3 protein zeta[Title] AND Mammalia[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]"
                ' NOT "Homo sapiens"[Organism] NOT "Mus musculus"[Organism]',
                "YWHAZ",
            ),
            (
                "ubiquitin C[Title] AND Mammalia[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]"
                ' NOT "Homo sapiens"[Organism] NOT "Mus musculus"[Organism]',
                "UBC",
            ),
        ],
    },

    "nematodes": {
        "file": DATA_DIR / "nematodes.fasta",
        "queries": [
            (
                "beta-tubulin[Title] AND Nematoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "beta-tubulin",
            ),
            (
                "ubiquitin[Title] AND Nematoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "ubiquitin",
            ),
            (
                "ribosomal protein L32[Title] AND Nematoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPL32",
            ),
            (
                "synaptobrevin[Title] AND Nematoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "synaptobrevin",
            ),
            (
                "ribosomal protein S18[Title] AND Nematoda[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPS18",
            ),
        ],
    },

    "plants": {
        "file": DATA_DIR / "plants.fasta",
        "queries": [
            (
                "glyceraldehyde-3-phosphate dehydrogenase[Title] AND Embryophyta[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "GAPDH",
            ),
            (
                "TIP41-like protein[Title] AND Embryophyta[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "TIP41",
            ),
            (
                "SAND family protein[Title] AND Embryophyta[Organism]"
                " AND biomol_mrna[Properties] AND 300:4000[SLEN]",
                "SAND",
            ),
            (
                "elongation factor 1-beta[Title] AND Embryophyta[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "EF1beta",
            ),
            (
                "ribosomal protein L13a[Title] AND Embryophyta[Organism]"
                " AND biomol_mrna[Properties] AND 300:3000[SLEN]",
                "RPL13a",
            ),
        ],
    },

    "human": {
        "file": DATA_DIR / "human.fasta",
        "queries": [
            ('ACTB[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "ACTB"),
            ('GAPDH[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "GAPDH"),
            ('HPRT1[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "HPRT1"),
            ('RPL13A[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "RPL13A"),
            ('B2M[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "B2M"),
            ('YWHAZ[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "YWHAZ"),
            ('TBP[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "TBP"),
            ('SDHA[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:5000[SLEN]', "SDHA"),
            ('PPIA[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "PPIA"),
            ('RPLP0[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "RPLP0"),
            ('UBC[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "UBC"),
            ('HMBS[Title] AND "Homo sapiens"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "HMBS"),
        ],
    },

    "mouse": {
        "file": DATA_DIR / "mouse.fasta",
        "queries": [
            ('Actb[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Actb"),
            ('Gapdh[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Gapdh"),
            ('Hprt1[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "Hprt1"),
            ('Rpl13a[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "Rpl13a"),
            ('B2m[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "B2m"),
            ('Ywhaz[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Ywhaz"),
            ('Tbp[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Tbp"),
            ('Sdha[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:5000[SLEN]', "Sdha"),
            ('Ppia[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:3000[SLEN]', "Ppia"),
            ('Rplp0[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Rplp0"),
            ('Ubc[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Ubc"),
            ('Eef1a1[Title] AND "Mus musculus"[Organism] AND biomol_mrna[Properties] AND 300:4000[SLEN]', "Eef1a1"),
        ],
    },
}

# ── NCBI helpers ──────────────────────────────────────────────────────────────

def _accession_rank(caption: str) -> int:
    """Lower is better. NM_ (curated RefSeq) > XM_ (predicted) > others."""
    if caption.startswith("NM_"):
        return 0
    if caption.startswith("XM_"):
        return 1
    return 2


def _search_ids(query: str, pool: int, delay: float) -> list[str]:
    """Run esearch and return up to pool UIDs."""
    time.sleep(delay)
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=pool)
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def _summarise(uid_list: list[str], delay: float) -> list[dict]:
    """Fetch esummary records for uid_list in batches of 500."""
    results: list[dict] = []
    for i in range(0, len(uid_list), 500):
        batch = uid_list[i:i + 500]
        time.sleep(delay)
        try:
            handle = Entrez.esummary(db="nucleotide", id=",".join(batch))
            records = Entrez.read(handle)
            handle.close()
            results.extend(records)
        except Exception as exc:
            logger.warning(f"    esummary batch {i//500 + 1} failed: {exc} — skipping batch")
    return results


def _pick_diverse_uids(
    uid_list: list[str],
    delay: float,
    per_species: int,
) -> tuple[list[str], int]:
    """
    From uid_list return at most per_species UIDs per unique TaxId,
    preferring NM_ > XM_ > GenBank. Returns (selected_uids, n_species).
    Falls back to uid_list[:per_species*500] if esummary fails entirely.
    """
    if not uid_list:
        return [], 0

    summaries = _summarise(uid_list, delay)
    if not summaries:
        logger.warning("    esummary returned nothing — using raw UID list (no diversity filter)")
        return uid_list, 0

    # Group by TaxId → [(rank, uid), ...]
    by_taxid: dict[str, list[tuple[int, str]]] = {}
    for rec in summaries:
        uid     = str(rec.get("Id", ""))
        taxid   = str(rec.get("TaxId", "0"))
        caption = str(rec.get("Caption", ""))
        rank    = _accession_rank(caption)
        if uid and taxid:
            by_taxid.setdefault(taxid, []).append((rank, uid))

    # For each species: sort by rank, keep best per_species
    selected: list[str] = []
    for taxid, entries in by_taxid.items():
        entries.sort(key=lambda x: x[0])          # best rank first
        for _, uid in entries[:per_species]:
            selected.append(uid)

    return selected, len(by_taxid)


def _fetch_fasta(uid_list: list[str], delay: float, retries: int = 3) -> list[SeqRecord]:
    """Fetch FASTA records for uid_list, with retry on failure."""
    if not uid_list:
        return []
    for attempt in range(1, retries + 1):
        try:
            time.sleep(delay)
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(uid_list),
                rettype="fasta",
                retmode="text",
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as exc:
            wait = delay * (2 ** attempt)
            logger.warning(f"    efetch attempt {attempt}/{retries} failed: {exc} — retry in {wait:.1f}s")
            time.sleep(wait)
    logger.error("    All fetch attempts failed — skipping.")
    return []


# ── Existing-accession loader ──────────────────────────────────────────────────

def _load_existing_accessions(fasta_path: Path) -> set[str]:
    if not fasta_path.exists():
        return set()
    accessions: set[str] = set()
    with open(fasta_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                accessions.add(line[1:].split()[0].strip())
    return accessions


# ── Core loop ─────────────────────────────────────────────────────────────────

def run(args: argparse.Namespace) -> None:
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    delay = 0.11 if args.api_key else 0.35   # NCBI rate limits

    kingdoms_to_run = args.kingdoms or list(KINGDOMS.keys())
    unknown = set(kingdoms_to_run) - set(KINGDOMS.keys())
    if unknown:
        logger.error(f"Unknown kingdom(s): {', '.join(sorted(unknown))}")
        logger.info(f"Available: {', '.join(KINGDOMS.keys())}")
        sys.exit(1)

    grand_total = 0

    for kingdom in kingdoms_to_run:
        cfg        = KINGDOMS[kingdom]
        fasta_path = cfg["file"]
        queries    = cfg["queries"]

        logger.info(f"\n{'─' * 62}")
        logger.info(f"Kingdom : {kingdom}  ({len(queries)} gene queries)")
        logger.info(f"File    : {fasta_path}")

        existing = _load_existing_accessions(fasta_path)
        logger.info(f"Existing: {len(existing)} sequences")

        new_records: list[SeqRecord] = []

        for query, gene_label in queries:

            # 1 ── Search: pull a large candidate pool
            logger.info(f"  [{gene_label}] searching (pool={args.pool})…")
            try:
                uid_pool = _search_ids(query, args.pool, delay)
            except Exception as exc:
                logger.error(f"  [{gene_label}] esearch failed: {exc}")
                time.sleep(2.0)
                continue

            if not uid_pool:
                logger.info(f"  [{gene_label}] no results in NCBI.")
                continue

            logger.info(f"  [{gene_label}] {len(uid_pool)} candidates — selecting diverse species…")

            # 2 ── Diversity filter: ≤ per_species per TaxId, best accession type
            selected_uids, n_species = _pick_diverse_uids(uid_pool, delay, args.per_species)
            if n_species:
                logger.info(f"  [{gene_label}] {n_species} unique species → {len(selected_uids)} selected UIDs")

            # 3 ── Remove UIDs whose accession is already in the FASTA
            #      We can't know the accession before fetching, but we skip
            #      after fetching based on rec.id (accession).
            # Fetch in one batch (already a manageable number after diversity filter)
            records = _fetch_fasta(selected_uids, delay)

            added = skipped = 0
            for rec in records:
                acc = rec.id
                if acc in existing:
                    skipped += 1
                else:
                    existing.add(acc)
                    new_records.append(rec)
                    added += 1

            logger.info(f"  [{gene_label}] +{added} new  |  {skipped} already present")

        if not new_records:
            logger.info(f"  → Nothing new for {kingdom}.")
            continue

        logger.info(f"  → {len(new_records)} new sequence(s) total for {kingdom}.")

        if args.dry_run:
            logger.info(f"  [dry-run] would write to {fasta_path.name}")
            for rec in new_records[:8]:
                logger.info(f"    {rec.id:<20} {rec.description[:55]}")
            if len(new_records) > 8:
                logger.info(f"    … and {len(new_records) - 8} more")
            continue

        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        mode = "a" if fasta_path.exists() else "w"
        with open(fasta_path, mode, encoding="utf-8") as fh:
            SeqIO.write(new_records, fh, "fasta")

        logger.success(f"  Written → {fasta_path.name}  (+{len(new_records)} sequences)")
        grand_total += len(new_records)

    logger.info(f"\n{'=' * 62}")
    if args.dry_run:
        logger.info("Dry-run complete — no files were modified.")
    else:
        logger.success(f"Done. {grand_total} new sequence(s) added across all kingdoms.")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--email", required=True,
        help="Your email address (required by NCBI Entrez policy).",
    )
    parser.add_argument(
        "--api-key",
        help="NCBI API key — allows 10 req/s instead of 3. Free at ncbi.nlm.nih.gov/account/",
    )
    parser.add_argument(
        "--pool", type=int, default=500,
        help="UIDs to retrieve per query before diversity filtering (default: 500).",
    )
    parser.add_argument(
        "--per-species", type=int, default=1, dest="per_species",
        help="Maximum sequences to keep per species per query (default: 1).",
    )
    parser.add_argument(
        "--kingdoms", nargs="+", metavar="KINGDOM",
        help=f"Kingdoms to process (default: all). Choices: {', '.join(KINGDOMS.keys())}",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be downloaded without writing any files.",
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
