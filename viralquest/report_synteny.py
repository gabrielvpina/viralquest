"""
report_synteny.py — General synteny across the samples of ``viralquest-report``.

Each sample's synteny step already picked a reference genome per sequence (its
best viral BLASTn hit). Here those references are pooled over every sample and
grouped by similarity — DIAMOND blastp of all their proteins against each other
— so different accessions of one virus (strains, isolates, partial records)
fall into one group. Each group gets one anchor reference (the one most
sequences already used), and every sequence of the group, from any sample, is
aligned again to that anchor with the per-sample machinery (synteny.py),
contig ↔ contig links included. The Synteny tab of the consolidated report
draws the result exactly like the per-sample tab.

Files (``<outdir>/features/synteny/``):
    refs/<accession>.gb   GenBank records, copied from the sample folders
                          (downloaded only when no sample has them)
    proteome_*            reference ↔ reference blastp (grouping)
    orfs.faa · refs.faa · blastp.tsv · contigs_*   ORF ↔ anchor alignments
    synteny.tsv           ORF ↔ CDS pairs on the anchors (sequence = gid)

Report data: ``general_synteny = {references: {acc: RefGenome}, clusters:
[{reference, members: [gid], samples, merged: [...]}]}`` and, on each
sequence, ``general_synteny`` (a SyntenyResult on the anchor of its group).
"""

from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path

from loguru import logger

from viralquest.blastn import acc_key
from viralquest.synteny import (
    GenBankFetcher,
    ParsedReference,
    SyntenyAligner,
    SyntenyAnalyzer,
    _seq_from_json,
    best_viral_hit,
    parse_genbank,
)

MIN_IDENTITY = 80.0   # bit-score weighted aa identity of two references' best protein hits (~ same species)
MIN_COVERAGE = 50.0   # % of the smaller-covered proteome those hits span

# Where a sample's synteny step leaves its GenBank cache (current, pre-3.0.14).
_SAMPLE_REF_DIRS = ("features/synteny/refs", "synteny/refs")


# ── References ──────────────────────────────────────────────────────────────

def _collect_genbank(samples, cache_dir: Path) -> None:
    """Copy the samples' cached GenBank records into *cache_dir* (no download)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    for s in samples:
        for rel in _SAMPLE_REF_DIRS:
            for gb in (Path(s.directory) / rel).glob("*.gb"):
                dst = cache_dir / gb.name
                if not dst.exists():
                    shutil.copyfile(gb, dst)


def _members_by_ref(sequences: list[dict], known: set[str]) -> dict[str, list[dict]]:
    """
    {reference key: sequences on it}: each sequence's per-sample synteny
    reference, else (sample without the synteny step) its best viral BLASTn
    hit when that is one of the *known* references.
    """
    out: dict[str, list[dict]] = defaultdict(list)
    for q in sequences:
        syn = q.get("synteny") or {}
        key = acc_key(syn["reference"]) if syn.get("reference") else None
        if key is None and q.get("blastn_hits"):
            hit = best_viral_hit(_seq_from_json(q))
            key = acc_key(hit.accession) if hit and acc_key(hit.accession) in known else None
        if key in known:
            out[key].append(q)
    return out


# ── Grouping ────────────────────────────────────────────────────────────────

def _similar_pairs(refs: dict[str, ParsedReference], aligner: SyntenyAligner,
                   min_identity: float, min_coverage: float) -> list[tuple[str, str]]:
    """
    Pairs of references whose proteomes match: for each protein of A, its best
    HSP on any protein of B; A ~ B when those hits span ≥ *min_coverage* % of A's
    proteome (or B's, the other way) at ≥ *min_identity* % weighted identity.
    """
    prots = {f"{k}|{i}": aa for k, r in refs.items() for i, aa in enumerate(r.proteins) if aa}
    if len({p.rpartition("|")[0] for p in prots}) < 2:
        return []
    rows = aligner.align(prots, prots, tag="proteome")
    total = defaultdict(int)
    for pid, aa in prots.items():
        total[pid.rpartition("|")[0]] += len(aa)

    best: dict[tuple[str, str], dict] = {}          # (query protein, subject ref) → best HSP
    for r in rows:
        qa, sb = r["qseqid"].rpartition("|")[0], r["sseqid"].rpartition("|")[0]
        if qa == sb:
            continue
        k = (r["qseqid"], sb)
        if k not in best or r["bitscore"] > best[k]["bitscore"]:
            best[k] = r
    span: dict[tuple[str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    for (qp, sb), r in best.items():
        acc = span[(qp.rpartition("|")[0], sb)]
        acc[0] += r["qend"] - r["qstart"] + 1          # aligned aa
        acc[1] += r["pident"] * r["bitscore"]
        acc[2] += r["bitscore"]

    pairs = []
    for a, b in sorted({tuple(sorted(k)) for k in span}):
        ab, ba = span.get((a, b)), span.get((b, a))
        cov = max(100 * ab[0] / max(1, total[a]) if ab else 0,
                  100 * ba[0] / max(1, total[b]) if ba else 0)
        bits = (ab[2] if ab else 0) + (ba[2] if ba else 0)
        ident = ((ab[1] if ab else 0) + (ba[1] if ba else 0)) / bits if bits else 0
        if cov >= min_coverage and ident >= min_identity:
            pairs.append((a, b))
    return pairs


def _groups(keys: list[str], pairs: list[tuple[str, str]]) -> list[list[str]]:
    """Connected components (union–find) of the reference keys."""
    parent = {k: k for k in keys}

    def find(k: str) -> str:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for a, b in pairs:
        parent[find(a)] = find(b)
    out: dict[str, list[str]] = defaultdict(list)
    for k in keys:
        out[find(k)].append(k)
    return list(out.values())


# ── Public entry point ──────────────────────────────────────────────────────

def build_general_synteny(samples, outdir: str | Path, threads: int = 2,
                          email: str | None = None,
                          min_identity: float = MIN_IDENTITY, min_coverage: float = MIN_COVERAGE,
                          fetcher: GenBankFetcher | None = None,
                          aligner: SyntenyAligner | None = None) -> dict | None:
    """
    Group the samples' synteny references, re-align every sequence of a group
    to its anchor and set ``general_synteny`` on each sequence dict (in place).
    Returns the report's ``general_synteny`` dict, or None when no sample has
    synteny references.
    """
    from viralquest.exporter import ReportExporter, _to_serializable

    outdir = Path(outdir)
    meta: dict[str, dict] = {}                           # reference key → RefGenome dict
    for s in samples:
        for acc, ref in ((s.report.get("synteny") or {}).get("references") or {}).items():
            meta.setdefault(acc_key(acc), ref)
    if not meta:
        logger.info("General synteny: no sample has synteny references — skipped.")
        return None

    sequences = [q for s in samples for q in s.sequences]
    members   = _members_by_ref(sequences, set(meta))
    keys      = [k for k in meta if members.get(k)]
    if not keys:
        return None

    fetcher = fetcher or GenBankFetcher(outdir / "refs", email=email)
    _collect_genbank(samples, fetcher.cache_dir)
    paths = fetcher.fetch([meta[k]["accession"] for k in keys])
    refs  = {k: p for k in keys if k in paths and (p := parse_genbank(paths[k]))}
    if not refs:
        logger.warning("General synteny: no reference GenBank record could be read — skipped.")
        return None

    aligner = aligner or SyntenyAligner(outdir, threads=threads)
    pairs   = _similar_pairs(refs, aligner, min_identity, min_coverage)
    groups  = _groups(list(refs), pairs)

    # Anchor: the reference most sequences already sit on, then the most CDS, the longest.
    assign: dict[str, str] = {}
    merged: dict[str, list[dict]] = {}
    for g in groups:
        anchor = max(g, key=lambda k: (len(members[k]), len(refs[k].genome.cds),
                                       refs[k].genome.length, k))
        acc = refs[anchor].genome.accession
        merged[acc] = [{"accession": refs[k].genome.accession, "organism": refs[k].genome.organism,
                        "length": refs[k].genome.length, "n_members": len(members[k])}
                       for k in sorted(g, key=lambda k: (k != anchor, -len(members[k]), k))]
        for k in g:
            for q in members[k]:
                assign[q["gid"]] = acc

    by_gid = {q["gid"]: q for q in sequences}
    seqs   = [_seq_from_json(by_gid[g], id_key="gid") for g in assign]
    used   = SyntenyAnalyzer(outdir, fetcher=fetcher, aligner=aligner).run(seqs, assign=assign)
    for s in seqs:
        by_gid[s.id]["general_synteny"] = _to_serializable(s.synteny) if s.synteny else None

    out = ReportExporter._build_synteny(seqs, used)
    for c in out["clusters"]:
        c["samples"] = sorted({by_gid[g].get("sample") for g in c["members"]})
        c["merged"]  = merged.get(c["reference"], [])
    # Groups shared by more samples first, then by size.
    out["clusters"].sort(key=lambda c: (-len(c["samples"]), -len(c["members"]), c["reference"]))
    n_multi = sum(1 for c in out["clusters"] if len(c["samples"]) > 1)
    logger.success(
        f"General synteny: {len(refs)} reference(s) in {len(out['clusters'])} group(s), "
        f"{n_multi} shared by several samples; {len(seqs)} sequence(s) on their group's anchor."
    )
    return out
