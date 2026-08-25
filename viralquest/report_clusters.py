"""
report_clusters.py — Cross-sample clustering for ``viralquest-report``.

Pools every confirmed sequence from all loaded samples and groups sequences
that are highly similar *across* samples into "general clusters".  This single
structure feeds both the network graph (Section 2) and the RNA quantification
boxplots (Section 4).

Method (similarity-first greedy clustering)
-------------------------------------------
1. Pool every confirmed sequence carrying a nucleotide sequence.
2. Greedily, longest-first: take the longest unassigned sequence as the
   representative, BLASTn-align all remaining sequences to it, and pull in
   those passing the identity/coverage floor.  Repeat on what's left.
3. A real cluster is sequence-similarity defined and gets a brand-new id
   ``VQR_CLU_####`` — two sequences that merely share a best-hit species label
   are NOT merged unless they actually align.
4. Keep only clusters whose members span ≥2 distinct samples.

Species is annotation, not a partition: the representative's best BLASTx
species (NR first, RefSeq fallback, then best BLASTn subject) labels the
cluster, while each member keeps its own best species so the report can flag
whether every member resolved to the same species or not.

The BLASTn floor (default 90% identity / 70% coverage) is the *lowest* the
report's client-side sliders can ever show; every surviving member carries its
identity + coverage so the UI filters in-browser without re-running BLASTn.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path

from loguru import logger

from viralquest.setup_env import resolve_tool

# Build-time floor; also the minimum the report sliders expose.
FLOOR_IDENTITY = 90.0
FLOOR_COVERAGE = 70.0

_OUTFMT = "6 qseqid pident length qlen qstart qend sstart send bitscore"


# ── Data structures ────────────────────────────────────────────────────────

@dataclass
class GeneralMember:
    gid:               str           # {sample}::{id}
    sample:            str
    seq_id:            str           # original id (for display / viewer jump)
    length:            int
    identity:          float         # vs representative (100.0 for the rep)
    coverage:          float         # query coverage vs representative
    species:           str | None    # this member's own best species label
    species_db:        str | None    # "nr" | "refseq" | "blastn" | None
    tpm:               float | None
    # Alignment coordinates vs the representative (for the alignment-bar view).
    aln_start:         int = 1       # subject (representative) start
    aln_end:           int = 1       # subject (representative) end
    query_start:       int = 1
    query_end:         int = 1
    is_representative: bool = False


@dataclass
class GeneralCluster:
    gid:            str              # VQR_CLU_####
    species:        str             # representative label (graph still hides it)
    species_db:     str | None
    representative: str             # representative gid
    members:        list[GeneralMember] = field(default_factory=list)
    samples:        list[str] = field(default_factory=list)
    homogeneous:    bool = True     # all members share the same NR/RefSeq species?

    def to_dict(self) -> dict:
        d = asdict(self)
        d["size"] = len(self.members)
        d["cluster_id"] = self.gid          # alias used by the alignment-bar view
        for m in d["members"]:
            m["query_coverage"] = m["coverage"]   # field name expected by the view
        return d


# ── Species / sequence extraction ──────────────────────────────────────────

def _best_species(seq: dict) -> tuple[str | None, str | None]:
    """(label, db) — NR BLASTx, then RefSeq BLASTx, then best BLASTn subject."""
    for key, db in (("blastx_nr_hits", "nr"), ("blastx_hits", "refseq")):
        hits = seq.get(key) or []
        if hits and hits[0].get("species"):
            return hits[0]["species"], db
    bn = seq.get("blastn_hits") or []
    if bn:
        label = bn[0].get("stitle") or bn[0].get("accession")
        if label:
            return label, "blastn"
    return None, None


# ── BLASTn (subject = representative) ──────────────────────────────────────

def _run_blastn(rep_seq: str, members: list[dict], blastn_bin: str, min_identity: float) -> dict[str, dict]:
    """Align each member's nucleotide sequence to the representative; best hit per member."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path  = Path(tmp)
        rep_fa    = tmp_path / "rep.fasta"
        qry_fa    = tmp_path / "members.fasta"
        rep_fa.write_text(f">rep\n{rep_seq}\n", encoding="utf-8")
        with open(qry_fa, "w", encoding="utf-8") as fh:
            for m in members:
                fh.write(f">{m['gid']}\n{m.get('sequence', '')}\n")

        cmd = [
            blastn_bin,
            "-query", str(qry_fa), "-subject", str(rep_fa),
            "-outfmt", _OUTFMT, "-task", "blastn",
            "-perc_identity", str(min_identity), "-dust", "no",
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            logger.error(f"blastn failed: {exc.stderr[:300]}")
            return {}
        except FileNotFoundError:
            raise RuntimeError(
                f"blastn binary not found: {blastn_bin!r}. "
                "Run 'viralquest-setup' to install BLAST+ via pixi, or pass "
                "--blastn-bin with an explicit path."
            )

    best: dict[str, dict] = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        qseqid, pident, _aln, qlen, qstart, qend, sstart, send, bitscore = parts
        qlen_i = int(qlen)
        qs, qe = int(qstart), int(qend)
        ss, se = int(sstart), int(send)
        aln_len = abs(qe - qs) + 1
        qcov = round(aln_len / qlen_i * 100, 2) if qlen_i > 0 else 0.0
        hit = {
            "pident":   round(float(pident), 2),
            "qcov":     qcov,
            "bitscore": float(bitscore),
            "qstart":   min(qs, qe), "qend": max(qs, qe),     # normalise strand
            "sstart":   min(ss, se), "send": max(ss, se),
        }
        if qseqid not in best or hit["bitscore"] > best[qseqid]["bitscore"]:
            best[qseqid] = hit
    return best


def _member(seq: dict, hit: dict | None, is_rep: bool) -> GeneralMember:
    """Build a GeneralMember. *hit* is None for the representative (full-length)."""
    label, db = _best_species(seq)
    length = seq.get("length", 0)
    if is_rep or hit is None:
        return GeneralMember(
            gid=seq["gid"], sample=seq.get("sample", ""), seq_id=seq.get("id", ""),
            length=length, identity=100.0, coverage=100.0,
            species=label, species_db=db, tpm=seq.get("tpm"),
            aln_start=1, aln_end=length, query_start=1, query_end=length,
            is_representative=True,
        )
    return GeneralMember(
        gid=seq["gid"], sample=seq.get("sample", ""), seq_id=seq.get("id", ""),
        length=length, identity=hit["pident"], coverage=hit["qcov"],
        species=label, species_db=db, tpm=seq.get("tpm"),
        aln_start=hit["sstart"], aln_end=hit["send"],
        query_start=hit["qstart"], query_end=hit["qend"],
        is_representative=False,
    )


# ── Main clustering ─────────────────────────────────────────────────────────

def build_general_clusters(
    sequences:    list[dict],
    blastn_bin:   str | None = None,
    min_identity: float = FLOOR_IDENTITY,
    min_coverage: float = FLOOR_COVERAGE,
) -> list[GeneralCluster]:
    """
    Build cross-sample clusters from stamped sequence dicts (each needs
    ``gid``, ``sample``, ``id``, ``sequence``, ``length`` and blast hits).

    *blastn_bin* is an optional explicit binary; by default blastn is resolved
    from the pixi environment (same installation the main pipeline uses).
    """
    # Resolve once: pixi env first, then PATH. Fall back to the bare name so
    # _run_blastn raises its actionable RuntimeError when nothing is installed.
    blastn_bin = resolve_tool("blastn", blastn_bin) or blastn_bin or "blastn"
    logger.debug(f"Cross-sample clustering will use blastn: {blastn_bin}")

    # Greedy similarity clustering over the whole pool, longest-first.
    remaining = sorted(
        (s for s in sequences if s.get("sequence")),
        key=lambda s: s.get("length", 0), reverse=True,
    )

    clusters: list[GeneralCluster] = []
    counter = 0

    while len(remaining) >= 2:
        rep = remaining[0]
        rest = remaining[1:]
        hits = _run_blastn(rep.get("sequence", ""), rest, blastn_bin, min_identity)

        aligned = [_member(rep, None, is_rep=True)]
        leftover: list[dict] = []
        for m in rest:
            h = hits.get(m["gid"])
            if h and h["pident"] >= min_identity and h["qcov"] >= min_coverage:
                aligned.append(_member(m, h, is_rep=False))
            else:
                leftover.append(m)

        samples = sorted({m.sample for m in aligned})
        if len(aligned) >= 2 and len(samples) >= 2:
            counter += 1
            rep_label, rep_db = _best_species(rep)
            # Homogeneous = every member resolved to the representative's species
            # (only members with an actual NR/RefSeq species name are considered).
            homogeneous = all(
                (m.species == rep_label) for m in aligned
                if m.species_db in ("nr", "refseq")
            )
            clusters.append(GeneralCluster(
                gid=f"VQR_CLU_{counter:04d}",
                species=rep_label or "Unknown",
                species_db=rep_db,
                representative=rep["gid"],
                members=aligned,
                samples=samples,
                homogeneous=homogeneous,
            ))
            logger.success(
                f"VQR_CLU_{counter:04d} | {rep_label or 'Unknown'!r} | "
                f"{len(aligned)} members across {len(samples)} samples"
            )

        # The representative (and any same-cluster matches) are consumed; the
        # rest continue to the next round.
        remaining = leftover

    logger.success(f"Cross-sample clustering: {len(clusters)} general cluster(s).")
    return clusters


# ── Artefact writers ────────────────────────────────────────────────────────

def write_artifacts(
    clusters:   list[GeneralCluster],
    seq_by_gid: dict[str, dict],
    out_dir:    Path,
) -> None:
    """Write per-cluster FASTA + quantification TSV and the combined BLASTn table."""
    clu_dir = out_dir / "clusters"
    qnt_dir = out_dir / "quantification"
    bln_dir = out_dir / "blastn"
    for d in (clu_dir, qnt_dir, bln_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Combined member→representative table.
    with open(bln_dir / "general_blastn.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster_id", "gid", "sample", "seq_id", "length",
                    "identity", "coverage", "is_representative", "representative_gid"])
        for c in clusters:
            for m in c.members:
                w.writerow([c.gid, m.gid, m.sample, m.seq_id, m.length,
                            m.identity, m.coverage, int(m.is_representative), c.representative])

    for c in clusters:
        # FASTA of all member sequences.
        with open(clu_dir / f"{c.gid}.fasta", "w", encoding="utf-8") as fh:
            for m in c.members:
                seq = seq_by_gid.get(m.gid, {})
                tag = " [representative]" if m.is_representative else ""
                fh.write(f">{m.gid} {c.species}{tag}\n{seq.get('sequence', '')}\n")

        # Quantification table (TPM may be absent for samples without salmon).
        with open(qnt_dir / f"{c.gid}.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gid", "sample", "seq_id", "tpm", "identity", "coverage"])
            for m in c.members:
                w.writerow([m.gid, m.sample, m.seq_id,
                            "" if m.tpm is None else m.tpm, m.identity, m.coverage])

    logger.info(f"Cluster artefacts written → {out_dir}/ (clusters, quantification, blastn).")
