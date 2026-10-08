"""
synteny.py — Synteny of each viral sequence with its BLASTn reference genome.

For every sequence whose best BLASTn hit has a viral title, the reference
(the hit's subject) is downloaded from NCBI as GenBank (efetch, nuccore,
gbwithparts) and its CDS — with their translations — and mat_peptides are read.
The ORFs predicted on the input sequences are then aligned to the reference
proteins with DIAMOND blastp. Each HSP is mapped back to nucleotide intervals on
both genomes, which is what the report's Synteny tab draws as ribbons. ORFs
with no hit are the "non-canonical" ones the tab hides by default.

Sequences are grouped by reference accession: one synteny cluster per
reference, holding every sequence whose best viral hit is that record.

Files (``<outdir>/synteny/``):
    refs/<accession>.gb   GenBank cache — reused on later runs
    refs.faa              reference proteins (DIAMOND database input)
    orfs.faa              input ORFs (queries)
    blastp.tsv            raw DIAMOND blastp output
    synteny.tsv           kept ORF ↔ CDS pairs, one per line

Network or DIAMOND failures never stop the run: affected sequences simply get
no synteny.
"""

from __future__ import annotations

import csv
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from loguru import logger

from viralquest.biodata import (
    BlastnResult,
    NucSequence,
    RefFeature,
    RefGenome,
    SyntenyHit,
    SyntenyResult,
)
from viralquest.blastn import acc_key, ncbi_efetch
from viralquest.blastn_tax import is_viral_title

# Larger than any known viral genome: a viral-titled hit longer than this is
# not downloaded (e.g. a host chromosome carrying an integrated virus).
MAX_REF_LEN = 3_000_000

_BLASTP_FIELDS = ("qseqid", "sseqid", "pident", "length", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore", "qlen", "slen")


# ── Reference selection ─────────────────────────────────────────────────────

def best_viral_hit(seq: NucSequence) -> BlastnResult | None:
    """Highest bit-score BLASTn hit whose title looks viral (with an accession)."""
    hits = [
        h for h in seq.blastn_hits
        if h.accession and is_viral_title(h.stitle)
        and not (h.slen and h.slen > MAX_REF_LEN)
    ]
    return max(hits, key=lambda h: h.bit_score) if hits else None


# ── GenBank download + cache ────────────────────────────────────────────────

class GenBankFetcher:
    """
    Downloads GenBank records (nuccore, gbwithparts) into ``cache_dir`` as
    ``<accession-without-version>.gb``; records already there are not fetched
    again. Requests are batched and spaced for NCBI's 3 requests/s limit.
    `http` and `sleep` are injectable for tests.
    """

    BATCH = 20
    GAP   = 0.4

    def __init__(self, cache_dir: str | Path, email: str | None = None,
                 tool: str = "viralquest", http=None, sleep=None):
        self.cache_dir = Path(cache_dir)
        self.email     = email
        self.tool      = tool
        self._http     = http or (lambda params: ncbi_efetch(params, self.tool))
        self._sleep    = sleep or (lambda s: time.sleep(s))

    def path(self, accession: str) -> Path:
        return self.cache_dir / f"{acc_key(accession)}.gb"

    def fetch(self, accessions: list[str]) -> dict[str, Path]:
        """{accession key: .gb path} for every record available (cached or fetched)."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        keys    = list(dict.fromkeys(acc_key(a) for a in accessions))
        missing = [k for k in keys if not self.path(k).exists()]
        for i in range(0, len(missing), self.BATCH):
            chunk  = missing[i : i + self.BATCH]
            params = {"db": "nuccore", "id": ",".join(chunk),
                      "rettype": "gbwithparts", "retmode": "text", "tool": self.tool}
            if self.email:
                params["email"] = self.email
            if i:
                self._sleep(self.GAP)
            try:
                self._store(self._http(params))
            except Exception as exc:
                logger.warning(f"Synteny: GenBank download of {len(chunk)} record(s) failed: {exc}")
        found = {k: self.path(k) for k in keys if self.path(k).exists()}
        if len(found) < len(keys):
            lost = [k for k in keys if k not in found]
            logger.warning(f"Synteny: no GenBank record for {len(lost)} reference(s): {', '.join(lost[:5])}")
        return found

    def _store(self, text: str) -> None:
        """Split a multi-record GenBank response and write one file per record."""
        for chunk in re.split(r"^//\s*$", text, flags=re.MULTILINE):
            m = re.search(r"^ACCESSION\s+(\S+)", chunk, flags=re.MULTILINE)
            if not m:
                continue
            self.path(m.group(1)).write_text(chunk.strip("\n") + "\n//\n", encoding="utf-8")


# ── GenBank parsing ─────────────────────────────────────────────────────────

def _q(feature, key: str) -> str:
    vals = feature.qualifiers.get(key)
    return vals[0] if vals else ""


def _parts(location) -> list[tuple[int, int, int]]:
    """(start, end, strand) of each part, in biological (5'→3') order."""
    return [(int(p.start), int(p.end), -1 if p.strand == -1 else 1) for p in location.parts]


def _offset_to_genome(parts: list[tuple[int, int, int]], offset: int) -> int:
    """0-based offset along a (possibly joined) CDS → 0-based genome position."""
    for start, end, strand in parts:
        n = end - start
        if offset < n:
            return start + offset if strand == 1 else end - 1 - offset
        offset -= n
    start, end, strand = parts[-1]
    return end - 1 if strand == 1 else start


def cds_aa_to_genome(parts, codon_start: int, aa_from: int, aa_to: int) -> list[int]:
    """1-based aa range on a CDS → [start, end) envelope on the genome."""
    shift = codon_start - 1
    a = _offset_to_genome(parts, shift + 3 * (aa_from - 1))
    b = _offset_to_genome(parts, shift + 3 * aa_to - 1)
    return [min(a, b), max(a, b) + 1]


def aligned_segments(parts, codon_start: int, aa_from: int, aa_to: int,
                     orf_nt: list[int], orf_strand: str) -> list[list[int]]:
    """
    Split an alignment over a (possibly joined) CDS into one
    [ref_start, ref_end, orf_start, orf_end] per CDS part it touches. The ORF
    side is walked in step (5'→3' of the ORF), assuming a gapless alignment,
    and clamped to the aligned ORF stretch.
    """
    lo_off = codon_start - 1 + 3 * (aa_from - 1)
    hi_off = codon_start - 1 + 3 * aa_to
    segs, base, used = [], 0, 0
    for start, end, strand in parts:
        n = end - start
        a, b = max(lo_off, base), min(hi_off, base + n)
        if a < b:
            ref = ([start + a - base, start + b - base] if strand == 1
                   else [end - (b - base), end - (a - base)])
            size = b - a
            if orf_strand == "-":
                orf = [max(orf_nt[0], orf_nt[1] - used - size), orf_nt[1] - used]
            else:
                orf = [orf_nt[0] + used, min(orf_nt[1], orf_nt[0] + used + size)]
            if orf[0] < orf[1]:
                segs.append(ref + orf)
            used += size
        base += n
    return segs


def orf_aa_to_contig(start: int, stop: int, strand: str, aa_from: int, aa_to: int) -> list[int]:
    """1-based aa range on an ORF ([start, stop) on the contig) → [start, end)."""
    if strand == "-":
        lo, hi = stop - 3 * aa_to, stop - 3 * (aa_from - 1)
    else:
        lo, hi = start + 3 * (aa_from - 1), start + 3 * aa_to
    return [max(start, lo), min(stop, hi)]


class ParsedReference:
    """A RefGenome plus what blastp needs and the JSON does not carry:
    the CDS translations and their location parts."""

    def __init__(self, genome: RefGenome, proteins: list[str],
                 parts: list[list[tuple[int, int, int]]], codon_starts: list[int]):
        self.genome       = genome
        self.proteins     = proteins
        self.parts        = parts
        self.codon_starts = codon_starts


def parse_genbank(path: str | Path) -> ParsedReference | None:
    """Read one GenBank file; None when it cannot be parsed or has no CDS."""
    from Bio import SeqIO

    try:
        record = SeqIO.read(str(path), "genbank")
    except Exception as exc:
        logger.warning(f"Synteny: cannot parse '{Path(path).name}': {exc}")
        return None

    genome = RefGenome(
        accession=record.id,
        title=record.description,
        organism=record.annotations.get("organism", ""),
        length=len(record.seq),
    )
    proteins: list[str] = []
    parts_l:  list[list[tuple[int, int, int]]] = []
    codons:   list[int] = []

    for f in record.features:
        if f.type != "CDS" or "pseudo" in f.qualifiers:
            continue
        aa = _q(f, "translation")
        if not aa:
            try:
                aa = str(f.extract(record.seq).translate(table=int(_q(f, "transl_table") or 1)))
            except Exception:
                aa = ""
        aa = aa.rstrip("*").replace("*", "X")
        if not aa:
            continue
        loc = f.location
        genome.cds.append(RefFeature(
            kind="CDS",
            product=_q(f, "product") or _q(f, "gene") or _q(f, "note") or f"CDS {len(genome.cds) + 1}",
            protein_id=_q(f, "protein_id") or None,
            start=int(loc.start), end=int(loc.end),
            strand="-" if loc.strand == -1 else "+",
            length_aa=len(aa),
        ))
        proteins.append(aa)
        parts_l.append(_parts(loc))
        codons.append(int(_q(f, "codon_start") or 1))

    seen: set[tuple[int, int, str]] = set()
    for f in record.features:
        if f.type != "mat_peptide":
            continue
        loc    = f.location
        strand = "-" if loc.strand == -1 else "+"
        start, end = int(loc.start), int(loc.end)
        # Records often list a peptide once per polyprotein (e.g. Gag and Gag-Pol).
        if (start, end, strand) in seen:
            continue
        seen.add((start, end, strand))
        # The longest CDS on the same strand that contains it is its polyprotein.
        parents = [i for i, c in enumerate(genome.cds)
                   if c.strand == strand and c.start <= start and end <= c.end]
        genome.peptides.append(RefFeature(
            kind="mat_peptide",
            product=_q(f, "product") or _q(f, "note") or "mature peptide",
            protein_id=_q(f, "protein_id") or None,
            start=start, end=end, strand=strand,
            length_aa=(end - start) // 3,
            parent=max(parents, key=lambda i: genome.cds[i].end - genome.cds[i].start) if parents else None,
        ))

    if not genome.cds:
        logger.warning(f"Synteny: '{record.id}' has no translatable CDS — skipped.")
        return None
    return ParsedReference(genome, proteins, parts_l, codons)


# ── DIAMOND blastp ──────────────────────────────────────────────────────────

class SyntenyAligner:
    """DIAMOND blastp of input ORFs against the reference proteins."""

    def __init__(self, workdir: str | Path, diamond_bin: str = "diamond",
                 threads: int = 2, e_value: float = 1e-5, run=None):
        self.workdir     = Path(workdir)
        self.diamond_bin = diamond_bin
        self.threads     = threads
        self.e_value     = e_value
        self._run        = run or subprocess.run

    @staticmethod
    def _write_faa(records: dict[str, str], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for rid, aa in records.items():
                fh.write(f">{rid}\n{aa}\n")

    def align(self, queries: dict[str, str], subjects: dict[str, str]) -> list[dict]:
        """Rows of the blastp table (dicts keyed by _BLASTP_FIELDS); [] on failure."""
        if not queries or not subjects:
            return []
        self.workdir.mkdir(parents=True, exist_ok=True)
        q, s = self.workdir / "orfs.faa", self.workdir / "refs.faa"
        db, out = self.workdir / "refs", self.workdir / "blastp.tsv"
        self._write_faa(queries, q)
        self._write_faa(subjects, s)

        cmds = [
            [self.diamond_bin, "makedb", "--in", str(s), "--db", str(db), "--quiet"],
            [self.diamond_bin, "blastp", "--db", str(db), "--query", str(q), "--out", str(out),
             "--outfmt", "6", *_BLASTP_FIELDS, "--evalue", str(self.e_value),
             "--more-sensitive", "--max-target-seqs", "0",
             "--threads", str(self.threads), "--quiet"],
        ]
        for cmd in cmds:
            try:
                res = self._run(cmd, capture_output=True, text=True)
            except Exception as exc:
                logger.error(f"Synteny: DIAMOND {cmd[1]} failed: {exc}")
                return []
            if res.returncode != 0:
                logger.error(f"Synteny: DIAMOND {cmd[1]} failed: {(res.stderr or '').strip()}")
                return []
        return self.parse(out)

    @staticmethod
    def parse(path: Path) -> list[dict]:
        rows: list[dict] = []
        if not path.exists():
            return rows
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) < len(_BLASTP_FIELDS):
                    continue
                try:
                    r = dict(zip(_BLASTP_FIELDS, row))
                    for k in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"):
                        r[k] = int(r[k])
                    for k in ("pident", "evalue", "bitscore"):
                        r[k] = float(r[k])
                    rows.append(r)
                except ValueError:
                    continue
        return rows


# ── Orientation / placement ─────────────────────────────────────────────────

def _weighted_median(values: list[float], weights: list[float]) -> float:
    pairs = sorted(zip(values, weights))
    half, acc = sum(weights) / 2, 0.0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return v
    return pairs[-1][0]


def placement(hits: list[SyntenyHit], orf_strand: dict[str, str],
              cds_strand: list[str]) -> tuple[bool, float]:
    """
    (flipped, offset) placing the contig under the reference.

    The contig is *flipped* (reverse-complemented relative to the reference)
    when most of the bit score comes from ORF/CDS pairs on opposite strands.
    Then reference ≈ offset − contig position; otherwise reference ≈ contig
    position + offset. The offset is the bit-score-weighted median over the
    hits that agree with that orientation.
    """
    if not hits:
        return False, 0.0
    same = lambda h: (orf_strand[h.orf_name] == "+") == (cds_strand[h.cds_index] == "+")
    w_same = sum(h.bit_score for h in hits if same(h))
    flipped = sum(h.bit_score for h in hits) - w_same > w_same
    agree = [h for h in hits if same(h) != flipped]
    mid = lambda iv: (iv[0] + iv[1]) / 2
    vals = [(mid(h.ref_nt) + mid(h.orf_nt)) if flipped else (mid(h.ref_nt) - mid(h.orf_nt))
            for h in agree]
    return flipped, _weighted_median(vals, [h.bit_score for h in agree])


# ── Orchestration ───────────────────────────────────────────────────────────

class SyntenyAnalyzer:
    """
    Runs the whole synteny step on the given sequences: sets ``seq.synteny``
    and returns the reference genomes used, keyed by accession.
    """

    def __init__(self, outdir: str | Path, email: str | None = None,
                 diamond_bin: str = "diamond", threads: int = 2, e_value: float = 1e-5,
                 fetcher: GenBankFetcher | None = None, aligner: SyntenyAligner | None = None):
        self.outdir  = Path(outdir)
        self.fetcher = fetcher or GenBankFetcher(self.outdir / "refs", email=email)
        self.aligner = aligner or SyntenyAligner(self.outdir, diamond_bin, threads, e_value)

    def run(self, seqs: list[NucSequence]) -> dict[str, RefGenome]:
        picked = {s.id: h for s in seqs if (h := best_viral_hit(s))}
        if not picked:
            logger.info("Synteny: no sequence has a viral BLASTn hit — skipped.")
            return {}

        paths = self.fetcher.fetch([h.accession for h in picked.values()])
        refs: dict[str, ParsedReference] = {}
        for key, path in paths.items():
            parsed = parse_genbank(path)
            if parsed:
                refs[key] = parsed

        # Queries: every ORF of every sequence with a parsed reference.
        by_id   = {s.id: s for s in seqs}
        queries: dict[str, str] = {}
        qmap:    dict[str, tuple[NucSequence, object]] = {}
        for sid, hit in picked.items():
            if acc_key(hit.accession) not in refs:
                continue
            for orf in by_id[sid].orfs:
                aa = (orf.aa_sequence or "").rstrip("*").replace("*", "X")
                if aa:
                    qid = f"q{len(queries)}"
                    queries[qid], qmap[qid] = aa, (by_id[sid], orf)
        subjects = {
            f"{key}|{i}": aa
            for key, ref in refs.items() for i, aa in enumerate(ref.proteins)
        }

        rows = self.aligner.align(queries, subjects)

        # Best HSP per (ORF, CDS), restricted to the sequence's own reference.
        best: dict[tuple[str, str], dict] = {}
        for r in rows:
            seq, _orf = qmap.get(r["qseqid"], (None, None))
            if seq is None:
                continue
            key, _, _ = r["sseqid"].rpartition("|")
            if key != acc_key(picked[seq.id].accession):
                continue
            k = (r["qseqid"], r["sseqid"])
            if k not in best or r["bitscore"] > best[k]["bitscore"]:
                best[k] = r

        hits_by_seq: dict[str, list[SyntenyHit]] = defaultdict(list)
        for (qid, sid), r in best.items():
            seq, orf = qmap[qid]
            key, _, idx = sid.rpartition("|")
            ref, ci = refs[key], int(idx)
            orf_nt = orf_aa_to_contig(orf.start_position, orf.stop_position, orf.strand,
                                      r["qstart"], r["qend"])
            hits_by_seq[seq.id].append(SyntenyHit(
                orf_name=orf.name, cds_index=ci,
                pident=r["pident"], evalue=r["evalue"], bit_score=r["bitscore"],
                q_start=r["qstart"], q_end=r["qend"], s_start=r["sstart"], s_end=r["send"],
                q_cov=round(100 * (r["qend"] - r["qstart"] + 1) / max(1, r["qlen"]), 1),
                s_cov=round(100 * (r["send"] - r["sstart"] + 1) / max(1, r["slen"]), 1),
                orf_nt=orf_nt,
                ref_nt=cds_aa_to_genome(ref.parts[ci], ref.codon_starts[ci],
                                        r["sstart"], r["send"]),
                segments=aligned_segments(ref.parts[ci], ref.codon_starts[ci],
                                          r["sstart"], r["send"], orf_nt, orf.strand),
            ))

        used: dict[str, RefGenome] = {}
        for sid, hit in picked.items():
            parsed = refs.get(acc_key(hit.accession))
            if parsed is None:
                continue
            seq  = by_id[sid]
            hits = sorted(hits_by_seq.get(sid, []), key=lambda h: h.orf_nt[0])
            flipped, offset = placement(
                hits, {o.name: o.strand for o in seq.orfs},
                [c.strand for c in parsed.genome.cds],
            )
            seq.synteny = SyntenyResult(reference=parsed.genome.accession,
                                        flipped=flipped, offset=round(offset, 1), hits=hits)
            used[parsed.genome.accession] = parsed.genome

        self._write_table(seqs, used)
        n_seq = sum(1 for s in seqs if s.synteny and s.synteny.hits)
        logger.success(
            f"Synteny: {len(used)} reference(s); {n_seq} sequence(s) with ORFs matching "
            f"their reference ({sum(len(v) for v in hits_by_seq.values())} ORF↔CDS pair(s))."
        )
        return used

    def _write_table(self, seqs: list[NucSequence], refs: dict[str, RefGenome]) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)
        with open(self.outdir / "synteny.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sequence", "orf", "reference", "cds_product", "protein_id",
                        "pident", "evalue", "bitscore", "orf_cov", "cds_cov",
                        "contig_start", "contig_end", "ref_start", "ref_end"])
            for s in seqs:
                if not s.synteny:
                    continue
                ref = refs[s.synteny.reference]
                for h in s.synteny.hits:
                    c = ref.cds[h.cds_index]
                    w.writerow([s.id, h.orf_name, ref.accession, c.product, c.protein_id or "",
                                h.pident, h.evalue, h.bit_score, h.q_cov, h.s_cov,
                                h.orf_nt[0] + 1, h.orf_nt[1], h.ref_nt[0] + 1, h.ref_nt[1]])
