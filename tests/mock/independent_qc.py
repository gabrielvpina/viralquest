#!/usr/bin/env python3
"""
independent_qc.py — independent coverage & base-quality oracle for tests/mock/.

Recomputes per-base **coverage** (total + sense/antisense) and **mean base
quality** for each mock contig directly from the FASTQ reads, using *only* the
Python standard library — no minimap2, no samtools, and crucially **no import of
viralquest**. It is therefore an independent source of truth to cross-check the
pipeline's `coverage.py` output (`CoveragePipeline`).

How it can be exact and aligner-free
------------------------------------
The mock reads (see generate.py) are error-free substrings of the contigs:
  * R1 = fragment[:150]           → forward strand, exact substring of the contig
  * R2 = revcomp(fragment[-150:]) → reverse strand; its revcomp is an exact
                                    substring of the contig
Every 150 bp read spans at least ~70 bp of unique (random) flank, so each read
occurs at exactly one position — a plain `str.find` recovers the true mapping.
R1 reads are counted as **sense** (forward), R2 as **antisense** (reverse),
mirroring the pipeline's strand split derived from the mpileup read-base case.

Outputs (written next to this script, under tests/mock/independent_qc/)
  * <contig>.tsv  — per-base: position, depth, sense_depth, antisense_depth, mean_quality
  * summary.tsv   — one row per contig: length, mean/max depth, breadth_1x, mean_quality
  * report.html   — self-contained coverage plots in the viralquest visual style
                    (two-sided sense/antisense track, bars coloured by base quality)
and a summary table printed to stdout.

Run from anywhere:  python tests/mock/independent_qc.py
"""

from __future__ import annotations

from pathlib import Path

HERE     = Path(__file__).parent
OUT_DIR  = HERE / "independent_qc"

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ── parsing ──────────────────────────────────────────────────────────────────

def read_fasta(path: Path) -> dict[str, str]:
    """Return {contig_id: sequence}. ID is the first whitespace-delimited token."""
    contigs: dict[str, str] = {}
    cid: str | None = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cid is not None:
                    contigs[cid] = "".join(chunks)
                cid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if cid is not None:
        contigs[cid] = "".join(chunks)
    return contigs


def read_fastq(path: Path):
    """Yield (name, sequence, quality) for each FASTQ record."""
    with open(path, encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq  = fh.readline().rstrip("\n")
            fh.readline()                      # '+' separator
            qual = fh.readline().rstrip("\n")
            name = header[1:].rstrip("\n").split()[0]
            yield name, seq, qual


def contig_of(read_name: str) -> str:
    """`contig_random_frag0007/1` → `contig_random` (strip mate + _fragNNNN)."""
    name = read_name.rsplit("/", 1)[0]         # drop /1 or /2 mate suffix
    return name.rsplit("_frag", 1)[0]


# ── core: map one read onto its contig ───────────────────────────────────────

def place_read(seq: str, qual: str, mate: int, contig: str) -> tuple[int, str, int] | None:
    """
    Return (start, forward_quality, strand) for an exact placement, or None.

    strand: +1 for sense (R1, forward), -1 for antisense (R2, reverse).
    forward_quality is the quality string oriented 5'->3' on the contig's
    forward strand (reversed for R2 so it aligns with the placed bases).
    """
    if mate == 1:
        pos = contig.find(seq)
        if pos < 0:
            return None
        return pos, qual, +1
    # mate 2: the read is the reverse complement of a forward substring
    fwd_seq = revcomp(seq)
    pos = contig.find(fwd_seq)
    if pos < 0:
        return None
    return pos, qual[::-1], -1                 # reverse qual to match forward bases


# ── build profiles ───────────────────────────────────────────────────────────

def build_profiles(contigs: dict[str, str], read_files):
    """
    Returns {cid: dict} with per-base arrays (plain lists):
        depth, sense, anti, qsum, qn   (qsum/qn → mean quality per base).
    """
    prof = {
        cid: {
            "seq":   seq,
            "depth": [0] * len(seq),
            "sense": [0] * len(seq),
            "anti":  [0] * len(seq),
            "qsum":  [0] * len(seq),   # summed Phred per position
            "qn":    [0] * len(seq),   # #bases stacked per position
        }
        for cid, seq in contigs.items()
    }
    unmapped = 0
    for path, mate in read_files:
        for name, seq, qual in read_fastq(path):
            cid = contig_of(name)
            p = prof.get(cid)
            if p is None:
                unmapped += 1
                continue
            placed = place_read(seq, qual, mate, p["seq"])
            if placed is None:
                unmapped += 1
                continue
            start, fqual, strand = placed
            for i in range(len(seq)):
                pos = start + i
                p["depth"][pos] += 1
                if strand > 0:
                    p["sense"][pos] += 1
                else:
                    p["anti"][pos] += 1
                p["qsum"][pos] += ord(fqual[i]) - 33
                p["qn"][pos]   += 1
    return prof, unmapped


# ── export ───────────────────────────────────────────────────────────────────

def export(prof) -> list[tuple]:
    """Write per-contig TSVs; return summary rows."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: list[tuple] = []
    for cid, p in prof.items():
        depth, sense, anti, qsum, qn = p["depth"], p["sense"], p["anti"], p["qsum"], p["qn"]
        n = len(depth)

        with open(OUT_DIR / f"{cid}.tsv", "w", encoding="utf-8") as fh:
            fh.write("position_nt\tdepth\tsense_depth\tantisense_depth\tmean_quality\n")
            for i in range(n):
                mq = round(qsum[i] / qn[i], 2) if qn[i] else 0.0
                fh.write(f"{i + 1}\t{depth[i]}\t{sense[i]}\t{anti[i]}\t{mq}\n")

        covered      = sum(1 for d in depth if d >= 1)
        mean_depth   = round(sum(depth) / n, 2) if n else 0.0
        max_depth    = max(depth) if depth else 0
        breadth_1x   = round(covered / n, 4) if n else 0.0
        cov_q        = [qsum[i] / qn[i] for i in range(n) if qn[i]]
        mean_quality = round(sum(cov_q) / len(cov_q), 2) if cov_q else 0.0
        summary.append((cid, n, mean_depth, max_depth, breadth_1x, mean_quality))

    with open(OUT_DIR / "summary.tsv", "w", encoding="utf-8") as fh:
        fh.write("contig\tlength\tmean_depth\tmax_depth\tbreadth_1x\tmean_quality\n")
        for row in summary:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return summary


# ── HTML plots (viralquest visual style) ─────────────────────────────────────

# Palette lifted from viralquest/components/base.css so the plots read as the
# same system. Quality colouring mirrors _qualColor() in section_viewer.js.
_PAL = {
    "success": "#1f9d6b",  # Q>=30
    "warning": "#e0921c",  # Q20-29
    "danger":  "#d24a3d",  # Q<20
    "accent":  "#2f86d6",
    "text3":   "#94a3b8",  # no reads
    "border":  "#c5d0dc",
}


def _qual_color(q: float) -> str:
    if q >= 30: return _PAL["success"]
    if q >= 20: return _PAL["warning"]
    if q > 0:   return _PAL["danger"]
    return _PAL["text3"]


def _coverage_svg(p: dict) -> str:
    """One two-sided coverage track (sense up / antisense down), coloured by Phred."""
    seq   = p["seq"]
    n     = len(seq)
    sense, anti, qsum, qn = p["sense"], p["anti"], p["qsum"], p["qn"]

    PAD_L, PAD_R = 40, 14
    DRAW_W, HALF = 820, 46          # draw width; half-height of the two-sided band
    top, mid, bottom = 8, 8 + HALF, 8 + 2 * HALF
    axis_y = bottom + 16
    total_h = axis_y + 6
    bw    = DRAW_W / n
    maxS  = max(1, max(sense), max(anti))
    xof   = lambda i: PAD_L + i * bw

    out = [f'<svg viewBox="0 0 {PAD_L + DRAW_W + PAD_R} {total_h}" '
           f'width="100%" preserveAspectRatio="xMidYMid meet" '
           f'font-family="Inter, system-ui, sans-serif">']

    # bars: sense grows up (opacity .85), antisense grows down (opacity .5)
    for i in range(n):
        mq = qsum[i] / qn[i] if qn[i] else 0.0
        col = _qual_color(mq)
        x = round(xof(i), 2)
        hu = sense[i] / maxS * HALF
        hd = anti[i]  / maxS * HALF
        if hu > 0:
            out.append(f'<rect x="{x:.2f}" y="{mid - hu:.2f}" width="{bw:.2f}" '
                       f'height="{hu:.2f}" fill="{col}" opacity="0.85"/>')
        if hd > 0:
            out.append(f'<rect x="{x:.2f}" y="{mid:.2f}" width="{bw:.2f}" '
                       f'height="{hd:.2f}" fill="{col}" opacity="0.5"/>')

    # centre line + gutter markers
    out.append(f'<line x1="{PAD_L}" x2="{PAD_L + DRAW_W}" y1="{mid}" y2="{mid}" '
               f'stroke="{_PAL["border"]}" stroke-width="1"/>')
    out.append(f'<text x="{PAD_L - 8}" y="{top + 7}" text-anchor="end" font-size="8" '
               f'font-weight="700" fill="#475569">＋</text>')
    out.append(f'<text x="{PAD_L - 8}" y="{mid + 3}" text-anchor="end" font-size="9" '
               f'font-weight="700" fill="#475569" font-family="monospace">Cov</text>')
    out.append(f'<text x="{PAD_L - 8}" y="{bottom - 1}" text-anchor="end" font-size="8" '
               f'font-weight="700" fill="#475569">－</text>')
    out.append(f'<text x="{PAD_L + DRAW_W}" y="{mid + 3}" text-anchor="end" font-size="8" '
               f'fill="#94a3b8" font-family="monospace">{maxS}× max</text>')

    # genome axis with ticks every 100 nt
    out.append(f'<line x1="{PAD_L}" x2="{PAD_L + DRAW_W}" y1="{axis_y}" y2="{axis_y}" '
               f'stroke="{_PAL["border"]}" stroke-width="1"/>')
    step = 100
    for pos in range(0, n + 1, step):
        x = PAD_L + (pos / n) * DRAW_W
        out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{axis_y}" y2="{axis_y + 3}" '
                   f'stroke="{_PAL["border"]}"/>')
        out.append(f'<text x="{x:.1f}" y="{axis_y + 12}" text-anchor="middle" '
                   f'font-size="8" fill="#94a3b8">{pos}</text>')
    out.append("</svg>")
    return "".join(out)


def render_html(prof: dict) -> Path:
    """Write independent_qc/report.html with one coverage panel per contig."""
    panels = []
    for cid, p in prof.items():
        n = len(p["seq"])
        mean_s = sum(p["sense"]) / n if n else 0.0
        mean_a = sum(p["anti"])  / n if n else 0.0
        breadth = sum(1 for d in p["depth"] if d) / n if n else 0.0
        cov_q = [p["qsum"][i] / p["qn"][i] for i in range(n) if p["qn"][i]]
        mean_q = sum(cov_q) / len(cov_q) if cov_q else 0.0
        cap = (f"sense {mean_s:.1f}× · anti {mean_a:.1f}× · "
               f"breadth {breadth * 100:.0f}% · Q̄ {mean_q:.0f}")
        panels.append(f"""
    <section class="panel">
      <div class="phead">
        <span class="cid">{cid}</span>
        <span class="len">{n} bp</span>
        <span class="cap">{cap}</span>
      </div>
      {_coverage_svg(p)}
    </section>""")

    legend = (f'<span class="lead">Read coverage</span>'
              f'<span class="sw" style="background:{_PAL["success"]}"></span>Q≥30'
              f'<span class="sw" style="background:{_PAL["warning"]}"></span>Q20–29'
              f'<span class="sw" style="background:{_PAL["danger"]}"></span>Q&lt;20')

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Independent QC — tests/mock</title>
<style>
  :root {{ --bg:#f3f6fa; --surface:#fff; --border:#e2e8f0; --text:#1e293b;
           --text2:#475569; --text3:#94a3b8; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg:#0f1620; --surface:#18212e; --border:#28323f; --text:#e2e8f0;
             --text2:#a9b6c6; --text3:#6b7a8d; }}
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; padding:28px; background:var(--bg); color:var(--text);
          font-family:Inter,'Segoe UI',system-ui,sans-serif; }}
  h1 {{ font-size:1.25rem; margin:0 0 4px; }}
  .sub {{ color:var(--text2); font-size:.85rem; margin:0 0 20px; max-width:70ch; }}
  .sub code {{ font-family:'JetBrains Mono',monospace; font-size:.8rem; }}
  .legend {{ display:flex; align-items:center; gap:6px; font-size:.72rem;
             color:var(--text3); margin:0 0 18px; }}
  .legend .lead {{ font-weight:600; color:var(--text2); margin-right:4px; }}
  .legend .sw {{ width:11px; height:11px; border-radius:2px; margin-left:8px; }}
  .panel {{ background:var(--surface); border:1px solid var(--border);
            border-radius:10px; padding:14px 16px; margin-bottom:16px; }}
  .phead {{ display:flex; align-items:baseline; gap:12px; margin-bottom:4px; }}
  .cid {{ font-family:'JetBrains Mono',monospace; font-weight:600; font-size:.92rem; }}
  .len {{ color:var(--text3); font-size:.75rem; }}
  .cap {{ margin-left:auto; color:var(--text3); font-size:.75rem; }}
</style></head>
<body>
  <h1>Independent coverage &amp; quality — <code>tests/mock</code></h1>
  <p class="sub">Computed directly from <code>reads_R1/R2.fastq</code> by exact
  substring placement (pure Python, no aligner, no viralquest import). An
  independent oracle to cross-check <code>CoveragePipeline</code>. Sense = R1
  (forward), antisense = R2 (reverse); bars coloured by mean base quality.</p>
  <div class="legend">{legend}</div>
  {"".join(panels)}
</body></html>"""

    path = OUT_DIR / "report.html"
    path.write_text(html, encoding="utf-8")
    return path


def main() -> None:
    contigs = read_fasta(HERE / "contigs.fasta")
    read_files = [(HERE / "reads_R1.fastq", 1), (HERE / "reads_R2.fastq", 2)]
    prof, unmapped = build_profiles(contigs, read_files)
    summary = export(prof)
    html_path = render_html(prof)

    print(f"independent QC oracle → {OUT_DIR.relative_to(HERE.parent.parent)}")
    print(f"  plots → {html_path.relative_to(HERE.parent.parent)}")
    if unmapped:
        print(f"  WARNING: {unmapped} read(s) could not be placed by exact match")
    print(f"{'contig':<24}{'len':>6}{'mean_dp':>9}{'max_dp':>8}{'breadth':>9}{'mean_q':>8}")
    for cid, n, md, mx, br, mq in summary:
        print(f"{cid:<24}{n:>6}{md:>9}{mx:>8}{br:>9}{mq:>8}")


if __name__ == "__main__":
    main()
