"""
seq_table.py — The aggregated results workbook of a ViralQuest run.

Written next to the JSON/HTML at the end of a run as ``<stem>_viralquest.xlsx``,
one sheet per kind of result:

    Sequences         one row per sequence: taxonomy, domains, best BLASTn /
                      BLASTx hit (accessions link to their NCBI record)
    Classification    heuristic / LLM scores, novelty tier, cluster
    Taxonomy          full lineage of the BLASTx hit + the BLASTn hit's taxonomy
    ORFs              every ORF of every sequence
    Domains           every HMM domain (RVDB, Vfam, EggNOG, Pfam)
    BLASTx hits       every RefSeq / NR hit, ranked by bit score
    BLASTn hits       every BLASTn hit, ranked by bit score
    Clusters          cluster membership with identity / coverage
    Synteny           ORF ↔ reference CDS blastp hits
    Contig quality    low complexity, internal repeats, k-mer repetitiveness
    Genome coverage   read depth, breadth and uniformity per sequence
    Quantification    every Salmon target (viral, reference HK, conserved HK, host)
    Host-viral hits   host transcripts similar to viral sequences (EVEs)

The consolidated report (viralquest-report) writes the same sheets over all
samples, with a Sample column, plus "Samples" and "General clusters".

Every row builder works on the exported report dict (the JSON schema). The
Sequences sheet mirrors the HTML viewer's "Table" export
(components*/export.js, vqSeqTableRows) — keep both in sync.

The XLSX is written with the standard library only (zipfile + SpreadsheetML),
so no spreadsheet dependency is needed.
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import escape

from loguru import logger

COLUMNS = [
    "Sample", "Sequence ID", "Order", "Family", "Genus", "Species", "Genome type",
    "Sequence (nt)", "Length", "Domains",
    "BLASTn hit", "BLASTn identity (%)", "BLASTn coverage (%)",
    "BLASTn hit ID", "BLASTn hit sequence", "BLASTn e-value",
    "BLASTx hit", "BLASTx identity (%)", "BLASTx coverage (%)",
    "BLASTx hit ID", "BLASTx hit sequence", "BLASTx e-value",
]

# Accession columns → the NCBI database their record lives in.
LINK_COLUMNS = {"BLASTn hit ID": "nuccore", "BLASTx hit ID": "protein"}

# Excel refuses cells longer than this; such values stay complete in the JSON.
XLSX_CELL_LIMIT = 32767

_ACC_RE = re.compile(r"^[A-Za-z0-9_]+(\.[0-9]+)?$")


def ncbi_url(accession: str | None, db: str) -> str | None:
    """NCBI record URL for a bare accession; None for pipe-packed / local ids."""
    if not accession or not _ACC_RE.match(accession):
        return None
    return f"https://www.ncbi.nlm.nih.gov/{db}/{accession}"


@dataclass
class Sheet:
    name:    str
    columns: list[str]
    rows:    list[dict]
    links:   dict[str, str] = field(default_factory=dict)   # column → NCBI db
    widths:  dict[str, int] = field(default_factory=dict)


# ── Sequences sheet ─────────────────────────────────────────────────────────

def _best(hits: list[dict]) -> dict | None:
    return max(hits, key=lambda h: h.get("bit_score") or 0) if hits else None


def _ranked(hits: list[dict]) -> list[tuple[int, dict]]:
    return list(enumerate(sorted(hits or [], key=lambda h: -(h.get("bit_score") or 0)), 1))


def _domains(seq: dict) -> str:
    seen: list[str] = []
    for orf in seq.get("orfs") or []:
        for d in orf.get("domains") or []:
            label = f"{d.get('target', '')} ({d.get('database', '')})"
            if label not in seen:
                seen.append(label)
    return "; ".join(seen)


def sample_name(seq: dict, report: dict) -> str:
    """Multi-sample stamp, then the per-sequence field, then the input file stem."""
    name = seq.get("sample") or seq.get("sample_name")
    if name:
        return name
    input_name = ((report.get("meta") or {}).get("input_file") or {}).get("name") or ""
    return Path(input_name).stem if input_name else ""


def seq_row(seq: dict, report: dict) -> dict:
    """One Sequences-sheet row (column name → value) for an exported sequence."""
    bt  = seq.get("blastn_taxonomy") or {}
    tax = seq.get("taxonomy") or bt.get("taxonomy") or {}
    bn  = _best(seq.get("blastn_hits") or [])
    bx  = _best(seq.get("blastx_nr_hits") or seq.get("blastx_hits") or [])
    return {
        "Sample":              sample_name(seq, report),
        "Sequence ID":         seq.get("id") or "",
        "Order":               tax.get("order") or "",
        "Family":              tax.get("family") or "",
        "Genus":               tax.get("genus") or "",
        "Species":             tax.get("species") or tax.get("scientific_name") or "",
        "Genome type":         tax.get("genome") or "",
        "Sequence (nt)":       seq.get("sequence") or "",
        "Length":              seq.get("length"),
        "Domains":             _domains(seq),
        "BLASTn hit":          bn.get("stitle") or "" if bn else "",
        "BLASTn identity (%)": bn.get("pident") if bn else None,
        "BLASTn coverage (%)": bn.get("qcovhsp") if bn else None,
        "BLASTn hit ID":       bn.get("accession") or "" if bn else "",
        "BLASTn hit sequence": bn.get("subject_seq") or "" if bn else "",
        "BLASTn e-value":      bn.get("evalue") if bn else None,
        "BLASTx hit":          bx.get("subject_title") or "" if bx else "",
        "BLASTx identity (%)": bx.get("pct_identity") if bx else None,
        "BLASTx coverage (%)": bx.get("query_coverage") if bx else None,
        "BLASTx hit ID":       bx.get("subject_id") or "" if bx else "",
        "BLASTx hit sequence": bx.get("subject_seq") or "" if bx else "",
        "BLASTx e-value":      bx.get("e_value") if bx else None,
    }


def table_rows(report: dict) -> list[dict]:
    return [seq_row(s, report) for s in report.get("sequences") or []]


# ── Other sheets ────────────────────────────────────────────────────────────

def _seqs(report: dict) -> list[dict]:
    return report.get("sequences") or []


def _cluster_of(report: dict) -> dict[str, str]:
    return {m.get("seq_id"): c.get("cluster_id")
            for c in report.get("clusters") or [] for m in c.get("members") or []}


def _classification(report: dict) -> Sheet:
    cl = _cluster_of(report)
    rows = []
    for s in _seqs(report):
        h, llm = s.get("heuristic_output") or {}, s.get("llm_output") or {}
        rows.append({
            "Sequence ID": s.get("id"), "Viral": "yes" if s.get("is_viral") else "no",
            "Heuristic score": h.get("vq_score"), "Classification": h.get("classification") or "",
            "Novelty": h.get("novelty") or "", "Novelty flags": "; ".join(h.get("novelty_flags") or []),
            "LLM score": llm.get("vq_score"), "LLM classification": llm.get("classification") or "",
            "LLM model": llm.get("model") or "", "Cluster": cl.get(s.get("id"), ""),
            "Heuristic analysis": h.get("analysis") or "", "LLM analysis": llm.get("analysis") or "",
        })
    return Sheet("Classification", list(rows[0]) if rows else [
        "Sequence ID", "Viral", "Heuristic score", "Classification", "Novelty", "Novelty flags",
        "LLM score", "LLM classification", "LLM model", "Cluster", "Heuristic analysis",
        "LLM analysis"], rows, widths={"Heuristic analysis": 60, "LLM analysis": 60})


_RANKS = [("Realm", "clade"), ("Kingdom", "kingdom"), ("Phylum", "phylum"), ("Class", "class"),
          ("Order", "order"), ("Family", "family"), ("Subfamily", "subfamily"),
          ("Genus", "genus"), ("Species", "species")]


def _taxonomy(report: dict) -> Sheet:
    cols = (["Sequence ID", "Tax ID", "Name"] + [r for r, _ in _RANKS] + ["Genome type",
            "BLASTn taxonomy status", "BLASTn matched name", "BLASTn tax ID",
            "BLASTn family", "BLASTn genus", "BLASTn species"])
    rows = []
    for s in _seqs(report):
        t  = s.get("taxonomy") or {}
        bt = s.get("blastn_taxonomy") or {}
        btt = bt.get("taxonomy") or {}
        rows.append({
            "Sequence ID": s.get("id"), "Tax ID": t.get("tax_id"), "Name": t.get("scientific_name") or "",
            **{r: t.get(k) or "" for r, k in _RANKS}, "Genome type": t.get("genome") or "",
            "BLASTn taxonomy status": bt.get("status") or "", "BLASTn matched name": bt.get("matched_name") or "",
            "BLASTn tax ID": bt.get("taxid"), "BLASTn family": btt.get("family") or "",
            "BLASTn genus": btt.get("genus") or "",
            "BLASTn species": btt.get("species") or btt.get("scientific_name") or "",
        })
    return Sheet("Taxonomy", cols, rows)


def _orfs(report: dict) -> Sheet:
    cols = ["Sequence ID", "ORF", "Start", "Stop", "Strand", "Frame", "Type",
            "Length (aa)", "Length (nt)", "Domains", "Protein sequence"]
    rows = [{
        "Sequence ID": s.get("id"), "ORF": o.get("name"), "Start": o.get("start_position"),
        "Stop": o.get("stop_position"), "Strand": o.get("strand"), "Frame": o.get("frame"),
        "Type": o.get("orf_type") or "", "Length (aa)": o.get("length_aa"),
        "Length (nt)": o.get("length_nt"), "Domains": len(o.get("domains") or []),
        "Protein sequence": o.get("aa_sequence") or "",
    } for s in _seqs(report) for o in s.get("orfs") or []]
    return Sheet("ORFs", cols, rows, widths={"ORF": 30, "Protein sequence": 30})


def _domain_sheet(report: dict) -> Sheet:
    cols = ["Sequence ID", "ORF", "Database", "Target", "Description", "Score",
            "E-value", "Start", "Stop"]
    rows = [{
        "Sequence ID": s.get("id"), "ORF": o.get("name"), "Database": d.get("database"),
        "Target": d.get("target"), "Description": d.get("description") if isinstance(d.get("description"), str) else "",
        "Score": d.get("score"), "E-value": d.get("e_value"), "Start": d.get("start"), "Stop": d.get("stop"),
    } for s in _seqs(report) for o in s.get("orfs") or [] for d in o.get("domains") or []]
    return Sheet("Domains", cols, rows, widths={"ORF": 30, "Target": 24, "Description": 50})


def _blastx(report: dict) -> Sheet:
    cols = ["Sequence ID", "Database", "Rank", "Subject ID", "Subject title", "Species",
            "Identity (%)", "Coverage (%)", "Alignment length", "E-value", "Bit score",
            "Query start", "Query end", "Subject start", "Subject end"]
    rows = []
    for s in _seqs(report):
        for db, key in (("RefSeq", "blastx_hits"), ("NR", "blastx_nr_hits")):
            for rank, x in _ranked(s.get(key)):
                rows.append({
                    "Sequence ID": s.get("id"), "Database": db, "Rank": rank,
                    "Subject ID": x.get("subject_id"), "Subject title": x.get("subject_title"),
                    "Species": x.get("species") or "", "Identity (%)": x.get("pct_identity"),
                    "Coverage (%)": x.get("query_coverage"), "Alignment length": x.get("aln_length"),
                    "E-value": x.get("e_value"), "Bit score": x.get("bit_score"),
                    "Query start": x.get("query_start"), "Query end": x.get("query_end"),
                    "Subject start": x.get("subject_start"), "Subject end": x.get("subject_end"),
                })
    return Sheet("BLASTx hits", cols, rows, links={"Subject ID": "protein"},
                 widths={"Subject ID": 16, "Subject title": 50, "Species": 30})


def _blastn(report: dict) -> Sheet:
    cols = ["Sequence ID", "Rank", "Accession", "Title", "Identity (%)", "Coverage (%)",
            "Query length", "Subject length", "E-value", "Bit score", "Query start", "Query end"]
    rows = [{
        "Sequence ID": s.get("id"), "Rank": rank, "Accession": x.get("accession"),
        "Title": x.get("stitle"), "Identity (%)": x.get("pident"), "Coverage (%)": x.get("qcovhsp"),
        "Query length": x.get("qlen"), "Subject length": x.get("slen"), "E-value": x.get("evalue"),
        "Bit score": x.get("bit_score"), "Query start": x.get("query_start"), "Query end": x.get("query_end"),
    } for s in _seqs(report) for rank, x in _ranked(s.get("blastn_hits"))]
    return Sheet("BLASTn hits", cols, rows, links={"Accession": "nuccore"},
                 widths={"Accession": 16, "Title": 60})


def _clusters(report: dict) -> Sheet:
    cols = ["Cluster", "Species", "Representative", "Size", "Sequence ID", "Is representative",
            "Identity (%)", "Coverage (%)", "Alignment start", "Alignment end"]
    rows = [{
        "Cluster": c.get("cluster_id"), "Species": c.get("species") or "",
        "Representative": c.get("representative_id"), "Size": c.get("size"),
        "Sequence ID": m.get("seq_id"), "Is representative": "yes" if m.get("is_representative") else "no",
        "Identity (%)": m.get("identity"), "Coverage (%)": m.get("query_coverage"),
        "Alignment start": m.get("aln_start"), "Alignment end": m.get("aln_end"),
    } for c in report.get("clusters") or [] for m in c.get("members") or []]
    return Sheet("Clusters", cols, rows, widths={"Species": 30})


def _synteny(report: dict) -> Sheet:
    cols = ["Sequence ID", "Reference", "Orientation", "ORF", "Reference CDS", "Protein ID",
            "Identity (%)", "E-value", "Bit score", "ORF aa start", "ORF aa end",
            "CDS aa start", "CDS aa end", "ORF coverage (%)", "CDS coverage (%)",
            "Reference start", "Reference end"]
    refs = (report.get("synteny") or {}).get("references") or {}
    rows = []
    for s in _seqs(report):
        syn = s.get("synteny") or {}
        cds = (refs.get(syn.get("reference")) or {}).get("cds") or []
        for h in syn.get("hits") or []:
            c = cds[h["cds_index"]] if 0 <= h.get("cds_index", -1) < len(cds) else {}
            ref_nt = h.get("ref_nt") or [None, None]
            rows.append({
                "Sequence ID": s.get("id"), "Reference": syn.get("reference"),
                "Orientation": "reverse complement" if syn.get("flipped") else "forward",
                "ORF": h.get("orf_name"), "Reference CDS": c.get("product") or "",
                "Protein ID": c.get("protein_id") or "", "Identity (%)": h.get("pident"),
                "E-value": h.get("evalue"), "Bit score": h.get("bit_score"),
                "ORF aa start": h.get("q_start"), "ORF aa end": h.get("q_end"),
                "CDS aa start": h.get("s_start"), "CDS aa end": h.get("s_end"),
                "ORF coverage (%)": h.get("q_cov"), "CDS coverage (%)": h.get("s_cov"),
                "Reference start": ref_nt[0] + 1 if ref_nt[0] is not None else None,
                "Reference end": ref_nt[1],
            })
    return Sheet("Synteny", cols, rows, links={"Reference": "nuccore", "Protein ID": "protein"},
                 widths={"ORF": 30, "Reference CDS": 40, "Reference": 16, "Protein ID": 16})


def _runs(ivs) -> str:
    return "; ".join(f"{a + 1}-{b + 1}" for a, b in ivs or [])


def _quality(report: dict) -> Sheet:
    cols = ["Sequence ID", "Length", "GC (%)", "N count", "Low-complexity fraction",
            "Low-complexity regions", "Direct repeats", "Inverted repeats", "k-mer size",
            "Distinct k-mers", "Total k-mers", "k-mer repeat score", "Low-complexity runs"]
    rows = []
    for s in _seqs(report):
        q = s.get("seq_quality") or {}
        reps = q.get("self_repeats") or []
        rows.append({
            "Sequence ID": s.get("id"), "Length": s.get("length"), "GC (%)": s.get("gc_content"),
            "N count": s.get("n_count"), "Low-complexity fraction": q.get("low_complexity_frac"),
            "Low-complexity regions": len(q.get("low_complexity_regions") or []) if q else None,
            "Direct repeats": sum(1 for r in reps if r.get("strand") == "plus") if q else None,
            "Inverted repeats": sum(1 for r in reps if r.get("strand") == "minus") if q else None,
            "k-mer size": q.get("kmer_size"), "Distinct k-mers": q.get("kmer_distinct"),
            "Total k-mers": q.get("kmer_total"), "k-mer repeat score": q.get("kmer_repeat_score"),
            "Low-complexity runs": _runs(q.get("low_complexity_regions")),
        })
    return Sheet("Contig quality", cols, rows, widths={"Low-complexity runs": 40})


def _coverage(report: dict) -> Sheet:
    cols = ["Sequence ID", "Length", "Mean depth", "Max depth", "CV", "Breadth ≥1x (%)",
            "Mean depth sense", "Mean depth antisense", "Mean base quality",
            "Low-coverage regions", "Low-coverage runs"]
    rows = []
    for s in _seqs(report):
        c = s.get("coverage")
        if not c:
            continue
        b = c.get("breadth_1x")
        rows.append({
            "Sequence ID": s.get("id"), "Length": c.get("length"), "Mean depth": c.get("mean_depth"),
            "Max depth": c.get("max_depth"), "CV": c.get("cv"),
            "Breadth ≥1x (%)": round(100 * b, 2) if b is not None else None,
            "Mean depth sense": c.get("mean_depth_fwd"), "Mean depth antisense": c.get("mean_depth_rev"),
            "Mean base quality": c.get("mean_quality"),
            "Low-coverage regions": len(c.get("low_cov_regions") or []),
            "Low-coverage runs": _runs(c.get("low_cov_regions")),
        })
    return Sheet("Genome coverage", cols, rows, widths={"Low-coverage runs": 40})


_PREFIXES = (("VQ_VIRAL_", "viral"), ("VQ_REFHK_", "reference HK"), ("VQ_CONS_", "conserved HK"))


def _quantification(report: dict) -> Sheet:
    cols = ["Target", "Type", "Kingdom", "Sequence ID", "Length", "Effective length",
            "TPM", "Reads"]
    sq = report.get("salmon_quant") or {}
    entries = (list(sq.get("viral_quant") or []) + list(sq.get("ref_hk_quant") or [])
               + list(sq.get("conserved_quant") or []) + list(sq.get("pfam_hk_quant") or [])
               + [r.get("transcript") or {} for r in sq.get("host_viral_hits") or []])
    rows = []
    for e in entries:
        name = e.get("name") or ""
        sid = name.removeprefix("VQ_VIRAL_") if name.startswith("VQ_VIRAL_") else ""
        rows.append({
            "Target": name, "Type": e.get("seq_type") or "", "Kingdom": (e.get("kingdom") or "").title(),
            "Sequence ID": sid, "Length": e.get("length"), "Effective length": e.get("eff_length"),
            "TPM": e.get("tpm"), "Reads": e.get("num_reads"),
        })
    return Sheet("Quantification", cols, rows, widths={"Target": 36})


def _host_viral(report: dict) -> Sheet:
    cols = ["Host transcript", "Sequence ID", "Identity (%)", "Coverage (%)", "E-value",
            "Bit score", "Transcript TPM", "Transcript reads"]
    rows = []
    for r in (report.get("salmon_quant") or {}).get("host_viral_hits") or []:
        t = r.get("transcript") or {}
        for h in r.get("blast_hits") or []:
            rows.append({
                "Host transcript": h.get("host_transcript_id"), "Sequence ID": h.get("viral_seq_id"),
                "Identity (%)": h.get("pident"), "Coverage (%)": h.get("qcovhsp"),
                "E-value": h.get("evalue"), "Bit score": h.get("bit_score"),
                "Transcript TPM": t.get("tpm"), "Transcript reads": t.get("num_reads"),
            })
    return Sheet("Host-viral hits", cols, rows, widths={"Host transcript": 28})


_WIDTHS = {"Sample": 18, "Sequence ID": 18, "Species": 30, "BLASTn hit": 40, "BLASTx hit": 40,
           "Domains": 30, "Sequence (nt)": 30, "BLASTn hit sequence": 30,
           "BLASTx hit sequence": 30, "BLASTn hit ID": 16, "BLASTx hit ID": 16}


def workbook_sheets(report: dict) -> list[Sheet]:
    """Every sheet of the results workbook, in order (empty sheets keep their header)."""
    return [
        Sheet("Sequences", COLUMNS, table_rows(report), LINK_COLUMNS, _WIDTHS),
        _classification(report), _taxonomy(report), _orfs(report), _domain_sheet(report),
        _blastx(report), _blastn(report), _clusters(report), _synteny(report),
        _quality(report), _coverage(report), _quantification(report), _host_viral(report),
    ]


# ── Multi-sample workbook ───────────────────────────────────────────────────

def _general_clusters(report_data: dict) -> Sheet:
    cols = ["Cluster", "Species", "Representative", "Samples", "Size", "Sample", "Sequence ID",
            "Is representative", "Identity (%)", "Coverage (%)", "TPM", "Reads"]
    rows = [{
        "Cluster": c.get("gid") or c.get("cluster_id"), "Species": c.get("species") or "",
        "Representative": c.get("representative"), "Samples": len(c.get("samples") or []),
        "Size": len(c.get("members") or []), "Sample": m.get("sample"), "Sequence ID": m.get("seq_id"),
        "Is representative": "yes" if m.get("is_representative") else "no",
        "Identity (%)": m.get("identity"), "Coverage (%)": m.get("coverage"),
        "TPM": m.get("tpm"), "Reads": m.get("reads"),
    } for c in report_data.get("clusters") or [] for m in c.get("members") or []]
    return Sheet("General clusters", cols, rows, widths={"Species": 30, "Representative": 26})


def _samples_sheet(report_data: dict) -> Sheet:
    cols = ["Sample", "Input sequences", "Viral flagged", "Confirmed viral", "Clusters",
            "Salmon pathway", "Mapping rate (%)", "Input reads", "Viral reads",
            "ViralQuest version", "Run started", "Run time (s)"]
    rows = []
    for s in report_data.get("samples") or []:
        sq, wf = s.get("salmon") or {}, s.get("workflow") or {}
        rows.append({
            "Sample": s.get("sample"), "Input sequences": s.get("n_input"),
            "Viral flagged": s.get("n_viral"), "Confirmed viral": s.get("n_confirmed"),
            "Clusters": s.get("n_clusters"), "Salmon pathway": sq.get("pathway") or "",
            "Mapping rate (%)": sq.get("mapping_rate"), "Input reads": sq.get("input_reads"),
            "Viral reads": sq.get("viral_reads"), "ViralQuest version": wf.get("version") or "",
            "Run started": wf.get("started_at") or "", "Run time (s)": wf.get("total_seconds"),
        })
    return Sheet("Samples", cols, rows, widths={"Sample": 20, "Run started": 26})


def merged_sample_sheets(samples) -> list[Sheet]:
    """
    Every per-sample sheet, concatenated over *samples* (SampleReport objects),
    with a leading Sample column (the Sequences sheet already has one). The
    per-sample cluster sheet is renamed "Sample clusters".
    """
    merged: dict[str, Sheet] = {}
    for smp in samples:
        for sh in workbook_sheets(smp.report):
            name = "Sample clusters" if sh.name == "Clusters" else sh.name
            has_sample = "Sample" in sh.columns
            if name not in merged:
                merged[name] = Sheet(name, sh.columns if has_sample else ["Sample"] + sh.columns,
                                     [], dict(sh.links), {"Sample": 18, **sh.widths})
            merged[name].rows += sh.rows if has_sample else [{"Sample": smp.sample, **r} for r in sh.rows]
    return list(merged.values())


def _general_synteny(report_data: dict) -> Sheet:
    """ORF ↔ CDS pairs of the general synteny (report_synteny.py): every sample's
    sequences on their group's anchor reference, with a Sample column."""
    gs = report_data.get("general_synteny") or {}
    rows = []
    for q in _seqs(report_data):
        if not q.get("general_synteny"):
            continue
        one = _synteny({"synteny": gs, "sequences": [{"id": q.get("id"),
                                                       "synteny": q["general_synteny"]}]})
        rows += [{"Sample": q.get("sample"), **r} for r in one.rows]
    base = _synteny({})
    return Sheet("General synteny", ["Sample"] + base.columns, rows,
                 links=base.links, widths=base.widths)


def write_multi_tables(samples, report_data: dict, path: str | Path) -> Path:
    """
    The consolidated report's workbook: Samples, every per-sample sheet with a
    Sample column, the cross-sample General clusters and, when it ran, the
    General synteny.
    """
    sheets = merged_sample_sheets(samples)
    at = next((i for i, sh in enumerate(sheets) if sh.name == "Sample clusters"), len(sheets))
    sheets = [_samples_sheet(report_data)] + sheets[:at] + [_general_clusters(report_data)] + sheets[at:]
    if report_data.get("general_synteny"):
        sheets.append(_general_synteny(report_data))
    xlsx = write_xlsx(sheets, path)
    logger.success(f"Results workbook written → '{Path(path).name}' ({len(sheets)} sheets).")
    return xlsx


def write_sheet_tsv(sheet: Sheet, path: str | Path) -> Path:
    """One sheet as a TSV (header + rows)."""
    import csv
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(sheet.columns)
        for r in sheet.rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in sheet.columns])
    return path


# ── XLSX (minimal SpreadsheetML) ────────────────────────────────────────────

_ILLEGAL_XML = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _col(i: int) -> str:
    """0-based column index → Excel letters (0 → A, 26 → AA)."""
    s = ""
    i += 1
    while i:
        i, rem = divmod(i - 1, 26)
        s = chr(65 + rem) + s
    return s


def _cell(ref: str, value, style: int = 0) -> str:
    st = f' s="{style}"' if style else ""
    if value is None or value == "":
        return f'<c r="{ref}"{st}/>' if style else ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if value != value or value in (float("inf"), float("-inf")):   # NaN / inf
            return ""
        return f'<c r="{ref}"{st}><v>{value!r}</v></c>'
    text = _ILLEGAL_XML.sub("", str(value))
    if len(text) > XLSX_CELL_LIMIT:
        text = f"[{len(text):,} characters: too long for an Excel cell, see the JSON]"
    return (f'<c r="{ref}"{st} t="inlineStr"><is><t xml:space="preserve">'
            f"{escape(text)}</t></is></c>")


_STYLES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    '<fonts count="3"><font><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><name val="Calibri"/></font>'
    '<font><u/><sz val="11"/><color rgb="FF0563C1"/><name val="Calibri"/></font></fonts>'
    '<fills count="2"><fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill></fills>'
    '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
    '<cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
    '</styleSheet>'
)

_NS  = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
_RNS = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def _sheet_xml(sh: Sheet) -> tuple[str, list[tuple[str, str]]]:
    """Worksheet XML plus its (cell ref, url) hyperlinks."""
    last = _col(len(sh.columns) - 1)
    xml_rows = ['<row r="1">' + "".join(
        _cell(f"{_col(i)}1", c, 1) for i, c in enumerate(sh.columns)) + "</row>"]
    links: list[tuple[str, str]] = []
    for n, r in enumerate(sh.rows, start=2):
        cells = []
        for i, c in enumerate(sh.columns):
            ref, v = f"{_col(i)}{n}", r.get(c)
            url = ncbi_url(v, sh.links[c]) if c in sh.links and isinstance(v, str) else None
            if url:
                links.append((ref, url))
            cells.append(_cell(ref, v, 2 if url else 0))
        xml_rows.append(f'<row r="{n}">' + "".join(cells) + "</row>")

    cols = "".join(
        f'<col min="{i + 1}" max="{i + 1}" width="{sh.widths.get(c, max(12, min(24, len(c) + 3)))}" customWidth="1"/>'
        for i, c in enumerate(sh.columns))
    hyper = ("<hyperlinks>" + "".join(
        f'<hyperlink ref="{ref}" r:id="rId{k}"/>' for k, (ref, _) in enumerate(links, 1))
        + "</hyperlinks>") if links else ""
    xml = (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet {_NS} {_RNS}>'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        f"<cols>{cols}</cols><sheetData>{''.join(xml_rows)}</sheetData>"
        f'<autoFilter ref="A1:{last}{max(1, len(sh.rows) + 1)}"/>{hyper}</worksheet>'
    )
    return xml, links


def _sheet_name(name: str) -> str:
    """Excel sheet names: ≤ 31 characters, none of []:*?/\\."""
    return re.sub(r"[\[\]:*?/\\]", "-", name)[:31]


def xlsx_parts(sheets: list[Sheet] | list[dict]) -> dict[str, str]:
    """Package part name → XML text. A plain row list is a one-sheet Sequences workbook."""
    if sheets and isinstance(sheets[0], dict) or not sheets:
        sheets = [Sheet("Sequences", COLUMNS, list(sheets), LINK_COLUMNS, _WIDTHS)]

    parts: dict[str, str] = {"xl/styles.xml": _STYLES}
    overrides, wb_rels, wb_sheets, names = [], [], [], []
    for k, sh in enumerate(sheets, 1):
        xml, links = _sheet_xml(sh)
        parts[f"xl/worksheets/sheet{k}.xml"] = xml
        overrides.append(f'<Override PartName="/xl/worksheets/sheet{k}.xml" '
                         'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
        wb_rels.append(f'<Relationship Id="rId{k}" Type="{_REL}/worksheet" Target="worksheets/sheet{k}.xml"/>')
        name = escape(_sheet_name(sh.name))
        wb_sheets.append(f'<sheet name="{name}" sheetId="{k}" r:id="rId{k}"/>')
        last = _col(len(sh.columns) - 1)
        names.append(f'<definedName name="_xlnm._FilterDatabase" localSheetId="{k - 1}" hidden="1">'
                     f"'{name}'!$A$1:${last}${max(1, len(sh.rows) + 1)}</definedName>")
        if links:
            parts[f"xl/worksheets/_rels/sheet{k}.xml.rels"] = (
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                + "".join(
                    f'<Relationship Id="rId{j}" Type="{_REL}/hyperlink" '
                    f'Target="{escape(url)}" TargetMode="External"/>'
                    for j, (_, url) in enumerate(links, 1))
                + "</Relationships>")
    n = len(sheets)
    wb_rels.append(f'<Relationship Id="rId{n + 1}" Type="{_REL}/styles" Target="styles.xml"/>')

    parts["[Content_Types].xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        + "".join(overrides) +
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        "</Types>")
    parts["_rels/.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'<Relationship Id="rId1" Type="{_REL}/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>")
    parts["xl/workbook.xml"] = (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook {_NS} {_RNS}>'
        f"<sheets>{''.join(wb_sheets)}</sheets>"
        f"<definedNames>{''.join(names)}</definedNames></workbook>")
    parts["xl/_rels/workbook.xml.rels"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(wb_rels) + "</Relationships>")
    return parts


def write_xlsx(sheets: list[Sheet] | list[dict], path: str | Path) -> Path:
    path = Path(path)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, xml in xlsx_parts(sheets).items():
            zf.writestr(name, xml)
    return path


def write_tables(report: dict, outdir: str | Path, stem: str) -> Path:
    """Write the results workbook <stem>_viralquest.xlsx into *outdir*."""
    sheets = workbook_sheets(report)
    xlsx = write_xlsx(sheets, Path(outdir) / f"{stem}_viralquest.xlsx")
    logger.success(f"Results workbook written → '{xlsx.name}' "
                   f"({len(sheets)} sheets, {len(sheets[0].rows)} sequences).")
    return xlsx
