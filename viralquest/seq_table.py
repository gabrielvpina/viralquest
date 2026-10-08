"""
seq_table.py — One-row-per-sequence summary table of a ViralQuest report.

Written next to the JSON/HTML at the end of a run:

    <stem>_sequences.tsv    every value in full (machine-friendly)
    <stem>_sequences.xlsx   same table; BLASTn / BLASTx accessions are links
                            to their NCBI record

The row builder works on the exported report dict (the JSON schema), so it
mirrors exactly what the HTML viewer shows — the viewer's "Table" export
(components*/export.js, vqSeqTableRows) builds the same columns in the browser.
Keep both in sync.

The XLSX is written with the standard library only (zipfile + SpreadsheetML),
so no spreadsheet dependency is needed.
"""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from loguru import logger

COLUMNS = [
    "Sample", "Order", "Family", "Genus", "Species", "Genome type",
    "Sequence (nt)", "Length", "Domains",
    "BLASTn hit", "BLASTn identity (%)", "BLASTn coverage (%)",
    "BLASTn hit ID", "BLASTn hit sequence", "BLASTn e-value",
    "BLASTx hit", "BLASTx identity (%)", "BLASTx coverage (%)",
    "BLASTx hit ID", "BLASTx hit sequence", "BLASTx e-value",
]

# Accession columns → the NCBI database their record lives in.
LINK_COLUMNS = {"BLASTn hit ID": "nuccore", "BLASTx hit ID": "protein"}

# Excel refuses cells longer than this; such values stay complete in the TSV.
XLSX_CELL_LIMIT = 32767

_ACC_RE = re.compile(r"^[A-Za-z0-9_]+(\.[0-9]+)?$")


def ncbi_url(accession: str | None, db: str) -> str | None:
    """NCBI record URL for a bare accession; None for pipe-packed / local ids."""
    if not accession or not _ACC_RE.match(accession):
        return None
    return f"https://www.ncbi.nlm.nih.gov/{db}/{accession}"


# ── Rows ────────────────────────────────────────────────────────────────────

def _best(hits: list[dict]) -> dict | None:
    return max(hits, key=lambda h: h.get("bit_score") or 0) if hits else None


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
    """One table row (column name → value) for an exported sequence record."""
    bt  = seq.get("blastn_taxonomy") or {}
    tax = seq.get("taxonomy") or bt.get("taxonomy") or {}
    bn  = _best(seq.get("blastn_hits") or [])
    bx  = _best(seq.get("blastx_nr_hits") or seq.get("blastx_hits") or [])
    return {
        "Sample":              sample_name(seq, report),
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


# ── TSV ─────────────────────────────────────────────────────────────────────

def _text(v) -> str:
    return "" if v is None else str(v)


def write_tsv(rows: list[dict], path: str | Path) -> Path:
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([_text(r.get(c)) for c in COLUMNS])
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
        return f'<c r="{ref}"{st}><v>{value!r}</v></c>'
    text = _ILLEGAL_XML.sub("", str(value))
    if len(text) > XLSX_CELL_LIMIT:
        text = f"[{len(text):,} characters: too long for an Excel cell, see the TSV]"
    return (f'<c r="{ref}"{st} t="inlineStr"><is><t xml:space="preserve">'
            f"{escape(text)}</t></is></c>")


_WIDTHS = {"Sample": 18, "Species": 30, "BLASTn hit": 40, "BLASTx hit": 40,
           "Domains": 30, "Sequence (nt)": 30, "BLASTn hit sequence": 30,
           "BLASTx hit sequence": 30, "BLASTn hit ID": 16, "BLASTx hit ID": 16}

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


def xlsx_parts(rows: list[dict], sheet: str = "Sequences") -> dict[str, str]:
    """Package part name → XML text for a one-sheet workbook."""
    ns  = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    last = _col(len(COLUMNS) - 1)

    xml_rows = ['<row r="1">' + "".join(
        _cell(f"{_col(i)}1", c, 1) for i, c in enumerate(COLUMNS)) + "</row>"]
    links: list[tuple[str, str]] = []
    for n, r in enumerate(rows, start=2):
        cells = []
        for i, c in enumerate(COLUMNS):
            ref, v = f"{_col(i)}{n}", r.get(c)
            url = ncbi_url(v, LINK_COLUMNS[c]) if c in LINK_COLUMNS else None
            if url:
                links.append((ref, url))
            cells.append(_cell(ref, v, 2 if url else 0))
        xml_rows.append(f'<row r="{n}">' + "".join(cells) + "</row>")

    cols = "".join(
        f'<col min="{i + 1}" max="{i + 1}" width="{_WIDTHS.get(c, 12)}" customWidth="1"/>'
        for i, c in enumerate(COLUMNS))
    hyper = ("<hyperlinks>" + "".join(
        f'<hyperlink ref="{ref}" r:id="rId{k}"/>' for k, (ref, _) in enumerate(links, 1))
        + "</hyperlinks>") if links else ""
    sheet_xml = (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet {ns} {rns}>'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        f"<cols>{cols}</cols><sheetData>{''.join(xml_rows)}</sheetData>"
        f'<autoFilter ref="A1:{last}{max(1, len(rows) + 1)}"/>{hyper}</worksheet>'
    )

    rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    parts = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            "</Types>"),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/>'
            "</Relationships>"),
        "xl/workbook.xml": (
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook {ns} {rns}>'
            f'<sheets><sheet name="{escape(sheet)}" sheetId="1" r:id="rId1"/></sheets>'
            '<definedNames><definedName name="_xlnm._FilterDatabase" localSheetId="0" hidden="1">'
            f"'{escape(sheet)}'!$A$1:${last}${max(1, len(rows) + 1)}</definedName></definedNames>"
            "</workbook>"),
        "xl/_rels/workbook.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/worksheet" Target="worksheets/sheet1.xml"/>'
            f'<Relationship Id="rId2" Type="{rel}/styles" Target="styles.xml"/>'
            "</Relationships>"),
        "xl/styles.xml": _STYLES,
        "xl/worksheets/sheet1.xml": sheet_xml,
    }
    if links:
        parts["xl/worksheets/_rels/sheet1.xml.rels"] = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join(
                f'<Relationship Id="rId{k}" Type="{rel}/hyperlink" '
                f'Target="{escape(url)}" TargetMode="External"/>'
                for k, (_, url) in enumerate(links, 1))
            + "</Relationships>")
    return parts


def write_xlsx(rows: list[dict], path: str | Path) -> Path:
    path = Path(path)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, xml in xlsx_parts(rows).items():
            zf.writestr(name, xml)
    return path


def write_tables(report: dict, outdir: str | Path, stem: str) -> tuple[Path, Path]:
    """Write <stem>_sequences.tsv and .xlsx into *outdir*; returns both paths."""
    outdir = Path(outdir)
    rows = table_rows(report)
    tsv  = write_tsv(rows, outdir / f"{stem}_sequences.tsv")
    xlsx = write_xlsx(rows, outdir / f"{stem}_sequences.xlsx")
    logger.success(f"Sequence table written → '{tsv.name}', '{xlsx.name}' ({len(rows)} rows).")
    return tsv, xlsx
