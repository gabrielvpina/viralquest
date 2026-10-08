"""Tests for viralquest.seq_table — the per-sequence TSV / XLSX summary table."""
import csv
import json
import zipfile
from pathlib import Path

import pytest

from viralquest.seq_table import (
    COLUMNS,
    XLSX_CELL_LIMIT,
    ncbi_url,
    seq_row,
    table_rows,
    write_tables,
    xlsx_parts,
)

EXAMPLE = Path(__file__).parent / "example.json"


def _seq(**kw):
    seq = {
        "id": "c1", "sequence": "ACGT" * 10, "length": 40,
        "taxonomy": {"order": "Amarillovirales", "family": "Flaviviridae",
                     "genus": "Orthoflavivirus", "species": "Orthoflavivirus denguei",
                     "genome": "ssRNA(+)"},
        "orfs": [{"domains": [{"target": "RdRp", "database": "Pfam"},
                              {"target": "RdRp", "database": "Pfam"},
                              {"target": "Flavi_NS5", "database": "RVDB"}]}],
        "blastn_hits": [
            {"stitle": "weak", "accession": "X1", "pident": 70.0, "qcovhsp": 10,
             "evalue": 1e-3, "bit_score": 50.0},
            {"stitle": "Dengue virus 2, complete genome", "accession": "NC_001474",
             "pident": 98.5, "qcovhsp": 95, "evalue": 0.0, "bit_score": 900.0,
             "subject_seq": "ACGTACGT"},
        ],
        "blastx_hits": [{"subject_title": "refseq", "subject_id": "YP_1.1",
                         "pct_identity": 80.0, "query_coverage": 50.0,
                         "e_value": 1e-20, "bit_score": 100.0}],
        "blastx_nr_hits": [{"subject_title": "polyprotein [Dengue virus 2]",
                            "subject_id": "AAA12345.1", "pct_identity": 99.0,
                            "query_coverage": 97.0, "e_value": 1e-150,
                            "bit_score": 600.0, "subject_seq": "MNNQRK"}],
    }
    seq.update(kw)
    return seq


REPORT = {"meta": {"input_file": {"name": "sampleA.fasta"}}}


class TestRow:

    def test_columns_and_values(self):
        r = seq_row(_seq(), REPORT)
        assert list(r) == COLUMNS
        assert r["Sample"] == "sampleA"
        assert r["Family"] == "Flaviviridae"
        assert r["Genome type"] == "ssRNA(+)"
        assert r["Domains"] == "RdRp (Pfam); Flavi_NS5 (RVDB)"

    def test_best_blastn_hit_by_bit_score(self):
        r = seq_row(_seq(), REPORT)
        assert r["BLASTn hit ID"] == "NC_001474"
        assert r["BLASTn hit sequence"] == "ACGTACGT"
        assert r["BLASTn identity (%)"] == 98.5

    def test_blastx_prefers_nr(self):
        r = seq_row(_seq(), REPORT)
        assert r["BLASTx hit ID"] == "AAA12345.1"
        assert r["BLASTx hit sequence"] == "MNNQRK"

    def test_blastx_falls_back_to_refseq(self):
        r = seq_row(_seq(blastx_nr_hits=[]), REPORT)
        assert r["BLASTx hit ID"] == "YP_1.1"
        assert r["BLASTx hit sequence"] == ""

    def test_no_hits_gives_empty_cells(self):
        r = seq_row(_seq(blastn_hits=[], blastx_hits=[], blastx_nr_hits=[]), REPORT)
        assert r["BLASTn hit"] == "" and r["BLASTn e-value"] is None
        assert r["BLASTx hit ID"] == ""

    def test_taxonomy_falls_back_to_blastn_taxonomy(self):
        r = seq_row(_seq(taxonomy=None,
                         blastn_taxonomy={"taxonomy": {"family": "Parvoviridae"}}), REPORT)
        assert r["Family"] == "Parvoviridae"

    def test_multi_sample_stamp_wins(self):
        assert seq_row(_seq(sample="S2"), REPORT)["Sample"] == "S2"


@pytest.mark.parametrize("acc, db, url", [
    ("NC_001474", "nuccore", "https://www.ncbi.nlm.nih.gov/nuccore/NC_001474"),
    ("YP_009305131.1", "protein", "https://www.ncbi.nlm.nih.gov/protein/YP_009305131.1"),
    ("gi|123|ref|NP_1.1|", "protein", None),
    ("", "nuccore", None),
])
def test_ncbi_url(acc, db, url):
    assert ncbi_url(acc, db) == url


class TestXlsx:

    def test_accessions_are_hyperlinks(self):
        parts = xlsx_parts([seq_row(_seq(), REPORT)])
        rels  = parts["xl/worksheets/_rels/sheet1.xml.rels"]
        assert "https://www.ncbi.nlm.nih.gov/nuccore/NC_001474" in rels
        assert "https://www.ncbi.nlm.nih.gov/protein/AAA12345.1" in rels
        assert parts["xl/worksheets/sheet1.xml"].count("<hyperlink ") == 2

    def test_overlong_cells_are_replaced(self):
        long_seq = "A" * (XLSX_CELL_LIMIT + 1)
        sheet = xlsx_parts([seq_row(_seq(sequence=long_seq), REPORT)])["xl/worksheets/sheet1.xml"]
        assert long_seq not in sheet
        assert "too long for an Excel cell" in sheet

    def test_text_is_xml_escaped(self):
        sheet = xlsx_parts([seq_row(_seq(id="x", sample="a<b&c"), REPORT)])["xl/worksheets/sheet1.xml"]
        assert "a&lt;b&amp;c" in sheet


def test_write_tables_from_example_report(tmp_path):
    report = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    tsv, xlsx = write_tables(report, tmp_path, "ex")
    assert tsv.name == "ex_sequences.tsv" and xlsx.name == "ex_sequences.xlsx"

    with open(tsv, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0] == COLUMNS
    assert len(rows) == len(report["sequences"]) + 1
    assert len(table_rows(report)) == len(report["sequences"])

    with zipfile.ZipFile(xlsx) as zf:
        assert zf.testzip() is None
        assert {"[Content_Types].xml", "xl/workbook.xml",
                "xl/worksheets/sheet1.xml"} <= set(zf.namelist())
