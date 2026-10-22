"""Tests for viralquest.seq_table — the aggregated results workbook."""
import json
import re
import zipfile
from pathlib import Path

import pytest

from viralquest.seq_table import (
    COLUMNS,
    XLSX_CELL_LIMIT,
    ncbi_url,
    seq_row,
    table_rows,
    workbook_sheets,
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
        assert r["Sequence ID"] == "c1"
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
    xlsx = write_tables(report, tmp_path, "ex")
    assert xlsx.name == "ex_viralquest.xlsx"
    assert not list(tmp_path.glob("*.tsv"))
    assert len(table_rows(report)) == len(report["sequences"])

    sheets = workbook_sheets(report)
    with zipfile.ZipFile(xlsx) as zf:
        assert zf.testzip() is None
        names = set(zf.namelist())
        assert {f"xl/worksheets/sheet{k}.xml" for k in range(1, len(sheets) + 1)} <= names
        wb = zf.read("xl/workbook.xml").decode()
    assert re.findall(r'<sheet name="([^"]+)"', wb) == [sh.name for sh in sheets]


class TestSheets:

    def _report(self):
        seq = _seq()
        seq["orfs"] = [{"name": "c1_ORF_1", "start_position": 0, "stop_position": 30, "strand": "+",
                        "frame": 1, "length_aa": 10, "length_nt": 30, "aa_sequence": "M" * 10,
                        "domains": [{"target": "RdRp", "database": "Pfam", "score": 50.0,
                                     "e_value": 1e-10, "start": 1, "stop": 9}]}]
        seq["heuristic_output"] = {"vq_score": 90, "classification": "viral-known",
                                   "novelty": "variant", "novelty_flags": ["low-coverage"]}
        seq["coverage"] = {"length": 40, "mean_depth": 12.5, "max_depth": 30, "cv": 0.2,
                           "breadth_1x": 0.975, "low_cov_regions": [[0, 4]]}
        seq["seq_quality"] = {"low_complexity_frac": 0.1, "low_complexity_regions": [[2, 5]],
                              "self_repeats": [{"strand": "plus"}, {"strand": "minus"},
                                               {"strand": "minus"}],
                              "kmer_size": 15, "kmer_repeat_score": 0.01}
        return {
            "meta": REPORT["meta"], "sequences": [seq],
            "clusters": [{"cluster_id": "VQ_CLU_0001", "species": "Dengue", "representative_id": "c1",
                          "size": 1, "members": [{"seq_id": "c1", "is_representative": True,
                                                  "identity": 100.0, "query_coverage": 100.0}]}],
            "salmon_quant": {
                "viral_quant": [{"name": "VQ_VIRAL_c1", "tpm": 10.0, "num_reads": 5.0, "seq_type": "viral"}],
                "ref_hk_quant": [{"name": "VQ_REFHK_ACTB", "tpm": 3.0, "num_reads": 2.0, "seq_type": "ref_hk"}],
                "conserved_quant": [], "host_viral_hits": [
                    {"transcript": {"name": "XM_1", "tpm": 1.0, "num_reads": 1.0, "seq_type": "host"},
                     "blast_hits": [{"host_transcript_id": "XM_1", "viral_seq_id": "c1", "pident": 90.0}]}],
            },
        }

    def test_every_sheet_has_its_rows(self):
        sh = {s.name: s for s in workbook_sheets(self._report())}
        assert sh["Classification"].rows[0]["Novelty"] == "variant"
        assert sh["Classification"].rows[0]["Cluster"] == "VQ_CLU_0001"
        assert sh["ORFs"].rows[0]["Protein sequence"] == "M" * 10
        assert sh["Domains"].rows[0]["Target"] == "RdRp"
        assert [r["Database"] for r in sh["BLASTx hits"].rows] == ["RefSeq", "NR"]
        assert [r["Accession"] for r in sh["BLASTn hits"].rows] == ["NC_001474", "X1"]
        assert sh["Genome coverage"].rows[0]["Breadth ≥1x (%)"] == 97.5
        assert sh["Genome coverage"].rows[0]["Low-coverage runs"] == "1-5"
        q = sh["Contig quality"].rows[0]
        assert (q["Direct repeats"], q["Inverted repeats"]) == (1, 2)
        assert [r["Type"] for r in sh["Quantification"].rows] == ["viral", "ref_hk", "host"]
        assert sh["Quantification"].rows[0]["Sequence ID"] == "c1"
        assert sh["Host-viral hits"].rows[0]["Sequence ID"] == "c1"

    def test_hit_accessions_link_in_every_sheet(self):
        parts = xlsx_parts(workbook_sheets(self._report()))
        rels = "".join(v for k, v in parts.items() if k.startswith("xl/worksheets/_rels/"))
        assert "nuccore/NC_001474" in rels and "protein/AAA12345.1" in rels and "protein/YP_1.1" in rels
