#!/usr/bin/env python3
"""
preview_report.py — Generate a mock ViralQuest HTML report for UI development.

Usage:
    python preview_report.py                  # writes preview.html
    python preview_report.py out.html         # custom output path
    python preview_report.py out.html --open  # open in browser automatically

The report contains realistic fake data covering all UI sections:
  - Stats, Clusters, All-Viruses viewer, Taxonomy tree, Salmon quant
"""

import json
import sys
import webbrowser
from pathlib import Path

_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT))

from viralquest.html_report import write_report

# ---------------------------------------------------------------------------
# Mock report data
# ---------------------------------------------------------------------------

REPORT: dict = {
    "meta": {
        "viralquest_version": "3.0.0",
        "timestamp":          "2026-05-26T14:30:00",
        "input_file": {
            "name":       "sample_metagenome.fasta",
            "size_bytes": 4_821_033,
            "n_sequences": 312,
        },
        "pipeline": {
            "cap3":        False,
            "blastn":      "local",
            "nr_db":       True,
            "salmon":      True,
            "ai_scoring":  False,
        },
    },

    # ── Sequences ──────────────────────────────────────────────────────────
    "sequences": [
        {
            "seq_id":   "contig_001",
            "length":   4820,
            "is_viral": True,
            "gc_content": 41.2,
            "taxonomy": {
                "phylum": "Negarnaviricota",
                "order":  "Amarillovirales",
                "family": "Flaviviridae",
                "genus":  "Flavivirus",
                "species": "Dengue virus 2",
            },
            "blastx_hits": [{
                "subject_title": "polyprotein [Dengue virus 2]",
                "pct_identity":  92.4,
                "e_value":       1e-120,
                "bit_score":     780.3,
                "query_coverage": 88.5,
            }],
            "blastx_nr_hits": [],
            "blastn_hits": [{
                "subject_title": "Dengue virus 2, complete genome",
                "pct_identity":  95.1,
                "e_value":       0.0,
                "bit_score":     8421.0,
                "query_coverage": 100.0,
            }],
            "orfs": [
                {
                    "name":    "contig_001|orf_1",
                    "start":   120,
                    "end":     4700,
                    "strand":  "+",
                    "length":  4580,
                    "domains": [
                        {
                            "database": "Pfam",
                            "target":   "PF00949",
                            "description": "Flavivirus polyprotein",
                            "score": 412.3,
                            "e_value": 1e-128,
                            "start": 10,
                            "stop":  1520,
                            "length": 1510,
                        },
                        {
                            "database": "RVDB",
                            "target":   "RVDB_00341",
                            "description": "RNA-dependent RNA polymerase",
                            "score": 298.1,
                            "e_value": 1e-90,
                            "start": 2800,
                            "stop":  4400,
                            "length": 1600,
                        },
                    ],
                }
            ],
            "salmon_tpm":      1842.3,
            "salmon_num_reads": 9241,
            "cluster_id": "CLU_001",
            "ai_score":   None,
        },
        {
            "seq_id":   "contig_002",
            "length":   2340,
            "is_viral": True,
            "gc_content": 38.9,
            "taxonomy": {
                "phylum": "Negarnaviricota",
                "order":  "Bunyavirales",
                "family": "Phenuiviridae",
                "genus":  "Phlebovirus",
                "species": "Uukuniemi virus",
            },
            "blastx_hits": [{
                "subject_title": "nucleocapsid protein [Uukuniemi virus]",
                "pct_identity":  78.3,
                "e_value":       3e-55,
                "bit_score":     214.8,
                "query_coverage": 72.1,
            }],
            "blastx_nr_hits": [],
            "blastn_hits": [],
            "orfs": [
                {
                    "name":   "contig_002|orf_1",
                    "start":  40,
                    "end":    2200,
                    "strand": "+",
                    "length": 2160,
                    "domains": [
                        {
                            "database": "Vfam",
                            "target":   "Vfam_00128",
                            "description": "Phlebovirus nucleocapsid",
                            "score": 188.5,
                            "e_value": 2e-58,
                            "start": 10,
                            "stop":  700,
                            "length": 690,
                        }
                    ],
                }
            ],
            "salmon_tpm":      312.7,
            "salmon_num_reads": 1564,
            "cluster_id": "CLU_002",
            "ai_score":   None,
        },
        {
            "seq_id":   "contig_003",
            "length":   1820,
            "is_viral": True,
            "gc_content": 44.5,
            "taxonomy": {
                "phylum": "Pisuviricota",
                "order":  "Picornavirales",
                "family": "Nodaviridae",
                "genus":  "Alphanodavirus",
                "species": "Flock House virus",
            },
            "blastx_hits": [{
                "subject_title": "coat protein [Flock House virus]",
                "pct_identity":  65.2,
                "e_value":       8e-38,
                "bit_score":     162.1,
                "query_coverage": 58.3,
            }],
            "blastx_nr_hits": [],
            "blastn_hits": [],
            "orfs": [
                {
                    "name":   "contig_003|orf_1",
                    "start":  80,
                    "end":    1750,
                    "strand": "+",
                    "length": 1670,
                    "domains": [
                        {
                            "database": "EggNOG",
                            "target":   "COG5412",
                            "description": "Nodavirus coat protein B2",
                            "score": 142.0,
                            "e_value": 1e-42,
                            "start": 5,
                            "stop":  550,
                            "length": 545,
                        }
                    ],
                }
            ],
            "salmon_tpm":      87.4,
            "salmon_num_reads": 437,
            "cluster_id": "CLU_003",
            "ai_score":   None,
        },
        {
            "seq_id":   "contig_014",
            "length":   980,
            "is_viral": True,
            "gc_content": 52.1,
            "taxonomy": {
                "phylum": "Negarnaviricota",
                "order":  "Amarillovirales",
                "family": "Flaviviridae",
                "genus":  "Flavivirus",
                "species": "West Nile virus",
            },
            "blastx_hits": [{
                "subject_title": "envelope protein [West Nile virus]",
                "pct_identity":  81.0,
                "e_value":       2e-70,
                "bit_score":     289.4,
                "query_coverage": 91.0,
            }],
            "blastx_nr_hits": [],
            "blastn_hits": [],
            "orfs": [
                {
                    "name":   "contig_014|orf_1",
                    "start":  30,
                    "end":    950,
                    "strand": "+",
                    "length": 920,
                    "domains": [
                        {
                            "database": "Pfam",
                            "target":   "PF00949",
                            "description": "Flavivirus envelope",
                            "score": 201.5,
                            "e_value": 5e-62,
                            "start": 5,
                            "stop":  300,
                            "length": 295,
                        }
                    ],
                }
            ],
            "salmon_tpm":      621.8,
            "salmon_num_reads": 3109,
            "cluster_id": "CLU_001",
            "ai_score":   None,
        },
    ],

    # ── Clusters ───────────────────────────────────────────────────────────
    "clusters": [
        {
            "cluster_id":     "CLU_001",
            "n_sequences":    2,
            "representative": "contig_001",
            "family":         "Flaviviridae",
            "genus":          "Flavivirus",
            "total_length":   5800,
        },
        {
            "cluster_id":     "CLU_002",
            "n_sequences":    1,
            "representative": "contig_002",
            "family":         "Phenuiviridae",
            "genus":          "Phlebovirus",
            "total_length":   2340,
        },
        {
            "cluster_id":     "CLU_003",
            "n_sequences":    1,
            "representative": "contig_003",
            "family":         "Nodaviridae",
            "genus":          "Alphanodavirus",
            "total_length":   1820,
        },
    ],

    # ── Salmon quant ───────────────────────────────────────────────────────
    "salmon_quant": {
        "mapping_rate":  72.4,
        "total_reads":   2_481_039,
        "viral_quant": [
            {"seq_id": "contig_001", "tpm": 1842.3, "num_reads": 9241, "length": 4820},
            {"seq_id": "contig_014", "tpm":  621.8, "num_reads": 3109, "length":  980},
            {"seq_id": "contig_002", "tpm":  312.7, "num_reads": 1564, "length": 2340},
            {"seq_id": "contig_003", "tpm":   87.4, "num_reads":  437, "length": 1820},
        ],
        "conserved_quant": [
            {"gene": "GAPDH",   "tpm":  4821.0, "num_reads": 24105},
            {"gene": "ACTB",    "tpm":  3240.5, "num_reads": 16202},
            {"gene": "RPL13A",  "tpm":  2180.2, "num_reads": 10901},
        ],
        "host_viral_hits": [
            {"seq_id": "contig_001", "host_gene": "GAPDH", "pct_identity": 28.4},
        ],
    },
}


# ---------------------------------------------------------------------------
# Build & open
# ---------------------------------------------------------------------------

def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("preview.html")
    open_browser = "--open" in sys.argv

    print(f"Building preview report → {out}")
    write_report(REPORT, out)
    print(f"Done. Open in browser: file://{out.resolve()}")

    if open_browser:
        webbrowser.open(f"file://{out.resolve()}")


if __name__ == "__main__":
    main()
