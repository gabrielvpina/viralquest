"""
report_db.py — Named, merge-ready exports of a ``viralquest-report`` run.

Every consolidated report is a *dataset* with a user-given name (a BioProject,
a host species, a study …).  Next to the HTML it is written as:

* ``<name>.viralquest.json``   — the full report data (meta, per-sample
  summaries, every sequence record, cross-sample clusters), with the dataset
  name and schema version at the top;
* ``<name>.viralquest.sqlite`` — the same data as relational tables.

Built to be merged later: every SQLite row carries the ``dataset`` column and
every primary key starts with it, so datasets can be combined with ATTACH +
INSERT (or UNION) without id clashes — sample names, contig ids and cluster ids
(VQR_CLU_0001 …) are only unique *within* one dataset.  ``dataset_info`` holds
one row per dataset with the schema version, so a merge tool can refuse
duplicate names or incompatible schemas.
"""

from __future__ import annotations

import json
import re
import sqlite3
from pathlib import Path

SCHEMA_VERSION = 1
DATASET_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")


def validate_dataset_name(name: str) -> str:
    """Return the name when it is a valid dataset name, else raise ValueError."""
    if not isinstance(name, str) or not DATASET_RE.match(name):
        raise ValueError(
            f"invalid dataset name {name!r}: use 1–64 letters, digits, '.', '_' or '-', "
            "starting with a letter or digit (e.g. PRJNA123456, Aedes_aegypti_2024)"
        )
    return name


def output_paths(out_dir: Path, name: str) -> tuple[Path, Path]:
    return out_dir / f"{name}.viralquest.json", out_dir / f"{name}.viralquest.sqlite"


# ── JSON ────────────────────────────────────────────────────────────────────

def write_json(report: dict, path: Path) -> Path:
    meta = report.get("meta", {})
    doc = {
        "schema":         "viralquest-report",
        "schema_version": SCHEMA_VERSION,
        "dataset":        meta.get("dataset"),
        "meta":           meta,
        "samples":        report.get("samples", []),
        "clusters":       report.get("clusters", []),
        "sequences":      report.get("sequences", []),
    }
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, separators=(",", ":"), default=str)
    tmp.replace(path)
    return path


# ── SQLite ──────────────────────────────────────────────────────────────────

SCHEMA = """
CREATE TABLE dataset_info (
    dataset            TEXT PRIMARY KEY,
    schema_version     INTEGER NOT NULL,
    viralquest_version TEXT,
    generated          TEXT,
    input_root         TEXT,
    n_samples          INTEGER,
    n_sequences        INTEGER,
    n_clusters         INTEGER
);
CREATE TABLE samples (
    dataset            TEXT NOT NULL,
    sample             TEXT NOT NULL,
    n_input            INTEGER,
    n_viral_flagged    INTEGER,
    n_confirmed        INTEGER,
    n_clusters         INTEGER,
    salmon_pathway     TEXT,
    salmon_mapping_rate REAL,
    salmon_input_reads INTEGER,
    salmon_viral_reads REAL,
    viralquest_version TEXT,
    run_started        TEXT,
    run_seconds        REAL,
    PRIMARY KEY (dataset, sample)
);
CREATE TABLE sequences (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,     -- sample::seq_id, unique within the dataset
    sample             TEXT NOT NULL,
    seq_id             TEXT NOT NULL,
    length             INTEGER,
    gc_content         REAL,
    n_count            INTEGER,
    is_viral           INTEGER,
    tax_id             INTEGER,           -- taxonomy of the best BLASTx hit
    realm              TEXT,
    kingdom            TEXT,
    phylum             TEXT,
    class              TEXT,
    "order"            TEXT,
    family             TEXT,
    subfamily          TEXT,
    genus              TEXT,
    species            TEXT,
    genome             TEXT,
    heuristic_score    REAL,
    classification     TEXT,
    novelty            TEXT,
    novelty_flags      TEXT,              -- ";"-separated
    llm_score          REAL,
    llm_classification TEXT,
    tpm                REAL,
    reads              REAL,
    rpm                REAL,              -- reads per million input reads (Salmon)
    sequence           TEXT,
    PRIMARY KEY (dataset, gid)
);
CREATE TABLE blastx_hits (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    source             TEXT NOT NULL,     -- refseq | nr
    rank               INTEGER NOT NULL,  -- 1 = best by bit score
    subject_id         TEXT,
    subject_title      TEXT,
    species            TEXT,
    pct_identity       REAL,
    query_coverage     REAL,
    e_value            REAL,
    bit_score          REAL,
    query_start        INTEGER,
    query_end          INTEGER,
    PRIMARY KEY (dataset, gid, source, rank)
);
CREATE TABLE blastn_hits (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    rank               INTEGER NOT NULL,
    accession          TEXT,
    title              TEXT,
    pct_identity       REAL,
    query_coverage     REAL,
    e_value            REAL,
    bit_score          REAL,
    query_start        INTEGER,
    query_end          INTEGER,
    PRIMARY KEY (dataset, gid, rank)
);
CREATE TABLE blastn_taxonomy (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    status             TEXT,
    method             TEXT,
    accession          TEXT,
    tax_id             INTEGER,
    matched_name       TEXT,
    synonym_of         TEXT,
    family             TEXT,
    genus              TEXT,
    species            TEXT,
    PRIMARY KEY (dataset, gid)
);
CREATE TABLE orfs (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    orf                TEXT NOT NULL,
    start              INTEGER,
    stop               INTEGER,
    strand             TEXT,
    frame              INTEGER,
    orf_type           TEXT,
    length_aa          INTEGER,
    PRIMARY KEY (dataset, gid, orf)
);
CREATE TABLE domains (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    orf                TEXT NOT NULL,
    database           TEXT,
    target             TEXT,
    description        TEXT,
    score              REAL,
    e_value            REAL,
    start              INTEGER,
    stop               INTEGER
);
CREATE TABLE clusters (
    dataset            TEXT NOT NULL,
    cluster_id         TEXT NOT NULL,     -- VQR_CLU_#### — unique within the dataset only
    species            TEXT,
    species_db         TEXT,
    representative     TEXT,              -- gid
    n_members          INTEGER,
    n_samples          INTEGER,
    homogeneous        INTEGER,
    PRIMARY KEY (dataset, cluster_id)
);
CREATE TABLE cluster_members (
    dataset            TEXT NOT NULL,
    cluster_id         TEXT NOT NULL,
    gid                TEXT NOT NULL,
    sample             TEXT,
    identity           REAL,
    coverage           REAL,
    is_representative  INTEGER,
    tpm                REAL,
    reads              REAL,
    PRIMARY KEY (dataset, cluster_id, gid)
);
CREATE INDEX idx_sequences_sample  ON sequences (dataset, sample);
CREATE INDEX idx_sequences_family  ON sequences (family);
CREATE INDEX idx_sequences_species ON sequences (species);
CREATE INDEX idx_domains_gid       ON domains (dataset, gid);
CREATE INDEX idx_members_gid       ON cluster_members (dataset, gid);
"""


def _num(v):
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def _text(v):
    return v if isinstance(v, str) else None


def _ranked(hits: list[dict], key: str) -> list[tuple[int, dict]]:
    return list(enumerate(sorted(hits or [], key=lambda h: -(h.get(key) or 0)), 1))


def write_sqlite(report: dict, path: Path) -> Path:
    meta = report.get("meta", {})
    ds = validate_dataset_name(meta.get("dataset"))
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    con = sqlite3.connect(tmp)
    try:
        con.executescript(SCHEMA)
        con.execute("PRAGMA user_version = %d" % SCHEMA_VERSION)

        samples = report.get("samples", []) or []
        input_reads = {s["sample"]: (s.get("salmon") or {}).get("input_reads") for s in samples}
        con.execute(
            "INSERT INTO dataset_info VALUES (?,?,?,?,?,?,?,?)",
            (ds, SCHEMA_VERSION, meta.get("viralquest_version"), meta.get("generated"),
             meta.get("input_root"), len(samples), len(report.get("sequences", []) or []),
             len(report.get("clusters", []) or [])))

        con.executemany(
            "INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [(ds, s["sample"], s.get("n_input"), s.get("n_viral"), s.get("n_confirmed"), s.get("n_clusters"),
              (s.get("salmon") or {}).get("pathway"), _num((s.get("salmon") or {}).get("mapping_rate")),
              (s.get("salmon") or {}).get("input_reads"), _num((s.get("salmon") or {}).get("viral_reads")),
              (s.get("workflow") or {}).get("version"), (s.get("workflow") or {}).get("started_at"),
              _num((s.get("workflow") or {}).get("total_seconds")))
             for s in samples])

        seq_rows, bx_rows, bn_rows, bnt_rows, orf_rows, dom_rows = [], [], [], [], [], []
        for q in report.get("sequences", []) or []:
            gid, sample = q.get("gid"), q.get("sample")
            t = q.get("taxonomy") or {}
            h = q.get("heuristic_output") or {}
            llm = q.get("llm_output") or {}
            reads, inp = q.get("reads"), input_reads.get(sample)
            seq_rows.append((
                ds, gid, sample, q.get("id"), q.get("length"), _num(q.get("gc_content")), q.get("n_count"),
                int(bool(q.get("is_viral"))), t.get("tax_id"), t.get("clade"), t.get("kingdom"), t.get("phylum"),
                t.get("class"), t.get("order"), t.get("family"), t.get("subfamily"), t.get("genus"),
                t.get("species") or t.get("scientific_name"), t.get("genome"),
                _num(h.get("vq_score")), h.get("classification"), h.get("novelty"),
                ";".join(h.get("novelty_flags") or []) or None,
                _num(llm.get("vq_score")), llm.get("classification"),
                _num(q.get("tpm")), _num(reads),
                (reads / inp * 1e6) if reads is not None and inp else None,
                q.get("sequence"),
            ))
            for source, key in (("refseq", "blastx_hits"), ("nr", "blastx_nr_hits")):
                for rank, x in _ranked(q.get(key), "bit_score"):
                    bx_rows.append((ds, gid, source, rank, x.get("subject_id"), x.get("subject_title"),
                                    x.get("species"), _num(x.get("pct_identity")), _num(x.get("query_coverage")),
                                    _num(x.get("e_value")), _num(x.get("bit_score")),
                                    x.get("query_start"), x.get("query_end")))
            for rank, x in _ranked(q.get("blastn_hits"), "bit_score"):
                bn_rows.append((ds, gid, rank, x.get("accession"), x.get("stitle"), _num(x.get("pident")),
                                _num(x.get("qcovhsp")), _num(x.get("evalue")), _num(x.get("bit_score")),
                                x.get("query_start"), x.get("query_end")))
            bt = q.get("blastn_taxonomy")
            if bt:
                bt_t = bt.get("taxonomy") or {}
                bnt_rows.append((ds, gid, bt.get("status"), bt.get("method"), bt.get("accession"), bt.get("taxid"),
                                 bt.get("matched_name"), bt.get("synonym_of"), bt_t.get("family"), bt_t.get("genus"),
                                 bt_t.get("species") or bt_t.get("scientific_name")))
            for o in q.get("orfs") or []:
                name = o.get("name")
                orf_rows.append((ds, gid, name, o.get("start_position"), o.get("stop_position"), o.get("strand"),
                                 o.get("frame"), o.get("orf_type"), o.get("length_aa")))
                for d in o.get("domains") or []:
                    dom_rows.append((ds, gid, name, d.get("database"), d.get("target"), _text(d.get("description")),
                                     _num(d.get("score")), _num(d.get("e_value")), d.get("start"), d.get("stop")))

        con.executemany(f"INSERT INTO sequences VALUES ({','.join('?' * 29)})", seq_rows)
        con.executemany("INSERT OR IGNORE INTO blastx_hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", bx_rows)
        con.executemany("INSERT OR IGNORE INTO blastn_hits VALUES (?,?,?,?,?,?,?,?,?,?,?)", bn_rows)
        con.executemany("INSERT INTO blastn_taxonomy VALUES (?,?,?,?,?,?,?,?,?,?,?)", bnt_rows)
        con.executemany("INSERT OR IGNORE INTO orfs VALUES (?,?,?,?,?,?,?,?,?)", orf_rows)
        con.executemany("INSERT INTO domains VALUES (?,?,?,?,?,?,?,?,?,?)", dom_rows)

        for c in report.get("clusters", []) or []:
            members = c.get("members", []) or []
            con.execute("INSERT INTO clusters VALUES (?,?,?,?,?,?,?,?)",
                        (ds, c.get("gid"), c.get("species"), c.get("species_db"), c.get("representative"),
                         len(members), len(c.get("samples") or []), int(bool(c.get("homogeneous")))))
            con.executemany(
                "INSERT OR IGNORE INTO cluster_members VALUES (?,?,?,?,?,?,?,?,?)",
                [(ds, c.get("gid"), m.get("gid"), m.get("sample"), _num(m.get("identity")),
                  _num(m.get("coverage")), int(bool(m.get("is_representative"))),
                  _num(m.get("tpm")), _num(m.get("reads"))) for m in members])
        con.commit()
    finally:
        con.close()
    tmp.replace(path)
    return path
