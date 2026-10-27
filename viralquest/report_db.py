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

SCHEMA_VERSION = 2   # 2: + coverage, seq_quality, synteny_hits, quantification, host_viral_hits
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
    if report.get("general_synteny"):          # report_synteny.py; sequences carry theirs
        doc["general_synteny"] = report["general_synteny"]
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
CREATE TABLE coverage (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    length             INTEGER,
    mean_depth         REAL,
    max_depth          REAL,
    cv                 REAL,
    breadth_1x         REAL,              -- fraction of bases with depth >= 1
    mean_depth_sense   REAL,
    mean_depth_antisense REAL,
    mean_quality       REAL,
    n_low_cov_regions  INTEGER,
    PRIMARY KEY (dataset, gid)
);
CREATE TABLE seq_quality (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    length             INTEGER,
    low_complexity_frac REAL,
    n_low_complexity   INTEGER,
    n_direct_repeats   INTEGER,
    n_inverted_repeats INTEGER,
    kmer_size          INTEGER,
    kmer_distinct      INTEGER,
    kmer_total         INTEGER,
    kmer_repeat_score  REAL,
    PRIMARY KEY (dataset, gid)
);
CREATE TABLE synteny_hits (
    dataset            TEXT NOT NULL,
    gid                TEXT NOT NULL,
    reference          TEXT,              -- accession of the best viral BLASTn hit
    flipped            INTEGER,           -- 1 = contig is the reverse complement
    orf                TEXT,
    cds_product        TEXT,
    protein_id         TEXT,
    pct_identity       REAL,
    e_value            REAL,
    bit_score          REAL,
    orf_aa_start       INTEGER,
    orf_aa_end         INTEGER,
    cds_aa_start       INTEGER,
    cds_aa_end         INTEGER,
    orf_coverage       REAL,
    cds_coverage       REAL,
    ref_start          INTEGER,           -- 1-based, on the reference genome
    ref_end            INTEGER
);
CREATE TABLE quantification (
    dataset            TEXT NOT NULL,
    sample             TEXT NOT NULL,
    target             TEXT NOT NULL,     -- Salmon target name (VQ_VIRAL_…, VQ_REFHK_…, VQ_CONS_…, host)
    seq_type           TEXT,              -- viral | ref_hk | conserved | pfam_hk | host
    kingdom            TEXT,
    seq_id             TEXT,              -- viral targets only
    length             INTEGER,
    eff_length         REAL,
    tpm                REAL,
    reads              REAL,
    PRIMARY KEY (dataset, sample, target)
);
CREATE TABLE host_viral_hits (
    dataset            TEXT NOT NULL,
    sample             TEXT NOT NULL,
    host_transcript    TEXT,
    seq_id             TEXT,
    pct_identity       REAL,
    coverage           REAL,
    e_value            REAL,
    bit_score          REAL,
    transcript_tpm     REAL,
    transcript_reads   REAL
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


def _sequence_fields(q: dict, tpm, reads, input_reads) -> tuple:
    """A sequence's columns after its key: length … sequence (shared by both schemas)."""
    t = q.get("taxonomy") or {}
    h = q.get("heuristic_output") or {}
    llm = q.get("llm_output") or {}
    return (
        q.get("length"), _num(q.get("gc_content")), q.get("n_count"),
        int(bool(q.get("is_viral"))), t.get("tax_id"), t.get("clade"), t.get("kingdom"), t.get("phylum"),
        t.get("class"), t.get("order"), t.get("family"), t.get("subfamily"), t.get("genus"),
        t.get("species") or t.get("scientific_name"), t.get("genome"),
        _num(h.get("vq_score")), h.get("classification"), h.get("novelty"),
        ";".join(h.get("novelty_flags") or []) or None,
        _num(llm.get("vq_score")), llm.get("classification"),
        _num(tpm), _num(reads),
        (reads / input_reads * 1e6) if reads is not None and input_reads else None,
        q.get("sequence"),
    )


def _child_rows(q: dict) -> dict[str, list[tuple]]:
    """Hit / taxonomy / ORF / domain rows of one sequence, without their key columns."""
    out: dict[str, list[tuple]] = {"blastx": [], "blastn": [], "blastn_taxonomy": [],
                                   "orfs": [], "domains": []}
    for source, key in (("refseq", "blastx_hits"), ("nr", "blastx_nr_hits")):
        for rank, x in _ranked(q.get(key), "bit_score"):
            out["blastx"].append((source, rank, x.get("subject_id"), x.get("subject_title"),
                                  x.get("species"), _num(x.get("pct_identity")), _num(x.get("query_coverage")),
                                  _num(x.get("e_value")), _num(x.get("bit_score")),
                                  x.get("query_start"), x.get("query_end")))
    for rank, x in _ranked(q.get("blastn_hits"), "bit_score"):
        out["blastn"].append((rank, x.get("accession"), x.get("stitle"), _num(x.get("pident")),
                              _num(x.get("qcovhsp")), _num(x.get("evalue")), _num(x.get("bit_score")),
                              x.get("query_start"), x.get("query_end")))
    bt = q.get("blastn_taxonomy")
    if bt:
        bt_t = bt.get("taxonomy") or {}
        out["blastn_taxonomy"].append((bt.get("status"), bt.get("method"), bt.get("accession"), bt.get("taxid"),
                                       bt.get("matched_name"), bt.get("synonym_of"), bt_t.get("family"),
                                       bt_t.get("genus"), bt_t.get("species") or bt_t.get("scientific_name")))
    for o in q.get("orfs") or []:
        name = o.get("name")
        out["orfs"].append((name, o.get("start_position"), o.get("stop_position"), o.get("strand"),
                            o.get("frame"), o.get("orf_type"), o.get("length_aa")))
        for d in o.get("domains") or []:
            out["domains"].append((name, d.get("database"), d.get("target"), _text(d.get("description")),
                                   _num(d.get("score")), _num(d.get("e_value")), d.get("start"), d.get("stop")))
    return out


def _analysis_rows(q: dict, refs: dict) -> dict[str, list[tuple]]:
    """Coverage / sequence-quality / synteny rows of one sequence, without key columns.
    *refs*: that sample's synteny references (accession → reference record)."""
    out: dict[str, list[tuple]] = {"coverage": [], "seq_quality": [], "synteny": []}
    c = q.get("coverage")
    if c:
        out["coverage"].append((c.get("length"), _num(c.get("mean_depth")), _num(c.get("max_depth")),
                                _num(c.get("cv")), _num(c.get("breadth_1x")), _num(c.get("mean_depth_fwd")),
                                _num(c.get("mean_depth_rev")), _num(c.get("mean_quality")),
                                len(c.get("low_cov_regions") or [])))
    sq = q.get("seq_quality")
    if sq:
        reps = sq.get("self_repeats") or []
        out["seq_quality"].append((sq.get("length"), _num(sq.get("low_complexity_frac")),
                                   len(sq.get("low_complexity_regions") or []),
                                   sum(1 for r in reps if r.get("strand") == "plus"),
                                   sum(1 for r in reps if r.get("strand") == "minus"),
                                   sq.get("kmer_size"), sq.get("kmer_distinct"), sq.get("kmer_total"),
                                   _num(sq.get("kmer_repeat_score"))))
    syn = q.get("synteny") or {}
    cds = (refs.get(syn.get("reference")) or {}).get("cds") or []
    for h in syn.get("hits") or []:
        ci = h.get("cds_index", -1)
        cd = cds[ci] if isinstance(ci, int) and 0 <= ci < len(cds) else {}
        ref_nt = h.get("ref_nt") or [None, None]
        out["synteny"].append((syn.get("reference"), int(bool(syn.get("flipped"))), h.get("orf_name"),
                               cd.get("product"), cd.get("protein_id"), _num(h.get("pident")),
                               _num(h.get("evalue")), _num(h.get("bit_score")), h.get("q_start"), h.get("q_end"),
                               h.get("s_start"), h.get("s_end"), _num(h.get("q_cov")), _num(h.get("s_cov")),
                               ref_nt[0] + 1 if ref_nt[0] is not None else None, ref_nt[1]))
    return out


def _salmon_rows(sq: dict | None) -> dict[str, list[tuple]]:
    """Every Salmon target and host-viral hit of one sample, without key columns."""
    out: dict[str, list[tuple]] = {"quantification": [], "host_viral": []}
    if not sq:
        return out
    entries = (list(sq.get("viral_quant") or []) + list(sq.get("ref_hk_quant") or [])
               + list(sq.get("conserved_quant") or []) + list(sq.get("pfam_hk_quant") or [])
               + [r.get("transcript") or {} for r in sq.get("host_viral_hits") or []])
    for e in entries:
        name = e.get("name") or ""
        out["quantification"].append((
            name, e.get("seq_type"), e.get("kingdom") or None,
            name.removeprefix("VQ_VIRAL_") if name.startswith("VQ_VIRAL_") else None,
            e.get("length"), _num(e.get("eff_length")), _num(e.get("tpm")), _num(e.get("num_reads"))))
    for r in sq.get("host_viral_hits") or []:
        t = r.get("transcript") or {}
        for h in r.get("blast_hits") or []:
            out["host_viral"].append((h.get("host_transcript_id"), h.get("viral_seq_id"), _num(h.get("pident")),
                                      _num(h.get("qcovhsp")), _num(h.get("evalue")), _num(h.get("bit_score")),
                                      _num(t.get("tpm")), _num(t.get("num_reads"))))
    return out


def write_sqlite(report: dict, path: Path, sample_reports=None) -> Path:
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

        # Per-sample data the consolidated report does not ship: synteny references
        # and the full Salmon tables (from the loaded SampleReport objects).
        refs_by_sample = {s.sample: ((s.report.get("synteny") or {}).get("references") or {})
                          for s in sample_reports or []}
        extra: dict[str, list[tuple]] = {"coverage": [], "seq_quality": [], "synteny": [],
                                         "quantification": [], "host_viral": []}
        for smp in sample_reports or []:
            for k, rs in _salmon_rows(smp.salmon_quant).items():
                extra[k] += [(ds, smp.sample) + r for r in rs]

        seq_rows, bx_rows, bn_rows, bnt_rows, orf_rows, dom_rows = [], [], [], [], [], []
        for q in report.get("sequences", []) or []:
            gid, sample = q.get("gid"), q.get("sample")
            reads, inp = q.get("reads"), input_reads.get(sample)
            seq_rows.append((ds, gid, sample, q.get("id"))
                            + _sequence_fields(q, q.get("tpm"), reads, inp))
            rows = _child_rows(q)
            key = (ds, gid)
            bx_rows  += [key + r for r in rows["blastx"]]
            bn_rows  += [key + r for r in rows["blastn"]]
            bnt_rows += [key + r for r in rows["blastn_taxonomy"]]
            orf_rows += [key + r for r in rows["orfs"]]
            dom_rows += [key + r for r in rows["domains"]]
            for k, rs in _analysis_rows(q, refs_by_sample.get(sample, {})).items():
                extra[k] += [key + r for r in rs]

        con.executemany(f"INSERT INTO sequences VALUES ({','.join('?' * 29)})", seq_rows)
        con.executemany("INSERT OR IGNORE INTO blastx_hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", bx_rows)
        con.executemany("INSERT OR IGNORE INTO blastn_hits VALUES (?,?,?,?,?,?,?,?,?,?,?)", bn_rows)
        con.executemany("INSERT INTO blastn_taxonomy VALUES (?,?,?,?,?,?,?,?,?,?,?)", bnt_rows)
        con.executemany("INSERT OR IGNORE INTO orfs VALUES (?,?,?,?,?,?,?,?,?)", orf_rows)
        con.executemany("INSERT INTO domains VALUES (?,?,?,?,?,?,?,?,?,?)", dom_rows)
        _insert_extra(con, extra, n_key=2, n_salmon_key=2)

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


def _insert_extra(con, extra: dict[str, list[tuple]], n_key: int, n_salmon_key: int) -> None:
    """Insert the coverage / quality / synteny / Salmon rows. *n_key*: key columns of a
    sequence row (dataset + gid, or seq_id); *n_salmon_key*: of a Salmon row
    (dataset + sample, or none)."""
    ph = lambda n: ",".join("?" * n)
    con.executemany(f"INSERT OR IGNORE INTO coverage VALUES ({ph(n_key + 9)})", extra["coverage"])
    con.executemany(f"INSERT OR IGNORE INTO seq_quality VALUES ({ph(n_key + 9)})", extra["seq_quality"])
    con.executemany(f"INSERT INTO synteny_hits VALUES ({ph(n_key + 16)})", extra["synteny"])
    con.executemany(f"INSERT OR IGNORE INTO quantification VALUES ({ph(n_salmon_key + 8)})", extra["quantification"])
    con.executemany(f"INSERT INTO host_viral_hits VALUES ({ph(n_salmon_key + 8)})", extra["host_viral"])


# ── Single-sample SQLite ────────────────────────────────────────────────────
#
# The tables of the consolidated report for one ``viralquest`` run: no
# dataset / sample tables or columns, sequences keyed by their own id.

SAMPLE_SCHEMA = """
CREATE TABLE sequences (
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
    PRIMARY KEY (seq_id)
);
CREATE TABLE blastx_hits (
    seq_id             TEXT NOT NULL,
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
    PRIMARY KEY (seq_id, source, rank)
);
CREATE TABLE blastn_hits (
    seq_id             TEXT NOT NULL,
    rank               INTEGER NOT NULL,
    accession          TEXT,
    title              TEXT,
    pct_identity       REAL,
    query_coverage     REAL,
    e_value            REAL,
    bit_score          REAL,
    query_start        INTEGER,
    query_end          INTEGER,
    PRIMARY KEY (seq_id, rank)
);
CREATE TABLE blastn_taxonomy (
    seq_id             TEXT NOT NULL,
    status             TEXT,
    method             TEXT,
    accession          TEXT,
    tax_id             INTEGER,
    matched_name       TEXT,
    synonym_of         TEXT,
    family             TEXT,
    genus              TEXT,
    species            TEXT,
    PRIMARY KEY (seq_id)
);
CREATE TABLE orfs (
    seq_id             TEXT NOT NULL,
    orf                TEXT NOT NULL,
    start              INTEGER,
    stop               INTEGER,
    strand             TEXT,
    frame              INTEGER,
    orf_type           TEXT,
    length_aa          INTEGER,
    PRIMARY KEY (seq_id, orf)
);
CREATE TABLE domains (
    seq_id             TEXT NOT NULL,
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
    cluster_id         TEXT NOT NULL,     -- VQ_CLU_####
    species            TEXT,
    species_db         TEXT,
    representative     TEXT,              -- seq_id
    n_members          INTEGER,
    homogeneous        INTEGER,
    PRIMARY KEY (cluster_id)
);
CREATE TABLE cluster_members (
    cluster_id         TEXT NOT NULL,
    seq_id             TEXT NOT NULL,
    identity           REAL,
    coverage           REAL,
    is_representative  INTEGER,
    tpm                REAL,
    reads              REAL,
    PRIMARY KEY (cluster_id, seq_id)
);
CREATE TABLE coverage (
    seq_id             TEXT NOT NULL,
    length             INTEGER,
    mean_depth         REAL,
    max_depth          REAL,
    cv                 REAL,
    breadth_1x         REAL,              -- fraction of bases with depth >= 1
    mean_depth_sense   REAL,
    mean_depth_antisense REAL,
    mean_quality       REAL,
    n_low_cov_regions  INTEGER,
    PRIMARY KEY (seq_id)
);
CREATE TABLE seq_quality (
    seq_id             TEXT NOT NULL,
    length             INTEGER,
    low_complexity_frac REAL,
    n_low_complexity   INTEGER,
    n_direct_repeats   INTEGER,
    n_inverted_repeats INTEGER,
    kmer_size          INTEGER,
    kmer_distinct      INTEGER,
    kmer_total         INTEGER,
    kmer_repeat_score  REAL,
    PRIMARY KEY (seq_id)
);
CREATE TABLE synteny_hits (
    seq_id             TEXT NOT NULL,
    reference          TEXT,              -- accession of the best viral BLASTn hit
    flipped            INTEGER,           -- 1 = contig is the reverse complement
    orf                TEXT,
    cds_product        TEXT,
    protein_id         TEXT,
    pct_identity       REAL,
    e_value            REAL,
    bit_score          REAL,
    orf_aa_start       INTEGER,
    orf_aa_end         INTEGER,
    cds_aa_start       INTEGER,
    cds_aa_end         INTEGER,
    orf_coverage       REAL,
    cds_coverage       REAL,
    ref_start          INTEGER,           -- 1-based, on the reference genome
    ref_end            INTEGER
);
CREATE TABLE quantification (
    target             TEXT NOT NULL,     -- Salmon target name (VQ_VIRAL_…, VQ_REFHK_…, VQ_CONS_…, host)
    seq_type           TEXT,              -- viral | ref_hk | conserved | pfam_hk | host
    kingdom            TEXT,
    seq_id             TEXT,              -- viral targets only
    length             INTEGER,
    eff_length         REAL,
    tpm                REAL,
    reads              REAL,
    PRIMARY KEY (target)
);
CREATE TABLE host_viral_hits (
    host_transcript    TEXT,
    seq_id             TEXT,
    pct_identity       REAL,
    coverage           REAL,
    e_value            REAL,
    bit_score          REAL,
    transcript_tpm     REAL,
    transcript_reads   REAL
);
CREATE INDEX idx_sequences_family  ON sequences (family);
CREATE INDEX idx_sequences_species ON sequences (species);
CREATE INDEX idx_domains_seq       ON domains (seq_id);
CREATE INDEX idx_members_seq       ON cluster_members (seq_id);
"""


def write_sample_sqlite(report: dict, path: Path) -> Path:
    """
    One ``viralquest`` run as ``<stem>_viralquest.sqlite``: the consolidated
    report's sequences / hits / ORFs / domains / clusters tables for a single
    sample (no dataset or sample columns; rows keyed by ``seq_id``).
    """
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    sq = report.get("salmon_quant") or {}
    quant = {(e.get("name") or "").removeprefix("VQ_VIRAL_"): (e.get("tpm"), e.get("num_reads"))
             for e in sq.get("viral_quant") or []}
    # Same RPM denominator as the consolidated report: input reads = mapped / mapping rate.
    from viralquest.report_data import _input_reads
    input_reads = _input_reads(sq)
    con = sqlite3.connect(tmp)
    try:
        con.executescript(SAMPLE_SCHEMA)
        con.execute("PRAGMA user_version = %d" % SCHEMA_VERSION)
        rows = {"blastx": [], "blastn": [], "blastn_taxonomy": [], "orfs": [], "domains": []}
        extra = {"coverage": [], "seq_quality": [], "synteny": [], **_salmon_rows(sq)}
        refs = (report.get("synteny") or {}).get("references") or {}
        seq_rows = []
        for q in report.get("sequences", []) or []:
            sid = q.get("id")
            tpm, reads = quant.get(sid, (None, None))
            seq_rows.append((sid,) + _sequence_fields(q, tpm, reads, input_reads))
            for k, rs in _child_rows(q).items():
                rows[k] += [(sid,) + r for r in rs]
            for k, rs in _analysis_rows(q, refs).items():
                extra[k] += [(sid,) + r for r in rs]
        con.executemany(f"INSERT OR IGNORE INTO sequences VALUES ({','.join('?' * 26)})", seq_rows)
        con.executemany("INSERT OR IGNORE INTO blastx_hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows["blastx"])
        con.executemany("INSERT OR IGNORE INTO blastn_hits VALUES (?,?,?,?,?,?,?,?,?,?)", rows["blastn"])
        con.executemany("INSERT OR IGNORE INTO blastn_taxonomy VALUES (?,?,?,?,?,?,?,?,?,?)",
                        rows["blastn_taxonomy"])
        con.executemany("INSERT OR IGNORE INTO orfs VALUES (?,?,?,?,?,?,?,?)", rows["orfs"])
        con.executemany("INSERT INTO domains VALUES (?,?,?,?,?,?,?,?,?)", rows["domains"])
        _insert_extra(con, extra, n_key=1, n_salmon_key=0)
        for c in report.get("clusters", []) or []:
            members = c.get("members", []) or []
            homogeneous = len({(m.get("species") or c.get("species")) for m in members}) <= 1
            con.execute("INSERT OR IGNORE INTO clusters VALUES (?,?,?,?,?,?)",
                        (c.get("cluster_id"), c.get("species"), None, c.get("representative_id"),
                         len(members), int(homogeneous)))
            con.executemany(
                "INSERT OR IGNORE INTO cluster_members VALUES (?,?,?,?,?,?,?)",
                [(c.get("cluster_id"), m.get("seq_id"), _num(m.get("identity")),
                  _num(m.get("query_coverage")), int(bool(m.get("is_representative"))),
                  _num(quant.get(m.get("seq_id"), (None, None))[0]),
                  _num(quant.get(m.get("seq_id"), (None, None))[1])) for m in members])
        con.commit()
    finally:
        con.close()
    tmp.replace(path)
    return path
