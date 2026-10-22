import csv
from pathlib import Path

from loguru import logger

from .biodata import NOVELTY_TIERS, BlastnResult, NucSequence


class OutputOrganizer:
    """
    Manages the structured output directory layout:

        outdir/
        ├── features/                ← every alignment / analysis table
        │   ├── diamond/             refseq.tsv (merged batches) · nr.tsv (--nr-db)
        │   ├── hmm/                 RVDB.tsv · Vfam.tsv · EggNOG.tsv · Pfam.tsv
        │   ├── blastn/              blastn.tsv
        │   ├── synteny/             synteny.tsv · refs/<accession>.gb
        │   ├── seq_quality/         <seq_id>.tsv · summary.tsv
        │   └── cap3/                (--cap3)
        ├── reads/                   ← read-based data (--reads)
        │   ├── coverage/            <seq_id>.tsv per-bin depth
        │   ├── mapping/             BAM of the reads on the viral sequences + mapped FASTQ
        │   ├── salmon/              Salmon index / quant output
        │   ├── coverage_summary.tsv
        │   └── quantification.tsv
        ├── fastas/
        │   ├── {stem}_viral_contigs.fasta      every confirmed viral sequence
        │   └── {stem}_<novelty tier>.fasta     one per tier present
        ├── {stem}_viralquest.json
        ├── {stem}_viralquest.html
        ├── {stem}_viralquest.xlsx   ← aggregated results, one sheet per table
        ├── {stem}_viralquest.sqlite ← the consolidated report's tables, one sample
        └── viralquest.log
    """

    _HMM_HEADERS    = ("hmm_profile", "orf_id", "score", "i_evalue", "env_from", "env_to")
    _BLASTN_HEADERS = ("qseqid", "qlen", "slen", "qcovhsp", "pident", "evalue", "bit_score", "stitle")

    def __init__(self, outdir: Path, stem: str) -> None:
        self.outdir          = outdir
        self.stem            = stem
        self.features_dir    = outdir / "features"
        self.diamond_dir     = self.features_dir / "diamond"
        self.blastn_dir      = self.features_dir / "blastn"
        self.hmm_dir         = self.features_dir / "hmm"
        self.synteny_dir     = self.features_dir / "synteny"
        self.seq_quality_dir = self.features_dir / "seq_quality"
        self.cap3_dir        = self.features_dir / "cap3"
        self.reads_dir       = outdir / "reads"
        self.coverage_dir    = self.reads_dir / "coverage"
        self.mapping_dir     = self.reads_dir / "mapping"
        self.salmon_dir      = self.reads_dir / "salmon"
        self.fastas_dir      = outdir / "fastas"
        # blastn/, synteny/, reads/ … appear only when their step runs.
        for d in (self.diamond_dir, self.hmm_dir, self.fastas_dir):
            d.mkdir(parents=True, exist_ok=True)

    def report_path(self, ext: str) -> Path:
        """Top-level report file: {stem}_viralquest.<ext> (json, html, xlsx, sqlite)."""
        return self.outdir / f"{self.stem}_viralquest.{ext}"

    # ── Diamond ───────────────────────────────────────────────────────────────

    def finalize_diamond_refseq(self, tsv_paths: list[Path]) -> Path:
        """Merge per-batch TSVs into diamond/refseq.tsv, then remove the batch files."""
        out = self.diamond_dir / "refseq.tsv"
        with open(out, "w", encoding="utf-8") as fout:
            for p in tsv_paths:
                if p.exists():
                    if p.stat().st_size > 0:
                        fout.write(p.read_text(encoding="utf-8"))
                    p.unlink()
        logger.info(f"Diamond RefSeq results → {out}")
        return out

    def finalize_diamond_nr(self, tsv_path: Path | None) -> Path | None:
        """Move the NR TSV to diamond/nr.tsv."""
        if tsv_path is None or not tsv_path.exists():
            return None
        out = self.diamond_dir / "nr.tsv"
        tsv_path.rename(out)
        logger.info(f"Diamond NR results → {out}")
        return out

    # ── HMM ───────────────────────────────────────────────────────────────────

    def save_hmm_table(self, hits: list[tuple], db_name: str) -> Path:
        """Write raw HMM hit tuples to hmm/{db_name}.tsv."""
        out = self.hmm_dir / f"{db_name}.tsv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(self._HMM_HEADERS)
            writer.writerows(hits)
        logger.info(f"HMM {db_name} ({len(hits)} hits) → {out}")
        return out

    # ── BLASTn ───────────────────────────────────────────────────────────────

    def save_blastn_table(self, hits: list[BlastnResult]) -> Path:
        """
        Write all BLASTn hits to blastn/blastn.tsv and remove any loose
        blastn_batch*.tsv files left in blastn_dir by BlastnRunner.
        """
        self.blastn_dir.mkdir(parents=True, exist_ok=True)
        out = self.blastn_dir / "blastn.tsv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(self._BLASTN_HEADERS)
            for h in hits:
                writer.writerow([
                    h.qseqid, h.qlen, h.slen, h.qcovhsp,
                    h.pident, h.evalue, h.bit_score, h.stitle,
                ])
        # clean up batch files that BlastnRunner left behind
        for batch_file in self.blastn_dir.glob("blastn_batch*.tsv"):
            batch_file.unlink(missing_ok=True)
        logger.info(f"BLASTn results ({len(hits)} hits) → {out}")
        return out

    # ── Viral contigs FASTA ───────────────────────────────────────────────────

    def save_viral_contigs(self, seqs: list[NucSequence]) -> Path:
        """Write confirmed viral sequences to fastas/{stem}_viral_contigs.fasta."""
        out = self.fastas_dir / f"{self.stem}_viral_contigs.fasta"
        self._write_fasta(seqs, out)
        logger.info(f"Viral contigs FASTA ({len(seqs)} seqs) → {out}")
        return out

    def save_novelty_fastas(self, seqs: list[NucSequence]) -> dict[str, Path]:
        """
        One FASTA per heuristic novelty tier (fastas/{stem}_<tier>.fasta), for
        the tiers present; sequences without a heuristic score go to
        {stem}_unscored.fasta. Returns {tier: path}.
        """
        by_tier: dict[str, list[NucSequence]] = {}
        for seq in seqs:
            h = seq.heuristic_output
            tier = (h.novelty if h and h.novelty else "unscored")
            by_tier.setdefault(tier, []).append(seq)
        order = [t for t in NOVELTY_TIERS if t in by_tier] + sorted(set(by_tier) - set(NOVELTY_TIERS))
        out: dict[str, Path] = {}
        for tier in order:
            path = self.fastas_dir / f"{self.stem}_{tier}.fasta"
            self._write_fasta(by_tier[tier], path)
            out[tier] = path
        if out:
            logger.info("Novelty FASTAs → " + ", ".join(f"{t} ({len(by_tier[t])})" for t in order))
        return out

    @staticmethod
    def _write_fasta(seqs: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in seqs:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    # ── Read tables ───────────────────────────────────────────────────────────

    def save_read_tables(self, report: dict) -> list[Path]:
        """
        reads/coverage_summary.tsv (per-sequence depth statistics) and
        reads/quantification.tsv (every Salmon target), from the exported
        report; each only when its data exists.
        """
        from .seq_table import _coverage, _host_viral, _quantification, write_sheet_tsv
        written = []
        for sheet, name in ((_coverage(report), "coverage_summary.tsv"),
                            (_quantification(report), "quantification.tsv"),
                            (_host_viral(report), "host_viral_hits.tsv")):
            if not sheet.rows:
                continue
            self.reads_dir.mkdir(parents=True, exist_ok=True)
            written.append(write_sheet_tsv(sheet, self.reads_dir / name))
        if written:
            logger.info("Read tables → " + ", ".join(p.name for p in written))
        return written

    # ── Log file ──────────────────────────────────────────────────────────────

    def add_log_sink(self) -> int:
        """Add a loguru file sink for realtime logging; returns the sink id."""
        log_path = self.outdir / "viralquest.log"
        return logger.add(
            str(log_path),
            format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
            rotation=None,
            encoding="utf-8",
        )

    @staticmethod
    def remove_log_sink(sink_id: int) -> None:
        logger.remove(sink_id)
