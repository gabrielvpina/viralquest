import csv
from pathlib import Path

from loguru import logger

from .biodata import BlastnResult, NucSequence


class OutputOrganizer:
    """
    Manages the structured output directory layout:

        outdir/
        ├── diamond/
        │   ├── refseq.tsv      ← merged from all RefSeq batches
        │   └── nr.tsv          ← single NR run (when --nr-db is used)
        ├── blastn/
        │   └── blastn.tsv      ← all BLASTn hits merged
        ├── hmm/
        │   ├── RVDB.tsv
        │   ├── Vfam.tsv
        │   ├── EggNOG.tsv
        │   └── Pfam.tsv
        ├── {stem}_viralquest.json
        ├── {stem}_viralquest.html
        ├── {stem}_viral_contigs.fasta
        └── viralquest.log
    """

    _HMM_HEADERS    = ("hmm_profile", "orf_id", "score", "i_evalue", "env_from", "env_to")
    _BLASTN_HEADERS = ("qseqid", "qlen", "slen", "qcovhsp", "pident", "evalue", "bit_score", "stitle")

    def __init__(self, outdir: Path, stem: str) -> None:
        self.outdir      = outdir
        self.stem        = stem
        self.diamond_dir = outdir / "diamond"
        self.blastn_dir  = outdir / "blastn"
        self.hmm_dir     = outdir / "hmm"
        self.diamond_dir.mkdir(parents=True, exist_ok=True)
        self.blastn_dir.mkdir(parents=True, exist_ok=True)
        self.hmm_dir.mkdir(parents=True, exist_ok=True)

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
        """Write confirmed viral sequences to {stem}_viral_contigs.fasta."""
        out   = self.outdir / f"{self.stem}_viral_contigs.fasta"
        viral = [s for s in seqs if s.is_viral]
        with open(out, "w", encoding="utf-8") as fh:
            for seq in viral:
                fh.write(f">{seq.id}\n{seq.sequence}\n")
        logger.info(f"Viral contigs FASTA ({len(viral)} seqs) → {out}")
        return out

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
