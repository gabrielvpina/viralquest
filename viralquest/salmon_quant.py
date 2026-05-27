import re
import shutil
import subprocess
from pathlib import Path

from loguru import logger

from viralquest.biodata import (
    HostViralHit,
    HostViralRecord,
    NucSequence,
    SalmonEntry,
    SalmonQuantReport,
)


# ---------------------------------------------------------------------------
# Kingdom FASTA paths
# ---------------------------------------------------------------------------

_KINGDOMS_DIR = Path(__file__).parent.parent / "data" / "salmon-quant"

_KINGDOM_FILES: dict[str, Path] = {
    name: _KINGDOMS_DIR / f"{name}.fasta"
    for name in (
        "mammals", "arthropods", "plants",
        "fish", "fungi", "bacteria", "nematodes",
    )
}


# ---------------------------------------------------------------------------
# ConservedFastaLoader
# ---------------------------------------------------------------------------

class ConservedFastaLoader:
    """
    Reads all kingdom FASTA files and prefixes each sequence ID with
    ``VQ_CONS_{KINGDOM}_`` so kingdom origin is trackable in quant.sf.

    Returns
    -------
    dict mapping upper-case kingdom name to a list of
    ``(prefixed_id, description, sequence)`` tuples.
    """

    @classmethod
    def load_all(cls) -> dict[str, list[tuple[str, str, str]]]:
        result: dict[str, list[tuple[str, str, str]]] = {}
        for kingdom, path in _KINGDOM_FILES.items():
            if not path.exists():
                logger.warning(f"Conserved FASTA not found, skipping: '{path}'")
                continue
            entries = cls._parse_fasta(path, kingdom.upper())
            result[kingdom.upper()] = entries
            logger.debug(f"Conserved [{kingdom.upper()}]: {len(entries)} sequence(s) loaded.")
        return result

    @staticmethod
    def _parse_fasta(path: Path, kingdom: str) -> list[tuple[str, str, str]]:
        entries:    list[tuple[str, str, str]] = []
        seq_id      = ""
        desc        = ""
        seq_lines:  list[str] = []

        def _flush() -> None:
            if seq_id:
                entries.append((f"VQ_CONS_{kingdom}_{seq_id}", desc, "".join(seq_lines)))

        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    _flush()
                    parts  = line[1:].split(None, 1)
                    seq_id = parts[0]
                    desc   = parts[1] if len(parts) > 1 else ""
                    seq_lines = []
                else:
                    seq_lines.append(line)
        _flush()
        return entries


# ---------------------------------------------------------------------------
# CombinedFastaWriter
# ---------------------------------------------------------------------------

class CombinedFastaWriter:
    """
    Writes a single merged reference FASTA for Salmon indexing:

      1. Viral sequences (prefix ``VQ_VIRAL_``)
      2. Conserved kingdom genes (prefix ``VQ_CONS_{KINGDOM}_``)
      3. User transcriptome (original IDs, unchanged)

    The prefix scheme lets the quant.sf parser split entries by origin
    without a lookup table.

    Returns
    -------
    tuple of:
      - ``set[str]``          viral prefixed IDs
      - ``dict[str, str]``    {conserved prefixed ID: kingdom}
    """

    @staticmethod
    def write(
        output_path:        Path,
        user_transcriptome: Path,
        viral_seqs:         list[NucSequence],
        conserved:          dict[str, list[tuple[str, str, str]]],
    ) -> tuple[set[str], dict[str, str]]:
        viral_ids:    set[str]       = set()
        conserved_map: dict[str, str] = {}     # {prefixed_id: kingdom}

        with open(output_path, "w", encoding="utf-8") as fh:

            # 1 — viral sequences
            for seq in viral_seqs:
                prefixed = f"VQ_VIRAL_{seq.id}"
                fh.write(f">{prefixed}\n{seq.sequence}\n")
                viral_ids.add(prefixed)

            # 2 — conserved genes (all kingdoms, with kingdom prefix)
            for kingdom, entries in conserved.items():
                for prefixed_id, desc, sequence in entries:
                    header = f">{prefixed_id}" + (f" {desc}" if desc else "")
                    fh.write(f"{header}\n{sequence}\n")
                    conserved_map[prefixed_id] = kingdom

            # 3 — user transcriptome (streamed line-by-line to avoid memory load)
            with open(user_transcriptome, encoding="utf-8") as txfh:
                for line in txfh:
                    fh.write(line if line.endswith("\n") else line + "\n")

        logger.info(
            f"Combined FASTA: {len(viral_ids)} viral + "
            f"{len(conserved_map)} conserved + host transcriptome "
            f"→ '{output_path.name}'"
        )
        return viral_ids, conserved_map


# ---------------------------------------------------------------------------
# SalmonIndexBuilder
# ---------------------------------------------------------------------------

class SalmonIndexBuilder:
    """Builds a Salmon quasi-mapping index from the combined reference FASTA."""

    def __init__(self, salmon_bin: str = "salmon", threads: int = 4):
        self.salmon_bin = salmon_bin
        self.threads    = threads

    def build(self, ref_fasta: Path, index_dir: Path) -> None:
        index_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.salmon_bin, "index",
            "-t", str(ref_fasta),
            "-i", str(index_dir),
            "-p", str(self.threads),
        ]
        logger.info(f"Salmon index: '{ref_fasta.name}' → '{index_dir.name}' ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            # show only the last 20 lines of stderr — the crash is always at the end
            tail = "\n".join(result.stderr.strip().splitlines()[-20:])
            raise RuntimeError(
                f"salmon index failed (exit {result.returncode}):\n{tail}\n\n"
                "Common causes: insufficient RAM for SSHash (try closing other "
                "processes), or a corrupted/empty combined reference FASTA."
            )
        logger.success(f"Salmon index built → '{index_dir}'")


# ---------------------------------------------------------------------------
# SalmonQuantRunner
# ---------------------------------------------------------------------------

class SalmonQuantRunner:
    """
    Runs ``salmon quant`` in mapping-based mode.

    Supports single-end (one reads file) and paired-end (two files).
    Library type is auto-detected (``-l A``).
    Mapping rate is extracted from salmon's stderr log.
    """

    _RATE_RE = re.compile(r"Mapping rate\s*=\s*([\d.]+)%", re.IGNORECASE)

    def __init__(self, salmon_bin: str = "salmon", threads: int = 4):
        self.salmon_bin = salmon_bin
        self.threads    = threads

    def run(
        self,
        index_dir: Path,
        reads:     list[str],
        out_dir:   Path,
    ) -> tuple[Path, float]:
        """
        Run salmon quant and return ``(quant_sf_path, mapping_rate_percent)``.
        Raises ``ValueError`` for wrong number of reads files.
        Raises ``RuntimeError`` on non-zero exit.
        """
        if len(reads) not in (1, 2):
            raise ValueError(f"Expected 1 or 2 reads files, got {len(reads)}.")
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.salmon_bin, "quant",
            "-i",  str(index_dir),
            "-l",  "A",
            "-p",  str(self.threads),
            "--validateMappings",
            "-o",  str(out_dir),
        ]
        if len(reads) == 1:
            cmd += ["-r", reads[0]]
        else:
            cmd += ["-1", reads[0], "-2", reads[1]]

        mode = "paired" if len(reads) == 2 else "single"
        logger.info(f"Salmon quant [{mode}-end, {self.threads} threads] → '{out_dir.name}' ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"salmon quant failed:\n{result.stderr.strip()}")

        mapping_rate = 0.0
        m = self._RATE_RE.search(result.stderr)
        if m:
            mapping_rate = float(m.group(1))

        quant_sf = out_dir / "quant.sf"
        logger.success(
            f"Salmon quant done — mapping rate {mapping_rate:.1f}% → '{quant_sf}'"
        )
        return quant_sf, mapping_rate


# ---------------------------------------------------------------------------
# QuantsfParser
# ---------------------------------------------------------------------------

class QuantsfParser:
    """
    Parses ``quant.sf`` (tab-separated, one header row) into
    ``SalmonEntry`` objects tagged by sequence origin.

    Tagging rules (based on the ID prefix written by CombinedFastaWriter):
      ``VQ_VIRAL_*``        → seq_type = "viral",     kingdom = ""
      ``VQ_CONS_{KG}_*``   → seq_type = "conserved",  kingdom = "{KG}"
      anything else        → seq_type = "host",       kingdom = ""
    """

    @staticmethod
    def parse(
        quant_sf:      Path,
        conserved_map: dict[str, str],   # {prefixed_id: kingdom}
    ) -> list[SalmonEntry]:
        entries: list[SalmonEntry] = []
        if not quant_sf.exists():
            logger.error(f"quant.sf not found: '{quant_sf}'")
            return entries

        with open(quant_sf, encoding="utf-8") as fh:
            next(fh)    # skip header: Name Length EffectiveLength TPM NumReads
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 5:
                    continue
                try:
                    name = parts[0]
                    if name.startswith("VQ_VIRAL_"):
                        seq_type = "viral"
                        kingdom  = ""
                    elif name.startswith("VQ_CONS_"):
                        seq_type = "conserved"
                        kingdom  = conserved_map.get(name, "")
                        if not kingdom:
                            # fallback: VQ_CONS_MAMMALS_accession → parts[2]
                            tokens  = name.split("_")
                            kingdom = tokens[2] if len(tokens) > 2 else ""
                    else:
                        seq_type = "host"
                        kingdom  = ""

                    entries.append(SalmonEntry(
                        name       = name,
                        length     = int(parts[1]),
                        eff_length = float(parts[2]),
                        tpm        = float(parts[3]),
                        num_reads  = float(parts[4]),
                        seq_type   = seq_type,
                        kingdom    = kingdom,
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(f"Skipping malformed quant.sf row: {exc}")

        viral_n     = sum(1 for e in entries if e.seq_type == "viral")
        conserved_n = sum(1 for e in entries if e.seq_type == "conserved")
        host_n      = sum(1 for e in entries if e.seq_type == "host")
        logger.info(
            f"quant.sf: {viral_n} viral | {conserved_n} conserved | {host_n} host entries."
        )
        return entries


# ---------------------------------------------------------------------------
# TranscriptomeViralAligner
# ---------------------------------------------------------------------------

class TranscriptomeViralAligner:
    """
    BLASTn of the user transcriptome against the viralquest viral sequences.

    Identifies host transcripts with nucleotide similarity to detected viral
    sequences — endogenous viral elements (EVEs), chimeric contigs, etc.

    The viral sequences form a small database (typically 5–100 sequences),
    so ``blastn -subject`` is used to avoid a ``makeblastdb`` dependency.

    Parameters
    ----------
    blastn_bin  : blastn executable name or full path
    e_value     : e-value cutoff
    threads     : ``-num_threads`` passed to blastn
    min_pident  : minimum % identity (applied in command via ``-perc_identity``)
    min_qcov    : minimum query coverage % (applied post-hoc)
    """

    _OUTFMT = "6 qseqid sseqid pident length qlen qcovhsp evalue bitscore"

    def __init__(
        self,
        blastn_bin:  str   = "blastn",
        e_value:     float = 1e-5,
        threads:     int   = 4,
        min_pident:  float = 70.0,
        min_qcov:    int   = 20,
    ):
        self.blastn_bin = blastn_bin
        self.e_value    = e_value
        self.threads    = threads
        self.min_pident = min_pident
        self.min_qcov   = min_qcov

    def align(
        self,
        transcriptome: Path,
        viral_seqs:    list[NucSequence],
        outdir:        Path,
    ) -> list[HostViralHit]:
        if not viral_seqs:
            logger.warning("TranscriptomeViralAligner: no viral sequences — skipping BLASTn.")
            return []

        viral_fasta = outdir / "vq_viral_subject.fa"
        self._write_fasta(viral_seqs, viral_fasta)

        out_tsv = outdir / "host_vs_viral.tsv"
        cmd = [
            self.blastn_bin,
            "-query",         str(transcriptome),
            "-subject",       str(viral_fasta),
            "-out",           str(out_tsv),
            "-outfmt",        self._OUTFMT,
            "-evalue",        str(self.e_value),
            "-num_threads",   str(self.threads),
            "-perc_identity", str(self.min_pident),
        ]

        logger.info(
            f"BLASTn [host→viral]: transcriptome vs {len(viral_seqs)} viral seq(s) ..."
        )
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"BLASTn [host→viral] failed: {result.stderr.strip()}")
            return []

        hits = self._parse_tsv(out_tsv)
        hits = [h for h in hits if h.qcovhsp >= self.min_qcov]
        logger.success(f"BLASTn [host→viral]: {len(hits)} hit(s) after filtering.")
        return hits

    @staticmethod
    def _write_fasta(seqs: list[NucSequence], path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for seq in seqs:
                fh.write(f">{seq.id}\n{seq.sequence}\n")

    def _parse_tsv(self, path: Path) -> list[HostViralHit]:
        hits: list[HostViralHit] = []
        if not path.exists() or path.stat().st_size == 0:
            return hits
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 8:
                    continue
                try:
                    hits.append(HostViralHit(
                        host_transcript_id = parts[0],
                        viral_seq_id       = parts[1],
                        pident             = float(parts[2]),
                        qcovhsp            = int(float(parts[5])),
                        evalue             = float(parts[6]),
                        bit_score          = float(parts[7]),
                    ))
                except (ValueError, IndexError) as exc:
                    logger.warning(f"Skipping malformed BLASTn row: {exc}")
        return hits


# ---------------------------------------------------------------------------
# SalmonQuantPipeline
# ---------------------------------------------------------------------------

class SalmonQuantPipeline:
    """
    Orchestrates the full Salmon quantification workflow.

    Steps
    -----
    1. Load all conserved kingdom FASTA files (VQ_CONS_{KINGDOM}_ prefix).
    2. Write combined reference FASTA (viral + conserved + user transcriptome).
    3. Build Salmon index from the combined FASTA.
    4. Run Salmon quant against user RNA-seq reads.
    5. Parse ``quant.sf`` into ``SalmonEntry`` objects tagged by origin.
    6. BLASTn: user transcriptome vs viral sequences → ``HostViralHit`` objects.
    7. Assemble and return a ``SalmonQuantReport``.
    8. Clean up the Salmon index and temporary files (if ``cleanup=True``).

    The quantification output directory (``outdir/vq_quant/``) is always
    preserved so the user can inspect the raw quant.sf and salmon logs.

    Parameters
    ----------
    salmon_bin  : salmon executable
    blastn_bin  : blastn executable
    threads     : threads for index building, quant, and BLASTn
    cleanup     : remove the Salmon index and temp files after the run
    e_value     : BLASTn e-value cutoff for host→viral alignment
    min_pident  : minimum % identity for host→viral BLASTn hits
    min_qcov    : minimum query coverage % for host→viral BLASTn hits
    """

    def __init__(
        self,
        salmon_bin:  str   = "salmon",
        blastn_bin:  str   = "blastn",
        threads:     int   = 4,
        cleanup:     bool  = True,
        e_value:     float = 1e-5,
        min_pident:  float = 70.0,
        min_qcov:    int   = 20,
    ):
        self._index_builder = SalmonIndexBuilder(salmon_bin, threads)
        self._quant_runner  = SalmonQuantRunner(salmon_bin, threads)
        self._aligner       = TranscriptomeViralAligner(
            blastn_bin, e_value, threads, min_pident, min_qcov,
        )
        self.cleanup = cleanup

    # --- public ---------------------------------------------------------------

    def run(
        self,
        user_transcriptome: Path,
        reads:              list[str],
        viral_seqs:         list[NucSequence],
        outdir:             Path,
    ) -> SalmonQuantReport:
        """
        Execute the full pipeline and return a ``SalmonQuantReport``.

        Parameters
        ----------
        user_transcriptome : path to the host transcriptome FASTA
        reads              : ``[single.fastq]`` or ``[r1.fastq, r2.fastq]``
        viral_seqs         : all ``NucSequence`` objects from the pipeline
                             (unconfirmed ones are silently filtered out)
        outdir             : base directory for all output (created if absent)
        """
        outdir.mkdir(parents=True, exist_ok=True)
        tmp_dir   = outdir / "vq_tmp"
        index_dir = outdir / "vq_salmon_idx"
        quant_dir = outdir / "vq_quant"
        tmp_dir.mkdir(exist_ok=True)

        confirmed = [s for s in viral_seqs if s.is_viral and s.blastx_nr_hits]
        if not confirmed:
            logger.warning(
                "SalmonQuantPipeline: no confirmed viral sequences — "
                "the viral_quant section will be empty."
            )

        try:
            conserved = ConservedFastaLoader.load_all()

            combined_fasta = tmp_dir / "combined_ref.fasta"
            _, conserved_map = CombinedFastaWriter.write(
                output_path        = combined_fasta,
                user_transcriptome = user_transcriptome,
                viral_seqs         = confirmed,
                conserved          = conserved,
            )

            self._index_builder.build(combined_fasta, index_dir)

            quant_sf, mapping_rate = self._quant_runner.run(index_dir, reads, quant_dir)

            all_entries = QuantsfParser.parse(quant_sf, conserved_map)

            blast_hits = self._aligner.align(user_transcriptome, confirmed, tmp_dir)

            report = self._assemble(
                reads        = reads,
                mapping_rate = mapping_rate,
                entries      = all_entries,
                blast_hits   = blast_hits,
            )

        finally:
            if self.cleanup:
                shutil.rmtree(index_dir, ignore_errors=True)
                shutil.rmtree(tmp_dir,   ignore_errors=True)
                logger.debug("SalmonQuantPipeline: temporary index and files removed.")

        return report

    # --- private --------------------------------------------------------------

    @staticmethod
    def _assemble(
        reads:        list[str],
        mapping_rate: float,
        entries:      list[SalmonEntry],
        blast_hits:   list[HostViralHit],
    ) -> SalmonQuantReport:
        viral_quant     = [e for e in entries if e.seq_type == "viral"]
        conserved_quant = [e for e in entries if e.seq_type == "conserved"]
        host_map        = {e.name: e for e in entries if e.seq_type == "host"}

        total_reads = int(sum(e.num_reads for e in entries))

        # group blast hits by host transcript ID
        hits_by_host: dict[str, list[HostViralHit]] = {}
        for hit in blast_hits:
            hits_by_host.setdefault(hit.host_transcript_id, []).append(hit)

        host_viral_hits: list[HostViralRecord] = []
        for host_id, hits in hits_by_host.items():
            hits_sorted = sorted(hits, key=lambda h: h.bit_score, reverse=True)
            entry = host_map.get(host_id, SalmonEntry(
                name       = host_id,
                length     = 0,
                eff_length = 0.0,
                tpm        = 0.0,
                num_reads  = 0.0,
                seq_type   = "host",
                kingdom    = "",
            ))
            host_viral_hits.append(HostViralRecord(
                transcript = entry,
                blast_hits = hits_sorted,
            ))

        # sort records: strongest blast evidence first
        host_viral_hits.sort(
            key=lambda r: r.blast_hits[0].bit_score if r.blast_hits else 0.0,
            reverse=True,
        )

        # sort conserved by kingdom then TPM descending — easier to read
        conserved_quant.sort(key=lambda e: (e.kingdom, -e.tpm))

        logger.info(
            f"SalmonQuantReport assembled: "
            f"{len(viral_quant)} viral | "
            f"{len(conserved_quant)} conserved | "
            f"{len(host_viral_hits)} host-viral record(s)."
        )
        return SalmonQuantReport(
            reads            = reads,
            mapping_rate     = mapping_rate,
            total_reads      = total_reads,
            viral_quant      = viral_quant,
            conserved_quant  = conserved_quant,
            host_viral_hits  = host_viral_hits,
        )
