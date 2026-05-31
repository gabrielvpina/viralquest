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
# Shared FASTA iterator
# ---------------------------------------------------------------------------

def _parse_fasta_iter(path: Path):
    """Yield (seq_id, description, sequence) for each record in a FASTA file."""
    seq_id  = ""
    desc    = ""
    buf: list[str] = []

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id:
                    yield seq_id, desc, "".join(buf)
                parts  = line[1:].split(None, 1)
                seq_id = parts[0]
                desc   = parts[1] if len(parts) > 1 else ""
                buf    = []
            else:
                buf.append(line)
    if seq_id:
        yield seq_id, desc, "".join(buf)


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
            entries = [
                (f"VQ_CONS_{kingdom.upper()}_{seq_id}", desc, seq)
                for seq_id, desc, seq in _parse_fasta_iter(path)
            ]
            result[kingdom.upper()] = entries
            logger.debug(f"Conserved [{kingdom.upper()}]: {len(entries)} sequence(s) loaded.")
        return result


# ---------------------------------------------------------------------------
# ReferenceHkLoader  (reference pathway — pathway 1)
# ---------------------------------------------------------------------------

class ReferenceHkLoader:
    """
    Extracts reference housekeeping gene sequences from the user transcriptome.

    Reads gene IDs from a plain-text file (one ID per line, ``#`` comments
    ignored) and finds matching sequences in the transcriptome FASTA.
    Matching is exact and case-sensitive against the first whitespace-delimited
    token of each header (same rule as the first word of a ``>`` line).

    Returns
    -------
    tuple of:
      - list[tuple[str, str, str]]   (``VQ_REFHK_``-prefixed entries)
      - set[str]                     original IDs (for transcriptome exclusion)
    """

    @classmethod
    def load(
        cls,
        hk_genes_file: Path,
        transcriptome:  Path,
    ) -> tuple[list[tuple[str, str, str]], set[str]]:
        wanted = cls._read_ids(hk_genes_file)
        if not wanted:
            logger.warning(f"ReferenceHkLoader: no IDs found in '{hk_genes_file}'.")
            return [], set()

        entries:   list[tuple[str, str, str]] = []
        found_ids: set[str] = set()

        for seq_id, desc, seq in _parse_fasta_iter(transcriptome):
            if seq_id in wanted:
                entries.append((f"VQ_REFHK_{seq_id}", desc, seq))
                found_ids.add(seq_id)

        missing = wanted - found_ids
        if missing:
            sample = ", ".join(sorted(missing)[:10])
            logger.warning(
                f"ReferenceHkLoader: {len(missing)} ID(s) not found in transcriptome: "
                f"{sample}" + ("…" if len(missing) > 10 else "")
            )

        logger.info(
            f"ReferenceHkLoader: {len(entries)} reference HK gene(s) loaded "
            f"from '{hk_genes_file.name}'."
        )
        return entries, found_ids

    @staticmethod
    def _read_ids(path: Path) -> set[str]:
        ids: set[str] = set()
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    ids.add(line)
        return ids


# ---------------------------------------------------------------------------
# ContigHkBlaster  (de-novo pathway — pathway 2)
# ---------------------------------------------------------------------------

class ContigHkBlaster:
    """
    BLASTs bundled housekeeping genes against assembled contigs to identify
    HK-matching contigs for de-novo pathway normalization.

    Returns ``{contig_id: (kingdom, hk_gene_id)}`` — best hit per contig.
    Kingdom is upper-case (e.g. ``"MAMMALS"``).
    """

    _OUTFMT = "6 qseqid sseqid pident qcovhsp evalue bitscore"

    def __init__(
        self,
        blastn_bin:  str   = "blastn",
        e_value:     float = 1e-5,
        threads:     int   = 4,
        min_pident:  float = 80.0,
        min_qcov:    int   = 50,
    ):
        self.blastn_bin = blastn_bin
        self.e_value    = e_value
        self.threads    = threads
        self.min_pident = min_pident
        self.min_qcov   = min_qcov

    def match(
        self,
        contigs_fasta: Path,
        outdir:        Path,
    ) -> dict[str, tuple[str, str]]:
        outdir.mkdir(parents=True, exist_ok=True)

        hk_query: Path = outdir / "hk_query.fasta"
        gene_kingdom: dict[str, str] = {}
        with open(hk_query, "w", encoding="utf-8") as fh:
            for kingdom, path in _KINGDOM_FILES.items():
                if not path.exists():
                    continue
                for seq_id, desc, seq in _parse_fasta_iter(path):
                    gene_kingdom[seq_id] = kingdom.upper()
                    header = f">{seq_id}" + (f" {desc}" if desc else "")
                    fh.write(f"{header}\n{seq}\n")

        if not gene_kingdom:
            logger.warning("ContigHkBlaster: no HK FASTA files found — skipping.")
            return {}

        out_tsv = outdir / "hk_vs_contigs.tsv"
        cmd = [
            self.blastn_bin,
            "-query",         str(hk_query),
            "-subject",       str(contigs_fasta),
            "-out",           str(out_tsv),
            "-outfmt",        self._OUTFMT,
            "-evalue",        str(self.e_value),
            "-num_threads",   str(self.threads),
            "-perc_identity", str(self.min_pident),
        ]

        logger.info(
            f"BLASTn [HK→contigs]: {len(gene_kingdom)} HK gene(s) vs assembled contigs ..."
        )
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"BLASTn [HK→contigs] failed: {result.stderr.strip()}")
            return {}

        hits = self._parse_best_hits(out_tsv, gene_kingdom)
        logger.info(f"BLASTn [HK→contigs]: {len(hits)} contig(s) matched HK genes.")
        return hits

    def _parse_best_hits(
        self,
        tsv:          Path,
        gene_kingdom: dict[str, str],
    ) -> dict[str, tuple[str, str]]:
        best: dict[str, tuple[str, str, float]] = {}
        if not tsv.exists() or tsv.stat().st_size == 0:
            return {}
        with open(tsv, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 6:
                    continue
                try:
                    gene_id   = parts[0]
                    contig_id = parts[1]
                    pident    = float(parts[2])
                    qcov      = int(float(parts[3]))
                    bitscore  = float(parts[5])
                    if pident < self.min_pident or qcov < self.min_qcov:
                        continue
                    kingdom = gene_kingdom.get(gene_id, "UNKNOWN")
                    if contig_id not in best or bitscore > best[contig_id][2]:
                        best[contig_id] = (kingdom, gene_id, bitscore)
                except (ValueError, IndexError):
                    continue
        return {cid: (k, gid) for cid, (k, gid, _) in best.items()}


# ---------------------------------------------------------------------------
# CombinedFastaWriter  (reference pathway — pathway 1)
# ---------------------------------------------------------------------------

class CombinedFastaWriter:
    """
    Writes the merged reference FASTA for reference-pathway Salmon indexing:

      1. Viral sequences                      (prefix ``VQ_VIRAL_``)
      2. Conserved bundled HK genes           (prefix ``VQ_CONS_{KINGDOM}_``)
      3. User reference HK genes (optional)   (prefix ``VQ_REFHK_``)
      4. User transcriptome (original IDs, skipping ref-HK IDs)

    Returns
    -------
    tuple of:
      - set[str]          viral prefixed IDs
      - dict[str, str]    {conserved prefixed ID: kingdom}
    """

    @staticmethod
    def write(
        output_path:        Path,
        user_transcriptome: Path,
        viral_seqs:         list[NucSequence],
        conserved:          dict[str, list[tuple[str, str, str]]],
        ref_hk:             list[tuple[str, str, str]] | None = None,
        hk_original_ids:    set[str] | None                  = None,
    ) -> tuple[set[str], dict[str, str]]:
        viral_ids:     set[str]       = set()
        conserved_map: dict[str, str] = {}
        skip_ids = hk_original_ids or set()

        with open(output_path, "w", encoding="utf-8") as fh:

            # 1 — viral sequences
            for seq in viral_seqs:
                prefixed = f"VQ_VIRAL_{seq.id}"
                fh.write(f">{prefixed}\n{seq.sequence}\n")
                viral_ids.add(prefixed)

            # 2 — conserved bundled HK genes (all kingdoms)
            for kingdom, entries in conserved.items():
                for prefixed_id, desc, sequence in entries:
                    header = f">{prefixed_id}" + (f" {desc}" if desc else "")
                    fh.write(f"{header}\n{sequence}\n")
                    conserved_map[prefixed_id] = kingdom

            # 3 — user reference HK genes (VQ_REFHK_)
            for prefixed_id, desc, sequence in (ref_hk or []):
                header = f">{prefixed_id}" + (f" {desc}" if desc else "")
                fh.write(f"{header}\n{sequence}\n")

            # 4 — user transcriptome, skipping ref-HK IDs
            for seq_id, desc, seq in _parse_fasta_iter(user_transcriptome):
                if seq_id in skip_ids:
                    continue
                header = f">{seq_id}" + (f" {desc}" if desc else "")
                fh.write(f"{header}\n{seq}\n")

        ref_hk_n = len(ref_hk) if ref_hk else 0
        logger.info(
            f"Combined FASTA: {len(viral_ids)} viral + "
            f"{len(conserved_map)} bundled HK + {ref_hk_n} ref HK + transcriptome "
            f"→ '{output_path.name}'"
        )
        return viral_ids, conserved_map


# ---------------------------------------------------------------------------
# DeNovoFastaWriter  (de-novo pathway — pathway 2)
# ---------------------------------------------------------------------------

class DeNovoFastaWriter:
    """
    Writes the reference FASTA for de-novo pathway Salmon indexing.

    Every contig from ``assembled_fasta`` is written with a tag prefix:
      - Confirmed viral contigs   → ``VQ_VIRAL_``
      - HK-matched contigs        → ``VQ_CONS_{KINGDOM}_``
      - All other contigs         → original ID (unchanged)

    Returns
    -------
    tuple of:
      - set[str]          viral prefixed IDs
      - dict[str, str]    {conserved prefixed ID: kingdom}
    """

    @staticmethod
    def write(
        output_path:     Path,
        assembled_fasta: Path,
        viral_seqs:      list[NucSequence],
        hk_contig_map:   dict[str, tuple[str, str]],
    ) -> tuple[set[str], dict[str, str]]:
        viral_id_set   = {s.id for s in viral_seqs if s.is_viral}
        viral_ids:     set[str]       = set()
        conserved_map: dict[str, str] = {}

        with open(output_path, "w", encoding="utf-8") as fh:
            for seq_id, desc, seq in _parse_fasta_iter(assembled_fasta):
                if seq_id in viral_id_set:
                    prefixed = f"VQ_VIRAL_{seq_id}"
                    viral_ids.add(prefixed)
                    header = f">{prefixed}" + (f" {desc}" if desc else "")
                elif seq_id in hk_contig_map:
                    kingdom, _ = hk_contig_map[seq_id]
                    prefixed   = f"VQ_CONS_{kingdom}_{seq_id}"
                    conserved_map[prefixed] = kingdom
                    header = f">{prefixed}" + (f" {desc}" if desc else "")
                else:
                    header = f">{seq_id}" + (f" {desc}" if desc else "")
                fh.write(f"{header}\n{seq}\n")

        logger.info(
            f"De-novo FASTA: {len(viral_ids)} viral + "
            f"{len(conserved_map)} HK-matched + rest unchanged "
            f"→ '{output_path.name}'"
        )
        return viral_ids, conserved_map


# ---------------------------------------------------------------------------
# SalmonIndexBuilder
# ---------------------------------------------------------------------------

class SalmonIndexBuilder:
    """Builds a Salmon quasi-mapping index from the combined reference FASTA."""

    def __init__(self, salmon_bin: str = "salmon", threads: int = 4, kmer_len: int = 31):
        self.salmon_bin = salmon_bin
        self.threads    = threads
        self.kmer_len   = kmer_len

    def build(self, ref_fasta: Path, index_dir: Path) -> None:
        index_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.salmon_bin, "index",
            "-t", str(ref_fasta),
            "-i", str(index_dir),
            "-p", str(self.threads),
            "-k", str(self.kmer_len),
        ]
        logger.info(f"Salmon index: '{ref_fasta.name}' → '{index_dir.name}' ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
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

    def __init__(self, salmon_bin: str = "salmon", threads: int = 4, low_memory: bool = False):
        self.salmon_bin = salmon_bin
        self.threads    = threads
        self.low_memory = low_memory

    def run(
        self,
        index_dir: Path,
        reads:     list[str],
        out_dir:   Path,
    ) -> tuple[Path, float]:
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
        if self.low_memory:
            cmd += ["--gcBias", "--reduceGCMemory"]
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

    Tagging rules (based on the ID prefix written by the FASTA writers):
      ``VQ_VIRAL_*``        → seq_type = "viral",     kingdom = ""
      ``VQ_CONS_{KG}_*``   → seq_type = "conserved",  kingdom = "{KG}"
      ``VQ_REFHK_*``        → seq_type = "ref_hk",   kingdom = ""
      anything else        → seq_type = "host",       kingdom = ""
    """

    @staticmethod
    def parse(
        quant_sf:      Path,
        conserved_map: dict[str, str],
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
                    elif name.startswith("VQ_REFHK_"):
                        seq_type = "ref_hk"
                        kingdom  = ""
                    elif name.startswith("VQ_CONS_"):
                        seq_type = "conserved"
                        kingdom  = conserved_map.get(name, "")
                        if not kingdom:
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
        ref_hk_n    = sum(1 for e in entries if e.seq_type == "ref_hk")
        host_n      = sum(1 for e in entries if e.seq_type == "host")
        logger.info(
            f"quant.sf: {viral_n} viral | {conserved_n} conserved | "
            f"{ref_hk_n} ref_hk | {host_n} host entries."
        )
        return entries


# ---------------------------------------------------------------------------
# TranscriptomeViralAligner  (reference pathway only)
# ---------------------------------------------------------------------------

class TranscriptomeViralAligner:
    """
    BLASTn of the user transcriptome against the viralquest viral sequences.

    Identifies host transcripts with nucleotide similarity to detected viral
    sequences — endogenous viral elements (EVEs), chimeric contigs, etc.

    The viral sequences form a small database (typically 5–100 sequences),
    so ``blastn -subject`` is used to avoid a ``makeblastdb`` dependency.
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

    Two pathways are supported:

    Reference (``user_transcriptome`` provided)
        Index = viral + bundled HK + ref-HK (optional) + transcriptome.
        BLASTn transcriptome vs viral sequences to find EVEs.
        ``ref_hk_quant`` populated when ``hk_genes_file`` is given.

    De-novo (no ``user_transcriptome``)
        Index = assembled contigs with viral/HK-matched ones prefixed.
        BLASTn bundled HK genes vs contigs to identify normalization anchors.
        No EVE detection (no reference transcriptome to compare against).
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
        low_memory:  bool  = False,
    ):
        self._index_builder = SalmonIndexBuilder(salmon_bin, threads, kmer_len=25 if low_memory else 31)
        self._quant_runner  = SalmonQuantRunner(salmon_bin, threads, low_memory=low_memory)
        self._aligner       = TranscriptomeViralAligner(
            blastn_bin, e_value, threads, min_pident, min_qcov,
        )
        self._hk_blaster = ContigHkBlaster(blastn_bin, e_value, threads)
        self.cleanup     = cleanup

    # --- public ---------------------------------------------------------------

    def run(
        self,
        reads:              list[str],
        viral_seqs:         list[NucSequence],
        outdir:             Path,
        assembled_fasta:    Path,
        user_transcriptome: Path | None = None,
        hk_genes_file:      Path | None = None,
    ) -> SalmonQuantReport:
        """
        Execute the full pipeline and return a ``SalmonQuantReport``.

        Parameters
        ----------
        reads              : ``[single.fastq]`` or ``[r1.fastq, r2.fastq]``
        viral_seqs         : all ``NucSequence`` objects from the pipeline
        outdir             : base directory for all output (created if absent)
        assembled_fasta    : source contigs (CAP3 output or original input)
        user_transcriptome : if given, run reference pathway; else de-novo
        hk_genes_file      : plain-text file of reference HK IDs (reference only)
        """
        if user_transcriptome is not None:
            return self._run_reference(
                reads, viral_seqs, outdir, assembled_fasta,
                user_transcriptome, hk_genes_file,
            )
        return self._run_de_novo(reads, viral_seqs, outdir, assembled_fasta)

    # --- reference pathway ----------------------------------------------------

    def _run_reference(
        self,
        reads:              list[str],
        viral_seqs:         list[NucSequence],
        outdir:             Path,
        assembled_fasta:    Path,
        user_transcriptome: Path,
        hk_genes_file:      Path | None,
    ) -> SalmonQuantReport:
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

            ref_hk: list[tuple[str, str, str]] = []
            hk_original_ids: set[str] = set()
            if hk_genes_file:
                ref_hk, hk_original_ids = ReferenceHkLoader.load(
                    hk_genes_file, user_transcriptome
                )

            combined_fasta = tmp_dir / "combined_ref.fasta"
            _, conserved_map = CombinedFastaWriter.write(
                output_path        = combined_fasta,
                user_transcriptome = user_transcriptome,
                viral_seqs         = confirmed,
                conserved          = conserved,
                ref_hk             = ref_hk or None,
                hk_original_ids    = hk_original_ids or None,
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
                pathway      = "reference",
            )

        finally:
            if self.cleanup:
                shutil.rmtree(index_dir, ignore_errors=True)
                shutil.rmtree(tmp_dir,   ignore_errors=True)
                logger.debug("SalmonQuantPipeline: temporary files removed.")

        return report

    # --- de-novo pathway ------------------------------------------------------

    def _run_de_novo(
        self,
        reads:           list[str],
        viral_seqs:      list[NucSequence],
        outdir:          Path,
        assembled_fasta: Path,
    ) -> SalmonQuantReport:
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
            hk_contig_map = self._hk_blaster.match(assembled_fasta, tmp_dir)

            denovo_fasta = tmp_dir / "denovo_ref.fasta"
            _, conserved_map = DeNovoFastaWriter.write(
                output_path     = denovo_fasta,
                assembled_fasta = assembled_fasta,
                viral_seqs      = confirmed,
                hk_contig_map   = hk_contig_map,
            )

            self._index_builder.build(denovo_fasta, index_dir)

            quant_sf, mapping_rate = self._quant_runner.run(index_dir, reads, quant_dir)

            all_entries = QuantsfParser.parse(quant_sf, conserved_map)

            report = self._assemble(
                reads        = reads,
                mapping_rate = mapping_rate,
                entries      = all_entries,
                blast_hits   = [],
                pathway      = "de_novo",
            )

        finally:
            if self.cleanup:
                shutil.rmtree(index_dir, ignore_errors=True)
                shutil.rmtree(tmp_dir,   ignore_errors=True)
                logger.debug("SalmonQuantPipeline: temporary files removed.")

        return report

    # --- assembly helper ------------------------------------------------------

    @staticmethod
    def _assemble(
        reads:        list[str],
        mapping_rate: float,
        entries:      list[SalmonEntry],
        blast_hits:   list[HostViralHit],
        pathway:      str,
    ) -> SalmonQuantReport:
        viral_quant     = [e for e in entries if e.seq_type == "viral"]
        conserved_quant = [e for e in entries if e.seq_type == "conserved"]
        ref_hk_quant    = [e for e in entries if e.seq_type == "ref_hk"]
        host_map        = {e.name: e for e in entries if e.seq_type == "host"}

        total_reads = int(sum(e.num_reads for e in entries))

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

        host_viral_hits.sort(
            key=lambda r: r.blast_hits[0].bit_score if r.blast_hits else 0.0,
            reverse=True,
        )

        conserved_quant.sort(key=lambda e: (e.kingdom, -e.tpm))

        logger.info(
            f"SalmonQuantReport assembled [{pathway}]: "
            f"{len(viral_quant)} viral | "
            f"{len(conserved_quant)} conserved | "
            f"{len(ref_hk_quant)} ref_hk | "
            f"{len(host_viral_hits)} host-viral record(s)."
        )
        return SalmonQuantReport(
            reads            = reads,
            mapping_rate     = mapping_rate,
            total_reads      = total_reads,
            pathway          = pathway,
            viral_quant      = viral_quant,
            conserved_quant  = conserved_quant,
            ref_hk_quant     = ref_hk_quant,
            host_viral_hits  = host_viral_hits,
        )
