import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from loguru import logger

from viralquest.biodata import (
    HostViralHit,
    HostViralRecord,
    NucSequence,
    PfamHkEntry,
    SalmonEntry,
    SalmonQuantReport,
)




# ---------------------------------------------------------------------------
# PfamHkFinder
# ---------------------------------------------------------------------------

class PfamHkFinder:
    """
    Identifies housekeeping-gene contigs by running a targeted pyhmmer search
    against 22 high-confidence Pfam HK domain profiles, using the pre-computed
    ORFs stored in each NucSequence, then links hits to Salmon quantification data.

    Only the de-novo pathway is supported.

    DISABLED — not wired into the pipeline (see SalmonQuantPipeline.__init__ and
    _run_de_novo).  ``pfam_hk_quant`` is therefore always ``[]`` and every
    downstream consumer already treats that as "panel absent".  Housekeeping
    normalization currently comes only from the predefined gene sets written
    into the Salmon index: the user's reference HK list (``VQ_REFHK_``) and the
    bundled per-kingdom conserved genes (``VQ_CONS_``).

    Known defects to settle before re-enabling:
      · ``_search`` reads ``top_hits.query_name``, removed in pyhmmer 0.10+;
        the attribute is ``top_hits.query.name`` (as in hmm.py's HmmSearcher).
        The AttributeError is swallowed by the broad except, so the step
        silently yields zero hits on every run.
      · ``_HK_FILE`` / ``_KINGDOMS_DIR`` resolve to ``<package>/../data``, which
        only exists in a source checkout — ``data/`` is not listed in
        pyproject's package-data.
      · ``_search``'s docstring claims GA inclusion thresholds, but hmmsearch is
        called without ``bit_cutoffs="gathering"``, so plain E-value cutoffs apply.
    """

    _HK_FILE = Path(__file__).parent.parent / "data" / "salmon-quant" / "hk_pfam_domains.json"

    def __init__(self, pfam_hmm_path: Path, threads: int = 4):
        self.pfam_hmm_path = pfam_hmm_path
        self.threads       = threads
        self._hk_meta: dict[str, dict] = self._load_hk_meta()

    # ── public ──────────────────────────────────────────────────────────────

    def find(self, seqs: list[NucSequence], quant_sf: Path) -> list[PfamHkEntry]:
        if not self._hk_meta:
            logger.debug("PfamHkFinder: HK domain file not found — skipping.")
            return []
        if not self.pfam_hmm_path.exists():
            logger.debug("PfamHkFinder: Pfam-A HMM not found — skipping.")
            return []

        tpm_map = self._parse_quant_sf(quant_sf)
        if not tpm_map:
            return []

        seq_block, name_map = self._build_seq_block(seqs)
        if seq_block is None:
            return []

        profiles = self._load_profiles()
        if not profiles:
            return []

        hits = self._search(profiles, seq_block, name_map)
        entries = self._build_entries(hits, tpm_map)
        logger.info(f"PfamHkFinder: {len(entries)} HK contig(s) identified via Pfam domains.")
        return entries

    # ── private ─────────────────────────────────────────────────────────────

    @classmethod
    def _load_hk_meta(cls) -> dict[str, dict]:
        if not cls._HK_FILE.exists():
            return {}
        with open(cls._HK_FILE, encoding="utf-8") as fh:
            entries = json.load(fh)
        return {e["pfam_target_id"]: e for e in entries}

    @staticmethod
    def _parse_quant_sf(quant_sf: Path) -> dict[str, tuple[float, float]]:
        """Returns {seq_id: (tpm, num_reads)}."""
        result: dict[str, tuple[float, float]] = {}
        if not quant_sf.exists():
            return result
        with open(quant_sf, encoding="utf-8") as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 5:
                    continue
                raw = parts[0]
                if raw.startswith("VQ_VIRAL_"):
                    name = raw[len("VQ_VIRAL_"):]
                elif raw.startswith("VQ_CONS_"):
                    name = raw.split("_", 3)[-1]  # drop VQ_CONS_{KINGDOM}_
                else:
                    name = raw
                try:
                    result[name] = (float(parts[3]), float(parts[4]))
                except ValueError:
                    continue
        return result

    def _build_seq_block(
        self, seqs: list[NucSequence]
    ) -> tuple["pyhmmer.easel.DigitalSequenceBlock | None", dict[str, str]]:
        """
        Build a pyhmmer DigitalSequenceBlock from the pre-computed ORFs stored
        in each NucSequence.  Returns the block and a name_map from orf.name
        back to the parent contig ID.
        """
        try:
            import pyhmmer.easel
        except ImportError:
            logger.warning("PfamHkFinder: pyhmmer not available — skipping.")
            return None, {}

        alphabet : pyhmmer.easel.Alphabet              = pyhmmer.easel.Alphabet.amino()
        digital  : list[pyhmmer.easel.DigitalSequence] = []
        name_map : dict[str, str]                      = {}

        for nuc_seq in seqs:
            for orf in nuc_seq.orfs:
                aa = orf.aa_sequence
                if not aa or len(aa) < 30:
                    continue
                name_map[orf.name] = nuc_seq.id
                digital.append(
                    pyhmmer.easel.TextSequence(
                        name=orf.name.encode(),
                        sequence=aa,
                    ).digitize(alphabet)
                )

        if not digital:
            return None, {}

        return pyhmmer.easel.DigitalSequenceBlock(alphabet, digital), name_map

    def _load_profiles(self) -> list:
        """Scan Pfam-A.hmm and return only the HK profiles."""
        try:
            import pyhmmer.plan7
        except ImportError:
            return []

        hk_ids   = set(self._hk_meta.keys())
        profiles = []
        try:
            with pyhmmer.plan7.HMMFile(str(self.pfam_hmm_path)) as hf:
                for hmm in hf:
                    raw = hmm.name
                    name = raw if isinstance(raw, str) else raw.decode("utf-8", errors="replace")
                    if name in hk_ids:
                        profiles.append(hmm)
        except Exception as exc:
            logger.warning(f"PfamHkFinder: error loading Pfam profiles: {exc}")
        logger.debug(f"PfamHkFinder: {len(profiles)}/{len(hk_ids)} HK profiles loaded.")
        return profiles

    def _search(
        self,
        profiles: list,
        seq_block: "pyhmmer.easel.DigitalSequenceBlock",
        name_map: dict[str, str],
    ) -> dict[str, dict[str, tuple[float, float]]]:
        """
        Returns {contig_id: {pfam_target: (score, evalue)}} keeping the best
        hit per contig per domain.  Only hits passing the GA inclusion threshold
        are kept (``hit.included`` from pyhmmer).
        """
        try:
            import pyhmmer
        except ImportError:
            return {}

        def _s(v) -> str:
            return v if isinstance(v, str) else v.decode("utf-8", errors="replace")

        results: dict[str, dict[str, tuple[float, float]]] = {}
        try:
            for top_hits in pyhmmer.hmmsearch(profiles, seq_block, cpus=self.threads):
                profile_name = _s(top_hits.query_name)
                for hit in top_hits:
                    if not hit.included:
                        continue
                    frag_name = _s(hit.name)
                    contig_id = name_map.get(frag_name)
                    if contig_id is None:
                        continue
                    contig_hits = results.setdefault(contig_id, {})
                    cur = contig_hits.get(profile_name)
                    if cur is None or hit.score > cur[0]:
                        contig_hits[profile_name] = (hit.score, hit.evalue)
        except Exception as exc:
            logger.warning(f"PfamHkFinder: search error: {exc}")
        return results

    def _build_entries(
        self,
        hits: dict[str, dict[str, tuple[float, float]]],
        tpm_map: dict[str, tuple[float, float]],
    ) -> list[PfamHkEntry]:
        entries: list[PfamHkEntry] = []
        for contig_id, domain_hits in hits.items():
            tpm_data = tpm_map.get(contig_id)
            if tpm_data is None:
                continue
            tpm, num_reads = tpm_data
            # Use the domain with the best score for this contig
            best_domain = max(domain_hits, key=lambda d: domain_hits[d][0])
            score, evalue = domain_hits[best_domain]
            meta = self._hk_meta.get(best_domain, {})
            entries.append(PfamHkEntry(
                seq_id       = contig_id,
                pfam_target  = best_domain,
                pfam_acc     = meta.get("pfam_accession", ""),
                pfam_desc    = meta.get("pfam_description", ""),
                pfam_details = meta.get("pfam_details", ""),
                tpm          = tpm,
                num_reads    = num_reads,
                score        = score,
                e_value      = evalue,
            ))
        entries.sort(key=lambda e: -e.tpm)
        return entries


# ---------------------------------------------------------------------------
# Kingdom FASTA paths
# ---------------------------------------------------------------------------

_KINGDOMS_DIR = Path(__file__).parent.parent / "data" / "salmon-quant"

_KINGDOM_FILES: dict[str, Path] = {
    name: _KINGDOMS_DIR / f"{name}.fasta"
    for name in (
        "mammals", "arthropods", "plants",
        "fish", "fungi", "bacteria", "nematodes",
        "human", "mouse",
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


def _parse_fasta_headers(path: Path):
    """Yield (seq_id, description) for each record — headers only, no sequence.

    Used when a first pass over a large transcriptome only needs to build an
    ID index; keeps peak memory flat regardless of transcriptome size.
    """
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip().split(None, 1)
            if not parts:
                continue
            yield parts[0], (parts[1] if len(parts) > 1 else "")


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

    Reads gene IDs from a plain-text file (one ID per line; ``#`` comments and
    blank lines ignored, ``>`` prefixes stripped, first column of a TSV/CSV
    line taken) and finds the matching records in the transcriptome FASTA.

    Matching is deliberately tolerant: HK lists are normally copied out of NCBI
    tables, spreadsheets or gene-symbol lists and rarely reproduce the FASTA
    header verbatim. For every wanted ID the tiers below are tried in order and
    the first hit wins. Tiers 2-5 only accept a key that resolves to exactly one
    transcript, so a loose ID can never silently grab the wrong sequence:

      1. exact first header token                  ``NM_001101.5``
      2. case-insensitive first token              ``nm_001101.5``
      3. accession without trailing version        ``NM_001101``
      4. gene symbol in parentheses in the header  ``... beta (ACTB), mRNA``
      5. identifier-looking token in the header    ``... ACTB ...``

    Returns
    -------
    tuple of:
      - list[tuple[str, str, str]]   (``VQ_REFHK_``-prefixed entries)
      - set[str]                     original IDs (for transcriptome exclusion)
    """

    _TIER_NAMES = (
        "exact ID",
        "case-insensitive ID",
        "unversioned accession",
        "gene symbol",
        "header token",
    )

    _SYMBOL_RE = re.compile(r"\(([A-Za-z0-9][A-Za-z0-9_.\-]{1,24})\)")
    _TOKEN_RE  = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.\-]{2,}")
    _VERSION_RE = re.compile(r"\.\d+$")

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

        logger.info(
            f"ReferenceHkLoader: {len(wanted)} HK ID(s) read from "
            f"'{hk_genes_file.name}'."
        )

        tiers = cls._build_header_index(transcriptome)

        chosen:     dict[str, str] = {}    # transcriptome seq_id -> wanted ID
        tier_stats: dict[str, int] = {}
        missing:    set[str]       = set()
        redundant:  set[str]       = set()

        for want in sorted(wanted):
            hit = cls._resolve(want, tiers, taken=set(chosen))
            if hit is None:
                # Distinguish "no such transcript" from "the transcript it points
                # at was already claimed by an earlier ID in the list".
                if cls._resolve(want, tiers, taken=set()) is None:
                    missing.add(want)
                else:
                    redundant.add(want)
                continue
            seq_id, tier = hit
            chosen[seq_id] = want
            tier_stats[tier] = tier_stats.get(tier, 0) + 1

        if missing:
            sample = ", ".join(sorted(missing)[:10])
            logger.warning(
                f"ReferenceHkLoader: {len(missing)}/{len(wanted)} ID(s) not found "
                f"in '{transcriptome.name}': {sample}"
                + ("…" if len(missing) > 10 else "")
            )
        if redundant:
            sample = ", ".join(sorted(redundant)[:10])
            logger.info(
                f"ReferenceHkLoader: {len(redundant)} ID(s) resolve to a transcript "
                f"already claimed by another entry in the list: {sample}"
                + ("…" if len(redundant) > 10 else "")
            )

        if not chosen:
            logger.warning(
                "ReferenceHkLoader: no HK ID matched the transcriptome headers — "
                "the 'Reference Housekeeping Genes' card will be missing from the "
                "HTML report. Check that the IDs correspond to the first word of "
                f"the '>' lines in '{transcriptome.name}'."
            )
            return [], set()

        entries:   list[tuple[str, str, str]] = []
        found_ids: set[str] = set()
        for seq_id, desc, seq in _parse_fasta_iter(transcriptome):
            if seq_id in chosen and seq_id not in found_ids:
                entries.append((f"VQ_REFHK_{seq_id}", desc, seq))
                found_ids.add(seq_id)

        detail = ", ".join(f"{n}× {t}" for t, n in tier_stats.items())
        logger.info(
            f"ReferenceHkLoader: {len(entries)} reference HK gene(s) loaded "
            f"from '{hk_genes_file.name}' ({detail})."
        )
        return entries, found_ids

    # --- internals -----------------------------------------------------------

    @classmethod
    def _build_header_index(
        cls,
        transcriptome: Path,
    ) -> list[dict[str, list[str]]]:
        """Build one lookup dict per matching tier: ``{key: [seq_id, ...]}``."""
        tiers: list[dict[str, list[str]]] = [{} for _ in cls._TIER_NAMES]

        def _add(tier: int, key: str, seq_id: str) -> None:
            if key:
                tiers[tier].setdefault(key, []).append(seq_id)

        for seq_id, desc in _parse_fasta_headers(transcriptome):
            low = seq_id.lower()
            _add(0, seq_id, seq_id)
            _add(1, low, seq_id)
            _add(2, cls._VERSION_RE.sub("", low), seq_id)
            for sym in cls._SYMBOL_RE.findall(desc):
                _add(3, sym.lower(), seq_id)
            for tok in cls._TOKEN_RE.findall(desc):
                _add(4, tok.lower(), seq_id)

        return tiers

    @classmethod
    def _resolve(
        cls,
        want:  str,
        tiers: list[dict[str, list[str]]],
        taken: set[str],
    ) -> tuple[str, str] | None:
        """Return ``(seq_id, tier_name)`` for the first tier that resolves."""
        low  = want.lower()
        keys = (want, low, cls._VERSION_RE.sub("", low), low, low)

        for i, key in enumerate(keys):
            candidates = [s for s in tiers[i].get(key, ()) if s not in taken]
            if not candidates:
                continue
            # The exact tier may legitimately repeat an ID; looser tiers must be
            # unambiguous or they are skipped.
            if i == 0 or len(set(candidates)) == 1:
                return candidates[0], cls._TIER_NAMES[i]
            logger.debug(
                f"ReferenceHkLoader: '{want}' is ambiguous at tier "
                f"'{cls._TIER_NAMES[i]}' ({len(set(candidates))} transcripts) — "
                "skipping this tier."
            )
        return None

    @staticmethod
    def _read_ids(path: Path) -> set[str]:
        ids: set[str] = set()
        with open(path, encoding="utf-8-sig") as fh:
            for raw in fh:
                line = raw.split("#", 1)[0].strip()
                if not line:
                    continue
                if line.startswith(">"):
                    line = line[1:].strip()
                token = re.split(r"[\s,;]+", line)[0].strip().strip("\"'")
                if token:
                    ids.add(token)
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
      2. User reference HK genes (optional)   (prefix ``VQ_REFHK_``)
      3. Conserved bundled HK genes           (prefix ``VQ_CONS_{KINGDOM}_``)
      4. User transcriptome (original IDs, skipping ref-HK IDs)

    Records are de-duplicated **by sequence content**, not only by ID. Salmon
    discards sequence-identical transcripts at index time unless told otherwise,
    and the discarded ones never reach ``quant.sf`` — so a user HK gene that is
    byte-identical to a bundled conserved gene (very common: both come from
    RefSeq) would silently vanish from the report. Doing the de-duplication here
    keeps it deterministic and lets the write order define the priority:
    viral > ref-HK > conserved > transcriptome.

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
        seen_ids:  set[str]            = set()
        seen_seqs: dict[bytes, str]    = {}
        dup_seq_counts: dict[str, int] = {}

        def _write(fh, record_id: str, sequence: str, desc: str = "",
                   category: str = "transcriptome") -> bool:
            if record_id in seen_ids:
                logger.warning(f"Duplicate FASTA ID '{record_id}' — skipping to avoid salmon index failure.")
                return False

            digest = hashlib.sha1(sequence.upper().encode()).digest() if sequence else None
            if digest is not None:
                owner = seen_seqs.get(digest)
                if owner is not None:
                    dup_seq_counts[category] = dup_seq_counts.get(category, 0) + 1
                    logger.debug(
                        f"Sequence-identical duplicate: '{record_id}' == '{owner}' "
                        "— skipping (salmon would discard it at index time anyway)."
                    )
                    return False
                seen_seqs[digest] = record_id

            seen_ids.add(record_id)
            header = f">{record_id}" + (f" {desc}" if desc else "")
            fh.write(f"{header}\n{sequence}\n")
            return True

        with open(output_path, "w", encoding="utf-8") as fh:

            # 1 — viral sequences (highest priority, never de-duplicated away)
            for seq in viral_seqs:
                prefixed = f"VQ_VIRAL_{seq.id}"
                if _write(fh, prefixed, seq.sequence, category="viral"):
                    viral_ids.add(prefixed)

            # 2 — user reference HK genes (VQ_REFHK_) — written before the
            #     bundled set so a shared RefSeq record is kept as ref-HK
            for prefixed_id, desc, sequence in (ref_hk or []):
                _write(fh, prefixed_id, sequence, desc, category="ref_hk")

            # 3 — conserved bundled HK genes (all kingdoms)
            for kingdom, entries in conserved.items():
                for prefixed_id, desc, sequence in entries:
                    if _write(fh, prefixed_id, sequence, desc, category="conserved"):
                        conserved_map[prefixed_id] = kingdom

            # 4 — user transcriptome, skipping ref-HK IDs
            for seq_id, desc, seq in _parse_fasta_iter(user_transcriptome):
                if seq_id in skip_ids:
                    continue
                _write(fh, seq_id, seq, desc, category="transcriptome")

        ref_hk_n = len(ref_hk) if ref_hk else 0
        logger.info(
            f"Combined FASTA: {len(viral_ids)} viral + {ref_hk_n} ref HK + "
            f"{len(conserved_map)} bundled HK + transcriptome "
            f"→ '{output_path.name}'"
        )
        if dup_seq_counts:
            detail = ", ".join(f"{n} {cat}" for cat, n in sorted(dup_seq_counts.items()))
            logger.info(
                f"Combined FASTA: {sum(dup_seq_counts.values())} sequence-identical "
                f"duplicate(s) removed ({detail}) — a lower-priority copy of a "
                "sequence already in the index."
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
        min_bg_len:      int = 0,
    ) -> tuple[set[str], dict[str, str]]:
        """
        min_bg_len: background contigs (not viral, not HK-matched) shorter than
        this are excluded from the index to reduce SSHash memory usage.
        Viral and HK-matched contigs are always included regardless of length.
        """
        viral_id_set   = {s.id for s in viral_seqs if s.is_viral}
        viral_ids:     set[str]       = set()
        conserved_map: dict[str, str] = {}
        seen:          set[str]       = set()
        skipped_bg = 0

        with open(output_path, "w", encoding="utf-8") as fh:
            for seq_id, desc, seq in _parse_fasta_iter(assembled_fasta):
                if seq_id in viral_id_set:
                    prefixed = f"VQ_VIRAL_{seq_id}"
                elif seq_id in hk_contig_map:
                    kingdom, _ = hk_contig_map[seq_id]
                    prefixed   = f"VQ_CONS_{kingdom}_{seq_id}"
                else:
                    if min_bg_len > 0 and len(seq) < min_bg_len:
                        skipped_bg += 1
                        continue
                    prefixed = seq_id

                if prefixed in seen:
                    logger.warning(
                        f"Duplicate FASTA ID '{prefixed}' in assembled contigs — "
                        "skipping to avoid salmon index failure."
                    )
                    continue
                seen.add(prefixed)

                if prefixed.startswith("VQ_VIRAL_"):
                    viral_ids.add(prefixed)
                elif prefixed.startswith("VQ_CONS_"):
                    conserved_map[prefixed] = kingdom  # type: ignore[possibly-undefined]

                header = f">{prefixed}" + (f" {desc}" if desc else "")
                fh.write(f"{header}\n{seq}\n")

        if skipped_bg:
            logger.info(f"De-novo FASTA: {skipped_bg} background contig(s) <{min_bg_len}bp excluded (low-memory).")
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
            # Without this, salmon silently drops sequence-identical transcripts
            # and they never appear in quant.sf. The FASTA writers already
            # de-duplicate with an explicit priority, so nothing is lost here —
            # the flag only stops salmon from making that decision on its own.
            "--keepDuplicates",
        ]
        logger.info(f"Salmon index: '{ref_fasta.name}' → '{index_dir.name}' ...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            tail = "\n".join(result.stderr.strip().splitlines()[-20:])
            raise RuntimeError(
                f"salmon index failed (exit {result.returncode}):\n{tail}\n\n"
                "Common causes: duplicate sequence IDs in the input FASTA (check "
                "the WARNING lines above), insufficient RAM for SSHash (try closing "
                "other processes), or a corrupted/empty combined reference FASTA."
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
        salmon_bin:    str        = "salmon",
        blastn_bin:    str        = "blastn",
        threads:       int        = 4,
        cleanup:       bool       = True,
        e_value:       float      = 1e-5,
        min_pident:    float      = 70.0,
        min_qcov:      int        = 20,
        low_memory:    bool       = False,
        pfam_hmm_path: Path | None = None,
    ):
        self._index_builder = SalmonIndexBuilder(salmon_bin, threads, kmer_len=21 if low_memory else 31)
        self._quant_runner  = SalmonQuantRunner(salmon_bin, threads, low_memory=low_memory)
        self._aligner       = TranscriptomeViralAligner(
            blastn_bin, e_value, threads, min_pident, min_qcov,
        )
        self._hk_blaster  = ContigHkBlaster(blastn_bin, e_value, threads)
        # Pfam-domain HK identification is parked — see PfamHkFinder's docstring.
        # ``pfam_hmm_path`` is still accepted so callers need no change when it
        # comes back.
        # self._pfam_hk   = PfamHkFinder(pfam_hmm_path, threads) if pfam_hmm_path else None
        self._pfam_hk     = None
        self.cleanup      = cleanup
        self.low_memory   = low_memory

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

            self._check_ref_hk(all_entries, ref_hk, hk_genes_file)

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

    @staticmethod
    def _check_ref_hk(
        entries:       list[SalmonEntry],
        ref_hk:        list[tuple[str, str, str]],
        hk_genes_file: Path | None,
    ) -> None:
        """Warn loudly when --hk-genes was given but produced nothing usable.

        The HTML report only renders the 'Reference Housekeeping Genes' card
        when ``ref_hk_quant`` is non-empty, so an empty list is otherwise
        indistinguishable from 'no --hk-genes given'.
        """
        if hk_genes_file is None:
            return

        got = {e.name for e in entries if e.seq_type == "ref_hk"}

        if not ref_hk:
            logger.warning(
                f"--hk-genes ('{hk_genes_file}') produced no reference HK gene: "
                "no ID matched the transcriptome headers. The 'Reference "
                "Housekeeping Genes' card will be absent from the report."
            )
            return

        missing = {pid for pid, _, _ in ref_hk} - got
        if missing:
            sample = ", ".join(sorted(missing)[:10])
            logger.warning(
                f"{len(missing)}/{len(ref_hk)} reference HK gene(s) were written to "
                f"the index but are absent from quant.sf: {sample}"
                + ("…" if len(missing) > 10 else "")
            )
        else:
            logger.success(
                f"{len(got)} reference HK gene(s) quantified — "
                "'Reference Housekeeping Genes' card will be rendered."
            )

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
                min_bg_len      = 500 if self.low_memory else 0,
            )

            self._index_builder.build(denovo_fasta, index_dir)

            quant_sf, mapping_rate = self._quant_runner.run(index_dir, reads, quant_dir)

            all_entries = QuantsfParser.parse(quant_sf, conserved_map)

            # Pfam-domain HK identification is parked; the de-novo run keeps the
            # HK-matched contigs that ContigHkBlaster already put in the index
            # (``VQ_CONS_*``), which is what conserved_quant reports.
            # pfam_hk: list[PfamHkEntry] = []
            # if self._pfam_hk is not None:
            #     pfam_hk = self._pfam_hk.find(viral_seqs, quant_sf)

            report = self._assemble(
                reads        = reads,
                mapping_rate = mapping_rate,
                entries      = all_entries,
                blast_hits   = [],
                pathway      = "de_novo",
                # pfam_hk    = pfam_hk,
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
        pfam_hk:      list[PfamHkEntry] | None = None,
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

        pfam_hk_quant = pfam_hk or []
        logger.info(
            f"SalmonQuantReport assembled [{pathway}]: "
            f"{len(viral_quant)} viral | "
            f"{len(conserved_quant)} conserved | "
            f"{len(ref_hk_quant)} ref_hk | "
            f"{len(host_viral_hits)} host-viral | "
            f"{len(pfam_hk_quant)} pfam-hk record(s)."
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
            pfam_hk_quant    = pfam_hk_quant,
        )
