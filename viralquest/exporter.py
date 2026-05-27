import dataclasses
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from loguru import logger

from viralquest.biodata import (
    InputFasta,
    NucSequence,
    SalmonQuantReport,
    Taxonomy,
    ViralCluster,
)


# ---------------------------------------------------------------------------
# Generic recursive converter
# ---------------------------------------------------------------------------

def _to_serializable(obj: Any) -> Any:
    """
    Recursively convert a dataclass tree to JSON-serializable Python types.

    Handles:
      - uuid.UUID  → str
      - pathlib.Path → str
      - dataclasses  → dict of field-name: converted-value
      - list / dict  → recursed element-wise
      - everything else → returned as-is (int, float, str, bool, None)
    """
    if isinstance(obj, uuid.UUID):
        return str(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, list):
        return [_to_serializable(item) for item in obj]
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {
            f.name: _to_serializable(getattr(obj, f.name))
            for f in dataclasses.fields(obj)
        }
    return obj


# ---------------------------------------------------------------------------
# ReportExporter
# ---------------------------------------------------------------------------

class ReportExporter:
    """
    Serialises pipeline output to a structured JSON report.

    The top-level JSON has three keys:

        meta       — run metadata (timestamp, version, input file info)
        sequences  — one record per NucSequence with all nested annotations
        clusters   — one record per ViralCluster with member alignments

    Notes
    -----
    - ViralFamilyInfo is intentionally excluded: it is an LLM input artefact
      and carries no additional information beyond what Taxonomy already holds.
    - Raw nucleotide and amino-acid sequences are always included so the HTML
      report is fully self-contained.
    - UUIDs are serialised as strings; None values are preserved.

    Parameters
    ----------
    include_llm : bool
        Include LlmOutput in each sequence record.  Set to False when the
        pipeline was run without LLM scoring so the field is simply absent
        from the output rather than appearing as null everywhere.
    force : bool
        When True, export every sequence that entered the pipeline regardless
        of viral confirmation status.  Equivalent to --force in the CLI.
        Useful when the user wants to inspect all sequences or perform their
        own downstream filtering.  Clusters are also exported in full.
    """

    def __init__(self, include_llm: bool = True, force: bool = False):
        self.include_llm = include_llm
        self.force       = force

    # --- public ---------------------------------------------------------------

    def export(
        self,
        nuc_seqs:      list[NucSequence],
        clusters:      list[ViralCluster],
        input_fasta:   InputFasta | None        = None,
        output_path:   str | Path | None        = None,
        version:       str                      = "unknown",
        salmon_report: SalmonQuantReport | None = None,
    ) -> dict:
        """
        Build the report dictionary, optionally write it to *output_path*,
        and return it for programmatic use.

        Parameters
        ----------
        nuc_seqs      : sequences produced by the pipeline
        clusters      : ViralCluster objects from SequenceTracker
        input_fasta   : InputFasta from FastaParser (None if not available)
        output_path   : if given, the report is written as indented JSON here
        version       : viralquest version string embedded in meta
        salmon_report : if given, a ``salmon_quant`` key is added to the JSON
        """
        if self.force:
            confirmed          = nuc_seqs
            confirmed_clusters = clusters
            logger.warning(
                f"Force mode: exporting all {len(nuc_seqs)} sequence(s) "
                f"without viral confirmation filters."
            )
        else:
            nr_was_run = any(s.blastx_nr_hits for s in nuc_seqs)
            if nr_was_run:
                confirmed = [s for s in nuc_seqs if s.is_viral and s.blastx_nr_hits]
                filter_label = "NR BLASTx viral confirmation"
            else:
                confirmed = [s for s in nuc_seqs if s.is_viral]
                filter_label = "RefSeq / HMM viral confirmation (NR not run)"
            dropped = len(nuc_seqs) - len(confirmed)
            if dropped:
                logger.info(
                    f"Exporter: {dropped} sequence(s) excluded — did not pass "
                    f"{filter_label}."
                )
            confirmed_ids      = {s.id for s in confirmed}
            confirmed_clusters = [
                c for c in clusters if c.representative_id in confirmed_ids
            ]

        report = {
            "meta":      self._build_meta(input_fasta, version),
            "sequences": [self._seq_to_dict(s) for s in confirmed],
            "clusters":  [self._cluster_to_dict(c) for c in confirmed_clusters],
        }

        if salmon_report is not None:
            report["salmon_quant"] = _to_serializable(salmon_report)

        if output_path is not None:
            path = Path(output_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(report, fh, indent=2, ensure_ascii=False)
            logger.success(
                f"Report written → '{path}'  "
                f"({len(confirmed)} sequences, {len(confirmed_clusters)} clusters)."
            )

        return report

    # --- private --------------------------------------------------------------

    def _build_meta(self, input_fasta: InputFasta | None, version: str) -> dict:
        meta: dict = {
            "viralquest_version": version,
            "timestamp":          datetime.now(timezone.utc).isoformat(),
        }
        if input_fasta is not None:
            meta["input_file"] = {
                "name":       input_fasta.name,
                "path":       str(input_fasta.path),
                "num_seqs":   input_fasta.num_seqs,
                "size_bytes": input_fasta.size,
            }
        return meta

    def _seq_to_dict(self, seq: NucSequence) -> dict:
        """
        Convert one NucSequence applying field-level rules:
          - viral_family_info is excluded (LLM input artefact)
          - llm_output is included only when self.include_llm is True
        """
        d: dict = {
            "id":             seq.id,
            "uid":            str(seq.uid),
            "sample_name":    seq.sample_name,
            "sequence":       seq.sequence,
            "length":         seq.length,
            "n_count":        seq.n_count,
            "gc_content":     seq.gc_content,
            "is_viral":       seq.is_viral,
            "cluster_id":     seq.cluster_id,
            "orfs":           [_to_serializable(o) for o in seq.orfs],
            "blastx_hits":    [_to_serializable(h) for h in seq.blastx_hits],
            "blastx_nr_hits": [_to_serializable(h) for h in seq.blastx_nr_hits],
            "blastn_hits":    [_to_serializable(h) for h in seq.blastn_hits],
            "taxonomy":       self._taxonomy_to_dict(seq.taxonomy) if seq.taxonomy else None,
        }
        if self.include_llm:
            d["llm_output"] = _to_serializable(seq.llm_output) if seq.llm_output else None
        return d

    @staticmethod
    def _cluster_to_dict(cluster: ViralCluster) -> dict:
        """
        Serialize a ViralCluster, explicitly adding the computed 'size'
        property which is not a dataclass field and would otherwise be missed.
        """
        return {
            "cluster_id":        cluster.cluster_id,
            "species":           cluster.species,
            "representative_id": cluster.representative_id,
            "size":              cluster.size,
            "members":           [_to_serializable(m) for m in cluster.members],
        }

    @staticmethod
    def _taxonomy_to_dict(tax: Taxonomy) -> dict:
        """
        Serialize Taxonomy, renaming the Python-safe 'class_' field to 'class'
        since JSON has no keyword restrictions.
        """
        d = _to_serializable(tax)
        d["class"] = d.pop("class_")
        return d
