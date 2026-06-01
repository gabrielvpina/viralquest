import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from loguru import logger

from viralquest.biodata import ClusterMember, NucSequence, ViralCluster


# ---------------------------------------------------------------------------
# ClusterBuilder
# ---------------------------------------------------------------------------

class ClusterBuilder:
    """
    Groups NucSequence objects by their best BLASTx species name.
    Sequences without a resolvable BLASTx species are silently skipped.
    """

    @staticmethod
    def build(nuc_seqs: list[NucSequence]) -> dict[str, list[NucSequence]]:
        """Returns {species: [NucSequence, ...]} sorted by species name."""
        groups: dict[str, list[NucSequence]] = defaultdict(list)
        skipped = 0
        for seq in nuc_seqs:
            best = seq.best_blastx
            if not best or not best.species:
                skipped += 1
                continue
            groups[best.species].append(seq)

        if skipped:
            logger.warning(
                f"ClusterBuilder: {skipped} sequence(s) skipped — no BLASTx species."
            )
        logger.success(
            f"ClusterBuilder: {len(groups)} cluster(s) from "
            f"{len(nuc_seqs) - skipped} sequence(s)."
        )
        return dict(sorted(groups.items()))


# ---------------------------------------------------------------------------
# RepresentativeSelector
# ---------------------------------------------------------------------------

class RepresentativeSelector:
    """Selects the longest NucSequence as the cluster representative."""

    @staticmethod
    def select(seqs: list[NucSequence]) -> NucSequence:
        return max(seqs, key=lambda s: s.length)


# ---------------------------------------------------------------------------
# BlastnAligner
# ---------------------------------------------------------------------------

class BlastnAligner:
    """
    Aligns every cluster member against the representative using BLASTn
    (-subject mode; no database needed).

    Parameters
    ----------
    blastn_bin    : path/command for blastn (default "blastn").
    min_identity  : minimum % identity passed to -perc_identity (default 90).
    """

    _OUTFMT = (
        "6 qseqid sseqid pident length qlen slen "
        "qstart qend sstart send evalue bitscore"
    )

    def __init__(
        self,
        blastn_bin:   str   = "blastn",
        min_identity: float = 90.0,
        min_coverage: float = 50.0,
    ):
        self.blastn_bin   = blastn_bin
        self.min_identity = min_identity
        self.min_coverage = min_coverage

    # --- public ---

    def align(
        self,
        representative: NucSequence,
        members: list[NucSequence],
    ) -> list[ClusterMember]:
        """
        Returns qualifying ClusterMembers (representative first).
        Members with query_coverage < min_coverage are silently dropped.
        """
        result: list[ClusterMember] = [
            ClusterMember(
                seq_id=representative.id,
                length=representative.length,
                is_representative=True,
                identity=100.0,
                query_coverage=100.0,
                aln_start=1,
                aln_end=representative.length,
                query_start=1,
                query_end=representative.length,
            )
        ]

        non_rep = [m for m in members if m.id != representative.id]
        if not non_rep:
            return result

        aln_map = self._run_blastn(representative, non_rep)

        for seq in non_rep:
            hit = aln_map.get(seq.id)
            if not hit:
                logger.debug(
                    f"No BLASTn alignment: {seq.id} vs representative {representative.id}"
                )
                continue
            if hit["qcov"] < self.min_coverage:
                logger.debug(
                    f"Coverage too low ({hit['qcov']:.1f}% < {self.min_coverage}%): "
                    f"{seq.id} — excluded from cluster."
                )
                continue
            result.append(ClusterMember(
                seq_id=seq.id,
                length=seq.length,
                is_representative=False,
                identity=hit["pident"],
                query_coverage=hit["qcov"],
                aln_start=hit["sstart"],
                aln_end=hit["send"],
                query_start=hit["qstart"],
                query_end=hit["qend"],
            ))

        return result

    # --- private ---

    def _run_blastn(
        self,
        representative: NucSequence,
        members: list[NucSequence],
    ) -> dict[str, dict]:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path   = Path(tmp)
            rep_fasta  = tmp_path / "representative.fasta"
            qry_fasta  = tmp_path / "members.fasta"

            rep_fasta.write_text(
                f">{representative.id}\n{representative.sequence}\n",
                encoding="utf-8",
            )
            with open(qry_fasta, "w", encoding="utf-8") as fh:
                for seq in members:
                    fh.write(f">{seq.id}\n{seq.sequence}\n")

            cmd = [
                self.blastn_bin,
                "-query",         str(qry_fasta),
                "-subject",       str(rep_fasta),
                "-outfmt",        self._OUTFMT,
                "-task",          "blastn",
                "-perc_identity", str(self.min_identity),
                "-dust",          "no",
            ]

            try:
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, check=True
                )
            except subprocess.CalledProcessError as exc:
                logger.error(f"blastn failed: {exc.stderr[:300]}")
                return {}
            except FileNotFoundError:
                logger.error(
                    f"blastn binary not found: {self.blastn_bin!r}. "
                    "Install BLAST+ or pass the full path to BlastnAligner."
                )
                return {}

        return self._parse_blastn(proc.stdout)

    @staticmethod
    def _parse_blastn(stdout: str) -> dict[str, dict]:
        """
        Parses outfmt-6 output and returns the best hit per query
        (highest bitscore). Subject coordinates are normalised so
        aln_start <= aln_end regardless of strand.
        """
        best: dict[str, dict] = {}
        for line in stdout.splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            (qseqid, _, pident, _aln_len, qlen, _slen,
             qstart, qend, sstart, send, _evalue, bitscore) = parts

            bitscore_f = float(bitscore)
            qlen_i     = int(qlen)
            qstart_i, qend_i   = int(qstart), int(qend)
            sstart_i, send_i   = int(sstart), int(send)
            aln_len    = abs(qend_i - qstart_i) + 1
            qcov       = round(aln_len / qlen_i * 100, 2) if qlen_i > 0 else 0.0

            hit = {
                "pident":   round(float(pident), 2),
                "qcov":     qcov,
                "qstart":   min(qstart_i, qend_i),   # normalise for minus-strand
                "qend":     max(qstart_i, qend_i),
                "sstart":   min(sstart_i, send_i),
                "send":     max(sstart_i, send_i),
                "bitscore": bitscore_f,
            }
            if qseqid not in best or bitscore_f > best[qseqid]["bitscore"]:
                best[qseqid] = hit

        return best


# ---------------------------------------------------------------------------
# SequenceTracker  (main entry point)
# ---------------------------------------------------------------------------

class SequenceTracker:
    """
    Builds ViralCluster objects from a list of NucSequence.

    Pipeline per cluster
    --------------------
    1. Group by best BLASTx species (ClusterBuilder).
    2. Pick the longest sequence as representative (RepresentativeSelector).
    3. Align all members against the representative with BLASTn (BlastnAligner).
    4. Assign cluster_id to every NucSequence.

    Parameters
    ----------
    blastn_bin   : blastn command or full path (default "blastn").
    min_identity : minimum % identity for BLASTn hits (default 90.0).
    """

    def __init__(
        self,
        blastn_bin:   str   = "blastn",
        min_identity: float = 90.0,
        min_coverage: float = 50.0,
    ):
        self._aligner     = BlastnAligner(blastn_bin, min_identity, min_coverage)
        self._min_coverage = min_coverage

    def track(self, nuc_seqs: list[NucSequence]) -> list[ViralCluster]:
        """
        Clusters sequences, runs alignments, annotates NucSequence.cluster_id,
        and returns the list of ViralCluster objects.

        A cluster is only formed when ≥2 members survive the identity and
        coverage thresholds (including the representative).
        """
        groups   = ClusterBuilder.build(nuc_seqs)
        clusters: list[ViralCluster] = []
        cluster_num = 0

        for species, seqs in groups.items():
            # Need at least 2 sequences before alignment
            if len(seqs) < 2:
                logger.debug(
                    f"ClusterBuilder: skipping '{species}' — only 1 sequence."
                )
                continue

            rep     = RepresentativeSelector.select(seqs)
            members = self._aligner.align(rep, seqs)

            # Members that didn't pass coverage/identity are already dropped by
            # the aligner; require ≥2 survivors (representative + ≥1 member).
            if len(members) < 2:
                logger.debug(
                    f"ClusterBuilder: skipping '{species}' — no members survived "
                    f"the coverage/identity filter."
                )
                continue

            cluster_num += 1
            cluster_id  = f"VQ_CLU_{cluster_num:04d}"

            cluster = ViralCluster(
                cluster_id=cluster_id,
                species=species,
                representative_id=rep.id,
                members=members,
            )
            clusters.append(cluster)

            valid_ids = {m.seq_id for m in members}
            for seq in seqs:
                if seq.id in valid_ids:
                    seq.cluster_id = cluster_id

            logger.success(
                f"{cluster_id} | {species!r} | "
                f"{len(members)} member(s) | representative: {rep.id!r} ({rep.length} nt)"
            )

        return clusters
