import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum

import pyhmmer
from loguru import logger

from viralquest.biodata import NucSequence, Orf, HmmDomain


class HmmRole(Enum):
    """
    FILTER      — RVDB, Vfam, EggNOG.
                  A hit here flags the sequence as a viral candidate.
    CHARACTERIZE — Pfam only.
                  Run after viral identity is confirmed; adds functional
                  domain annotation to the JSON report.
    """
    FILTER      = "filter"
    CHARACTERIZE = "characterize"


# Which database name belongs to which role.
HMM_ROLES: dict[str, HmmRole] = {
    "RVDB":   HmmRole.FILTER,
    "Vfam":   HmmRole.FILTER,
    "EggNOG": HmmRole.FILTER,
    "Pfam":   HmmRole.CHARACTERIZE,
}


class HmmMetadataLoader:
    """Loads and normalises JSON metadata for a given HMM database."""

    _DB_KEYS: dict[str, tuple] = {
        "EggNOG": ("EggNOG_TargetID", "EggNOG_Description", "", ""),
        "Pfam":   ("Pfam_TargetID",   "Pfam_Description",   "Pfam_Type", "Pfam_Details"),
        "RVDB":   ("RVDB_TargetID",   "RVDB_Description",   "", ""),
        "Vfam":   ("Vfam_TargetID",   "Vfam_Description",   "", ""),
    }

    @classmethod
    def load(cls, json_path: str, db_name: str) -> dict[str, dict]:
        keys = cls._DB_KEYS.get(db_name)
        if not keys:
            logger.warning(f"Unknown HMM database '{db_name}'; returning empty metadata.")
            return {}
        id_key, desc_key, type_key, details_key = keys
        try:
            with open(json_path, encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception as exc:
            logger.error(f"Cannot read {json_path}: {exc}")
            return {}

        result = {
            item[id_key]: {
                "description": item.get(desc_key, ""),
                "type":        item.get(type_key, "")    if type_key    else "",
                "details":     item.get(details_key, "") if details_key else "",
            }
            for item in data
            if item.get(id_key)
        }
        logger.success(f"{len(result)} metadata records loaded for {db_name}.")
        return result


class HmmSequencePreparer:
    """Converts ORF amino-acid sequences into a pyhmmer DigitalSequenceBlock."""

    @staticmethod
    def prepare(
        nuc_seqs: list[NucSequence],
    ) -> tuple[pyhmmer.easel.DigitalSequenceBlock | None, dict[str, Orf]]:
        alphabet  = pyhmmer.easel.Alphabet.amino()
        text_seqs: list[pyhmmer.easel.TextSequence] = []
        orf_map:   dict[str, Orf] = {}

        for seq in nuc_seqs:
            for orf in seq.orfs:
                if not orf.aa_sequence:
                    continue
                orf_map[orf.name] = orf
                text_seqs.append(
                    pyhmmer.easel.TextSequence(
                        name=orf.name.encode(),
                        sequence=orf.aa_sequence,
                    )
                )

        if not text_seqs:
            return None, orf_map

        digital_seqs = [s.digitize(alphabet) for s in text_seqs]
        return pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_seqs), orf_map


class HmmSearcher:
    """
    Runs hmmsearch against one HMM file.
    Stateless — safe to call from multiple threads (pyhmmer releases the GIL).
    """

    def __init__(self, cpus: int = 0, score_threshold: float = 50.0):
        self.cpus = cpus
        self.score_threshold = score_threshold

    def search(
        self,
        seq_block: pyhmmer.easel.DigitalSequenceBlock,
        hmm_path: str,
    ) -> list[tuple]:
        """
        Returns raw hit tuples:
            (query_name: str, target_name: str, score, i_evalue, env_from, env_to)
        target_name is always str (normalised from bytes if needed).
        """
        hits: list[tuple] = []
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            for top_hits in pyhmmer.hmmsearch(hmm_file, seq_block, cpus=self.cpus):
                raw_name   = top_hits.query.name
                query_name = raw_name.decode("utf-8") if isinstance(raw_name, bytes) else raw_name

                for hit in top_hits:
                    if not hit.included:
                        continue
                    raw_target  = hit.name
                    target_name = raw_target.decode("utf-8") if isinstance(raw_target, bytes) else raw_target

                    for domain in hit.domains:
                        if domain.score < self.score_threshold:
                            continue
                        hits.append((
                            query_name,
                            target_name,
                            round(domain.score, 2),
                            domain.i_evalue,
                            domain.env_from,
                            domain.env_to,
                        ))
        return hits

    def search_parallel(
        self,
        seq_block: pyhmmer.easel.DigitalSequenceBlock,
        hmm_paths: list[str],
        max_workers: int = 4,
    ) -> dict[str, list[tuple]]:
        """Searches multiple HMM files concurrently. Returns {hmm_path: hits}."""
        results: dict[str, list[tuple]] = {}
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(self.search, seq_block, path): path
                for path in hmm_paths
            }
            for future in as_completed(futures):
                path = futures[future]
                try:
                    results[path] = future.result()
                    logger.success(f"HMM done: {path} ({len(results[path])} hits)")
                except Exception as exc:
                    logger.error(f"HMM failed: {path}: {exc}")
                    results[path] = []
        return results


class HmmResultAttacher:
    """Attaches raw hit tuples onto ORF objects and optionally flags viral sequences."""

    @staticmethod
    def attach(
        hits: list[tuple],
        orf_map: dict[str, Orf],
        metadata: dict[str, dict],
        db_name: str,
    ) -> int:
        attached = 0
        for query_name, target_name, score, e_value, env_from, env_to in hits:
            orf = orf_map.get(target_name)
            if not orf:
                continue
            meta = metadata.get(query_name, {})
            orf.domains.append(HmmDomain(
                database=db_name,
                target=query_name,
                score=score,
                e_value=e_value,
                start=env_from,
                stop=env_to,
                length=env_to - env_from,
                description=meta.get("description", ""),
                type=meta.get("type", ""),
                details=meta.get("details", ""),
            ))
            attached += 1
        return attached


class HmmViralFlagSetter:
    """
    After FILTER-role HMM searches, walks all sequences and marks those
    whose ORFs have at least one domain from a filter database as is_viral=True.
    """

    FILTER_DATABASES = {"RVDB", "Vfam", "EggNOG"}

    @classmethod
    def flag(cls, nuc_seqs: list[NucSequence]) -> int:
        """Returns count of newly flagged sequences"""
        newly_flagged = 0
        for seq in nuc_seqs:
            if seq.is_viral:
                continue   # already confirmed by BLASTx
            has_viral_domain = any(
                d.database in cls.FILTER_DATABASES
                for orf in seq.orfs
                for d in orf.domains
            )
            if has_viral_domain:
                seq.is_viral = True
                newly_flagged += 1
        return newly_flagged