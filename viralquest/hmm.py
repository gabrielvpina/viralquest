# viralquest/hmm.py
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import pyhmmer
from loguru import logger
from .biodata import NucSequence, Orf, HmmDomain


class HmmMetadataLoader:
    """Loads and normalises JSON metadata for a given HMM database."""

    _DB_KEYS = {
        "EggNOG": ("EggNOG_TargetID", "EggNOG_Description", "", ""),
        "Pfam":   ("Pfam_TargetID",   "Pfam_Description",   "Pfam_Type", "Pfam_Details"),
        "RVDB":   ("RVDB_TargetID",   "RVDB_Description",   "", ""),
        "Vfam":   ("Vfam_TargetID",   "Vfam_Description",   "", ""),
    }

    @classmethod
    def load(cls, json_path: str, db_name: str) -> dict[str, dict]:
        keys = cls._DB_KEYS.get(db_name)
        if not keys:
            logger.warning(f"Unknown database '{db_name}'; returning empty metadata.")
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
                "type":        item.get(type_key, "") if type_key else "",
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
    ) -> tuple[pyhmmer.easel.DigitalSequenceBlock | None, dict[bytes, Orf]]:
        alphabet = pyhmmer.easel.Alphabet.amino()
        text_seqs, orf_map = [], {}

        for seq in nuc_seqs:
            for orf in seq.orfs:
                if not orf.aa_sequence:
                    continue
                # name_bytes = orf.name.encode()
                orf_map[orf.name] = orf
                text_seqs.append(
                    pyhmmer.easel.TextSequence(name=orf.name.encode(), sequence=orf.aa_sequence)
                )

        if not text_seqs:
            return None, orf_map

        digital_seqs = [seq.digitize(alphabet) for seq in text_seqs]
        return pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_seqs), orf_map


class HmmSearcher:
    """
    Runs hmmsearch against one HMM file.

    This class is intentionally stateless so its search() method
    can be called safely from a worker process.
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
        Returns raw hit tuples (query_name, target_name, score, i_evalue,
        env_from, env_to) for every passing domain.
        No ORF objects touched here — safe to pickle for multiprocessing.
        """
        hits = []
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            for top_hits in pyhmmer.hmmsearch(hmm_file, seq_block, cpus=self.cpus):
                name = top_hits.query.name
                query_name = name.decode('utf-8') if isinstance(name, bytes) else name
                for hit in top_hits:
                    if not hit.included:
                        continue
                    for domain in hit.domains:
                        if domain.score < self.score_threshold:
                            continue
                        hits.append((
                            query_name,
                            hit.name,          # bytes
                            round(domain.score, 2),
                            domain.i_evalue,
                            domain.env_from,
                            domain.env_to,
                        ))
        return hits

    def search_parallel(self, seq_block, hmm_paths, max_workers=4):
        results = {}
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(self.search, seq_block, path): path
                for path in hmm_paths
            }
            for future in as_completed(futures):
                path = futures[future]
                try:
                    results[path] = future.result()
                    logger.success(f"Done: {path} ({len(results[path])} hits)")
                except Exception as exc:
                    logger.error(f"Search failed for {path}: {exc}")
                    results[path] = []
        return results


class HmmResultAttacher:
    """Attaches raw hit tuples back onto ORF objects using a metadata dict."""

    @staticmethod
    def attach(
        hits: list[tuple],
        orf_map: dict[bytes, Orf],
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