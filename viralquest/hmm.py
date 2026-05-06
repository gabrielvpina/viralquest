import json
import pyhmmer
from loguru import logger
from typing import Dict, Any, Tuple
from .biodata import NucSequence, Orf, HmmDomain

class HmmMetadataLoader:
    """
    class that reads JSON hmm metadata
    """

    @staticmethod
    def load(json_path: str, db_name: str) -> Dict[str, Dict[str, str]]:
        """
        read json file and returns a dict 
        """
        logger.info(f"Loading metadata to {db_name}")
        
        try:
            with open(json_path, 'r', encoding='utf-8') as file:
                data = json.load(file)
        except Exception as e:
            logger.error(f"Error to read JSON file: {e}")
            return {}

        metadata_dict = {}

        # mapping to get a standart dictionary
        for item in data:
            target_id = ""
            description = ""
            type = ""
            details = ""

            if db_name == "EggNOG":
                target_id = item.get("EggNOG_TargetID", "")
                description = item.get("EggNOG_Description", "")
                
            elif db_name == "Pfam":
                target_id = item.get("Pfam_TargetID", "")
                description = item.get("Pfam_Description", "")
                type = item.get("Pfam_Type", "")
                details = item.get("Pfam_Details", "")
                
            elif db_name == "RVDB":
                target_id = item.get("RVDB_TargetID", "")
                description = item.get("RVDB_Description", "")
                
            elif db_name == "Vfam":
                target_id = item.get("Vfam_TargetID", "")
                description = item.get("Vfam_Description", "")

            # if find a valid, saves it
            if target_id:
                metadata_dict[target_id] = {
                    "description": description,
                    "details": details,
                    "type": type
                }

        logger.success(f"{len(metadata_dict)} metadata records saved for {db_name}.")
        return metadata_dict


class HmmSearcher:
    """
    starts all hmm searchs using pyhhmer.
    read the memory saved ORFs 
    """

    def __init__(self, cpus: int = 0, score_threshold: float = 50.0):
        self.cpus = cpus
        self.score_threshold = score_threshold
        logger.info(f"Hmmsearche started (CPUs: {self.cpus}, Minimal Score: {self.score_threshold})")

    def _prepare_sequences(self, nuc_seqs: list[NucSequence]) -> Tuple[pyhmmer.easel.DigitalSequenceBlock, dict]:
        """
        convert aa sequences of ORFs in 'Easel Digital Sequences' to pyhmmer.
        returns digital block and a dict mapping the ids back to the ORF object.
        """
        alphabet = pyhmmer.easel.Alphabet.amino()
        text_sequences = []
        orf_map = {} 

        for seq in nuc_seqs:
            for i, orf in enumerate(seq.orfs):
                if not orf.aa_sequence:
                    continue 

                # creates uique name in bytes required by pyhmmer
                unique_name = f"{seq.id}_ORF_{i}".encode('utf-8')
                orf_map[unique_name] = orf 
                
                easel_seq = pyhmmer.easel.TextSequence(name=unique_name, sequence=orf.aa_sequence)
                text_sequences.append(easel_seq)

        if text_sequences:
            seq_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, text_sequences)
            return seq_block, orf_map
        
        return None, orf_map

    def search_domains(self, nuc_seqs: list[NucSequence], hmm_path: str, db_name: str, metadata_dict: dict):
        """
        Roda o hmmsearch contra o perfil fornecido e salva os hits nos objetos Orf.
        
        Args:
            nuc_seqs: Lista de objetos NucSequence.
            hmm_path: Caminho para o arquivo .hmm.
            db_name: Nome do banco ('Pfam', 'RVDB', etc.).
            metadata_dict: Dicionário padronizado gerado pelo HmmMetadataLoader.
        """
        logger.info(f"Preparando sequências para busca no {db_name}...")
        digital_seqs, orf_map = self._prepare_sequences(nuc_seqs)

        if not digital_seqs or len(digital_seqs) == 0:
            logger.warning("Nenhuma sequência de aminoácidos válida encontrada para a busca.")
            return

        logger.info(f"Iniciando busca pyhmmer usando {hmm_path}...")
        hits_count = 0
        
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            # Busca cada modelo HMM contra as nossas ORFs
            for top_hits in pyhmmer.hmmsearch(hmm_file, digital_seqs, cpus=self.cpus):
                
                # Nome do modelo HMM (ex: 'PF10417' ou '1-cysPrx_C')
                hmm_query_name = top_hits.query.name.decode('utf-8') 
                
                # Resgata os metadados do dicionário padronizado
                meta = metadata_dict.get(hmm_query_name, {})
                desc = meta.get("description", "Descrição não encontrada")
                details = meta.get("details", "")

                for hit in top_hits:
                    if not hit.included:
                        continue
                        
                    target_name = hit.name # Recupera o unique_name da ORF
                    target_orf = orf_map.get(target_name)

                    if not target_orf:
                        continue
                    
                    for domain in hit.domains:
                        if domain.score < self.score_threshold:
                            continue
                            
                        start_pos = domain.env_from
                        stop_pos = domain.env_to
                        length = stop_pos - start_pos
                        
                        # Instancia o objeto HmmDomain com todas as informações
                        new_domain = HmmDomain(
                            database=db_name,
                            target=hmm_query_name,
                            score=round(domain.score, 2),
                            e_value=domain.i_evalue,
                            start=start_pos,
                            stop=stop_pos,
                            length=length,
                            description=desc,
                            details=details
                        )

                        # Salva o domínio diretamente na lista da ORF
                        target_orf.domains.append(new_domain)
                        hits_count += 1

        logger.success(f"Busca finalizada no {db_name}. Encontrados {hits_count} domínios válidos.")

    def filter_overlapping_domains(self, nuc_seqs: list[NucSequence]):
        """
        Remove domínios sobrepostos na mesma ORF, mantendo apenas aquele com o maior score.
        Substitui toda a lógica complexa de Pandas que você usava antes.
        """
        logger.info("Filtrando domínios sobrepostos (overlaps)...")
        filtered_count = 0

        for seq in nuc_seqs:
            for orf in seq.orfs:
                if not orf.domains:
                    continue

                # Ordena os domínios da ORF com base na posição inicial (start)
                sorted_domains = sorted(orf.domains, key=lambda d: d.start)
                best_domains = []

                for current_domain in sorted_domains:
                    if not best_domains:
                        best_domains.append(current_domain)
                        continue

                    last_domain = best_domains[-1]

                    # Checa se o início do atual está antes do fim do anterior (sobreposição)
                    if current_domain.start <= last_domain.stop:
                        # Em caso de conflito, mantém o que tem maior score
                        if current_domain.score > last_domain.score:
                            best_domains[-1] = current_domain
                            filtered_count += 1
                        else:
                            filtered_count += 1 # Descarta o atual
                    else:
                        best_domains.append(current_domain)

                # Atualiza a lista da ORF apenas com os melhores domínios resolvidos
                orf.domains = best_domains
                
        logger.info(f"Removidos {filtered_count} domínios redundantes ou sobrepostos.")