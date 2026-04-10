import orfipy_core
from Bio.Seq import Seq
from loguru import logger
from .biodata import NucSequence, Orf

class OrfAnalyzer:
    """Find, filter and manipulate ORFs in valid sequences."""
    

    def __init__(self, min_len_nt: int = 150):
        self.min_len_nt = min_len_nt
        logger.info(f"Starting ORF serch (Minimal length: {self.min_len_nt} nt)")



    def _parse_orfipy_description(self, description: str) -> dict:

        desc_dict = {}
        for item in description.split(';'):
            if '=' in item:
                k, v = item.split('=')
                desc_dict[k] = v
            elif ':' in item:
                k, v = item.split(':')
                desc_dict[k] = v
        return desc_dict



    def find_n_save_orfs(self, nuc_seq: NucSequence) -> None:
        """uses orfipy to find and index the ORFs in NucSequence"""
        
        logger.debug(f"Searching ORFs in sequence {nuc_seq.id}...")
        
        seq_str = nuc_seq.sequence
        
        results = orfipy_core.orfs(seq_str, minlen=self.min_len_nt)
        
        orfs_number = 0
        
        for start, stop, strand, description in results:

            metadata = self._parse_orfipy_description(description)
            
            bio_seq = Seq(seq_str[start:stop])
            
            # reverse complement to negative strands 
            if strand == '-':
                bio_seq = bio_seq.reverse_complement()
                
            nuc_sequence = str(bio_seq)
            
            # translate to aminoacids  
            try:
                aa_sequence = str(bio_seq.translate())
            except Exception as e:
                logger.warning(f"Translate ORF from {start}-{stop}: {e}")
                aa_sequence = ""
            
            # Orf dataclass
            nova_orf = Orf(
                start_codon=metadata.get('Start', ''),
                stop_codon=metadata.get('Stop', ''),
                start_position=start,
                stop_position=stop,
                strand=strand,
                frame=int(metadata.get('ORF_frame', 0)),
                length_nt=len(nuc_sequence),
                length_aa=len(aa_sequence),
                bigger_than_50=len(aa_sequence) > 50,
                aa_sequence=aa_sequence,
                nuc_sequence=nuc_sequence
            )
            

            nuc_seq.orfs.append(nova_orf)
            orfs_number += 1
            
        logger.info(f"Success: {orfs_number} ORFs found in {nuc_seq.id}.")

    # some future methods
    #
    #
    #
    
    def get_longest_orf(self, nuc_seq: NucSequence) -> Orf | None:
        """returns biggest ORF"""
        if not nuc_seq.orfs:
            return None
        return max(nuc_seq.orfs, key=lambda orf: orf.length_aa)
        
    def filter_complete_orfs(self, nuc_seq: NucSequence) -> list[Orf]:
        """returns ORFs with valid start and stop codons"""
        return [orf for orf in nuc_seq.orfs if orf.start_codon and orf.stop_codon]
