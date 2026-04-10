
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
from .biodata import NucSequence

# This script will handle all data input and parsing of fasta files.
# I'm still learning some concepts about OOP, but my plan is aplly some of this concepts here.

class FastaParser:

    def __init__(self, file_path: str):

        self.file_path: str = file_path
        self.sequences: list[NucSequence] = []
        self.allowed_nucleotides: set[str] = set("ATCGUN-")

        logger.info("Starting fasta parser")


    def _is_nuc_only(self, record: SeqRecord) -> bool: # ensure that the sequences are composed by nucleotides

        seq_text = str(record.seq).upper()
        seq_chars = set(seq_text)

        return seq_chars.issubset(self.allowed_nucleotides)


    def read_input_file(self): # raw processing of input files

        if not os.path.exists(self.file_path):
            logger.error(f"The file {self.file_path} was not found.")
            raise ValueError("File not found.")

        records_iterator = SeqIO.parse(handle=self.file_path, format="fasta")
        
        for record in records_iterator: # read each record and check if they're valid
            if not self._is_nuc_only(record):
                logger.error(f"The sequence {record.id} contains invalid characters.")
                raise ValueError("The file contains sequences that are invalid to the analysis.")
            
            # save sequence info to Biosequence dataclass
            sequence_string = str(record.seq).upper()
            new_sequence = NucSequence(
                id=record.id,
                sequence=sequence_string,
            )

            self.sequences.append(new_sequence)

        logger.success(f"Found {len(self.sequences)} valid sequences in {self.file_path}")


    # some extra informations of the sequences to send to Biosequence
    # Length and GC
    # def get_info_seqs(self):

    





