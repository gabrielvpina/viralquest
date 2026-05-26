import os
import shutil
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
from .biodata import NucSequence, InputFasta, Cap3

# This script will handle all data input and parsing of fasta files.
# I'm still learning some concepts about OOP, but my plan is aplly some of this concepts here.

class FastaParser:

    def __init__(self, file_path: str):

        self.file_path: Path = Path(file_path)
        self.sequences: list[NucSequence] = []
        self.allowed_nucleotides: set[str] = set("ATCGUN-")
        self.input_fasta: InputFasta | None = None

        logger.info("Starting fasta parser")


    def _is_nuc_only(self, record: SeqRecord) -> bool: # ensure that the sequences are composed by nucleotides

        seq_text = str(record.seq).upper()
        seq_chars = set(seq_text)

        return seq_chars.issubset(self.allowed_nucleotides)


    def read_input_file(self): # raw processing of input files

        if not os.path.exists(self.file_path):
            logger.error(f"The file {self.file_path} was not found.")
            raise ValueError("File not found.")


        if not os.access(self.file_path, os.R_OK):
            logger.error(f"Permission denied to read the file: {self.file_path}")
            raise PermissionError("Permission denied.")

        try:
            records_iterator = SeqIO.parse(handle=self.file_path, format="fasta")
            
            for record in records_iterator: # read each record and check if they're valid
                if not self._is_nuc_only(record):
                    logger.error(f"The sequence {record.id} contains invalid characters.")
                    raise ValueError("The file contains sequences that are invalid to the analysis.")
                
                # save sequence info to NucSequence dataclass
                sequence_string = str(record.seq).upper()
                new_sequence = NucSequence(
                    id=record.id,
                    sequence=sequence_string,
                )

                self.sequences.append(new_sequence) # sequences stores all information from fasta

            if not self.sequences:
                logger.warning(f"No valid FASTA sequences found in {self.file_path}. Is the file empty or in the wrong format?")
            else:
                logger.success(f"Found {len(self.sequences)} valid sequences in {self.file_path}")

            self.input_fasta = InputFasta(
                name=self.file_path.name,
                num_seqs=len(self.sequences),
                path=self.file_path,
                size=os.path.getsize(self.file_path)
            )

        except Exception as e:
            logger.exception(f"An unexpected error occurred while reading the file: {e}")
            raise


# ========================================================================
# USING CAP3 to improve some sequence length
# ========================================================================

class Cap3Runner:
    """RECOMMENDED QUERY: cap3 file.fasta -p 98 -o 100"""

    def __init__(
        self,
        file_path: str | Path,
        cap3_bin: str = "cap3",
        pident: int = 98,
        overlap: int = 100,
        outdir: str | None = None,
    ):
        self.original_fasta = Path(file_path).resolve()
        self.cap3_bin = cap3_bin

        self.pident = str(pident)
        self.overlap = str(overlap)

        self.outdir = Path(outdir).resolve() if outdir else Path(tempfile.gettempdir()) / "cap3_run"
        self.outdir.mkdir(parents=True, exist_ok=True)

        self.working_fasta = self.outdir / self.original_fasta.name



    def _setup_files(self):
        shutil.copy(self.original_fasta, self.working_fasta)



    def _run_cap3(self) -> Path:
        """execute and return log file"""
        cmd = [
            self.cap3_bin,
            str(self.working_fasta),
            "-p", self.pident,
            "-o", self.overlap
        ]
        log_path = self.working_fasta.with_suffix(".cap3.log")

        try:
            logger.info(f"Running CAP3 on raw file: '{self.original_fasta.name}'")
            logger.debug(f"Working directory copy: '{self.working_fasta}'")
            logger.debug(f"Full command executed: {' '.join(cmd)}")


            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True
            )
            log_path.write_text(result.stdout)

            logger.info(f"CAP3 executed successfully for '{self.original_fasta.name}'.")
            logger.info(f"Internal CAP3 execution log saved to: '{log_path}'")
            if "warning" in result.stdout.lower():
                logger.warning("CAP3 finished with internal warnings. Check the log file for details.")

            return log_path
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Fatal error executing CAP3 on '{self.original_fasta.name}'.")
            logger.error(f"Command line exit code: {e.returncode}")            
            if e.stderr:
                logger.error(f"CAP3 Stderr: {e.stderr.strip()}")
            if e.stdout:
                logger.debug(f"CAP3 Stdout before crash: {e.stdout.strip()}")
            raise e        





    def _process_cap3_out(self, log_path: Path) -> Cap3:
        """map result data to cap3 dataclass"""
        base_name = str(self.working_fasta)
        
        return Cap3(
            input_fasta=self.working_fasta,
            contigs=Path(f"{base_name}.cap.contigs"),
            singlets=Path(f"{base_name}.cap.singlets"),
            info=Path(f"{base_name}.cap.info"),
            ace=Path(f"{base_name}.cap.ace"),
            log=log_path
        )




    def cap3_runner(self) -> Cap3:
        """main method to invoke cap3"""
        self._setup_files()
        log_path = self._run_cap3()
        return self._process_cap3_out(log_path)











    





