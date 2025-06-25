import os, lzma
import glob
import subprocess
import pandas as pd
from Bio import SeqIO

def filterBLAST(vvFolder, database, CPU, diamond_path=None, log_file=None):
    """
    Filter BLAST using diamond binary with subprocess
    """
    # diamond binary path
    if diamond_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        diamond_path = os.path.join(project_root, "viralquest", "bin", "diamond")
        diamond_compressed = os.path.join(project_root, "viralquest", "bin", "diamond.xz")
        
        # Check if compressed version exists and extract it
        if os.path.exists(diamond_compressed) and not os.path.exists(diamond_path):
            try:
                with lzma.open(diamond_compressed, 'rb') as compressed_file:
                    with open(diamond_path, 'wb') as extracted_file:
                        extracted_file.write(compressed_file.read())
            except Exception as e:
                raise Exception(f"Failed to extract {diamond_compressed}: {e}")
        
        # Check if binary exists (either was already there or just extracted)
        if not os.path.exists(diamond_path):
            raise FileNotFoundError(f"Diamond binary not found at: {diamond_path}")
        
        # Make binary executable
        if not os.access(diamond_path, os.X_OK):
            os.chmod(diamond_path, 0o755)
    
    fileFasta = os.path.join(vvFolder, "*_vq.fasta")
    viral_files = glob.glob(fileFasta)
    
    for viral_file in viral_files:
        infile = viral_file
        sample = os.path.basename(viral_file).replace("_vq.fasta", "")
        outfile = os.path.join(vvFolder, sample)
        
        Dblastx_input = [
            diamond_path, "blastx",
            "--db", database,
            "--query", infile,
            "--threads", str(CPU),
            "--outfmt", "6", "qseqid", "qlen", "slen", "qcovhsp", "pident", "evalue", "stitle",
            "--max-target-seqs", "1",
            "--out", f"{outfile}_ref.tsv"
        ]

        try:
            # Executa o comando com subprocess.run
            subprocess.run(Dblastx_input, check=True, capture_output=True)

        except subprocess.CalledProcessError as e:
            # Captura e exibe informações sobre o erro
            print(f"Error: Diamond BLASTx failed with return code {e.returncode}.")
            print(f"Command output: {e.stdout}")
            print(f"Command error: {e.stderr}")
        
    # Process results (same as original)
    fasta_files = glob.glob(os.path.join(vvFolder, "*_vq.fasta"))
    if not fasta_files:
        raise FileNotFoundError("No file with '_vq.fasta' suffix.")
    fasta_file = fasta_files[0]
    
    hmm_table_file = f"{outfile}_ref.tsv"
    
    # Read and process results
    hmmTable = pd.read_csv(hmm_table_file, sep="\t", low_memory=False)
    hmmTable.columns = ["QueryID", "Querylength", "Subjlength", "Cover", "Identity", "Evalue", "SubjTitle"]
    
    query_names = set(hmmTable['QueryID'])
    output1 = f"{outfile}_ref.csv"
    hmmTable.to_csv(output1, index=False)
    
    output_file = f"{outfile}_filtered.fasta"
    
    def filter_seqs(input_fasta, output_fasta, queries):
        with open(output_fasta, "w") as out:
            for record in SeqIO.parse(input_fasta, "fasta"):
                if record.id in queries:
                    SeqIO.write(record, out, "fasta")
    
    filter_seqs(fasta_file, output_file, query_names)