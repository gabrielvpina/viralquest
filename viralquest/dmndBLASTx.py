import os, glob, re, subprocess
from Bio import SeqIO
import pandas as pd           

def diamond_blastx(vvFolder, database, CPU, diamond_path=None, log_file=None):

    if diamond_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        # bin/ directory
        diamond_path = os.path.join(project_root, "viralquest", "bin", "diamond")

    # get basenames
    fileFasta = os.path.join(vvFolder, "*_viralSeq.fasta")
    viral_files = glob.glob(fileFasta)

    for viral_file in viral_files:
        infile = viral_file
        sample = os.path.basename(viral_file).replace("_viralSeq.fasta", "")
        outfile = os.path.join(vvFolder,sample)

        Dblastx_input = [
            diamond_path, "blastx",
            "--db", database,
            "--query", infile,
            "--threads", str(CPU),
            "--outfmt", "6", "qseqid", "qlen", "slen", "qcovhsp", "pident", "evalue", "stitle",
            "--max-target-seqs", "1",
            "--out", f"{outfile}_blastx.tsv"
        ]

        try:
            # Executa o comando com subprocess.run
            subprocess.run(Dblastx_input, check=True, capture_output=True)

        except subprocess.CalledProcessError as e:
            # Captura e exibe informações sobre o erro
            print(f"Error: Diamond BLASTx failed with return code {e.returncode}.")
            print(f"Command output: {e.stdout}")
            print(f"Command error: {e.stderr}")

    os.remove(diamond_path) if os.path.exists(diamond_path) else None


# extract species name
def extract_organism_name(stitle):
    if stitle == "no_hits":
        return ""
    match = re.search(r'\[([^\]]+)\]', stitle)
    if match:
        return match.group(1)
    return ""


def generateFasta_blastn(vvFolder):

  # get name
  name_files = glob.glob(os.path.join(vvFolder,  "*_pfam_ORFs.fasta"))
  if not name_files:
        raise FileNotFoundError("No file with '_pfam_ORFs.fasta' suffix.")
  name_file = name_files[0]

  output_file = name_file.replace("_pfam_ORFs.fasta","_viral.fa")

  # get table
  blastx_file = glob.glob(os.path.join(vvFolder, "*_blastx.tsv"))
  if not blastx_file:
    raise FileNotFoundError("No file with '_blastx.tsv' suffix.")

  blastx_table_file = blastx_file[0]

  fasta_files = glob.glob(os.path.join(vvFolder, "*_viralSeq.fasta"))
  if not fasta_files:
    raise FileNotFoundError("No file with '_vq.fasta' suffix.")

  fasta_file = fasta_files[0]

  bxTable = pd.read_csv(blastx_table_file, sep="\t", header=None)
  bxTable_viral = bxTable[bxTable[6].str.contains(r'\[.*(?:virus[-\s]?|phage|riboviria).*\]', case=False, na=False)]
  query_names = set(bxTable_viral[0])

  def filter_seqs(input_fasta, output_fasta, queries):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in queries:
                SeqIO.write(record, out, "fasta")

  filter_seqs(fasta_file, output_file, query_names)