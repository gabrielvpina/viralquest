import os, shutil, shutil
from Bio import SeqIO
import pandas as pd

def copyfasta(vvFolder):

    input_fasta = None
    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_vq.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
            break

    # copy original fasta
    oldFasta = input_fasta.replace("_vq.fasta","_orig.fasta")
    shutil.copy(input_fasta, oldFasta)    

    
def countOriginalFasta(vvFolder):

    input_fasta = None
    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_orig.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
            break

    number_originalSeqs = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))

    return number_originalSeqs

def filterfasta(vvFolder):

    min_length = 500
    input_fasta = None
    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_vq.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
            break
    
    if input_fasta is None:
        print("Fasta file not found")
        return None, None

    temp_fasta = input_fasta + ".temp"

    filtered_seqs = [seq for seq in SeqIO.parse(input_fasta, "fasta") if len(seq.seq) >= min_length]
    number_filteredSeqs = len(filtered_seqs)

    # if any of the sequences is bigger than 500 - the function filter by length
    if number_filteredSeqs >= 0:
        with open(temp_fasta, "w") as output_handle:
            SeqIO.write(filtered_seqs, output_handle, "fasta")
        # replace original fasta
        os.replace(temp_fasta, input_fasta)

    return number_filteredSeqs

    

def renameFasta(vvFolder):
    # it's an important step when all the JSON files are merged to track down the source of each sequence

    fastaIndex_table = []

    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_vq.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
        
    name = os.path.basename(input_fasta.replace("_vq.fasta",""))

    mod_headers = []

    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
        original_id = record.id
        new_id = f"{name}_seq{i:02d}"
        record.id = new_id
        record.description = new_id

        mod_headers.append(record)
        fastaIndex_table.append([original_id, new_id])
    # fasta file
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(mod_headers, output_handle, "fasta")

    # original names table
    fastaIndex_table = pd.DataFrame(fastaIndex_table, columns=["Source", "QueryID"])

    return fastaIndex_table