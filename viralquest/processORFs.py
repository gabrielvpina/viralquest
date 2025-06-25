import os, glob, subprocess
from collections import defaultdict
from Bio import SeqIO


def findorf(vvFolder):
  """
  run orfipy to find all possible ORFs in the sequences
  """

  fileFasta = os.path.join(vvFolder, "*_vq.fasta")
  viral_files = glob.glob(fileFasta)
  viral_file = viral_files[0]

  entrada = viral_file
  nome = os.path.basename(viral_file).replace("_vq.fasta", "")
  saida = vvFolder

  # Executa orfipy - genetic code 01
  orfipy = ["orfipy","--partial-3","--partial-5","--table", "1", "--outdir", saida, entrada, "--pep", f"{nome}_all_ORFs.fasta"]

  subprocess.run(orfipy, check=True, capture_output=True)

  # remove log files
  log_files = [os.path.join(vvFolder, file) for file in os.listdir(vvFolder) if file.endswith(".log")]
  for log_file in log_files:
    os.remove(log_file)


def ORFs(vvFolder, numORF):
    """
    select the biggest non-overlapping ORFs 
    """
    orfs_dict = defaultdict(list)

    # get fasta
    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_all_ORFs.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)

    name_path = input_fasta
    name = os.path.basename(name_path).replace("_all_ORFs.fasta","")

    # roder y contig
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.description
        contig = header.split(".")[0]  # name of contig
        orfs_dict[contig].append(record)

    def orfs_do_not_overlap(orf1, orf2):
        start1, end1 = sorted([int(orf1.description.split("[")[1].split("]")[0].split("-")[0]),
                               int(orf1.description.split("[")[1].split("]")[0].split("-")[1])])
        start2, end2 = sorted([int(orf2.description.split("[")[1].split("]")[0].split("-")[0]),
                               int(orf2.description.split("[")[1].split("]")[0].split("-")[1])])
        return end1 <= start2 or end2 <= start1

    biggest_orfs = []

    # contigs
    for contig, orfs in orfs_dict.items():
        # order by length
        orfs.sort(key=lambda x: len(x.seq), reverse=True)
        selected_orfs = []

        # select biggest orfs
        for orf in orfs:
            if len(selected_orfs) >= int(numORF):
                break
            if len(orf.seq) > 50 and all(orfs_do_not_overlap(orf, selected) for selected in selected_orfs):
                selected_orfs.append(orf)
            # select only ORFs w/ more than 50 AA
        biggest_orfs.extend(selected_orfs)

    # Gravar novas ORFs em um arquivo FASTA
    output_fasta = os.path.join(vvFolder, f"{name}_biggest_ORFs.fasta")
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(biggest_orfs, output_handle, "fasta")

    numberORFs = len(biggest_orfs)

    return numberORFs



def countBiggestORFsFasta(vvFolder):

    input_fasta = None
    for fasta in os.listdir(vvFolder):
        if fasta.endswith("_biggest_ORFs.fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
            break

    number_ORFsSeqs = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))

    return number_ORFsSeqs 