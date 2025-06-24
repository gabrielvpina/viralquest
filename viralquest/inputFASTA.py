import os, subprocess, shutil, shutil
from Bio import SeqIO
import pkg_resources
import os

def get_binary_path(binary_name):
    """Get the path to a binary file included in the package"""
    try:
        # Try to get the path from the installed package
        return pkg_resources.resource_filename('viralquest', f'bin/{binary_name}')
    except:
        # Fallback to relative path for development
        script_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(os.path.dirname(script_dir), 'bin', binary_name)

def validateFasta(inputContig):
    """
    check the fasta file - only full nucleotide FASTA 
    """
    try:

        records = SeqIO.parse(inputContig, "fasta")

        # empty file
        first_record = next(records, None)
        if first_record is None:
            print("\nError: The Fasta file is empty or don't have valid content.")
            return False

        # validate first record
        if not all(c in "ACGTN" for c in first_record.seq.upper()):
            print(f"\nErro: The sequence '{first_record.id}' have invalid characters.")
            return False

        # validate all records
        for record in records:
            if not all(c in "ACGTN" for c in record.seq.upper()):
                print(f"\nError: The sequence '{record.id}' have invalid characters.")
                return False
        return True

    except FileNotFoundError:
        print("Error: Fasta file not found.")
        return False
    except Exception as e:
        print(f"Error to process FASTA file: {e}")
        return False



def cap3(inputContig, vvFolder):

    # directory w/ results
    if not os.path.exists(vvFolder):
      os.makedirs(vvFolder)

    shutil.copy(os.path.normpath(inputContig), os.path.normpath(vvFolder))

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    cap3_path = os.path.join(project_root, "viralquest", "bin", "cap3")

    # check if exists
    if not os.path.exists(cap3_path):
        raise FileNotFoundError(f"cap3 binary not found at: {cap3_path}")
    
    # binary is executable?
    if not os.access(cap3_path, os.X_OK):
        os.chmod(cap3_path, 0o755)

    for cpContig in os.listdir(vvFolder):
      if cpContig.endswith(".fasta" or ".fa"):
        # run CAP3
        cap3_command = [cap3_path, os.path.join(vvFolder,cpContig)]
        subprocess.run(cap3_command, check=True, capture_output=True)
      if not os.path.exists(vvFolder):
        print("####### Fasta file not found")

    for contigs in os.listdir(vvFolder):
        if contigs.endswith(".cap.contigs"):
            sample = contigs.replace(".fasta.cap.contigs","")
            #print(f"####### Assembling sample {sample} with CAP3")

            for singlets in os.listdir(vvFolder):
                if singlets.endswith(".cap.singlets"):

                    cat_command = f"cat {os.path.join(vvFolder,contigs)} {os.path.join(vvFolder,singlets)} >> {os.path.join(vvFolder, f'{sample}_vq.fasta')}"
                    #print(cat_command)
                    subprocess.run(cat_command, shell=True, check=True, capture_output=True)


    log_files = [os.path.join(vvFolder, file) for file in os.listdir(vvFolder)
                 if file.endswith((".links",".qual",".info",".contigs",".singlets",".ace"))]
    for log_file in log_files:
        os.remove(log_file)
    # remove original fasta file
    os.remove(os.path.join(vvFolder,os.path.basename(inputContig)))



def noCAP3(inputContig, vvFolder):
    # directory w/ results
    if not os.path.exists(vvFolder):
      os.makedirs(vvFolder)

    # copy original fasta in result directory
    #copy_command = ["cp", os.path.normpath(inputContig), os.path.normpath(vvFolder)]
    #subprocess.run(copy_command, check=True, capture_output=True)
    infile = os.path.normpath(inputContig)
    namebase = os.path.basename(inputContig).replace(".fasta","")
    outfile = os.path.join(vvFolder,f"{namebase}_vq.fasta")

    shutil.copy(infile, outfile)