import os, glob, subprocess, time, lzma
from Bio import SeqIO


def blastn(vvFolder, database, CPU, blastn_path=None):

    if blastn_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        # bin/ directory
        blastn_path = os.path.join(project_root, "viralquest", "bin", "blastn")
        blastn_compressed = os.path.join(project_root, "viralquest", "bin", "blastn.xz")

        # Check if compressed version exists and extract it
        if os.path.exists(blastn_compressed) and not os.path.exists(blastn_path):
            try:
                with lzma.open(blastn_compressed, 'rb') as compressed_file:
                    with open(blastn_path, 'wb') as extracted_file:
                        extracted_file.write(compressed_file.read())
            except Exception as e:
                raise Exception(f"Failed to extract {blastn_compressed}: {e}")
        
        # Check if binary exists (either was already there or just extracted)
        if not os.path.exists(blastn_path):
            raise FileNotFoundError(f"BLASTn binary not found at: {blastn_path}")

    # run as executable
    if not os.access(blastn_path, os.X_OK):
        os.chmod(blastn_path, 0o755)


    fileFasta = os.path.join(vvFolder, "*_viral.fa")
    viral_files = glob.glob(fileFasta)

    for viral_file in viral_files:
        infile = viral_file
        sample = os.path.basename(viral_file).replace("_viral.fa", "")
        outfile = os.path.join(vvFolder,sample)

        blastn_input = [blastn_path,"-query",infile,"-db",database,"-out",f"{outfile}_blastn.tsv","-outfmt",
        "6 qseqid qlen slen qcovs pident evalue stitle","-num_threads",str(CPU),"-max_target_seqs","1"]

        subprocess.run(blastn_input, check=True, capture_output=True)

    os.remove(blastn_path) if os.path.exists(blastn_path) else None
  


import warnings
from Bio import BiopythonWarning

def blastn_online(vvFolder, database, email):
    from Bio.Blast import NCBIXML, NCBIWWW
    from Bio import SeqIO
    import os
    import glob
    import time

    # define NCBI e-mail
    NCBIWWW.email = email

    fileFasta = os.path.join(vvFolder, "*_viral.fa")
    viral_files = glob.glob(fileFasta)

    for viral_file in viral_files:
        infile = viral_file
        sample = os.path.basename(viral_file).replace("_viral.fa", "")
        outfile = os.path.join(vvFolder, f"{sample}_blastn.tsv")

        if not os.path.exists(infile):
            continue

        try:
            records = list(SeqIO.parse(infile, "fasta"))
        except Exception as e:
            continue

        if not records:
            continue

        with open(outfile, "w") as output_file:
            for record in records:
                try:
                    # supress BiopythonWarning for qblast
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", BiopythonWarning)
                        result_handle = NCBIWWW.qblast(
                            program="blastn",
                            database=database,
                            sequence=str(record.seq),
                            format_type="XML",
                            hitlist_size=1
                        )

                    blast_records = NCBIXML.read(result_handle)
                    result_handle.close()

                    if blast_records.alignments:
                        alignment = blast_records.alignments[0]
                        hsp = alignment.hsps[0]

                        query_coverage = (hsp.align_length / blast_records.query_length) * 100
                        percent_identity = (hsp.identities / hsp.align_length) * 100

                        output_line = (
                            f"{record.id}\t"
                            f"{blast_records.query_length}\t"
                            f"{alignment.length}\t"
                            f"{query_coverage:.2f}\t"
                            f"{percent_identity:.2f}\t"
                            f"{hsp.expect}\t"
                            f"{alignment.title}\n"
                        )

                        output_file.write(output_line)
                    else:
                        # line w/ NaN
                        output_line = (
                            f"{record.id}\t"
                            f"{len(record.seq)}\t"
                            f"NaN\t"
                            f"NaN\t"
                            f"NaN\t"
                            f"NaN\t"
                            f"NaN\n"
                        )
                        output_file.write(output_line)

                    time.sleep(1) # NCBI limit rate

                except Exception as e:
                    # write row with errors
                    output_file.write(f"{record.id}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")
                    continue
