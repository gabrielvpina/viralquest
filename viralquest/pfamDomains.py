import os, glob
from Bio import SeqIO
import pandas as pd
import collections
import pyhmmer

from .data import pfam_metadata
from .data import pfam_complete_desc


ResultHMM = collections.namedtuple("Result", ["query", "subjID", "bitscore", "start", "end"])


def pfamFasta(vvFolder):

  rvdb_files = glob.glob(os.path.join(vvFolder, "*_RVDB.csv"))
  vfam_files = glob.glob(os.path.join(vvFolder, "*_Vfam.csv"))
  egg_files = glob.glob(os.path.join(vvFolder, "*_EggNOG.csv"))


  common_key = "QueryID"  

  rvdb_table = pd.read_csv(rvdb_files[0]) if rvdb_files else None
  vfam_table = pd.read_csv(vfam_files[0]) if vfam_files else None
  egg_table = pd.read_csv(egg_files[0]) if egg_files else None

  # Full Join
  if rvdb_table is not None and vfam_table is not None and egg_table is not None:
      merged_table = rvdb_table.merge(vfam_table, on=common_key, how="outer").merge(egg_table, on=common_key, how="outer")

  elif rvdb_table is not None and vfam_table is not None:
      merged_table = rvdb_table.merge(vfam_table, on=common_key, how="outer")

  elif rvdb_table is not None and egg_table is not None:
      merged_table = rvdb_table.merge(egg_table, on=common_key, how="outer")

  elif vfam_table is not None and egg_table is not None:
      merged_table = vfam_table.merge(egg_table, on=common_key, how="outer")

  elif rvdb_table is not None:
      merged_table = rvdb_table

  elif vfam_table is not None:
      merged_table = vfam_table

  elif egg_table is not None:
      merged_table = egg_table

  else:
      raise FileNotFoundError("Tables *_RVDB.csv, *_Vfam.csv or *_EggNOG.csv not found.")


  # Unify columns
  merged_table["Query_name"] = merged_table.filter(like="Query_name").bfill(axis=1).iloc[:, 0]
  merged_table["Sample_name"] = merged_table.filter(like="Sample_name").bfill(axis=1).iloc[:, 0]
  # Remove duplicated columns
  merged_table = merged_table.loc[:, ~merged_table.columns.str.contains("_x|_y|_z")]
  # relocate
  cols = ["Sample_name", "Query_name"] + [col for col in merged_table.columns if col not in ["Query_name", "Sample_name"]]
  merged_table = merged_table.reindex(columns=cols)
  merged_table = merged_table.drop_duplicates()


  orf_files = glob.glob(os.path.join(vvFolder,  "*_biggest_ORFs.fasta"))
  orf_file = orf_files[0]


  output_path = orf_file.replace("_biggest_ORFs.fasta","_hmm.csv")
  merged_table.to_csv(output_path, index=False)

  output_file = output_path.replace("_hmm.csv","_pfam_ORFs.fasta")

  fasta_files = glob.glob(os.path.join(vvFolder, "*_biggest_ORFs.fasta"))
  if not fasta_files:
    raise FileNotFoundError("No file with '_biggest_ORFs.fasta' suffix.")

  fasta_file = fasta_files[0]


  hmmTable = pd.read_csv(output_path)
  query_names = set(hmmTable['Query_name'])

  def filter_seqs(input_fasta, output_fasta, queries):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in queries:
                SeqIO.write(record, out, "fasta")

  filter_seqs(fasta_file, output_file, query_names)


def process_Pfam(vvFolder, hmmProfile, CPU):
    orf_files = glob.glob(os.path.join(vvFolder,  "*_pfam_ORFs.fasta"))
    if not orf_files:
        raise FileNotFoundError("No file with '_pfam_ORFs.fasta' suffix.")

    orf_file = orf_files[0]
    output_tsv = orf_file.replace("_pfam_ORFs.fasta", "_hmmsearch.tsv")

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(orf_file, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)

    results = []
    with pyhmmer.plan7.HMMFile(hmmProfile) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=CPU):
            hmm_name = hits.query.name.decode()

            for hit in hits:
                if hit.included:
                    for domain in hit.domains:
                        results.append(ResultHMM(hit.name.decode(), hmm_name, hit.score, domain.env_from, domain.env_to))



    # Write in tsv
    with open(output_tsv, "w") as out:
        out.write(f"Query_name\tPfam_TargetID\tPfam_Score\tPfam_Start\tPfam_End\n")
        for result in results:
            out.write(f"{result.query}\t{result.subjID}\t{result.bitscore}\t{result.start}\t{result.end}\n")

    hmm_files = glob.glob(os.path.join(vvFolder, "*_hmmsearch.tsv"))
    if not hmm_files:
        raise FileNotFoundError("No file with '_hmmsearch.tsv' suffix.")

    hmm_table_file = hmm_files[0]

    data_parser = pd.read_csv(hmm_table_file, sep="\t", low_memory=False)

    # Extract the QueryID from the Query_name
    data_parser["QueryID"] = data_parser["Query_name"].str.extract(r"(.*?)_ORF")
    data_parser = data_parser[["QueryID", "Query_name", "Pfam_TargetID", "Pfam_Score", "Pfam_Start", "Pfam_End"]]
    data_parser["Pfam_length"] = (data_parser["Pfam_End"] - data_parser["Pfam_Start"]).round(2)

    # merge with pfam metadata
    data_parser = pd.merge(data_parser, pfam_metadata, on="Pfam_TargetID", how="left")

    data_parser = data_parser[["QueryID", "Query_name", "Pfam_TargetID", "Pfam_Description", "Pfam_Type", "Pfam_Score", "Pfam_Start", "Pfam_End", "Pfam_length"]]

    # merge with complete descriptions
    data_parser = data_parser.merge(pfam_complete_desc, left_on="Pfam_Description", right_on="name", how="left")
    data_parser = data_parser[["QueryID", "Query_name", "accession", "Pfam_TargetID", "Pfam_Description", "description", "Pfam_Type", "Pfam_Score", "Pfam_Start", "Pfam_End", "Pfam_length"]]
    data_parser.columns = ["QueryID", "Query_name", "Pfam_Accession", "Pfam_TargetID", "Pfam_Description", "Pfam_Info", "Pfam_Type", "Pfam_Score", "Pfam_Start", "Pfam_End", "Pfam_length"]

    # sample name from folder
    sample_name = os.path.basename(vvFolder)
    data_parser["Sample_name"] = sample_name

    # handle overlapping domains
    result_df = pd.DataFrame()

    # group by query_name to process each ORF
    for query_name, group in data_parser.groupby("Query_name"):
        # if only one domain for a ORF, add it directly
        if len(group) == 1:
            result_df = pd.concat([result_df, group])
            continue

        # sort by start position
        group = group.sort_values("Pfam_Start")

        # remove duplicates
        group = group.drop_duplicates()

        # check for overlaps
        overlapping_groups = []
        current_group = [0]  # first domain

        for i in range(1, len(group)):
            current_row = group.iloc[i]
            previous_rows = [group.iloc[j] for j in current_group]

            # check if domain overlaps with any domain in the current group
            overlaps = False
            for prev_row in previous_rows:
                # check for overlap
                if (current_row["Pfam_Start"] <= prev_row["Pfam_End"] and
                    current_row["Pfam_End"] >= prev_row["Pfam_Start"]):
                    overlaps = True
                    break

            if overlaps:
                # add to current overlap group
                current_group.append(i)
            else:
                # this domain doesn't overlap with current group
                # save current group and start a new one
                overlapping_groups.append(current_group)
                current_group = [i]

        # add the last group
        if current_group:
            overlapping_groups.append(current_group)

        # group of overlapping domains
        combined_rows = []
        for overlap_group in overlapping_groups:
            if len(overlap_group) == 1:
                # no overlap, just one domain
                combined_rows.append(group.iloc[overlap_group[0]])
            else:
                # for overlapping domains, keep the one with high score
                overlapping_domains = group.iloc[overlap_group]
                best_domain = overlapping_domains.loc[overlapping_domains["Pfam_Score"].idxmax()]
                combined_rows.append(best_domain)

        # if there is multiple non-overlapping domains, combine them into one row
        if len(combined_rows) > 1:
            base_row = combined_rows[0].copy()

            # new row with all domains combined
            for i, row in enumerate(combined_rows[1:], start=2):
                for col in row.index:
                    if col not in ["Sample_name", "Query_name", "QueryID"]:
                        base_row[f"{col}_{i}"] = row[col]

            result_df = pd.concat([result_df, pd.DataFrame([base_row])])
        else:
            # only one domain of overlapping domains
            result_df = pd.concat([result_df, pd.DataFrame([combined_rows[0]])])

    # reorder columns 
    base_columns = ["Sample_name", "Query_name", "QueryID", "Pfam_Accession", "Pfam_Description",
                    "Pfam_Info", "Pfam_Type", "Pfam_Score", "Pfam_Start", "Pfam_End", "Pfam_length"]

    # all columns with numeric suffixes
    all_columns = []
    for col in result_df.columns:
        if col in base_columns or any(col.startswith(base_col + "_") for base_col in base_columns):
            all_columns.append(col)

    data_parser = result_df[all_columns]


    hmmFinal = glob.glob(os.path.join(vvFolder,  "*_hmm.csv"))
    hmmPath = hmmFinal[0]
    hmmTable = pd.read_csv(hmmPath, sep=",", low_memory=False)

    # merge with others hmm table
    hmmTable = pd.merge(hmmTable, data_parser, on=["QueryID", "Query_name"], how="left")
    hmmTable["Query_name"] = hmmTable.filter(like="Query_name").bfill(axis=1).iloc[:, 0]
    hmmTable["Sample_name"] = hmmTable.filter(like="Sample_name").bfill(axis=1).iloc[:, 0]
    hmmTable = hmmTable.loc[:, ~hmmTable.columns.str.contains("_x|_y|_z")]
    hmmTable["Pfam_Score"] = hmmTable["Pfam_Score"].round(2).fillna(1)
    cols = ["Sample_name", "Query_name"] + [col for col in hmmTable.columns if col not in ["Query_name", "Sample_name"]]
    hmmTable = hmmTable.reindex(columns=cols)

    # column fullseq
    fasta_pfam = hmmPath.replace("_hmm.csv","_pfam_ORFs.fasta")
    seq_dict = {} # store sequences

    for record in SeqIO.parse(fasta_pfam, "fasta"):
      seq_dict[record.id] = str(record.seq)

    def get_sequence(query_name):
      return seq_dict.get(query_name, None)

    hmmTable['FullSequence'] = hmmTable['Query_name'].apply(get_sequence)

    hmmTable = hmmTable.sort_values(by="Pfam_Score", ascending=False).drop_duplicates(subset=["Query_name"], keep="first")

    # Write table
    csv_output_path = hmmPath
    hmmTable = hmmTable.drop_duplicates()
    hmmTable.to_csv(csv_output_path, index=False)

    os.remove(hmm_table_file)




def generateFasta(vvFolder):

  name_files = glob.glob(os.path.join(vvFolder,  "*_pfam_ORFs.fasta"))
  if not name_files:
        raise FileNotFoundError("No file with '_pfam_ORFs.fasta' suffix.")
  name_file = name_files[0]

  output_file = name_file.replace("_pfam_ORFs.fasta","_viralHMM.fasta")

  hmm_files = glob.glob(os.path.join(vvFolder, "*_hmm.csv"))
  if not hmm_files:
    raise FileNotFoundError("No file with '_hmm.csv' suffix.")

  hmm_table_file = hmm_files[0]

  fasta_files = glob.glob(os.path.join(vvFolder, "*_vq.fasta"))
  if not fasta_files:
    raise FileNotFoundError("No file with '_vq.fasta' suffix.")

  fasta_file = fasta_files[0]

  hmmTable = pd.read_csv(hmm_table_file, low_memory=False)
  query_names = set(hmmTable['QueryID'])

  def filter_seqs(input_fasta, output_fasta, queries):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in queries:
                SeqIO.write(record, out, "fasta")

  filter_seqs(fasta_file, output_file, query_names)



def mergeFASTA(vvFolder):

    name_files = glob.glob(os.path.join(vvFolder,  "*_vq.fasta"))
    name_file = name_files[0]
    name = os.path.basename(name_file).replace("_vq.fasta","")

    # input fasta
    fasta1 = os.path.join(vvFolder,f"{name}_viralHMM.fasta")
    fasta2 = os.path.join(vvFolder,f"{name}_filtered.fasta")

    # out
    outfile = os.path.join(vvFolder,f"{name}_viralSeq.fasta")

    # dict to store sequences
    unique_seqs = {}

    # read and store uniq sequences
    def add_fasta(file):
        for record in SeqIO.parse(file, "fasta"):
            seq_str = str(record.seq)
            if seq_str not in unique_seqs:
                unique_seqs[seq_str] = record

    # process input files
    add_fasta(fasta1)
    add_fasta(fasta2)

    # write out
    SeqIO.write(unique_seqs.values(), outfile, "fasta")