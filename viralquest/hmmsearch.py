import os, glob
import pandas as pd
import collections
import pyhmmer

from .data import EGG_metadata
from .data import VFAM_metadata
from .data import RVDB_metadata

RVDB_metadata = RVDB_metadata.rename(columns={'RVDB_Family': 'RVDB_Description'})


ResultHMM = collections.namedtuple("Result", ["query", "subjID", "bitscore", "start", "end"])

def process_HMM(vvFolder, hmmProfile, CPU, database_type, score_threshold=50):
    """
    Process HMM search for EggNOG, VFAM, or RVDB databases.
    """
    # Validate database type
    if database_type not in ["EggNOG", "VFAM", "RVDB"]:
        raise ValueError("database_type must be either 'EggNOG', 'VFAM', or 'RVDB'")
    
    # Set up database-specific parameters
    if database_type == "EggNOG":
        target_col = "EggNOG_TargetID"
        score_col = "EggNOG_Score"
        start_col = "EggNOG_Start"
        end_col = "EggNOG_End"
        length_col = "EggNOG_length"
        desc_col = "EggNOG_Description"
        metadata_df = EGG_metadata
        suffix = "EggNOG"
    elif database_type == "VFAM":
        target_col = "Vfam_TargetID"
        score_col = "Vfam_Score"
        start_col = "Vfam_Start"
        end_col = "Vfam_End"
        length_col = "Vfam_length"
        desc_col = "Vfam_Description"
        metadata_df = VFAM_metadata
        suffix = "Vfam"
    else:  # RVDB
        target_col = "RVDB_TargetID"
        score_col = "RVDB_Score"
        start_col = "RVDB_Start"
        end_col = "RVDB_End"
        length_col = "RVDB_length"
        desc_col = "RVDB_Description"
        metadata_df = RVDB_metadata
        suffix = "RVDB"
    
    # Encontrar arquivos de ORFs
    orf_files = glob.glob(os.path.join(vvFolder, "*_biggest_ORFs.fasta"))
    if not orf_files:
        raise FileNotFoundError("No file with '_biggest_ORFs.fasta' suffix.")

    orf_file = orf_files[0]
    name = os.path.basename(orf_file).replace("_biggest_ORFs.fasta", "")
    output_tsv = os.path.join(vvFolder, f"{name}_hmmsearch.tsv")

    # Realizar a busca HMM
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

    # Escrever resultados em TSV
    with open(output_tsv, "w") as out:
        out.write(f"Query_name\t{target_col}\t{score_col}\t{start_col}\t{end_col}\n")
        for result in results:
            out.write(f"{result.query}\t{result.subjID}\t{result.bitscore}\t{result.start}\t{result.end}\n")

    # Filtragem dos resultados
    hmm_files = glob.glob(os.path.join(vvFolder, "*_hmmsearch.tsv"))
    if not hmm_files:
        raise FileNotFoundError("No file with '_hmmsearch.tsv' suffix.")

    hmm_table_file = hmm_files[0]
    data_parser = pd.read_csv(hmm_table_file, sep="\t", low_memory=False)
    data_parser["Sample_name"] = name
    data_parser["QueryID"] = data_parser["Query_name"].str.extract(r"(.*?)_ORF")
    data_parser[score_col] = pd.to_numeric(data_parser[score_col]).round(2).fillna(1)
    data_parser[length_col] = (data_parser[end_col] - data_parser[start_col]).round(2)

    data_parser = data_parser[data_parser[score_col] >= score_threshold]
    data_parser = (
        data_parser.loc[data_parser.groupby(["Query_name"])[score_col].idxmax()]
        .reset_index(drop=True)
    )

    data_parser = pd.merge(data_parser, metadata_df, on=target_col, how="left")

    data_parser = data_parser[["Sample_name", "QueryID", "Query_name", target_col, desc_col, score_col, start_col, end_col, length_col]]

    data_parser = data_parser.loc[data_parser.groupby(["Query_name"])[score_col].idxmax()]

    csv_output_path = os.path.join(vvFolder, f"{name}_{suffix}.csv")
    data_parser.to_csv(csv_output_path, index=False)

    os.remove(hmm_table_file)


# maintains speciticity of models throught new functions
def process_EggNOG(vvFolder, hmmProfile, CPU, score_threshold=50):
    """Wrapper for EggNOG processing"""
    return process_HMM(vvFolder, hmmProfile, CPU, "EggNOG", score_threshold)

def process_VFAM(vvFolder, hmmProfile, CPU, score_threshold=50):
    """Wrapper for VFAM processing"""
    return process_HMM(vvFolder, hmmProfile, CPU, "VFAM", score_threshold)

def process_RVDB(vvFolder, hmmProfile, CPU, score_threshold=50):
    """Wrapper for RVDB processing"""
    return process_HMM(vvFolder, hmmProfile, CPU, "RVDB", score_threshold)