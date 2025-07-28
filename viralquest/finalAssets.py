import os, glob, csv, json, shutil, shutil, warnings, sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

from .data import viral_metadata
from .dmndBLASTx import extract_organism_name

def finalTable(vvFolder, args):

    tables = glob.glob(os.path.join(vvFolder, "*_vq.fasta"))
    first_name = tables[0]
    name = os.path.basename(first_name).replace("_vq.fasta","")

    ## table BLASTn
    BlastnTable = os.path.join(vvFolder, f"{name}_blastn.tsv")
    # Verifica se o arquivo existe e se não está vazio
    if os.path.exists(BlastnTable) and os.path.getsize(BlastnTable) != 0:
        inputblastn = pd.read_csv(BlastnTable, sep='\t')
    else:
        print(f"The file {BlastnTable} doesn't exist or is empty.")


    ## table BLASTx
    BlastxTable = os.path.join(vvFolder, f"{name}_blastx.tsv")
    # Verifica se o arquivo existe e se não está vazio
    if os.path.exists(BlastxTable) and os.path.getsize(BlastxTable) != 0:
        inputblastx = pd.read_csv(BlastxTable, sep='\t')
    else:
        print(f"The file {BlastxTable} doesn't exist or is empty.")

    if os.path.exists(BlastnTable) and os.path.getsize(BlastnTable) != 0:
        inputblastn.columns = ['QueryID','BLASTn_Qlength','BLASTn_Slength','BLASTn_Cover','BLASTn_Ident','BLASTn_evalue','BLASTn_Subject_Title']
        inputblastn = inputblastn.sort_values(by=["QueryID", "BLASTn_evalue", "BLASTn_Ident"], ascending=[True, True, False])
        inputblastn = inputblastn.drop_duplicates(subset="QueryID", keep="first")
        inputblastn = inputblastn.drop_duplicates()
        inputfile = inputblastn
        inputfile["BLASTn_Subject_Title"] = inputfile["BLASTn_Subject_Title"].fillna("no_hits_BLASTn")

    if os.path.exists(BlastxTable) and os.path.getsize(BlastxTable) != 0 and os.path.exists(BlastnTable) and os.path.getsize(BlastnTable) != 0:
        inputblastx.columns = ['QueryID','BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx_Subject_Title']
        inputblastx = inputblastx.sort_values(by=["QueryID", "BLASTx_evalue", "BLASTx_Ident"], ascending=[True, True, False])
        inputblastx = inputblastx.drop_duplicates(subset="QueryID", keep="first")
        inputblastx = inputblastx.drop_duplicates()
        inputblastx['Organism_Name'] = inputblastx['BLASTx_Subject_Title'].apply(extract_organism_name)
        inputblastx = pd.merge(inputblastx, viral_metadata, left_on="Organism_Name",right_on="ScientificName", how="left")
        inputblastx.columns = ['QueryID','BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx_Subject_Title', 'BLASTx_Organism_Name', 'TaxId','ScientificName','NoRank','Clade','Kingdom','Phylum','Class','Order','Family','Subfamily','Genus','Species','Genome']
        inputfile = inputfile.merge(inputblastx, on='QueryID', how='outer')
        inputfile["BLASTx_Subject_Title"] = inputfile["BLASTx_Subject_Title"].fillna("no_hits_BLASTx")

    # if the blastn table exists, but is empty
    if os.path.exists(BlastxTable) and os.path.getsize(BlastxTable) != 0 and not os.path.exists(BlastnTable):
        inputblastx.columns = ['QueryID','BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx_Subject_Title']
        inputblastx = inputblastx.sort_values(by=["QueryID", "BLASTx_evalue", "BLASTx_Ident"], ascending=[True, True, False])
        inputblastx = inputblastx.drop_duplicates(subset="QueryID", keep="first")
        inputblastx['Organism_Name'] = inputblastx['BLASTx_Subject_Title'].apply(extract_organism_name)
        inputblastx = pd.merge(inputblastx, viral_metadata, left_on="Organism_Name", right_on="ScientificName", how="left")
        inputfile = inputblastx
        inputfile["BLASTx_Subject_Title"] = inputfile["BLASTx_Subject_Title"].fillna("no_hits_BLASTx")

    # if the blastn table not exists
    if os.path.exists(BlastxTable) and os.path.getsize(BlastxTable) != 0 and not os.path.exists(BlastnTable):
        inputblastx.columns = ['QueryID','BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover','BLASTx_Ident','BLASTx_evalue','BLASTx_Subject_Title']
        inputblastx = inputblastx.sort_values(by=["QueryID", "BLASTx_evalue", "BLASTx_Ident"], ascending=[True, True, False])
        inputblastx = inputblastx.drop_duplicates(subset="QueryID", keep="first")
        inputblastx['Organism_Name'] = inputblastx['BLASTx_Subject_Title'].apply(extract_organism_name)
        inputblastx = pd.merge(inputblastx, viral_metadata, left_on="Organism_Name", right_on="ScientificName", how="left")
        inputfile = inputblastx
        inputfile["BLASTx_Subject_Title"] = inputfile["BLASTx_Subject_Title"].fillna("no_hits_BLASTx")


    # Insert sample name column
    sample_name = f"{name}"
    inputfile['Sample_name'] = sample_name

    # Order cols by
    if os.path.exists(BlastxTable) and os.path.getsize(BlastxTable) != 0 and os.path.exists(BlastnTable) and os.path.getsize(BlastnTable) != 0:
        inputfile = inputfile[['Sample_name', 'QueryID', 'BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover', 'BLASTx_Ident', 'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name','BLASTn_Qlength', 'BLASTn_Slength', 'BLASTn_Cover', 'BLASTn_Ident', 'BLASTn_evalue', 'BLASTn_Subject_Title','TaxId','ScientificName','NoRank','Clade','Kingdom','Phylum','Class','Order','Family','Subfamily','Genus','Species','Genome']]

    # if the blastn table exists, but is empty
    if os.path.getsize(BlastxTable) != 0 and os.path.getsize(BlastnTable) == 0 and not os.path.exists(BlastnTable) :
        inputfile = inputfile[['Sample_name', 'QueryID', 'BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover', 'BLASTx_Ident', 'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name','TaxId','ScientificName','NoRank','Clade','Kingdom','Phylum','Class','Order','Family','Subfamily','Genus','Species','Genome']]
    # if the blastn table not exists
    if os.path.getsize(BlastxTable) != 0 and not os.path.exists(BlastnTable) :
        inputfile = inputfile[['Sample_name', 'QueryID', 'BLASTx_Qlength','BLASTx_Slength','BLASTx_Cover', 'BLASTx_Ident', 'BLASTx_evalue', 'BLASTx_Subject_Title', 'BLASTx_Organism_Name','TaxId','ScientificName','NoRank','Clade','Kingdom','Phylum','Class','Order','Family','Subfamily','Genus','Species','Genome']]

    if os.path.getsize(BlastxTable) == 0 and os.path.getsize(BlastnTable) != 0:
        raise ValueError("BLASTx file is empty.")

    #### column fullseq
    fasta_cap3 = os.path.join(vvFolder,f"{name}_vq.fasta")
    seq_dict = {} # store sequences

    for record in SeqIO.parse(fasta_cap3, "fasta"):
        seq_dict[record.id] = str(record.seq)

    def get_sequence(query_name):
        return seq_dict.get(query_name, None)

    inputfile['FullSeq'] = inputfile['QueryID'].apply(get_sequence)

    csv_output_path = os.path.join(vvFolder, f"{name}_all-BLAST.csv")
    inputfile = inputfile.drop_duplicates()
    inputfile.to_csv(csv_output_path, index=False)

    
    
    
    #####
    # best viral hits table
    viralSeqs = inputfile[inputfile['BLASTx_Organism_Name'].str.contains(r'virus[-\s]?|phage|riboviria', case=False, na=False)]

    # CALCULATE ViralQuest SCORE
    def calculate_vq_score(row):
        
        score = 0 #star points
        # score for BLASTx Ident and cover
        for col in ['BLASTx_Ident']:
            value = row[col]
            if value > 80:
                score += 25
            elif 70 < value <= 80:
                score += 20
            elif 60 < value <= 70:
                score += 15
            elif 50 < value <= 60:
                score += 10
            elif 40 < value <= 50:
                score += 5
            else: # less than 40
                score += 1

        for col in ['BLASTx_Cover']:
            value = row[col]
            if value > 80:
                score += 25
            elif 70 < value <= 80:
                score += 20
            elif 60 < value <= 70:
                score += 15
            elif 50 < value <= 60:
                score += 10
            elif 40 < value <= 50:
                score += 5
            else: # less than 40
                score += 1

            subject_title = str(row['BLASTn_Subject_Title']).strip()
        if pd.isna(row['BLASTn_Subject_Title']) or subject_title == '':
            # if NaN, 0 points
            pass
        else:
            # if not, 10 pts
            score += 10
            # if not, and have a viral subject, 20pts
            if pd.Series([subject_title]).str.contains(r'virus[-\s]?|phage|riboviria', case=False, na=False).iloc[0]:
                score += 10 
                # BLASTn for viral sequences values - same than BLASTx
                for col in ['BLASTn_Ident']:
                    value = row[col]
                    if value > 80:
                        score += 15
                    elif 70 < value <= 80:
                        score += 12
                    elif 60 < value <= 70:
                        score += 10
                    elif 50 < value <= 60:
                        score += 8
                    elif 40 < value <= 50:
                        score += 5
                    else: # less 40
                        score += 1
                        
                for col in ['BLASTn_Cover']:
                    value = row[col]
                    if value > 80:
                        score += 15
                    elif 70 < value <= 80:
                        score += 12
                    elif 60 < value <= 70:
                        score += 10
                    elif 50 < value <= 60:
                        score += 8
                    elif 40 < value <= 50:
                        score += 5
                    else: # less 40
                        score += 1
            
        return score

    # best vq_score

    #viralSeqs['vq_score'] = viralSeqs.apply(calculate_vq_score, axis=1)
    warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)
    viralSeqs.loc[:, 'vq_score'] = viralSeqs.apply(calculate_vq_score, axis=1)
    # order by biggest score
    viralSeqs = viralSeqs.sort_values('vq_score', ascending=False)
    number_viralSeqs = viralSeqs['QueryID'].nunique()
    viralSeqs_path = os.path.join(vvFolder, f"{name}_viral-BLAST.csv")
    viralSeqs.to_csv(viralSeqs_path, index=False)





    # move tables
    mydir = "hit_tables"
    os.makedirs(os.path.join(vvFolder, mydir))
    mydir_path = os.path.join(vvFolder, mydir)
    shutil.move(BlastnTable, mydir_path)
    shutil.move(BlastxTable, mydir_path)

    if os.path.exists(os.path.join(vvFolder,f"{name}_RVDB.csv")):
        shutil.move(os.path.join(vvFolder,f"{name}_RVDB.csv"),mydir_path)
    if os.path.exists(os.path.join(vvFolder,f"{name}_Vfam.csv")):
        shutil.move(os.path.join(vvFolder,f"{name}_Vfam.csv"),mydir_path)
    if os.path.exists(os.path.join(vvFolder,f"{name}_EggNOG.csv")):
        shutil.move(os.path.join(vvFolder,f"{name}_EggNOG.csv"),mydir_path)

    # remove unused table
    os.remove(os.path.join(vvFolder,f"{name}_ref.tsv"))






    ################## CREATE JSON #####################
    # hmm hits 
    hmmTable_path = os.path.join(vvFolder, f"{name}_hmm.csv")
    hmmTable = pd.read_csv(hmmTable_path, low_memory=False)

    if not viralSeqs.empty:
        viralHmmTable = hmmTable[hmmTable["QueryID"].isin(viralSeqs["QueryID"])]
        viralSeqs = viralSeqs.sort_values('BLASTx_Qlength', ascending=False)
        viralSeqs = viralSeqs.dropna(axis=1, how='all')
        viralHmmTable.to_csv(hmmTable_path, index=False)
    else:
        if args.rvdb_hmm:
                hmmTable = hmmTable[hmmTable["RVDB_Score"] > 800]
        if args.vfam_hmm:
                hmmTable = hmmTable[hmmTable["Vfam_Score"] > 800]
        if args.eggnog_hmm:
                hmmTable = hmmTable[hmmTable["EggNOG_Score"] > 800]

        hmmTable.to_csv(hmmTable_path, index=False)
    

    # blastx hits
    refTable_path = os.path.join(vvFolder, f"{name}_ref.csv")
    refTable = pd.read_csv(refTable_path)
    
    # join and remove reps
    if not viralSeqs.empty:
        hmm_viral_ids = set(hmmTable[hmmTable["QueryID"].isin(viralSeqs["QueryID"])]["QueryID"])
        ref_ids = set(refTable["QueryID"].unique())
        valid_query_ids = hmm_viral_ids.union(ref_ids)
    else:
        valid_query_ids = set(refTable["QueryID"].unique())



    def csv_to_json(csv_file_path1, csv_file_path2):
        # Create JSON file
        # Dealing with high values in FullSequence field (for large genomes)
        max_int = sys.maxsize
        while True:
            try:
                csv.field_size_limit(max_int)
                break
            except OverflowError:
                max_int = int(max_int / 10)
        ####################################################################

        with open(csv_file_path1, mode='r', encoding='utf-8') as file:
            csv_reader = csv.DictReader(file)
            viral_hits = [row for row in csv_reader]

        with open(csv_file_path2, mode='r', encoding='utf-8') as file:
            csv_reader = csv.DictReader(file)
            hmm_hits = [row for row in csv_reader if row['QueryID'] in valid_query_ids]

        return {"Viral_Hits": viral_hits, "HMM_hits": hmm_hits}
    

    def parse_orf_fastas(file_path):
        # Process FASTA file
        orf_data_by_code = defaultdict(list)

        for record in SeqIO.parse(file_path, "fasta"):
            header = record.description
            parts = header.split()

            contig_full = parts[0]
            if '_ORF.' in contig_full:
                contig = contig_full[:contig_full.index('_ORF.')]
            else:
                contig = contig_full

            if contig not in valid_query_ids:
                continue

            coordinates_strand = parts[1].replace('[', '').replace(']', '')
            if '(' in coordinates_strand:
                coordinates, strand = coordinates_strand.split('(')
                coordinates = coordinates.strip()
                strand = strand.replace(')', '').strip()
                start, end = map(int, coordinates.split('-'))
            else:
                continue

            additional_info = parts[2:]
            additional_info_dict = {info.split(':')[0]: info.split(':')[1] for info in additional_info}

            orf = {
                'Query_name': contig_full,
                'QueryID': contig,
                'start': start,
                'end': end,
                'strand': strand,
                'sequence': str(record.seq),
                'type': additional_info_dict.get('type'),
                'length': int(additional_info_dict.get('length', 0)),
                'frame': additional_info_dict.get('frame'),
                'start_codon': additional_info_dict.get('start'),
                'stop_codon': additional_info_dict.get('stop')
            }

            orf_data_by_code[contig].append(orf)

        return dict(orf_data_by_code)

    def combined_json(csv_file_path1, csv_file_path2, fasta_file_path, output_json_path):
        csv_data = csv_to_json(csv_file_path1, csv_file_path2)
        orf_data = parse_orf_fastas(fasta_file_path)

        # Merge results
        combined_data = {**csv_data, "ORF_Data": orf_data}

        with open(output_json_path, mode='w', encoding='utf-8') as json_file:
            json.dump(combined_data, json_file, indent=4)



    combined_json(
        viralSeqs_path,
        hmmTable_path,
        os.path.join(vvFolder,f'{name}_biggest_ORFs.fasta'),
        os.path.join(vvFolder,f'{name}_bestSeqs.json')
    )

    #####################################################

    # move fasta files
    fastaDir = "fasta-files"
    os.makedirs(os.path.join(vvFolder, fastaDir))
    fastaDir_path = os.path.join(vvFolder, fastaDir)

    for fasta_text in os.listdir(vvFolder):
        if fasta_text.endswith(".fasta"):
            shutil.move(os.path.join(vvFolder, fasta_text), fastaDir_path)

    return number_viralSeqs


def getViralSeqsNumber(vvFolder):
    fileFasta = os.path.join(vvFolder, "*_viral-BLAST.csv")
    viralSeqs_path = glob.glob(fileFasta)
    table = viralSeqs_path[0]
    viralSeqs = pd.read_csv(table)
    number_viralSeqs = viralSeqs['QueryID'].nunique()

    return number_viralSeqs


### CLEAN JSON ###

def clean_hmm_hits_data(data):
    
    # Make a copy to avoid modifying the original data
    cleaned_data = data.copy()
    
    if 'HMM_hits' not in cleaned_data:
        return cleaned_data
    
    # Process each hit in HMM_hits
    cleaned_hits = []
    
    for hit in cleaned_data['HMM_hits']:
        cleaned_hit = clean_single_hit(hit)
        cleaned_hits.append(cleaned_hit)
    
    cleaned_data['HMM_hits'] = cleaned_hits
    return cleaned_data

def clean_single_hit(hit):

    # group fields by suffix (_2, _3, etc.)
    field_groups = defaultdict(list)
    base_fields = []
    
    for field_name in hit.keys():
        if '_' in field_name and field_name.split('_')[-1].isdigit():
            # Extract the base name and suffix
            parts = field_name.split('_')
            suffix = parts[-1]
            field_groups[suffix].append(field_name)
        else:
            base_fields.append(field_name)
    
    # Create cleaned hit starting with base fields
    cleaned_hit = {}
    for field in base_fields:
        cleaned_hit[field] = hit[field]
    
    # Check each suffix group and only keep if not all empty
    for suffix, fields in field_groups.items():
        # Check if all fields in this suffix group are empty
        all_empty = all(
            hit[field] == "" or hit[field] is None or hit[field] == "0.0" 
            for field in fields
        )
        
        # If not all empty, keep these fields
        if not all_empty:
            for field in fields:
                cleaned_hit[field] = hit[field]
    
    return cleaned_hit

def load_and_clean_json(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        cleaned_data = clean_hmm_hits_data(data)
        return cleaned_data
    
    except FileNotFoundError:
        print(f"File {file_path} not found")
        return None
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}")
        return None

def save_cleaned_json(cleaned_data, output_path):
    try:
        with open(output_path, 'w') as f:
            json.dump(cleaned_data, f, indent=2)
    except Exception as e:
        print(f"Error saving file: {e}")


def exec_cleaning(vvFolder):
    fileJSON = os.path.join(vvFolder, "*.json")
    viralJSON_path = glob.glob(fileJSON)
    json_file = viralJSON_path[0]

    cleaned_data = load_and_clean_json(json_file)
    if cleaned_data:
        save_cleaned_json(cleaned_data, json_file)



