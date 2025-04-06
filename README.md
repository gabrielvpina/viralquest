<br>

<div align="center">

<img src="https://github.com/gabrielvpina/viralquest/blob/main/images/headerLogo.png" width="530" height="180">
  
  <p align="center">
    <strong>A pipeline for viral diversity analysis</strong>
  </p>
</div>

- [Introduction](#introduction)
- [Setup](#setup)
- [Install  Databases](#install-databases)
- [Viral HMM Models](#viral-hmm-models)
- [Output Files](#output-files)


## Introduction
ViralQuest is a Python-based bioinformatics pipeline designed to detect, identify, and characterize viral sequences from assembled contig datasets. It streamlines the analysis of metagenomic or transcriptomic data by integrating multiple steps—such as sequence alignment, taxonomic classification, and annotation—into a cohesive and automated workflow. ViralQuest is particularly useful for virome studies, enabling researchers to uncover viral diversity, assess potential host-virus interactions, and explore the ecological or clinical significance of detected viruses.



<img src="https://github.com/gabrielvpina/viralquest/blob/main/images/VQscheme.png" width="850" height="1220">



## Setup
### Python Enviroment and Packages
Create conda enviroment
```
conda create -n viralquest python=3.12.6
```
Activate conda enviroment
```
conda activate viralquest
```
Install required packages (conda):
```
conda install -c bioconda cap3 diamond blast
```
Install required packages (pip):
```
pip install orfipy pandas Bio more-itertools pyfiglet pyhmmer
```

## Install Databases

### RefSeq Viral release
The RefSeq viral release is a curated collection of viral genome and protein sequences provided by the NCBI Reference Sequence (RefSeq) database. It includes high-quality, non-redundant, and well-annotated reference sequences for viruses, maintained and updated regularly by NCBI. The required file is `viral.1.protein.faa.gz`, download via [this link](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz).

### BLAST nr/nt Databases
The BLAST nr (non-redundant protein) and nt (nucleotide) databases are essential resources for viral identification. The nt database is useful for identifying viral genomes or transcripts using nucleotide similarity, while nr is especially powerful for detecting and annotating viral proteins, even in divergent or novel viruses, through translated searches like blastx.
Download the nr/nt databases in fasta format via [this link](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
- The file `nr.gz` is the nr database in FASTA
- The `nt.gz` file correspond to nt.fasta

## Viral HMM Models
The `Vfam` and `eggNOG` models needs to be 

### Vfam HMM
The VFam HMM models are profile Hidden Markov Models (HMMs) specifically designed for the identification of viral proteins. 

**Steps to Install**
1) Download `vfam.hmm.tar.gz` via [this link](https://fileshare.lisc.univie.ac.at/vog/vog228/vfam.hmm.tar.gz):
```
wget https://fileshare.lisc.univie.ac.at/vog/vog228/vfam.hmm.tar.gz
```
2) Extract the file:
```
tar -xzvf vfam.hmm.tar.gz
```
3) Unify all `.hmm` models in one model:
```
cat *.hmm >> vfam228.hmm
```
Now it's possible to use the `vfam228.hmm` file in the **ViralQuest** pipeline.

### eggNOG Viral HMM
The eggNOG viral OGs HMM models are part of the eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) resource and are designed to identify and annotate viral genes and proteins based on orthologous groups (OGs).

**Steps to Install**
1) Download the viral OGs in the eggNOG Database via [this link](http://eggnog45.embl.de/#/app/viruses). The HMM models download are in the last column.

2) Or download the data via this BASH script:
```
#!/bin/bash

mkdir eggNOG
cd eggNOG

wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/ssRNA/ssRNA.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Retrotranscribing/Retrotranscribing.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/dsDNA/dsDNA.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Viruses/Viruses.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Herpesvirales/Herpesvirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/ssDNA/ssDNA.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/ssRNA_positive/ssRNA_positive.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Retroviridae/Retroviridae.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Ligamenvirales/Ligamenvirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Caudovirales/Caudovirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Mononegavirales/Mononegavirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Tymovirales/Tymovirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Nidovirales/Nidovirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/Picornavirales/Picornavirales.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/dsRNA/dsRNA.hmm.tar.gz
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/viruses/ssRNA_negative/ssRNA_negative.hmm.tar.gz

for i in *.tar.gz; do tar -zxvf "$i" ;done
```
Save as `download_eggNOG.sh`. Now let's execute:
```
chmod +x download_eggNOG.sh && ./download_eggNOG.sh
```
3) Now join all result files:
```
cat eggNOG/hmm_files/*.hmm >> eggNOG.hmm
```

Now it's possible to use the `eggNOG.hmm` file in the **ViralQuest** pipeline.

### RVDB Viral HMM
The Reference Viral Database (RVDB) is a curated collection of viral sequences, and its protein HMM models—RVDB-prot and RVDB-prot-HMM—are designed to enhance the detection and annotation of viral proteins.

**Download RVDB hmm model**
1) Visit the RVDB Protein database via [this link](https://rvdb-prot.pasteur.fr/) and download the hmm model version 29.0.
2) Or download directly via linux termnial:
```
wget https://rvdb-prot.pasteur.fr/files/U-RVDBv29.0-prot.hmm.xz
```
3) Decompress the model:
```
unxz -v U-RVDBv29.0-prot.hmm.xz
```
Now it's possible to use the `U-RVDBv29.0-prot.hmm` file in the **ViralQuest** pipeline.

## Output Files
This is the output files:
```
├── fasta-files
│   ├── SRR1234_all_ORFs.fasta
│   ├── SRR1234_biggest_ORFs.fasta
│   ├── SRR1234_filtered.fasta
│   ├── SRR1234_orig.fasta
│   ├── SRR1234_pfam_ORFs.fasta
│   ├── SRR1234_viralHMM.fasta
│   ├── SRR1234_viralSeq.fasta
│   └── SRR1234_vq.fasta
├── hit_tables
│   ├── SRR1234_blastn.tsv
│   ├── SRR1234_blastx.tsv
│   ├── SRR1234_EggNOG.csv
│   ├── SRR1234_ref.tsv
│   ├── SRR1234_RVDB.csv
│   └── SRR1234_Vfam.csv
├── SRR1234_all-BLAST.csv
├── SRR1234_bestSeqs.json
├── SRR1234_hmm.csv
├── SRR1234.log
├── SRR1234_ref.csv
├── SRR1234_viral-BLAST.csv
└── SRR1234_visualization.html
```


[Example of HTML Viral Report Output (Click Here)](https://chocolate-yetta-73.tiiny.site)
* Suport to Reference Viral Database (RVDB) HMM Profile version 29.0
* Suport to Protein Family Database (Pfam) HMM Profile version 37.2
* Suport to Vfam HMM Profile version 228
* Suport to EggNOG Virus HMM Profile version 4.5
