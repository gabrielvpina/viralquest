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
- [Output Files](#output-files)


## Introduction
ViralQuest is a Python-based bioinformatics pipeline designed to detect, identify, and characterize viral sequences from assembled contig datasets. It streamlines the analysis of metagenomic or transcriptomic data by integrating multiple steps—such as sequence alignment, taxonomic classification, and annotation—into a cohesive and automated workflow. ViralQuest is particularly useful for virome studies, enabling researchers to uncover viral diversity, assess potential host-virus interactions, and explore the ecological or clinical significance of detected viruses.



<img src="https://github.com/gabrielvpina/viralquest/blob/main/images/VQscheme.png" width="850" height="1220">



## Setup
### Create conda enviroment
```
conda create -n viralquest python=3.12.6
```
### Install required packages (conda):

```
conda install -c bioconda cap3 diamond blast
```
### Pip modules:
```
pip install orfipy pandas Bio more-itertools pyfiglet pyhmmer
```
## Install databases






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
