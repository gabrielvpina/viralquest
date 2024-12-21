<br>

<div align="center">

<img src="https://github.com/gabrielvpina/my_images/blob/main/vz_blueBG.png" width="490" height="180">
  
  <p align="center">
    <strong>A pipeline for viral diversity analysis</strong>
  </p>
</div>
<br>

# Setup
## Install required packages:
'''
conda install -c bioconda cap3 diamond hmmer blast
'''
### Pip modules:
'''
pip install orfipy pandas Bio more-itertools pyfiglet
'''

# Usage
__     ___           _ _____                
\ \   / (_)_ __ __ _| |__  /___  _ __   ___ 
 \ \ / /| | '__/ _` | | / // _ \| '_ \ / _ \
  \ V / | | | | (_| | |/ /| (_) | | | |  __/
   \_/  |_|_|  \__,_|_/____\___/|_| |_|\___|
                                            

usage: ViralZone.py [-h] -in INPUT -out OUTDIR -vir VIRALDB [-N BLASTN] [-X BLASTX] [-dX DIAMOND_BLASTX]
                    [-rvdb RVDB_HMM] [-norf NUMORFS] [-cpu CPU] [-dmnd_path DIAMOND_PATH] [-v]

None A tool for viral diversity analysis. More info in: https://github.com/gabrielvpina/ Dependencies: CAP3; Diamond;
BLAST; HMMER and ORFiPy.

options:
  -h, --help            show this help message and exit
  -in INPUT, --input INPUT
                        Fasta file containing non-host contigs to be analyzed.
  -out OUTDIR, --outdir OUTDIR
                        Directory where the output files will be saved.
  -vir VIRALDB, --viralDB VIRALDB
                        Path to the Diamond database (.dmnd) for viral sequences.
  -N BLASTN, --blastn BLASTN
                        Path to the BLASTn database for nucleotide sequence comparison.
  -X BLASTX, --blastx BLASTX
                        Path to the BLASTx database for protein sequence comparison.
  -dX DIAMOND_BLASTX, --diamond_blastx DIAMOND_BLASTX
                        Path to the Diamond BLASTx database for protein sequence comparison.
  -rvdb RVDB_HMM, --rvdb_hmm RVDB_HMM
                        Path to the RVDB hmm for conserved domain analysis.
  -norf NUMORFS, --numORFs NUMORFS
                        Number of largest ORFs to select from the input sequences.
  -cpu CPU, --cpu CPU   Number of CPU cores to be used for computation.
  -dmnd_path DIAMOND_PATH, --diamond_path DIAMOND_PATH
                        OPTIONAL - Diamond bin application path for BLAST databases: path/to/./diamond
  -v, --version         show program's version number and exit


* Suport to RVDB HMM Profile version 29.0
