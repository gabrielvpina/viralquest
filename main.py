import time
from viralquest.input_parser import FastaParser

parser = FastaParser("data/COV.fasta")
parser.read_input_file()

time.sleep(5)

for seq in parser.sequences:

    print(f"ID: {seq.id}")
    print(f"Length: {seq.length} bp")
    print(f"N Bases: {seq.n_count}")
    print(f"GC Content: {seq.gc_content:.2f}%\n")

print(f"\n\nTotal sequences processed: {len(parser.sequences)}\n")
