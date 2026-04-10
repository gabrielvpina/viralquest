import time
from viralquest.parser import FastaParser

Parser = FastaParser("data/empty.fasta")
Parser.read_input_file()

time.sleep(5)

for seq in Parser.sequences:

    print(f"ID: {seq.id}")
    print(f"Length: {seq.length} bp")
    print(f"N Bases: {seq.n_count}")
    print(f"GC Content: {seq.gc_content:.2f}%\n")

print(f"\n\nTotal sequences processed: {len(Parser.sequences)}\n")
