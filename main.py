from viralquest.input_parser import FastaParser

meu_parser = FastaParser("data/meu_dna.fasta")
meu_parser.read_input_file()

primeira_seq = meu_parser.sequences[0]

print(f"ID: {primeira_seq.id}")
print(f"Tamanho: {primeira_seq.length} bp")
print(f"Bases N: {primeira_seq.n_count}")
print(f"Conteúdo GC: {primeira_seq.gc_content:.2f}%")
