import time
from viralquest.parser import FastaParser
from viralquest.orfs import OrfAnalyzer

Parser = FastaParser("data/test.fasta")
Parser.read_input_file()

orf_anlyzer = OrfAnalyzer(min_len_nt=150)


for seq in Parser.sequences:

    orf_anlyzer.find_n_save_orfs(seq)
    biggest_orf = orf_anlyzer.get_longest_orf(seq)

    print(f"\n# ===================================")
    print("-> SEQ INFO")
    print(f"ID: {seq.id}")
    print(f"Length: {seq.length} bp")
    print(f"N Bases: {seq.n_count}")
    print(f"GC Content: {seq.gc_content:.2f}%")
    print(f"UUID: {seq.uid}\n")

    if biggest_orf is not None:
        print(f"-> ORF INFO")
        print(f"ORF type: {biggest_orf.orf_type}")
        print(f"Length: {biggest_orf.length_aa} AA")
        print(f"Frame: {biggest_orf.frame} | Strand: {biggest_orf.strand}")
        print(f"Sequence (AA): {biggest_orf.aa_sequence}")
        print(f"# ===================================\n")

    else:
        print(f"-> ORF INFO")
        print(f"No valid ORFs :(\n")

print(f"\n\nTotal sequences processed: {len(Parser.sequences)}\n")
