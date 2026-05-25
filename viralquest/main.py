# ==================================================================
# USING CAP3
# ==================================================================
# 1. Instancia o runner
runner = Cap3Runner(file_path="./meus_dados/contigs.fasta", outdir="./resultados_cap3")
# 2. Executa e recebe a dataclass preenchida
resultados = runner.cap3_runner()
# 3. Interage com os dados de forma elegante
if resultados.is_successful:
    resultados.get_combined_fasta(Path("./resultados_cap3/montagem_final.fasta"))



