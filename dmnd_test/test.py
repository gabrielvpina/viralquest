import os
import math
import subprocess
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm
from loguru import logger

class DiamondPipeline:
    """
    Classe para gerenciamento e execução otimizada de alinhamentos com Diamond,
    desenhada para lidar com arquivos massivos e o banco de dados NR do NCBI.
    """

    def __init__(self, fasta_input: str, db_path: str, temp_dir: str, max_ram_gb: float, max_threads: int):
        self.fasta_input = Path(fasta_input)
        self.db_path = Path(db_path)
        self.temp_dir = Path(temp_dir)
        self.max_ram_gb = float(max_ram_gb)
        self.max_threads = int(max_threads)
        
        # Garante que o diretório temporário exista
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Calcula os parâmetros otimizados na inicialização da classe
        self.parametros = self._calcular_parametros_hardware()

    def _calcular_parametros_hardware(self) -> Dict[str, float]:
        """
        Calcula os parâmetros dinâmicos para o Diamond baseado em QUALQUER valor
        de RAM e CPU fornecido, usando uma função matemática contínua.
        """
        # 10% de margem de segurança para o Sistema Operacional (evita OOM Killer)
        ram_segura = self.max_ram_gb * 0.90
        
        # O banco NR (346GB) exige aproximadamente 64GB de RAM para carregar o índice completo
        RAM_NECESSARIA_INDICE = 64.0 
        
        if ram_segura >= RAM_NECESSARIA_INDICE:
            # Tem RAM sobrando: Carrega o índice inteiro de uma vez
            index_chunks = 1
            ram_para_blocos = ram_segura - RAM_NECESSARIA_INDICE
        else:
            # RAM limitada: O índice precisará ser fatiado.
            # Destinamos 60% da RAM segura para carregar pedaços do índice 
            # e deixamos 40% para o processamento das sequências (blocos).
            ram_para_indice = ram_segura * 0.60
            index_chunks = math.ceil(RAM_NECESSARIA_INDICE / ram_para_indice)
            
            # Recalcula a RAM exata que sobrou após fatiar o índice
            ram_para_blocos = ram_segura - (RAM_NECESSARIA_INDICE / index_chunks)
        
        # A fórmula de consumo do Diamond é: ~6GB de RAM por cada 1.0 de block_size
        block_size = ram_para_blocos / 6.0
        
        # Limitações de segurança:
        # Nunca menor que 0.5 (para não quebrar a ferramenta)
        # Nunca maior que 20.0 (ganhos marginais acima disso, aumenta tempo de Garbage Collection)
        block_size = max(0.5, min(block_size, 20.0))
        
        params = {
            "threads": self.max_threads,
            "block_size": round(block_size, 2),
            "index_chunks": int(index_chunks),
            "memory_limit": int(ram_segura)
        }
        
        return params

    def _dividir_fasta_por_lote(self, seqs_por_lote: int = 10000) -> List[Path]:
        """
        Lê o FASTA de forma eficiente (O(1) em RAM) e o divide em vários chunks.
        """
        arquivos_gerados = []
        lote_atual = 1
        contador_seqs = 0
        
        caminho_chunk = self.temp_dir / f"chunk_{lote_atual}.fasta"
        out_file = open(caminho_chunk, "w")
        arquivos_gerados.append(caminho_chunk)
        
        with open(self.fasta_input, "r") as f_in:
            for linha in f_in:
                if linha.startswith(">"):
                    if contador_seqs == seqs_por_lote:
                        out_file.close()
                        lote_atual += 1
                        contador_seqs = 0
                        caminho_chunk = self.temp_dir / f"chunk_{lote_atual}.fasta"
                        out_file = open(caminho_chunk, "w")
                        arquivos_gerados.append(caminho_chunk)
                    contador_seqs += 1
                out_file.write(linha)
                
        out_file.close()
        return arquivos_gerados

    def _executar_subprocesso(self, chunk_input: Path, chunk_output: Path):
        """
        Isola a execução do subprocesso do Diamond.
        """
        comando = [
            "diamond", "blastx",
            "-d", str(self.db_path),
            "-q", str(chunk_input),
            "-o", str(chunk_output),
            "--threads", str(self.parametros["threads"]),
            "--block-size", str(self.parametros["block_size"]),
            "--index-chunks", str(self.parametros["index_chunks"]),
            # "--memory-limit", str(self.parametros["memory_limit"]),
            "--tmpdir", "/dev/shm", # Utiliza disco em RAM para temporários
            "--quiet"
        ]
        subprocess.run(comando, check=True, text=True, capture_output=True)

    def run_com_barra_progresso(self, seqs_por_lote: int = 10000):
        """
        Executa a pipeline exibindo uma barra de progresso no terminal via tqdm.
        """
        print("\n=== Inicializando Pipeline Diamond ===")
        print(f"Hardware Dinâmico -> Threads: {self.parametros['threads']} | Limite RAM: {self.parametros['memory_limit']}GB")
        print(f"Parâmetros Otimizados -> Block Size: {self.parametros['block_size']} | Index Chunks: {self.parametros['index_chunks']}")
        
        print(f"\nFatiando FASTA de entrada (Lotes de {seqs_por_lote} sequências)...")
        chunks = self._dividir_fasta_por_lote(seqs_por_lote)
        
        for chunk in tqdm(chunks, desc="Progresso do Alinhamento", unit="lote", colour="cyan"):
            saida_chunk = chunk.with_name(chunk.name.replace(".fasta", "_resultado.tsv"))
            try:
                self._executar_subprocesso(chunk, saida_chunk)
            except subprocess.CalledProcessError as e:
                tqdm.write(f"[ERRO CRÍTICO] Falha no {chunk.name}: {e.stderr}")
                
        print("\nProcessamento UI concluído com sucesso.")

    def run_com_logs(self, seqs_por_lote: int = 10000):
        """
        Executa a pipeline em background registrando o status no loguru.
        """
        logger.add(self.temp_dir / "pipeline_diamond_{time}.log", rotation="10 MB", level="INFO")
        logger.info(f"Iniciando Pipeline Diamond. Parâmetros: {self.parametros}")
        
        chunks = self._dividir_fasta_por_lote(seqs_por_lote)
        total_chunks = len(chunks)
        logger.info(f"FASTA fatiado em {total_chunks} lotes de até {seqs_por_lote} sequências.")
        
        for indice, chunk in enumerate(chunks, start=1):
            saida_chunk = chunk.with_name(chunk.name.replace(".fasta", "_resultado.tsv"))
            logger.info(f"Processando lote {indice}/{total_chunks} [{chunk.name}]...")
            try:
                self._executar_subprocesso(chunk, saida_chunk)
                logger.success(f"Lote {indice} processado. Resultados em: {saida_chunk.name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Falha no lote {indice}: {e.stderr}")
                
        logger.info("Execução de background finalizada.")


# ==========================================
# EXEMPLO DE USO NO SEU SCRIPT PRINCIPAL
# ==========================================
if __name__ == "__main__":
    # Aqui você pode pedir o input do usuário ou ler os recursos reais da máquina usando os.cpu_count()
    try:
        ram_usuario = float(input("Digite a quantidade de RAM disponível (em GB, ex: 16, 64.5, 128): "))
        cpu_usuario = int(input("Digite o número de threads/núcleos a serem utilizados (ex: 8, 32): "))
    except ValueError:
        print("Entrada inválida. Usando valores padrão: 16GB RAM, 4 Threads.")
        ram_usuario = 16.0
        cpu_usuario = 4

    # Instancia a classe do Pipeline
    pipeline = DiamondPipeline(
        fasta_input="test.fasta",
        db_path="/home/gabriel/Downloads/viralDB.dmnd",
        temp_dir="./resultados_temporarios",
        max_ram_gb=ram_usuario,
        max_threads=cpu_usuario
    )

    # Inicia a execução (Escolha o método desejado)
    #pipeline.run_com_barra_progresso(seqs_por_lote=10000)
    pipeline.run_com_logs(seqs_por_lote=1000)
