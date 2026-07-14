# Read QC & Sequence Quality module

Um módulo acionado por `--reads` que agrega dois sinais complementares de qualidade:

1. **Qualidade das reads** (fastp) + **coloração da track de cobertura** pela qualidade de base, reaproveitando o BAM já produzido pelo `coverage.py`.
2. **Qualidade estrutural da sequência** — regiões de baixa complexidade (dustmasker), repetições internas (self-BLASTn), repetitividade de k-mers (jellyfish) e um **dot plot em SVG nativo**.

Ferramentas fixadas: `fastp` · `dustmasker` · `blastn` (self-BLAST) · `jellyfish` · dot plot SVG nativo · coloração da cobertura.

---

## Status

| Fase | Escopo | Situação |
|------|--------|----------|
| **0** | Dependências (`pixi.toml`) | ✅ **Implementada** |
| **1** | Modelo de dados (`biodata.py`) | ✅ **Implementada** |
| **2a** | Reaproveitar BAM p/ profundidade + qualidade (`coverage.py`) | ✅ **Implementada** |
| **2b** | Módulos `read_qc.py` / `seq_quality.py` | ✅ **Implementada** |
| **3** | Fiação no CLI + export JSON (`cli.py`, `exporter.py`) | ✅ **Implementada** |
| **4** | Enrich / stats (`html_report.py`) | ✅ **Implementada** |
| **5** | Componentes do relatório (`components/`) | ✅ **Implementada** |
| **6** | Testes | ⏳ Planejada (validado inline) |

---

## Fase 1 — Modelo de dados (implementado)

Arquivo: `viralquest/biodata.py`.

### `CoverageProfile` — dois campos novos
```python
quality_bins: list[float] = field(default_factory=list)  # qualidade média (Phred) por bin, alinhada a `bins`
mean_quality: float = 0.0                                 # qualidade média de base (Phred) na sequência
```
`quality_bins` é **1:1** com `bins` (mesma binagem/`_N_BINS`), então cada bin da cobertura tem sua cor de qualidade. Ambos os campos têm default, portanto construtores existentes (e testes) não quebram.

### Dataclasses novas

- **`SelfRepeat`** — um HSP interno do self-BLASTn (a diagonal self-completa é filtrada a montante):
  `q_start, q_end, s_start, s_end` (0-based, inclusivos), `strand` (`plus` = repetição direta, `minus` = repetição invertida), `pct_id`, `length`.

- **`DotPlot`** — dot plot pré-computado para virar `<svg>` nativo:
  `word` (tamanho do k-mer), `length` (lado, em bp), `segments: list[[x1,y1,x2,y2]]`. Segmentos na diagonal principal = a própria sequência; fora da diagonal = repetição direta; anti-diagonal = repetição invertida.

- **`SequenceQuality`** — sinais estruturais por sequência:
  `seq_id, length, low_complexity_regions: list[[start,end]]` (dustmasker), `self_repeats: list[SelfRepeat]`, `kmer_size, kmer_distinct, kmer_total, kmer_repeat_score` (fração de k-mers com multiplicidade > 1, 0–1), `low_complexity_frac`, `dotplot: DotPlot | None`.

- **`ReadQcReport`** — QC global das reads (fastp), com vetores por-ciclo já downsampleados:
  `reads, read_type, reads_before/after, bases_before/after, q20_rate, q30_rate, gc_pct, mean_length, dup_rate, adapter_rate, per_base_quality, per_base_gc`.

### `NucSequence` — campo novo
```python
seq_quality: "SequenceQuality | None" = field(default=None, init=False)
```
Fica ao lado de `coverage`, populado quando `--reads` é usado.

> Serialização: `exporter._to_serializable` usa `dataclasses.fields()` genericamente, logo todas as novas dataclasses serializam automaticamente assim que forem plugadas no export (Fase 4).

---

## Fase 2a — Reaproveitar o BAM (implementado)

Arquivo: `viralquest/coverage.py`.

O passo antigo `samtools depth -a` (só profundidade) foi **substituído** por um único passe de `samtools mpileup` sobre o **mesmo BAM** já gerado pelo alinhamento minimap2 — sem novo alinhamento, sem I/O extra além de reler o BAM.

### Método `_compute_depth_and_quality`
```
samtools mpileup -a -A -B -Q 0 -q 0 -d 0 coverage.bam
```
- Colunas usadas: `4` = profundidade (idêntica à contagem anterior) · `6` = qualidades de base (ASCII, Phred+33), 1 char por base empilhada.
- Filtros todos desligados (`-A -B -Q 0 -q 0`) para a profundidade bater com o `depth -a` anterior; `-a` emite toda posição; `-d 0` remove o teto de profundidade.
- Retorna `({seq_id: depth}, {seq_id: mean_quality})`, ambos `np.zeros(length)`; qualidade média por posição = média de `ord(c) - 33` sobre a coluna 6.

### `_build_profiles`
- Passa a receber `quality` e preenche `quality_bins` + `mean_quality`.
- **Qualidade só onde há reads:** a média é feita sobre posições com `depth > 0` (via máscara `covered`), para gaps não diluírem a cor.

### `_downsample_quality` (novo)
Binagem alinhada 1:1 com `_downsample`, mas mediando apenas posições cobertas dentro de cada bin (bin sem cobertura → `0.0`).

**Validação:** módulos importam e instanciam OK; a suíte `test_exporter.py` mantém exatamente o mesmo resultado de antes das mudanças (`12 failed, 32 passed` — falhas **pré-existentes**, não introduzidas aqui). Não há outros construtores de `CoverageProfile` fora do `coverage.py`.

---

## Fase 0 — Dependências (implementado)

Arquivo: `pixi.toml`. Adicionados ao bloco `[dependencies]`:
```toml
fastp = ">=0.24,<1"
kmer-jellyfish = ">=2.3,<3"
```
`dustmasker`, `blastn` e `makeblastdb` já vêm do `blast` (`>=2.16,<3`). Todas as ferramentas degradam com aviso se ausentes do PATH (padrão do `coverage.py`: `shutil.which` → warning → segue sem o sinal).

---

## Fase 2b — Módulos de análise (implementado)

### `viralquest/read_qc.py` → `ReadQcPipeline`

`run(reads, outdir) -> ReadQcReport | None`. Roda fastp uma vez sobre as reads e parseia `fastp.json`.

- **Comando:** `fastp --thread N -j fastp.json -h fastp.html -i R1 [-I R2] -o clean… [-O clean…]`. Single/paired escolhidos pelo nº de arquivos (mesma convenção do `CoveragePipeline`).
- **Long reads** (`ont/pb/hifi`): adiciona `--disable_adapter_trimming` e força `adapter_rate = 0`.
- **QC-only por padrão:** as reads limpas são escritas em `outdir` e o caminho fica em `self.cleaned_reads` (para o `--trim-reads` da Fase 3 optar por usá-las a jusante); o resto do pipeline continua com as reads originais.
- **Parse do JSON:** `reads/bases before/after` (`summary`), `q20/q30_rate`, `gc_content×100`, `read1_mean_length`, `duplication.rate`, `adapter_cutting.adapter_trimmed_reads / reads_before`, e as curvas por-ciclo `read1_after_filtering.quality_curves.mean` e `content_curves.GC`.
- **`_downsample`:** curvas por-ciclo reduzidas a ≤ 300 pontos (médias de bloco) — evita SVG gigante em long reads.
- Falha graciosa: fastp ausente, nº de arquivos ≠ 1–2, ou JSON ilegível → `None`.

### `viralquest/seq_quality.py` → `SequenceQualityPipeline`

`run(viral_seqs, outdir) -> dict[seq_id, SequenceQuality]`. Escreve um `viral_ref.fasta` multi-record e roda os 4 sinais; cada sequência sempre recebe uma entrada (listas vazias se um sinal faltar).

- **dustmasker** (`_dustmask`): `dustmasker -in ref -outfmt interval` → parse de `>seqid` + linhas `a - b` (0-based, inclusivas) → `low_complexity_regions`; `low_complexity_frac` = bases mascaradas / comprimento.
- **self-BLASTn** (`_self_blast`): uma chamada `blastn -query ref -subject ref -dust yes -evalue 1e-5 -word_size 11 -outfmt "6 qseqid sseqid qstart qend sstart send pident length sstrand"`. Filtra: hits cross-sequência (`qseqid≠sseqid`), a diagonal self-completa (plus + coords idênticas) e duplicatas simétricas (A→B/B→A) via chave de triângulo superior. Coords `minus` normalizadas (`s_lo/s_hi`). → `SelfRepeat` (strand `plus`=direta, `minus`=invertida).
- **jellyfish** (`_jellyfish`, **por sequência**): `jellyfish count -m K -s 10M -C -o out.jf seq.fa` + `jellyfish histo`. `kmer_repeat_score` = `(total − singletons) / total` (fração de *instâncias* de k-mer que recorrem, multiplicidade > 1); expõe também `kmer_distinct`/`kmer_total`. `K` default = 15.
- **dot plot SVG nativo** (`_dotplot`, numpy/puro-Python): amostra a sequência a cada `stride = L // 500` bases; índice de k-mers forward (`word` default 12). Matches forward agrupados por offset diagonal `d = y − x` (triângulo superior incl. diagonal identidade `d=0`); matches contra o **reverse-complement** agrupados por anti-diagonal `s = x + y` marcam repetições invertidas. `_emit_segments` funde âncoras colineares (gap ≤ 2·stride) em segmentos `[x1,y1,x2,y2]` (bp), limitados a 4000 por sequência. → `DotPlot`.

**Validação inline:** parse do fastp (incl. adapter=0 p/ long read) e curva downsampleada OK; dot plot detecta corretamente repetição direta (diagonal) e invertida (anti-diagonal) numa sequência sintética; parsers do dustmasker e do self-BLASTn conferidos com saída mockada (diagonal trivial, dup simétrica e cross-seq descartadas; direta e invertida mantidas com coords 0-based normalizadas). Ambos os módulos importam limpos. Testes formais ficam na Fase 6.

---

## Fase 3 — CLI + export JSON (implementado)

### Argumentos novos (`cli.py`, grupo *"read / sequence quality (optional, needs --reads)"*)
- `--skip-read-qc` — pula o passo de read QC (fastp); os sinais de qualidade de sequência ainda rodam.
- `--trim-reads` — alimenta as reads limpas do fastp ao Salmon/coverage no lugar das brutas (implica read QC; ignorado com `--skip-read-qc`).
- `--kmer K` — tamanho de k-mer do score do jellyfish (default 15).

Validações (`_validate_args`): `--trim-reads`/`--skip-read-qc` exigem `--reads`; `--kmer >= 2`.

### Passos e ordem (`_build_steps` ⟷ `_run_pipeline`)
Helper novo `_read_qc_enabled(args)` = `bool(args.reads) and not args.skip_read_qc`. O bloco de reads passa a ter **quatro** passos, na mesma ordem em ambas as funções:

1. **Read QC — fastp** (`_read_qc_enabled`) — roda `ReadQcPipeline` em `outdir/read_qc` → `read_qc_report`. Com `--trim-reads`, se as reads limpas existem e casam em número, `args.reads` é trocado por elas (afeta Salmon **e** coverage a jusante).
2. **Salmon quantification** (`_salmon_enabled`) — inalterado.
3. **Sequence quality — dustmask · self-BLAST · jellyfish** (`if args.reads`) — roda `SequenceQualityPipeline(kmer_size=args.kmer)` no mesmo conjunto `select_confirmed_sequences(seqs, force)` do coverage; anexa `seq.seq_quality`.
4. **Read coverage profiling** (`if args.reads`) — inalterado.

Contagem de passos validada: sem reads = 10; `sr` = 14 (read QC + Salmon + seq-quality + coverage); `ont` = 13 (sem Salmon); `sr --skip-read-qc` = 13 (sem read QC).

### Export (`exporter.py`)
- `export(...)` ganhou o parâmetro `read_qc_report: ReadQcReport | None` → adiciona `report["read_qc"]` (via `_to_serializable`).
- `_seq_to_dict` ganhou `"seq_quality"` por sequência.
- `cli.py` passa `read_qc_report=read_qc_report` ao `exporter.export(...)`.

**Validação inline:** `_build_steps` conferido nos 4 combos de flags; ordem de execução idêntica à lista de passos. Smoke test end-to-end de serialização com `NucSequence` populado (coverage `quality_bins`/`mean_quality`, `seq_quality` com dot plot + self-repeats) + `ReadQcReport`: `report["read_qc"]` e `sequences[].seq_quality` presentes e com round-trip JSON OK. Suíte `test_exporter.py` mantém o baseline (`12 failed, 32 passed` — falhas pré-existentes). Testes formais dos módulos ficam na Fase 6.

---

## Fase 4 — Enrich / stats (implementado)

Arquivo: `viralquest/html_report.py`. Duas funções novas espelhando `_build_salmon_stats`, plugadas no `_enrich_report` ao lado da linha do salmon:

- **`_build_readqc_stats(read_qc)`** → `report["readqc_stats"]`. Card-resumo do fastp: `{present: False}` sem dados; senão pass-through de `read_type`, `reads_before/after`, `bases_after`, `q20/q30_rate`, `gc_pct`, `mean_length`, `dup_rate`, `adapter_rate` + derivado `pass_rate = reads_after/reads_before` (guarda divisão por zero). As curvas `per_base_quality`/`per_base_gc` **não** são reexpostas aqui — o componente as lê direto de `report["read_qc"]` para plotar.
- **`_build_seq_quality_stats(seqs)`** → `report["seq_quality_stats"]`. Agrega `sequences[].seq_quality` (só entradas não-nulas): `{present: False}` se nenhuma; senão `seqs_analyzed`, `with_repeats`, `with_low_complexity`, `total_repeats`, `mean_repeat_score`, `max_repeat_score` (sobre `kmer_repeat_score`), `kmer_size`.

**Decisão registrada:** o placeholder `{{VQ_READQC}}` foi **deliberadamente adiado** para a Fase 5. Registrá-lo no `_render_template` agora dispararia `(_COMPONENTS / "section_readqc.js").read_text()` num arquivo inexistente (quebra de build + warning de placeholder não-resolvido). Ele entra junto com o componente.

**Validação inline:** `_build_readqc_stats(None)` e `_build_seq_quality_stats([])` → `{present: False}`; com report populado (via `ReportExporter.export`) o `_enrich_report` produz `readqc_stats` (com `pass_rate=0.95`) e `seq_quality_stats` (`with_repeats=1`, `total_repeats=1`, `kmer_size=15`) corretos; `write_report(...)` end-to-end **sem** warning de placeholder — prova de que `{{VQ_READQC}}` não foi introduzido. `test_report.py`/`test_html_report.py` mantêm o baseline (3 falhas pré-existentes, idênticas com o arquivo revertido).

---

## Fase 5 — Relatório HTML / componentes (implementado)

### `components/section_readqc.js` (novo) — aba "Read QC"
`vqInitReadQc(readQc, readqc_stats, seq_quality_stats)` renderiza em `#section-readqc`:
KPI chips do fastp (reads after/pass rate, Q20/Q30, GC, mean length, duplication, adapter),
duas curvas D3 por-ciclo (qualidade média e GC) via `_lineChart`, e um painel-resumo dos
sinais de qualidade de sequência (com link para o Sequence Viewer). Usa `VQ.esc` e chips
`.vq-stat-chip` existentes; containers via flex inline (sem novas classes CSS).

### `components/section_viewer.js` (estendido)
- **Cobertura colorida por qualidade** — `_drawCoverageTrack` agora, quando
  `cov.quality_bins` casa com `cov.bins`, desenha uma barra por bin colorida por
  `_qualColor(q)` (verde ≥Q30 · âmbar ≥Q20 · vermelho <Q20; cinza sem reads) em vez da
  área monocromática; legenda ganha `Q̄` e o tooltip mostra `Q` por posição. Sem
  `quality_bins` (relatórios antigos) mantém a área original — retrocompatível.
- **Bandas de baixa complexidade (dustmasker)** — `_drawLowComplexityBands` desenha
  faixas verticais translúcidas (âmbar) atrás das lanes de ORF, com tooltip da região.
- **Dot plot SVG nativo** — `_dotPlotSVG(seq_quality)` gera `<svg>` puro (eixos query ×
  subject): diagonal identidade (cinza), segmentos slope +1 = repetição direta (azul),
  slope −1 = invertida (vermelho). Renderizado *lazy* no `_toggleCard`, ao lado de um
  resumo textual (`_repeatSummary`: nº diretas/invertidas, `kmer_repeat_score`,
  `low_complexity_frac`).

### `components/shell.html`
Aba `data-section="readqc"` + `<section id="section-readqc">` + `<script>{{VQ_READQC}}</script>`;
visibilidade (esconde a aba quando `!VQ_REPORT.read_qc`, espelhando o Salmon); init
`vqInitReadQc(VQ_REPORT.read_qc, VQ_REPORT.readqc_stats, VQ_REPORT.seq_quality_stats)`.

### `viralquest/html_report.py`
Placeholder `{{VQ_READQC}}` → `section_readqc.js` registrado no `_render_template` (agora
que o componente existe — cf. decisão da Fase 4).

**Validação inline:** balanceamento de chaves/parênteses/colchetes OK nos dois JS
(node indisponível no ambiente); `write_report(...)` end-to-end **sem** warning de
placeholder, com `tab-readqc`/`section-readqc`/`vqInitReadQc`/`_dotPlotSVG`/`_qualColor`/
`vq-cov-bar`/`_drawLowComplexityBands`/`quality_bins` presentes no HTML; `test_html_report.py`/
`test_report.py` mantêm o baseline (3 falhas pré-existentes). Uma demo self-contained
(`readqc_demo.html`) foi gerada com cobertura colorida (dip de qualidade), banda de baixa
complexidade e dot plot (repetição direta + invertida) para inspeção visual.

---

## Roteiro restante (planejado)

### Fase 6 — Testes
`tests/test_read_qc.py` e `tests/test_seq_quality.py` mockando subprocess (padrão de `tests/test_salmon_quant.py`): parse do JSON do fastp, intervalos do dustmasker, filtragem da diagonal no self-BLAST, score do jellyfish, geração de segmentos do dot plot.

---

## Decisões default (reversíveis)
1. **fastp = QC-only** por padrão; reads limpas a jusante só com `--trim-reads`.
2. **Feature 2 gated em `--reads`** conforme pedido (tecnicamente roda só nos contigs; pode virar standalone depois).
3. **Long reads:** fastp degrada métricas de adapter; `fastplong` plugável se for preciso paridade total com `ont/pb/hifi`.
