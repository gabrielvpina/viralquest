# Sequence Quality module

Um módulo acionado por `--reads` que agrega sinais de qualidade das contigs montadas:

1. **Coloração da track de cobertura** pela qualidade de base, reaproveitando o BAM já produzido pelo `coverage.py`.
2. **Qualidade estrutural da sequência** — regiões de baixa complexidade (dustmasker), repetições internas (self-BLASTn), repetitividade de k-mers (jellyfish) e um **dot plot em SVG nativo**.

Ferramentas fixadas: `dustmasker` · `blastn` (self-BLAST) · `jellyfish` · dot plot SVG nativo · coloração da cobertura (`samtools mpileup`).

> **Histórico:** uma seção de QC de reads da biblioteca (fastp → aba "Read QC") foi
> implementada e depois **removida** — o fastp não tinha função na visualização da
> qualidade das contigs (que vem do `coverage.py` via `samtools mpileup`, não do fastp).
> O resumo de qualidade de sequência que ficava naquela aba virou um **card condicional
> na aba General Stats**. As menções a `fastp` / `read_qc.py` / aba "Read QC" abaixo
> descrevem o estado anterior e ficam como registro do histórico de implementação.

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
- **QC-only:** as reads filtradas são **descartadas** (`-o os.devnull`, e `-O` no paired) — o objetivo é informar a qualidade, não trimar. Só os relatórios `fastp.json`/`fastp.html` vão para `outdir/read_qc`; nenhum FASTQ limpo é gerado. O pipeline sempre consome as reads originais.
- **Parse do JSON:** `reads/bases before/after` (`summary`), `q20/q30_rate`, `gc_content×100`, `read1_mean_length`, `duplication.rate`, `adapter_cutting.adapter_trimmed_reads / reads_before`, e as curvas por-ciclo `read1_after_filtering.quality_curves.mean` e `content_curves.GC`.
- **`_downsample`:** curvas por-ciclo reduzidas a ≤ 300 pontos (médias de bloco) — evita SVG gigante em long reads.
- Falha graciosa: fastp ausente, nº de arquivos ≠ 1–2, ou JSON ilegível → `None`.

### `viralquest/seq_quality.py` → `SequenceQualityPipeline`

`run(viral_seqs, outdir) -> dict[seq_id, SequenceQuality]`. Escreve um `viral_ref.fasta` multi-record e roda os 4 sinais; cada sequência sempre recebe uma entrada (listas vazias se um sinal faltar).

- **dustmasker** (`_dustmask`): `dustmasker -in ref -outfmt interval` → parse de `>seqid` + linhas `a - b` (0-based, inclusivas) → `low_complexity_regions`; `low_complexity_frac` = bases mascaradas / comprimento.
- **self-BLASTn** (`_self_blast`): uma chamada `blastn -query ref -subject ref -dust yes -evalue 1e-5 -word_size 11 -outfmt "6 qseqid sseqid qstart qend sstart send pident length sstrand"`. Filtra: hits cross-sequência (`qseqid≠sseqid`), a diagonal self-completa (plus + coords idênticas) e duplicatas simétricas (A→B/B→A) via chave de triângulo superior. Coords `minus` normalizadas (`s_lo/s_hi`). → `SelfRepeat` (strand `plus`=direta, `minus`=invertida).
- **jellyfish** (`_jellyfish`, **por sequência**): `jellyfish count -m K -s 10M -C -o out.jf seq.fa` + `jellyfish histo`. `kmer_repeat_score` = `(total − singletons) / total` (fração de *instâncias* de k-mer que recorrem, multiplicidade > 1); expõe também `kmer_distinct`/`kmer_total`. `K` default = 15.
- **dot plot SVG nativo** (`_dotplot`, puro-Python): constrói um índice de k-mers forward sobre **todas** as posições (`word` default 12) e o consulta a partir de sementes amostradas a cada `stride = L // 500` bases. Indexar todas as posições (e não só as amostradas) é essencial: garante que uma repetição seja detectada **qualquer que seja o offset** — a versão anterior amostrava também o índice e perdia repetições (diretas e invertidas) cujas duas cópias não caíam na mesma grade de amostragem quando `stride > 1`. Matches forward → offset diagonal `d = j − i ≥ 0` (incl. diagonal identidade `d=0`); matches contra o **reverse-complement** → anti-diagonal `s = i + j` (repetições invertidas). k-mers ultra-repetitivos (> `_MAX_KMER_HITS = 200` ocorrências) são pulados. `_emit_segments` funde âncoras colineares (gap ≤ 2·stride) em segmentos `[x1,y1,x2,y2]` (bp), limitados a 4000. → `DotPlot`. **No relatório**, a diagonal identidade é tracejada/recuada e as repetições (azul = direta, vermelho = invertida) são desenhadas por cima, grossas e opacas, para não se confundirem com a diagonal.

**Validação inline:** parse do fastp (incl. adapter=0 p/ long read) e curva downsampleada OK; dot plot detecta corretamente repetição direta (diagonal) e invertida (anti-diagonal) numa sequência sintética; parsers do dustmasker e do self-BLASTn conferidos com saída mockada (diagonal trivial, dup simétrica e cross-seq descartadas; direta e invertida mantidas com coords 0-based normalizadas). Ambos os módulos importam limpos. Testes formais ficam na Fase 6.

---

## Fase 3 — CLI + export JSON (implementado)

### Argumentos novos (`cli.py`, grupo *"read / sequence quality (optional, needs --reads)"*)
- `--skip-read-qc` — pula o passo de read QC (fastp); os sinais de qualidade de sequência ainda rodam.
- `--kmer K` — tamanho de k-mer do score do jellyfish (default 15).

Validações (`_validate_args`): `--skip-read-qc` exige `--reads`; `--kmer >= 2`.

### Passos e ordem (`_build_steps` ⟷ `_run_pipeline`)
Helper novo `_read_qc_enabled(args)` = `bool(args.reads) and not args.skip_read_qc`. O bloco de reads passa a ter **quatro** passos, na mesma ordem em ambas as funções:

1. **Read QC — fastp** (`_read_qc_enabled`) — roda `ReadQcPipeline` em `outdir/read_qc` → `read_qc_report`. QC-only: fastp descarta as reads filtradas; o pipeline sempre usa as reads originais.
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
  faixas verticais translúcidas (âmbar) atrás das lanes de ORF. São **puramente visuais**
  (`pointer-events: none`) para **não bloquear** o hover da track de cobertura por baixo
  (o detalhe da região fica no realce do FASTA preview e no dot plot).
- **Dot plot SVG nativo** — `_dotPlotSVG(seq_quality)` gera `<svg>` puro (eixos query ×
  subject). **Fonte única de verdade:** um classificador compartilhado (`_segClass`)
  rotula cada segmento como `diag` (diagonal identidade, cinza tracejado), `lowcx`
  (auto-match dentro de região dustmask → âmbar), `direct` (slope +1 → azul) ou
  `inverted` (slope −1 → vermelho). A **legenda** (`_repeatSummary` via `_dotplotRepeats`)
  e o **realce no FASTA** (`_fastaHighlightedHTML`) contam/pintam exatamente os segmentos
  `direct`/`inverted` desenhados — então *o que aparece = o que é contado = o que é
  realçado*. Auto-matches de baixa complexidade aparecem em âmbar no dot plot e são
  representados pelo `low-complexity %`, **não** inflam a contagem de repetições
  estruturais. Renderizado *lazy* no `_toggleCard`.

  > **Nota de padronização:** antes, a legenda contava os HSPs do self-BLASTn
  > (`self_repeats`, corte por e-value) enquanto o dot plot vinha do hashing de k-mers —
  > por isso repetições pequenas apareciam no plot mas não na contagem. Agora ambos
  > derivam dos segmentos do dot plot; `self_repeats` (BLAST) segue no JSON como dado
  > suplementar.

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

## Parâmetros de qualidade — referência completa

Consolida **todos** os parâmetros e mostradores de qualidade do módulo: o que cada um
significa, como é calculado e como aparece no relatório. Inclui as últimas modificações
no dot plot e nos mostradores.

### A. Read QC da biblioteca (fastp) — **REMOVIDO**

Esta seção (fastp → `report["read_qc"]` → `readqc_stats` → aba "Read QC") foi **removida**.
O fastp não participava da visualização da qualidade das contigs; a coloração de qualidade
no genome viewer vem do `coverage.py` (`samtools mpileup`), descrito em **B**. Removidos:
`read_qc.py`, `section_readqc.js`, a aba/seção no `shell.html`, o placeholder `{{VQ_READQC}}`,
`_build_readqc_stats`, a dataclass `ReadQcReport`, a dependência `fastp` (`pixi.toml` /
`setup_env.py`) e a flag `--skip-read-qc`.

### B. Cobertura de reads — por sequência

Fonte: minimap2 → BAM → `samtools mpileup` (uma passada dá profundidade **e** qualidade).
Campos em `sequences[].coverage` (`CoverageProfile`).

| Parâmetro | Significado |
|-----------|-------------|
| `bins` | profundidade média por bin (≤ 600 pontos), alinhada a `quality_bins` |
| `mean_depth` / `max_depth` | profundidade média / máxima |
| `breadth_1x` | fração de bases com profundidade ≥ 1 |
| `cv` | `stdev/mean` da profundidade — uniformidade (menor = mais uniforme; degrau ⇒ possível quimera) |
| `quality_bins` | qualidade de base **média** (Phred) por bin — média só sobre posições cobertas |
| `mean_quality` | qualidade de base média (Phred) na sequência |

**Mostrador (genome viewer):** a track de cobertura tem **altura = profundidade** e
**cor = qualidade** (uma barra por bin), via `_qualColor`:

| Qualidade média do bin | Cor |
|------------------------|-----|
| Q ≥ 30 | verde (`--vq-success`) |
| 20 ≤ Q < 30 | âmbar (`--vq-warning`) |
| 0 < Q < 20 | vermelho (`--vq-danger`) |
| sem reads | cinza (`--vq-text-3`) |

Hover mostra posição, profundidade, **Q** naquele ponto, média e máx. Legenda mostra
`Q̄` (qualidade média). Relatórios antigos (sem `quality_bins`) mantêm a área monocromática.

### C. Qualidade estrutural da sequência — por sequência

Campos em `sequences[].seq_quality` (`SequenceQuality`).

| Parâmetro | Fonte | Significado |
|-----------|-------|-------------|
| `low_complexity_regions` | dustmasker (SDUST) | intervalos `[início,fim]` (0-based, incl.) de baixa complexidade |
| `low_complexity_frac` | derivado | fração de bases dentro dessas regiões |
| `kmer_size` | jellyfish (`--kmer`, default 15) | tamanho do k-mer do score |
| `kmer_distinct` / `kmer_total` | jellyfish `histo` | nº de k-mers distintos / total de instâncias |
| `kmer_repeat_score` | jellyfish | fração de instâncias de k-mer que recorrem: `(total − singletons) / total` |
| `self_repeats` | self-BLASTn | HSPs internos (direta/invertida) — **dado suplementar no JSON**, ver nota de padronização |
| `dotplot` | k-mer hashing (`_dotplot`) | segmentos `[x1,y1,x2,y2]` (bp) para o dot plot SVG |

#### Dot plot — detecção e classificação (últimas modificações)

- **Detecção independente de offset:** o índice de k-mers cobre **todas** as posições da
  sequência; as sementes são amostradas a cada `stride = L // 500` bases, mas a busca usa o
  índice completo. Assim uma repetição é detectada **qualquer que seja o offset** — a versão
  anterior amostrava também o índice e perdia repetições (diretas e invertidas) quando as
  duas cópias não caíam na mesma grade (`stride > 1`).
- **Teto anti-ruído:** k-mers com mais de `_MAX_KMER_HITS = 200` ocorrências são ignorados;
  no máximo 4000 segmentos por sequência.
- **Classificador único (`_segClass`)** — usado pelo dot plot, pela contagem da legenda e
  pelo realce do FASTA, garantindo que **os três concordem**:

  | Classe | Geometria | Cor no dot plot | Conta como repetição? |
  |--------|-----------|-----------------|-----------------------|
  | `diag` | diagonal identidade (`y = x`) | cinza tracejado | não |
  | `lowcx` | segmento com ponto médio dentro de região dustmask | âmbar | não (vai p/ `low-complexity %`) |
  | `direct` | slope +1 (fora de baixa complexidade) | azul, grosso | **sim** |
  | `inverted` | slope −1 (fora de baixa complexidade) | vermelho, grosso | **sim** |

> **Nota de padronização.** A contagem da legenda e o realce colorido no FASTA derivam
> dos **mesmos segmentos desenhados** (`_dotplotRepeats`), não mais do self-BLASTn — por
> isso *o que aparece no dot plot = o que é contado = o que é realçado*. Auto-similaridade
> de baixa complexidade é `lowcx` (âmbar, não contada), já representada por `low-complexity %`.

**Mostradores (genome viewer, card da sequência):**
- **Genome viewer:** bandas âmbar de baixa complexidade atrás das ORFs (visuais, sem
  bloquear o hover da cobertura).
- **Coluna esquerda:** dot plot SVG (300×300) + legenda empilhada (nº diretas, invertidas,
  `low-complexity %`, `k{n} repeat score`) com um **info-hover “?”** descrevendo cada fator.
- **Coluna central:** bloco **Top hit & taxonomy** — melhor hit BLASTx (NR preferido, senão
  RefSeq): título, accession, identidade, cobertura, e-value, bit score; + linhagem
  família/gênero/espécie de `sequences[].taxonomy`.
- **Coluna direita:** **FASTA preview** rolável (altura 300px, simétrica ao dot plot), com as
  regiões realçadas nas mesmas cores (âmbar/azul/vermelho); header contém **apenas o nome
  original** do contig (sem sufixo de espécie); botão **Copy FASTA** (header + sequência).

### D. Stats agregados (`html_report.py`)

- `seq_quality_stats` — agrega `seq_quality` das sequências: `seqs_analyzed`, `with_repeats`,
  `with_low_complexity`, `total_repeats`, `mean_repeat_score`, `max_repeat_score`, `kmer_size`.
  **`with_repeats`/`total_repeats` derivam dos segmentos do dot plot** (helper
  `_structural_repeats`, mesma classificação `_segClass` do viewer — exclui diagonal e
  auto-matches de baixa complexidade), então batem com a contagem da legenda do viewer.

**Mostrador:** card **condicional** "Sequence Quality" na aba **General Stats**
(`section_stats.js`, `#stats-seqqual-card`, aparece só quando `seq_quality_stats.present`):
KPI grande = `seqs_analyzed`; mini-rows = `with_repeats`, `with_low_complexity`,
`total_repeats`, `repeat score (k{n})` = `mean_repeat_score` (com `max`). Substitui o antigo
resumo que ficava na aba "Read QC" (removida).

---

## Roteiro restante (planejado)

### Fase 6 — Testes
`tests/test_read_qc.py` e `tests/test_seq_quality.py` mockando subprocess (padrão de `tests/test_salmon_quant.py`): parse do JSON do fastp, intervalos do dustmasker, filtragem da diagonal no self-BLAST, score do jellyfish, geração de segmentos do dot plot.

---

## Decisões default (reversíveis)
1. **fastp = QC-only**; as reads filtradas são descartadas (sem FASTQ limpo no output) — o módulo informa a qualidade, não trima.
2. **Feature 2 gated em `--reads`** conforme pedido (tecnicamente roda só nos contigs; pode virar standalone depois).
3. **Long reads:** fastp degrada métricas de adapter; `fastplong` plugável se for preciso paridade total com `ont/pb/hifi`.
