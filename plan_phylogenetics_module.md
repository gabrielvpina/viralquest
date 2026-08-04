# Plano — Módulo de Filogenia ancorada em perfil HMM

> Status: rascunho / em discussão. Rooting definido como **MAD** (ver §7).

## Contexto

Hoje o viralquest confirma identidade viral por hits de perfis HMM (RVDB, Vfam,
EggNOG) sobre as ORFs traduzidas, mas **para por aí**: não há contexto
evolutivo. A proposta é, quando uma sequência tem hit de modelo viral, colocá-la
em uma **árvore filogenética** junto a homólogos de referência do próprio banco
viral, para dar interpretação evolutiva/taxonômica ao achado.

A estratégia é um *phylogenetic placement* leve, ancorado no perfil: o mesmo HMM
que gerou o melhor hit define as colunas homólogas do alinhamento, e os
homólogos de referência vêm do `viralDB.dmnd` (mesma base já usada no BLASTx).

**Resultado esperado:** por sequência viral elegível, um alinhamento da região de
match + uma árvore (Newick) com rótulos taxonômicos, salva no JSON e plotada no
HTML final.

## Avaliação metodológica (resumo das decisões já tomadas)

- **Alinhar com `hmmalign` (pyhmmer), não MAFFT.** `pyhmmer.hmmalign(hmm, seqs,
  trim=True)` alinha a ORF-alvo + os homólogos **diretamente contra o mesmo
  perfil HMM**. Os *match states* do perfil *são* as colunas homólogas, e
  `trim=True` descarta os resíduos dos estados flanqueadores N/C — ou seja,
  entrega nativamente **só a região de match do modelo**, que é exatamente o
  requisito. Fica in-process (pyhmmer já é dependência), determinístico e
  consistente com o modelo que gerou o hit. MAFFT é descartado.
- **Nível de proteína.** Perfis são de aminoácido → árvore de proteína. É o
  correto para vírus divergentes; deve ficar explícito no relatório que a árvore
  reflete o domínio proteico, não o genoma.
- **Reportar suporte de ramo** (ver §7) — sem isso a topologia engana,
  especialmente em domínios curtos com poucos sítios variáveis.

## Etapas do módulo (fluxo)

Novo módulo: `viralquest/phylo.py`. **Roda após o scoring** (heurístico, passo
11 em `cli.py`; LLM opcional, passo 12) — porque o gate de escopo depende da
`classification`, que só existe depois do scoring. Nesse ponto `seq.is_viral`,
`orf.domains` e a taxonomia já estão resolvidos.

**Escopo — só `viral-unknown` (classificador heurístico).** O módulo processa
apenas sequências com `seq.heuristic_output.classification == "viral-unknown"`.
Fazer filogenia de vírus já conhecidos (`viral-known`) é redundante e de
`non-viral` é desperdício. **Sempre o classificador heurístico** — ele é
determinístico e sempre presente; o score do LLM depende fortemente do modelo
usado e pode ser instável, então **não** entra no gate. Sequência fora de
`viral-unknown` → pulada.

### 1. Selecionar o modelo — melhor score absoluto

- Para cada sequência viral, varrer todos os `orf.domains` cujo `database ∈
  {RVDB, Vfam, EggNOG}` e escolher o domínio de **maior `score` absoluto**,
  independentemente do banco de origem (`HmmDomain.score`, já disponível em
  `biodata.py`).
- Guardar: nome do modelo (`domain.target`), banco (`domain.database`), a ORF
  correspondente e a **região de match** na ORF (`domain.start` / `domain.stop`
  = `env_from` / `env_to`).
- Isolar a subsequência de aminoácidos da ORF nessa região de match — é a
  sequência-consulta que entra no alinhamento.

### 2. Recuperar o perfil HMM selecionado do arquivo mesclado

- Reusar o padrão já existente em `salmon_quant.py::PfamHkFinder._load_profiles()`:
  abrir `pyhmmer.plan7.HMMFile(<banco>.hmm)`, iterar, e ficar com o `hmm` cujo
  `hmm.name` == modelo selecionado.
- O caminho do banco vem de `cli.py::_resolve_db_paths()` (`rvdb`/`vfam`/`eggnog`).

### 3. Preparar o `viralDB.dmnd` como alvo pesquisável

- `viralDB.dmnd` está em formato DIAMOND binário; converter para FASTA
  **temporário** com `diamond getseq -d viralDB.dmnd` (subcomando confirmado no
  ambiente pixi). Arquivo temporário no scratch/outdir, removido no `finally`
  (mesma disciplina de cleanup de `coverage.py`).
- Digitalizar as sequências FASTA em um `DigitalSequenceBlock` (alfabeto amino),
  como em `HmmSequencePreparer.prepare()`.

### 4. Buscar homólogos: perfil selecionado × viralDB → top 25

- Rodar `pyhmmer.hmmsearch([hmm_selecionado], viraldb_block)` — um único perfil
  contra todo o viralDB.
- Ordenar os hits por **score** (mesma métrica da seleção) e pegar o **top 25**.
- Para cada hit, guardar: id/accession da sequência do viralDB, score, e a
  **região de match** (env_from/env_to) — só essa fatia entra no alinhamento,
  garantindo que todos os táxons cubram a mesma região do modelo.

### 5. Anexar taxonomia (família / filo) — reusar o que já existe

- Os headers do FASTA do viralDB carregam `[Species]` no final (padrão RefSeq).
  Extrair species com o mesmo regex de `BlastxResult.__post_init__`
  (`r'\[([^\[\]]+)\]\s*$'`).
- Resolver **família e filo** via `tax.py::TaxonomyMatcher.match(species)` →
  objeto `Taxonomy` (carregado de `viralTax.json.xz`, já usado no pipeline).
- Montar rótulos de tip informativos: `accession | species | family` (accession
  puro é inútil na árvore). Táxons sem resolução → `Unclassified`.

### 6. Alinhamento da região de match — `hmmalign`

- Conjunto de sequências = [ region-da-ORF-consulta ] + [ region de cada um dos
  top-25 ].
- `msa = pyhmmer.hmmalign(hmm_selecionado, digital_seqs, trim=True)`.
- Passo opcional de limpeza: descartar colunas com fração de gap acima de um
  limiar (ex.: >50%) antes de inferir a árvore, para reduzir artefato. (Sem
  dependência externa — filtro simples sobre a MSA.)
- Exportar a MSA (FASTA alinhado) como arquivo temporário para o inferidor de
  árvore.

### 7. Inferência da árvore + suporte de ramo (bootstrap)

Este é o ponto ainda em discussão. Como funciona o suporte:

- **Bootstrap clássico:** reamostrar as colunas do alinhamento com reposição N
  vezes (tipicamente 1000), reinferir a árvore em cada réplica, e medir em que %
  das réplicas cada clado aparece → valor de suporte (0–100) anexado aos nós
  internos no Newick. Caro (N inferências).
- **Ultrafast bootstrap (UFBoot, IQ-TREE):** aproximação muito mais rápida do
  bootstrap, apropriada para rodar por-sequência no pipeline.
- **SH-aLRT:** teste de razão de verossimilhança por ramo, complementar ao
  UFBoot (recomenda-se reportar os dois).
- **FastTree:** dá suporte "SH-like" local aproximado, **não** é bootstrap real,
  porém é o mais rápido.

Duas opções de ferramenta (nova dependência de binário — nenhum tree-builder
existe no ambiente pixi hoje; MAFFT **não** é necessário graças ao hmmalign):

| Ferramenta | Suporte | Velocidade | Observação |
|-----------|---------|-----------|-----------|
| **IQ-TREE 2** (DEFINIDO) | UFBoot + SH-aLRT (bootstrap real aproximado) | rápido para ~26 táxons/domínio curto | seleção automática de modelo (`-m MFP`) |
| ~~FastTree~~ | SH-like local (não é bootstrap) | mais rápido | descartado |

**Definido: IQ-TREE 2** com `-B 1000 -alrt 1000` (26 táxons numa região curta
rodam em segundos). **ModelFinder escolhe o modelo** (`-m MFP`) — e o **modelo
de substituição escolhido deve ser extraído e salvo no JSON** (parseável do
arquivo `.iqtree`, campo "Best-fit model") para constar no report final. Os
valores de suporte já saem embutidos como labels dos nós internos no Newick.

**Rooting — MAD (Minimal Ancestor Deviation):**
- Definido: usar **MAD rooting** em vez de midpoint. MAD **não assume relógio
  molecular** (ao contrário do midpoint) e tolera heterotaquia — avalia cada
  ramo como raiz candidata e escolhe o que minimiza os desvios relativos das
  distâncias ancestral-descendente entre todos os pares de folhas. Melhor para
  vírus divergentes, onde o midpoint enraíza viciado em táxons de ramo longo.
- É **pós-processamento**: entra a Newick não-enraizada com comprimentos de ramo
  (o IQ-TREE já produz), sai a Newick enraizada. Não toca na inferência.
- **Ambiguity index (AI):** o MAD devolve uma medida de confiança do
  enraizamento (quão melhor a melhor raiz é vs. a segunda). Guardar no JSON e
  exibir no HTML como sinal de confiabilidade da raiz.
- **Implementação:** usar o script de referência `mad.py` (Tria et al. 2017,
  lab Dagan) — Python self-contained, **só numpy** (já é dependência), lê Newick
  → escreve Newick enraizada + AI. Vendorizar (checar licença) ou chamar via
  subprocess. Nenhum binário novo. Alternativas descartadas: MADroot (C++, só
  vale p/ árvores grandes — não é o caso, ≤26 táxons) e PyMAD (repo marcado como
  "not working yet"). Referência: Tria, Landan & Dagan, *Nat Ecol Evol* 2017.

### 8. Saída: Newick → JSON → HTML

- **Dados brutos em disco (sub-diretórios por sequência):** dois diretórios de
  saída de nível superior, cada um com um sub-diretório por sequência
  (nome sanitizado como em `coverage._export_tsv`):
  - `outdir/alignment/<seq_id>/` — arquivos brutos do alinhamento: a MSA do
    `hmmalign` (FASTA alinhado da região de match).
  - `outdir/phylogeny/<seq_id>/` — todos os arquivos brutos da inferência:
    saída do IQ-TREE (`.treefile`, `.iqtree`, `.log`, `.bionj`, etc.) **e** a
    Newick enraizada por MAD (`<seq_id>.rooted.nwk` — o "arquivo plano"
    canônico). Nada é limpo (`cleanup=False` para esses); só os temporários
    intermediários (ex.: dump FASTA do dmnd) são removidos.
  - Ambos criados no `OutputOrganizer` (`output.py`), no mesmo padrão de
    `self.hmm_dir = outdir / "hmm"`.
- **JSON (Newick cru + árvore aninhada — DEFINIDO):** adicionar um bloco por
  sequência ao relatório via `exporter.py` (padrão de
  `salmon_report`/`CoverageProfile`), contendo **as duas** representações:
  - `newick`: string Newick crua enraizada (MAD) — o "arquivo plano", canônico,
    parseável por qualquer ferramenta externa.
  - `tree`: dict **aninhado** parseado em Python, no formato que o D3
    (`d3.hierarchy`) consome direto — sem parser de Newick em JS. É onde ficam
    os metadados por tip (taxonomia família/filo, accession, flag `is_query`) e
    o suporte por nó, que o Newick não carrega de forma limpa. Espelha o padrão
    já existente em `html_report.py::_build_taxonomy_tree` (linha 220).
  - Metadados do bloco (para o report final): modelo HMM usado + banco de
    origem, região de match, **modelo de substituição escolhido pelo ModelFinder
    do IQ-TREE**, valores de suporte, ambiguity index do MAD, e a lista de tips
    (top-25 + consulta) com taxonomia. Também os caminhos relativos dos dados
    brutos em `alignment/<seq_id>/` e `phylogeny/<seq_id>/`.
- **HTML:** novo componente em `viralquest/components/` (ex.:
  `section_phylo.js`), plotado com JS inline (sem CDN — o HTML é self-contained,
  CSP estrita). Registrar o placeholder em `html_report.py` (`{{VQ_PHYLO}}`),
  como os demais componentes.

## Arquivos a criar / modificar

- **Novo:** `viralquest/phylo.py` — pipeline do módulo (seleção, dmnd→fasta,
  hmmsearch top-25, taxonomia, hmmalign, árvore, export).
- **Novo:** `viralquest/components/section_phylo.js` — render da árvore no HTML.
- `viralquest/biodata.py` — novo dataclass (ex.: `PhyloResult`) anexado à
  `NucSequence`.
- `viralquest/cli.py` — nova flag (ex.: `--phylo` / `--phylo-topn`), chamada do
  módulo no `_run_pipeline`, e passo em `_build_steps`.
- `viralquest/exporter.py` — serializar `PhyloResult` no JSON.
- `viralquest/output.py` — criar `self.alignment_dir = outdir / "alignment"` e
  `self.phylogeny_dir = outdir / "phylogeny"` (padrão de `self.hmm_dir`), com
  helpers para os sub-diretórios por sequência.
- `viralquest/html_report.py` — placeholder + injeção do componente.
- `viralquest/setup_env.py` + config pixi — adicionar o binário de árvore
  (IQ-TREE 2 ou FastTree).

## Reaproveitamento (não reinventar)

- Filtro de modelo por nome: `salmon_quant.py::PfamHkFinder._load_profiles()`.
- Preparo de bloco digital: `hmm.py::HmmSequencePreparer.prepare()`.
- Resolução família/filo: `tax.py::TaxonomyMatcher` + `viralTax.json.xz`.
- Regex de species do header: `biodata.py::BlastxResult.__post_init__`.
- Disciplina de tmp/cleanup e export TSV/plano: `coverage.py`.
- Caminhos dos bancos: `cli.py::_resolve_db_paths()`.

## Dependências novas

- Um inferidor de árvore (IQ-TREE 2 recomendado; FastTree alternativa) — via
  bioconda no ambiente pixi.
- **Nenhuma** para alinhamento (hmmalign já vem no pyhmmer) e **nenhuma** para
  dump do dmnd (`diamond getseq` já disponível).
- **MAD rooting:** nenhum binário — script `mad.py` (numpy puro) vendorizado ou
  chamado via subprocess.

## Verificação (end-to-end)

1. Rodar o pipeline com `--reads`/entrada de teste que contenha ao menos uma
   sequência com hit viral forte (usar `tests/mock/`), com a flag do módulo.
2. Conferir: `outdir/alignment/<seq_id>/` tem a MSA do hmmalign e
   `outdir/phylogeny/<seq_id>/` tem a saída bruta do IQ-TREE + a Newick MAD
   (`<seq_id>.rooted.nwk`), que é Newick válida (parsear com um parser leve no
   teste).
3. Conferir no JSON: bloco phylo com o **modelo HMM** usado, o **modelo de
   substituição do IQ-TREE** (ex.: `LG+G4`), top-25 com taxonomia resolvida
   (família/filo não todos `Unclassified`), valores de suporte nos nós e o
   ambiguity index do MAD.
4. Confirmar que só sequências `viral-unknown` (heurístico) geraram filogenia —
   `viral-known` e `non-viral` devem ser puladas.
5. Abrir o HTML e confirmar que a árvore renderiza com rótulos legíveis.
6. Caso-limite: sequência `viral-unknown` sem nenhum hit RVDB/Vfam/EggNOG →
   módulo pula graciosamente (sem árvore, sem erro), como os demais passos
   opcionais.

## Decisões (fechadas)

1. **Rooting:** MAD rooting (§7).
2. **Ferramenta de árvore:** IQ-TREE 2 (`-B 1000 -alrt 1000 -m MFP`); o modelo
   escolhido pelo ModelFinder é extraído do `.iqtree` e salvo no JSON.
3. **Escopo:** apenas `seq.heuristic_output.classification == "viral-unknown"`
   (classificador **heurístico**, nunca o LLM), rodando após o scoring.
4. **JSON:** Newick cru **+** árvore aninhada parseada (ambos) + modelo IQ-TREE
   + metadados de filogenia para o report final.
5. **Dados brutos em disco:** `outdir/alignment/<seq_id>/` (MSA do hmmalign) e
   `outdir/phylogeny/<seq_id>/` (saída bruta do IQ-TREE + Newick MAD), um
   sub-diretório por sequência.
