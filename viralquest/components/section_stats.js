/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_stats.js — the two overview tabs
   Run    : how the pipeline ran — workflow, KPIs, then grouped cards
            (Detection · Scoring · Data & assembly).
   Virome : what was found — KPIs, the detected-viruses table, then grouped
            cards (Taxonomy · Similarity to known viruses · Sequences & proteins).
   Cards sit in titled groups on a 3-column grid (.vq-group / .vq-grid).
   ============================================================ */

function vqInitStats(report) {
  _initRunTab(report);
  _initViromeTab(report);
}

/* Re-render a panel's charts when it is first shown or resized: a hidden tab
   has zero width, so its charts are drawn at a fallback size until revealed. */
function _redrawOnResize(el, redraw) {
  let lastW = 0;
  new ResizeObserver(() => {
    const w = el.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    redraw();
  }).observe(el);
}

/* Titled group of cards on a grid. `slots` = grid columns the cards fill
   (a span-2 card counts as 2): 1–2 slots use a 2-column grid so a short group
   does not leave an empty third column. Empty groups vanish. */
function _group(title, hint, cardsHtml, slots = 3) {
  if (!cardsHtml.trim() || !slots) return '';
  return `
    <section class="vq-group">
      <div class="vq-group__head">
        <h3 class="vq-group__title">${VQ.esc(title)}</h3>
        ${hint ? `<span class="vq-group__hint">${VQ.esc(hint)}</span>` : ''}
      </div>
      <div class="vq-grid${slots <= 2 ? ' vq-grid--2' : ''}">${cardsHtml}</div>
    </section>`;
}

const _fmtNum = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();

/* Chart height that fills its card body. Grid rows stretch a card to its
   tallest neighbour; drawing at a fixed height would leave a blank band above
   the chart. Call after clearing the wrap, so the body's height reflects the
   row rather than the previous drawing. */
function _fillHeight(wrap, minH) {
  const body = wrap.parentElement;
  const h = body ? body.clientHeight - 8 : 0;
  return Math.max(minH, Math.floor(h));
}

// ════════════════════════════════════════════════════════════════════════
//  RUN — how the pipeline ran: workflow, filtering, databases, scoring, QC
// ════════════════════════════════════════════════════════════════════════

function _initRunTab(report) {
  const el = document.getElementById('section-run');
  if (!el) return;

  const esc     = VQ.esc;
  const meta    = report.meta             || {};
  const sum     = report.summary          || {};
  const blast   = report.blast_stats      || {};
  const hmm     = report.hmm_stats        || {};
  const llm     = report.llm_stats        || {};
  const heur    = report.heuristic_stats  || {};
  const salmon  = report.salmon_stats     || {};
  const seqQual = report.seq_quality_stats || {};
  const seqs    = report.sequences        || [];
  const wf      = report.workflow || (report.pipeline_stats || {}).workflow || null;
  const hasWf   = !!(wf && (wf.steps || []).length);

  const totalSeqs      = sum.total_sequences ?? seqs.length;
  const confirmedViral = sum.confirmed_viral ?? seqs.filter(s => s.is_viral).length;
  const totalOrfs      = sum.total_orfs      ?? seqs.reduce((a, s) => a + (s.orfs || []).length, 0);
  const totalClusters  = sum.total_clusters  ?? (report.clusters || []).length;
  const hmmTotal = (hmm.rvdb_hits || 0) + (hmm.vfam_hits || 0)
                 + (hmm.eggnog_hits ?? hmm.eggnong_hits ?? 0) + (hmm.pfam_hits || 0);

  const detection = `
    <div class="vq-chart-card" id="stats-funnel-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Detection Pipeline</div>
          <div class="vq-chart-card__sub">sequences at each filtering step</div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="padding:6px 14px 10px">
        <div id="stats-funnel-svg" style="width:100%"></div>
      </div>
    </div>
    <div class="vq-chart-card" id="stats-hmm-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">HMM Database Hits</div>
          <div class="vq-chart-card__sub">threshold-passing hits per database</div>
        </div>
        <div class="vq-chart-card__big" id="stats-hmm-total">${_fmtNum(hmmTotal)}</div>
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-hmm-svg" style="width:100%"></div>
      </div>
    </div>`;

  const scoring = `
    ${heur.present ? `
    <div class="vq-chart-card" id="stats-heur-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Heuristic Scoring</div>
          <div class="vq-chart-card__sub">rule-based · deterministic</div>
        </div>
        <div class="vq-chart-card__big">${heur.avg_score != null ? Number(heur.avg_score).toFixed(1) : '—'}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start;gap:0">
        <div id="stats-heur-gauge" style="width:100%;margin-bottom:10px"></div>
        ${_miniRow('Scored',        _fmtNum(heur.scored))}
        ${_miniRow('Viral known',   _fmtNum(heur.viral_known))}
        ${_miniRow('Viral unknown', _fmtNum(heur.viral_unknown))}
        ${heur.non_viral ? _miniRow('Non-viral', _fmtNum(heur.non_viral)) : ''}
      </div>
    </div>` : ''}
    ${llm.present ? `
    <div class="vq-chart-card" id="stats-llm-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">LLM Scoring</div>
          <div class="vq-chart-card__sub">${esc(llm.model || '—')}${llm.mode ? ' · ' + esc(llm.mode) + '-token mode' : ''}</div>
        </div>
        <div class="vq-chart-card__big">${llm.avg_score != null ? Number(llm.avg_score).toFixed(1) : '—'}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('Sent to LLM',   _fmtNum(llm.scored))}
        ${_miniRow('Viral known',   _fmtNum(llm.viral_known))}
        ${_miniRow('Viral unknown', _fmtNum(llm.viral_unknown))}
        ${llm.api_error   ? _miniRow('API errors',        _fmtNum(llm.api_error))   : ''}
        ${llm.parse_error ? _miniRow('Invalid responses', _fmtNum(llm.parse_error)) : ''}
      </div>
    </div>` : ''}`;

  const dataQc = `
    ${seqQual.present ? `
    <div class="vq-chart-card" id="stats-seqqual-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Sequence Quality</div>
          <div class="vq-chart-card__sub">structural signals across sequences</div>
        </div>
        <div class="vq-chart-card__big">${_fmtNum(seqQual.seqs_analyzed)}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('With repeats',        _fmtNum(seqQual.with_repeats))}
        ${_miniRow('With low complexity', _fmtNum(seqQual.with_low_complexity))}
        ${_miniRow('Total repeats',       _fmtNum(seqQual.total_repeats))}
        ${_miniRow(`Repeat score (k${seqQual.kmer_size ?? '?'})`,
                   seqQual.mean_repeat_score != null
                     ? (seqQual.mean_repeat_score * 100).toFixed(1) + '%'
                     + (seqQual.max_repeat_score != null
                         ? ' · max ' + (seqQual.max_repeat_score * 100).toFixed(1) + '%' : '')
                     : '—')}
      </div>
    </div>` : ''}
    ${sum.cap3_used ? `
    <div class="vq-chart-card" id="stats-cap3-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">CAP3 Assembly</div>
          <div class="vq-chart-card__sub">input re-assembled before analysis</div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('Contigs',  _fmtNum(sum.cap3_contigs))}
        ${_miniRow('Singlets', _fmtNum(sum.cap3_singlets))}
      </div>
    </div>` : ''}
    ${salmon.present ? `
    <div class="vq-chart-card" id="stats-salmon-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Read Mapping · Salmon</div>
          <div class="vq-chart-card__sub">
            ${_fmtNum(salmon.total_reads)} reads
            ${salmon.pathway ? ' · ' + esc(salmon.pathway) + ' pathway' : ''}
          </div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="flex-direction:column;justify-content:flex-start;gap:0">
        <div id="stats-salmon-svg" style="width:100%"></div>
        <div id="stats-salmon-rows" style="width:100%;padding:0 4px 4px"></div>
      </div>
    </div>` : ''}`;

  el.innerHTML = `
    ${hasWf ? '' : `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Run</div>
        <div class="vq-section-sub">
          ${esc(meta.input_file?.name || 'Sample report')}
          &nbsp;·&nbsp; run ${esc(meta.timestamp ? new Date(meta.timestamp).toLocaleString() : '—')}
          &nbsp;·&nbsp; ViralQuest v${esc(meta.viralquest_version || '?')}
        </div>
      </div>
    </div>`}

    <div class="vq-stats-page">
      ${hasWf ? _workflowCardHtml(wf, meta.viralquest_version) : ''}

      <div class="vq-stats-row vq-stats-row--top">
        ${_chip('Input Sequences', _fmtNum(totalSeqs), 'accent',
                meta.input_file?.name ? esc(meta.input_file.name) : '')}
        ${_chip('ORFs Predicted', _fmtNum(totalOrfs), '', 'across confirmed sequences')}
        ${_chip('Confirmed Viral', _fmtNum(confirmedViral), 'success',
                _pct(confirmedViral, totalSeqs) + ' of input')}
        ${hasWf
          ? _chip('Runtime', _fmtDur(wf.total_seconds), '',
                  `${(wf.steps || []).filter(x => x.status !== 'skipped').length} steps run`)
          : _chip('Clusters', _fmtNum(totalClusters), 'accent', 'species clusters')}
      </div>

      ${_group('Detection & scoring', 'how sequences were filtered, confirmed and scored',
               detection + scoring, 2 + !!heur.present + !!llm.present)}
      ${_group('Data & assembly', 'input processing and read-level quality',
               dataQc, !!seqQual.present + !!sum.cap3_used + !!salmon.present)}
    </div>`;

  if (hasWf) _renderWorkflow(wf);
  const draw = () => {
    _renderFunnel(blast);
    _renderHMMBars(hmm);
    if (heur.present)   _renderHeuristicGauge(heur);
    if (salmon.present) _renderSalmonChart(salmon);
  };
  draw();
  _redrawOnResize(el, () => { draw(); if (hasWf) _renderWorkflow(wf); });
}

// ════════════════════════════════════════════════════════════════════════
//  VIROME — what was found: viruses, taxonomy, similarity, sequences
// ════════════════════════════════════════════════════════════════════════

function _initViromeTab(report) {
  const el = document.getElementById('section-virome');
  if (!el) return;

  const esc       = VQ.esc;
  const seqs      = report.sequences || [];
  const viral     = seqs.filter(s => s.is_viral);
  const hasBlastn = viral.some(s => (s.blastn_hits || []).length > 0);
  const domainDbs = _domainDatabases(seqs);
  const viruses   = _virusRows(viral);
  const families  = new Set(viral.map(s => s.taxonomy?.family).filter(Boolean));
  const known     = viruses.filter(v => v.status === 'known').length;

  const table = `
    <div class="vq-chart-card vq-span-3" id="stats-viruses-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Detected Viruses</div>
          <div class="vq-chart-card__sub" id="stats-viruses-sub"></div>
        </div>
        <div class="vq-vt-tools">
          <input class="vq-input vq-input--sm" type="search" id="stats-viruses-search"
                 placeholder="Filter species, family…" aria-label="Filter detected viruses">
          <div class="vq-toggle" id="stats-viruses-status" role="tablist" aria-label="Status">
            <button class="vq-toggle__btn active" type="button" data-status="all">All</button>
            <button class="vq-toggle__btn" type="button" data-status="known">Known</button>
            <button class="vq-toggle__btn" type="button" data-status="related">Related</button>
            <button class="vq-toggle__btn" type="button" data-status="divergent">Divergent</button>
          </div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="justify-content:flex-start">
        <div class="vq-vt-wrap" id="stats-viruses-table"></div>
      </div>
    </div>`;

  const taxonomy = `
    <div class="vq-chart-card" id="stats-tax-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Viral Family Distribution</div>
          <div class="vq-chart-card__sub" id="stats-tax-sub">across confirmed sequences</div>
        </div>
        <div class="vq-toggle" role="tablist" aria-label="Taxonomy rank">
          <button class="vq-toggle__btn active" type="button" data-rank="family">Family</button>
          <button class="vq-toggle__btn"        type="button" data-rank="phylum">Phylum</button>
        </div>
      </div>
      <div class="vq-chart-card__body">
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;align-items:center">
          <div id="stats-tax-svg" style="display:flex;justify-content:center;align-items:center"></div>
          <div id="stats-tax-legend" class="vq-chart-legend"></div>
        </div>
      </div>
    </div>
    <div class="vq-chart-card" id="stats-genome-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Genome Types</div>
          <div class="vq-chart-card__sub" id="stats-genome-sub">sequences per genome composition</div>
        </div>
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-genome-svg" style="width:100%"></div>
      </div>
    </div>
    <div class="vq-chart-card" id="stats-nrclass-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title" id="stats-nrclass-title">NR Classification</div>
          <div class="vq-chart-card__sub" id="stats-nrclass-sub">virus · phage breakdown</div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="justify-content:flex-start">
        <div id="stats-nrclass-legend" style="width:100%"></div>
      </div>
    </div>`;

  const similarity = `
    <div class="vq-chart-card vq-span-2" id="stats-scatter-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Identity × Coverage</div>
          <div class="vq-chart-card__sub" id="stats-scatter-sub">best hit per sequence · click a point to open it</div>
        </div>
        ${hasBlastn ? `
        <div class="vq-toggle" id="stats-scatter-toggle" role="tablist" aria-label="Search">
          <button class="vq-toggle__btn active" type="button" data-scatter="blastx">BLASTx</button>
          <button class="vq-toggle__btn"        type="button" data-scatter="blastn">BLASTn</button>
        </div>` : ''}
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-scatter-svg" style="width:100%"></div>
      </div>
    </div>
    <div class="vq-stack">
      <div class="vq-chart-card" id="stats-identity-card">
        <div class="vq-chart-card__head">
          <div class="vq-chart-card__title" id="stats-blastx-title">BLASTx Identity</div>
          <div class="vq-toggle" id="blast-toggle-x" role="tablist">
            <button class="vq-toggle__btn active" type="button" data-blast-metric="identity">Identity</button>
            <button class="vq-toggle__btn"        type="button" data-blast-metric="coverage">Coverage</button>
          </div>
        </div>
        <div class="vq-chart-card__body" style="padding-top:4px">
          <div id="stats-identity-svg" style="width:100%"></div>
        </div>
      </div>
      ${hasBlastn ? `
      <div class="vq-chart-card" id="stats-blastn-card">
        <div class="vq-chart-card__head">
          <div class="vq-chart-card__title" id="stats-blastn-title">BLASTn Identity</div>
          <div class="vq-toggle" id="blast-toggle-n" role="tablist">
            <button class="vq-toggle__btn active" type="button" data-blast-metric="identity">Identity</button>
            <button class="vq-toggle__btn"        type="button" data-blast-metric="coverage">Coverage</button>
          </div>
        </div>
        <div class="vq-chart-card__body" style="padding-top:4px">
          <div id="stats-blastn-svg" style="width:100%"></div>
        </div>
      </div>` : ''}
    </div>`;

  const sequences = `
    <div class="vq-chart-card${domainDbs.length ? ' vq-wide-narrow' : ''}" id="stats-length-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Viral Length Distribution</div>
          <div class="vq-chart-card__sub" id="stats-length-sub">viral sequences</div>
        </div>
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-length-svg" style="width:100%"></div>
      </div>
    </div>
    ${domainDbs.length ? `
    <div class="vq-chart-card vq-span-2" id="stats-domains-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Top Domains</div>
          <div class="vq-chart-card__sub" id="stats-domains-sub">sequences carrying each model</div>
        </div>
        <div class="vq-toggle" id="stats-domains-toggle" role="tablist" aria-label="HMM database">
          ${domainDbs.map((db, i) => `
            <button class="vq-toggle__btn${i === 0 ? ' active' : ''}" type="button"
                    data-domain-db="${esc(db)}">${esc(db)}</button>`).join('')}
        </div>
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-domains-svg" style="width:100%"></div>
      </div>
    </div>` : ''}`;

  el.innerHTML = `
    <div class="vq-stats-page">
      <div class="vq-stats-row vq-stats-row--top">
        ${_chip('Viral Sequences', _fmtNum(viral.length), 'success', 'confirmed contigs')}
        ${_chip('Species', _fmtNum(viruses.length), 'accent', 'distinct best-hit species')}
        ${_chip('Families', _fmtNum(families.size), '', 'with resolved taxonomy')}
        ${_chip('Known Viruses', _fmtNum(known), 'success',
                `≥${_KNOWN_ID}% identity · ≥${_KNOWN_COV}% coverage`)}
      </div>

      ${_group('Viruses', 'one row per species · best BLASTx hit', table)}
      ${_group('Taxonomy', 'what kinds of viruses were found', taxonomy)}
      ${_group('Similarity to known viruses', 'how close each sequence is to its best reference', similarity)}
      ${_group('Sequences & proteins', 'contig sizes and protein domains',
               sequences, 1 + (domainDbs.length ? 2 : 0))}
    </div>`;

  // ── Detected viruses table (search + status filter + sorting) ──
  const tableState = { q: '', status: 'all', sort: 'contigs', dir: -1 };
  const drawTable = () => _renderVirusTable(viruses, tableState);
  document.getElementById('stats-viruses-search')?.addEventListener('input', e => {
    tableState.q = e.target.value.trim().toLowerCase(); drawTable();
  });
  const statusBtns = document.querySelectorAll('#stats-viruses-status [data-status]');
  statusBtns.forEach(btn => btn.addEventListener('click', () => {
    statusBtns.forEach(b => b.classList.toggle('active', b === btn));
    tableState.status = btn.dataset.status; drawTable();
  }));
  document.getElementById('stats-viruses-table')?.addEventListener('click', e => {
    const th = e.target.closest('th[data-sort]');
    if (th) {
      const key = th.dataset.sort;
      tableState.dir = tableState.sort === key ? -tableState.dir : (key === 'species' || key === 'family' ? 1 : -1);
      tableState.sort = key;
      drawTable();
      return;
    }
    if (e.target.closest('a')) return;                      // accession link
    const row = e.target.closest('tr[data-seq]');
    if (row) VQ.jumpToViewer(row.dataset.seq);
  });
  drawTable();

  // ── Taxonomy donut ──
  let currentRank = 'family';
  const renderTax = () => _renderTaxDonut(seqs, currentRank);
  document.querySelectorAll('#stats-tax-card [data-rank]').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('#stats-tax-card [data-rank]').forEach(b =>
        b.classList.toggle('active', b === btn));
      currentRank = btn.dataset.rank;
      renderTax();
    });
  });

  // ── Identity × coverage scatter ──
  let scatterSrc = 'blastx';
  const scatterBtns = document.querySelectorAll('#stats-scatter-toggle [data-scatter]');
  scatterBtns.forEach(btn => btn.addEventListener('click', () => {
    scatterBtns.forEach(b => b.classList.toggle('active', b === btn));
    scatterSrc = btn.dataset.scatter;
    _renderSimilarityScatter(viral, scatterSrc);
  }));

  // ── BLAST histograms: each card has its own Identity / Coverage toggle ──
  const blastMetric = { x: 'identity', n: 'identity' };
  const blastCards = {
    x: { toggle: 'blast-toggle-x', title: 'stats-blastx-title', label: 'BLASTx',
         draw: m => _renderIdentityHistogram(seqs, m) },
    n: { toggle: 'blast-toggle-n', title: 'stats-blastn-title', label: 'BLASTn',
         draw: m => _renderBLAstnHistogram(seqs, m) },
  };
  const drawBlast = () => {
    blastCards.x.draw(blastMetric.x);
    blastCards.n.draw(blastMetric.n);
  };
  Object.entries(blastCards).forEach(([key, card]) => {
    const btns = document.querySelectorAll(`#${card.toggle} [data-blast-metric]`);
    btns.forEach(btn => btn.addEventListener('click', () => {
      blastMetric[key] = btn.dataset.blastMetric;
      btns.forEach(b => b.classList.toggle('active', b === btn));
      const title = document.getElementById(card.title);
      if (title) title.textContent =
        `${card.label} ${blastMetric[key] === 'coverage' ? 'Coverage' : 'Identity'}`;
      card.draw(blastMetric[key]);
    }));
  });

  // ── Top domains ──
  let domainDb = domainDbs[0];
  const domainBtns = document.querySelectorAll('#stats-domains-toggle [data-domain-db]');
  domainBtns.forEach(btn => btn.addEventListener('click', () => {
    domainBtns.forEach(b => b.classList.toggle('active', b === btn));
    domainDb = btn.dataset.domainDb;
    _renderTopDomains(seqs, domainDb);
  }));

  // Fixed-height charts first; the scatter and length histogram then fill
  // the row height their neighbours set.
  const draw = () => {
    renderTax();
    _renderGenomeTypes(viral);
    _renderNRClassification(seqs);
    drawBlast();
    if (domainDb) _renderTopDomains(seqs, domainDb);
    _renderSimilarityScatter(viral, scatterSrc);
    _renderLengthHistogram(seqs);
  };
  draw();
  _redrawOnResize(el, draw);
}

// ────────────────────────────────────────────────────────────────────────
//  Virome helpers — per-species rows, status, genome types, scatter
// ────────────────────────────────────────────────────────────────────────

// "Known" mirrors the heuristic/LLM rule: >=90% identity AND >=70% coverage.
const _KNOWN_ID  = 90;
const _KNOWN_COV = 70;
const _RELATED_ID = 70;
const _STATUS_LABEL = { known: 'Known', related: 'Related', divergent: 'Divergent' };
const _STATUS_COLOR = {
  known:     'var(--vq-success)',
  related:   'var(--vq-accent)',
  divergent: 'var(--vq-warning)',
};

function _hitStatus(identity, coverage) {
  if (identity >= _KNOWN_ID && coverage >= _KNOWN_COV) return 'known';
  if (identity >= _RELATED_ID) return 'related';
  return 'divergent';
}

/* Best BLASTx hit of a sequence: NR first, RefSeq as fallback. */
function _bestBlastx(s) {
  const hits = (s.blastx_nr_hits && s.blastx_nr_hits.length) ? s.blastx_nr_hits : (s.blastx_hits || []);
  return hits.length ? hits.reduce((a, h) => (h.bit_score > a.bit_score ? h : a)) : null;
}

function _bestBlastn(s) {
  const hits = s.blastn_hits || [];
  return hits.length ? hits.reduce((a, h) => ((h.bit_score ?? 0) > (a.bit_score ?? 0) ? h : a)) : null;
}

/* Species of a hit: parsed `species`, else the bracketed organism in the title. */
function _hitSpecies(h) {
  if (h.species) return h.species;
  const m = /\[([^\]]+)\]\s*$/.exec(h.subject_title || '');
  return m ? m[1] : (h.subject_title || 'Unknown');
}

function _mostCommon(values) {
  const c = new Map();
  values.filter(Boolean).forEach(v => c.set(v, (c.get(v) || 0) + 1));
  let best = null, n = 0;
  c.forEach((k, v) => { if (k > n) { best = v; n = k; } });
  return best;
}

const _PROTEIN_ACC = /^[A-Za-z0-9_]+(\.[0-9]+)?$/;
function _proteinUrl(acc) {
  return acc && _PROTEIN_ACC.test(acc)
    ? `https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(acc)}` : null;
}

/* One row per best-hit species, summarising its contigs. */
function _virusRows(viral) {
  const bySpecies = new Map();
  viral.forEach(s => {
    const hit = _bestBlastx(s);
    if (!hit) return;
    const sp = _hitSpecies(hit);
    if (!bySpecies.has(sp)) bySpecies.set(sp, []);
    bySpecies.get(sp).push({ seq: s, hit, bn: _bestBlastn(s) });
  });
  return [...bySpecies.entries()].map(([species, items]) => {
    const best = items.reduce((a, b) => (b.hit.pct_identity > a.hit.pct_identity ? b : a));
    const bnBest = items.map(i => i.bn).filter(Boolean)
      .reduce((a, b) => (!a || b.pident > a.pident ? b : a), null);
    const identity = best.hit.pct_identity;
    const coverage = best.hit.query_coverage ?? 0;
    return {
      species,
      family:   _mostCommon(items.map(i => i.seq.taxonomy?.family)) || '',
      genome:   _mostCommon(items.map(i => i.seq.taxonomy?.genome)) || '',
      contigs:  items.length,
      length:   items.reduce((a, i) => a + (i.seq.length || 0), 0),
      identity, coverage,
      blastn:   bnBest ? bnBest.pident : null,
      accession: best.hit.subject_id,
      seqId:    best.seq.id,
      status:   _hitStatus(identity, coverage),
    };
  });
}

function _pctCell(v) {
  if (v == null) return '<span class="vq-vt-na">—</span>';
  const w = Math.max(0, Math.min(100, v));
  return `
    <span class="vq-vt-pct">
      <span class="vq-vt-pct__bar"><span style="width:${w}%"></span></span>
      <span class="vq-vt-pct__val">${v.toFixed(1)}</span>
    </span>`;
}

function _renderVirusTable(rows, st) {
  const wrap = document.getElementById('stats-viruses-table');
  if (!wrap) return;
  const esc = VQ.esc;

  const counts = { all: rows.length, known: 0, related: 0, divergent: 0 };
  rows.forEach(r => { counts[r.status]++; });
  document.querySelectorAll('#stats-viruses-status [data-status]').forEach(b => {
    const k = b.dataset.status;
    b.textContent = `${k === 'all' ? 'All' : _STATUS_LABEL[k]} · ${counts[k]}`;
  });

  let view = rows.filter(r =>
    (st.status === 'all' || r.status === st.status) &&
    (!st.q || `${r.species} ${r.family} ${r.genome}`.toLowerCase().includes(st.q)));
  const key = st.sort;
  view = view.sort((a, b) => {
    const x = a[key], y = b[key];
    if (typeof x === 'string' || typeof y === 'string')
      return st.dir * String(x || '').localeCompare(String(y || ''));
    return st.dir * ((x ?? -1) - (y ?? -1)) || a.species.localeCompare(b.species);
  });

  const sub = document.getElementById('stats-viruses-sub');
  if (sub) sub.textContent =
    `${rows.length} species · status from the best BLASTx hit (known ≥${_KNOWN_ID}% id & ≥${_KNOWN_COV}% cov, related ≥${_RELATED_ID}% id) · click a row to open its best contig`;

  if (!rows.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px;font-size:12px">No BLASTx hits.</div>';
    return;
  }

  const th = (k, label, cls = '') => {
    const on = st.sort === k;
    return `<th data-sort="${k}" class="${cls}${on ? ' vq-vt-sorted' : ''}" scope="col">
      ${label}<span class="vq-vt-arrow">${on ? (st.dir > 0 ? '▲' : '▼') : ''}</span></th>`;
  };
  const hasBn = rows.some(r => r.blastn != null);

  wrap.innerHTML = view.length ? `
    <table class="vq-table vq-vt">
      <thead><tr>
        ${th('species',  'Species')}
        ${th('status',   'Status')}
        ${th('family',   'Family')}
        ${th('genome',   'Genome')}
        ${th('contigs',  'Contigs',  'vq-vt-num')}
        ${th('length',   'Total bp', 'vq-vt-num')}
        ${th('identity', 'BLASTx id %')}
        ${th('coverage', 'Coverage %')}
        ${hasBn ? th('blastn', 'BLASTn id %') : ''}
        <th scope="col">Best hit</th>
      </tr></thead>
      <tbody>
        ${view.map(r => {
          const url = _proteinUrl(r.accession);
          return `
          <tr data-seq="${esc(r.seqId)}" title="Open ${esc(r.seqId)} in the Sequence Viewer">
            <td class="vq-vt-species">${esc(r.species)}</td>
            <td><span class="vq-vt-status vq-vt-status--${r.status}">${_STATUS_LABEL[r.status]}</span></td>
            <td>${r.family ? esc(r.family) : '<span class="vq-vt-na">—</span>'}</td>
            <td>${r.genome ? esc(r.genome) : '<span class="vq-vt-na">—</span>'}</td>
            <td class="vq-vt-num">${r.contigs.toLocaleString()}</td>
            <td class="vq-vt-num">${r.length.toLocaleString()}</td>
            <td>${_pctCell(r.identity)}</td>
            <td>${_pctCell(r.coverage)}</td>
            ${hasBn ? `<td>${_pctCell(r.blastn)}</td>` : ''}
            <td class="vq-td--mono">${url
              ? `<a class="vq-acc-link" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(r.accession)}</a>`
              : esc(r.accession || '—')}</td>
          </tr>`;
        }).join('')}
      </tbody>
    </table>` : '<div class="vq-empty" style="padding:24px;font-size:12px">No species match the filter.</div>';
}

const _GENOME_COLORS = [
  'var(--vq-accent)', 'var(--vq-primary-light)', 'var(--vq-success)',
  'var(--vq-accent-dark)', 'var(--vq-warning)', 'var(--vq-dom-5)', 'var(--vq-dom-6)',
];

function _renderGenomeTypes(viral) {
  const wrap = document.getElementById('stats-genome-svg');
  if (!wrap) return;
  VQ.tooltipHide();
  wrap.innerHTML = '';

  const counts = new Map();
  viral.forEach(s => {
    const g = s.taxonomy?.genome || 'Unresolved';
    counts.set(g, (counts.get(g) || 0) + 1);
  });
  // Resolved types by size, "Unresolved" always last.
  const data = [...counts.entries()]
    .sort((a, b) => (a[0] === 'Unresolved') - (b[0] === 'Unresolved') || b[1] - a[1])
    .map(([name, value], i) => ({
      name, value,
      color: name === 'Unresolved' ? 'var(--vq-border-dark)' : _GENOME_COLORS[i % _GENOME_COLORS.length],
    }));

  const sub = document.getElementById('stats-genome-sub');
  const resolved = viral.length - (counts.get('Unresolved') || 0);
  if (sub) sub.textContent = `${resolved} of ${viral.length} sequences with a resolved genome type`;

  if (!data.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px;font-size:12px">No sequences.</div>';
    return;
  }

  const W = Math.max(wrap.clientWidth || 0, 220);
  const PAD_L = 84, PAD_R = 44, ROW_H = 30, BAR_H = 14;
  const H = data.length * ROW_H + 8;
  const drawW = W - PAD_L - PAD_R;
  const xScale = d3.scaleLinear([0, d3.max(data, d => d.value) || 1], [0, drawW]);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  data.forEach((d, i) => {
    const y = 6 + i * ROW_H;
    svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + BAR_H / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', 11.5).attr('font-weight', 500)
      .attr('fill', d.name === 'Unresolved' ? 'var(--vq-text-3)' : 'var(--vq-text-2)')
      .text(d.name);
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y).attr('width', drawW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2).attr('fill', 'var(--vq-bg)');
    const bw = Math.max(xScale(d.value), 4);
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y).attr('width', bw).attr('height', BAR_H)
      .attr('rx', BAR_H / 2).attr('fill', d.color)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(d.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${d.value.toLocaleString()}</span>
          <span class="vq-tooltip__key">Share</span><span>${_pct(d.value, viral.length)}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
    svg.append('text')
      .attr('x', PAD_L + bw + 8).attr('y', y + BAR_H / 2 + 4)
      .attr('font-size', 11.5).attr('font-weight', 600)
      .attr('fill', 'var(--vq-text)').attr('font-variant-numeric', 'tabular-nums')
      .text(d.value.toLocaleString());
  });

  wrap.appendChild(svg.node());
}

function _renderSimilarityScatter(viral, source) {
  const wrap = document.getElementById('stats-scatter-svg');
  if (!wrap) return;
  VQ.tooltipHide();
  wrap.innerHTML = '';
  const esc = VQ.esc;
  const isX = source !== 'blastn';

  const pts = viral.map(s => {
    const h = isX ? _bestBlastx(s) : _bestBlastn(s);
    if (!h) return null;
    const id  = isX ? h.pct_identity : h.pident;
    const cov = isX ? h.query_coverage : h.qcovhsp;
    if (id == null || cov == null) return null;
    return {
      s, id, cov, status: _hitStatus(id, cov),
      species: isX ? _hitSpecies(h) : (h.stitle || ''),
    };
  }).filter(Boolean);

  const sub = document.getElementById('stats-scatter-sub');
  if (sub) sub.textContent =
    `${pts.length} sequences · best ${isX ? 'BLASTx' : 'BLASTn'} hit · shaded area = known (≥${_KNOWN_ID}% id, ≥${_KNOWN_COV}% cov) · click a point to open it`;

  if (!pts.length) {
    wrap.innerHTML = `<div class="vq-empty" style="padding:24px;font-size:12px">No ${isX ? 'BLASTx' : 'BLASTn'} hits.</div>`;
    return;
  }

  const W = Math.max(wrap.clientWidth || 0, 300);
  const H = _fillHeight(wrap, 260);
  const PAD_L = 40, PAD_R = 12, PAD_T = 26, PAD_B = 34;
  const drawW = W - PAD_L - PAD_R, drawH = H - PAD_T - PAD_B;

  // Axes start at the lowest observed value (rounded down to 5) and end at 100%.
  const floor5 = v => Math.max(0, Math.floor(v / 5) * 5);
  const xLo = Math.min(floor5(d3.min(pts, p => p.id)),  95);
  const yLo = Math.min(floor5(d3.min(pts, p => p.cov)), 95);
  const x = d3.scaleLinear([xLo, 100], [0, drawW]);
  const y = d3.scaleLinear([yLo, 100], [drawH, 0]);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');
  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // "Known" region
  const kx = x(Math.max(_KNOWN_ID, xLo)), ky = y(Math.max(_KNOWN_COV, yLo));
  g.append('rect').attr('x', kx).attr('y', 0).attr('width', drawW - kx).attr('height', ky)
    .attr('fill', 'var(--vq-success)').attr('opacity', 0.07);

  // Grid + ticks
  x.ticks(6).forEach(t => {
    g.append('line').attr('x1', x(t)).attr('x2', x(t)).attr('y1', 0).attr('y2', drawH)
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text').attr('x', x(t)).attr('y', drawH + 13)
      .attr('text-anchor', _tickAnchor(t, xLo, 100)).attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(t + '%');
  });
  y.ticks(5).forEach(t => {
    g.append('line').attr('x1', 0).attr('x2', drawW).attr('y1', y(t)).attr('y2', y(t))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text').attr('x', -5).attr('y', y(t) + 3)
      .attr('text-anchor', 'end').attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(t + '%');
  });
  g.append('line').attr('x1', 0).attr('x2', drawW).attr('y1', drawH).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)');
  g.append('line').attr('x1', 0).attr('x2', 0).attr('y1', 0).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)');
  svg.append('text').attr('x', PAD_L + drawW / 2).attr('y', H - 3)
    .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('fill', 'var(--vq-c-tick)')
    .text('identity (%)');
  svg.append('text')
    .attr('transform', `translate(10,${PAD_T + drawH / 2}) rotate(-90)`)
    .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('fill', 'var(--vq-c-tick)')
    .text('query coverage (%)');

  // Points — divergent drawn last so they stay visible over the dense known cluster
  const order = { known: 0, related: 1, divergent: 2 };
  pts.sort((a, b) => order[a.status] - order[b.status]).forEach(p => {
    g.append('circle')
      .attr('cx', x(p.id)).attr('cy', y(p.cov)).attr('r', 4.2)
      .attr('fill', _STATUS_COLOR[p.status]).attr('fill-opacity', 0.72)
      .attr('stroke', 'var(--vq-surface)').attr('stroke-width', 1)
      .style('cursor', 'pointer')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(p.s.id)}</div>
        <div style="margin-bottom:4px;font-style:italic">${esc(p.species)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Identity</span><span>${p.id.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Coverage</span><span>${Number(p.cov).toFixed(1)}%</span>
          <span class="vq-tooltip__key">Family</span><span>${esc(p.s.taxonomy?.family || '—')}</span>
          <span class="vq-tooltip__key">Status</span><span>${_STATUS_LABEL[p.status]}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide)
      .on('click', () => { VQ.tooltipHide(); VQ.jumpToViewer(p.s.id); });
  });

  // Legend
  const counts = { known: 0, related: 0, divergent: 0 };
  pts.forEach(p => counts[p.status]++);
  const lg = svg.append('g').attr('transform', `translate(${PAD_L + 4},10)`);
  let lx = 0;
  ['known', 'related', 'divergent'].forEach(k => {
    const item = lg.append('g').attr('transform', `translate(${lx},0)`);
    item.append('circle').attr('r', 4).attr('cy', 0).attr('fill', _STATUS_COLOR[k]);
    const t = item.append('text').attr('x', 8).attr('y', 3.5).attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)').text(`${_STATUS_LABEL[k]} ${counts[k]}`);
    lx += 8 + (t.node().getComputedTextLength?.() || 60) + 16;
  });

  wrap.appendChild(svg.node());
}

// ────────────────────────────────────────────────────────────────────────
//  Chip / detail helpers
// ────────────────────────────────────────────────────────────────────────

function _chip(label, value, mod = '', sub = '') {
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${VQ.esc(label)}</div>
      <div class="vq-stat-chip__value" title="${VQ.esc(value)}">${value}</div>
      ${sub ? `<div class="vq-stat-chip__sub" title="${VQ.esc(sub)}">${sub}</div>` : ''}
    </div>`;
}

function _blobChip(label, value, sub = '') {
  /* Flower-blob KPI chip: 6 overlapping semi-transparent circles form the
     petal ring; a white disc sits on top carrying the number and label. */
  const petals = [
    [104, 85, 'var(--vq-accent)'],
    [95,  68, 'var(--vq-accent-dark)'],
    [75,  68, 'var(--vq-primary-light)'],
    [66,  85, 'var(--vq-accent)'],
    [75, 102, 'var(--vq-accent-dark)'],
    [95, 102, 'var(--vq-primary-light)'],
  ];
  const rings = petals.map(([px, py, col]) =>
    `<circle cx="${px}" cy="${py}" r="54" fill="${col}" opacity="0.28"/>`
  ).join('');
  const subText = sub
    ? `<text x="85" y="108" text-anchor="middle" font-size="8.5"
             font-family="var(--vq-font)" fill="var(--vq-text-3)">${VQ.esc(sub)}</text>`
    : '';
  return `
    <div class="vq-stat-chip vq-stat-chip--blob">
      <svg viewBox="0 0 170 170" width="170" height="170"
           style="overflow:visible;display:block;margin:auto">
        ${rings}
        <circle cx="85" cy="85" r="47" fill="var(--vq-surface)"/>
        <text x="85" y="79" text-anchor="middle"
              font-size="26" font-weight="700" font-family="var(--vq-font)"
              fill="var(--vq-primary)">${VQ.esc(value)}</text>
        <text x="85" y="95" text-anchor="middle"
              font-size="8.5" font-weight="600" letter-spacing=".08em"
              font-family="var(--vq-font)" fill="var(--vq-text-3)"
        >${VQ.esc(label.toUpperCase())}</text>
        ${subText}
      </svg>
    </div>`;
}

function _detailCard(title, bodyHtml) {
  return `
    <div class="vq-chart-card" style="min-height:auto;padding:14px 16px">
      <div class="vq-chart-card__title" style="margin-bottom:10px">${VQ.esc(title)}</div>
      ${bodyHtml}
    </div>`;
}

function _miniRow(label, value) {
  return `
    <div style="display:flex;justify-content:space-between;align-items:baseline;
                padding:6px 0;border-bottom:1px dashed var(--vq-border);font-size:12px;
                gap:12px;min-width:0">
      <span style="color:var(--vq-text-3);flex-shrink:0">${VQ.esc(label)}</span>
      <span style="color:var(--vq-text);font-weight:600;font-variant-numeric:tabular-nums;
                   white-space:nowrap;overflow:hidden;text-overflow:ellipsis">${value}</span>
    </div>`;
}

// ────────────────────────────────────────────────────────────────────────
//  Pipeline workflow card — flowchart of every step (run / skipped),
//  success path between completed steps, per-step time and details.
// ────────────────────────────────────────────────────────────────────────

const _WF_SHORT = {
  parse: 'Parse FASTA', orfs: 'ORFs', refseq: 'RefSeq', hmm: 'HMM filter',
  nr: 'Diamond NR', blastn: 'BLASTn', pfam: 'Pfam', taxonomy: 'Taxonomy',
  clusters: 'Clusters', salmon: 'Salmon', seq_quality: 'Seq quality',
  coverage: 'Coverage', heuristic: 'Heuristic', llm: 'LLM', export: 'Export',
};
const _WF_STATUS = {
  done:    { label: 'Completed', color: 'var(--vq-success)', glyph: '✓' },
  partial: { label: 'Partial',   color: 'var(--vq-warning)', glyph: '!' },
  error:   { label: 'Error',     color: 'var(--vq-danger)',  glyph: '✕' },
  skipped: { label: 'Skipped',   color: 'var(--vq-text-3)',  glyph: ''  },
};

function _fmtDur(sec) {
  if (sec == null || isNaN(sec)) return '—';
  if (sec < 1)    return (sec * 1000).toFixed(0) + ' ms';
  if (sec < 60)   return sec.toFixed(1) + ' s';
  const m = Math.floor(sec / 60), s = Math.round(sec % 60);
  if (sec < 3600) return `${m}m ${String(s).padStart(2, '0')}s`;
  return `${Math.floor(m / 60)}h ${String(m % 60).padStart(2, '0')}m`;
}

function _wfOverall(steps) {
  if (steps.some(s => s.status === 'error'))   return 'error';
  if (steps.some(s => s.status === 'partial')) return 'partial';
  return 'done';
}

function _wfPill(status, text) {
  return `<span class="vq-wf-pill vq-wf-pill--${status}">${VQ.esc(text ?? _WF_STATUS[status].label)}</span>`;
}

function _workflowCardHtml(wf, version) {
  const esc     = VQ.esc;
  const steps   = wf.steps || [];
  const ran     = steps.filter(s => s.status !== 'skipped');
  const skipped = steps.length - ran.length;
  const overall = _wfOverall(steps);
  const issues  = steps.filter(s => s.status === 'error' || s.status === 'partial');
  const o       = wf.options || {};

  const opt = (label, value, on = true) =>
    `<span class="vq-wf-opt${on ? '' : ' vq-wf-opt--off'}">
       <span class="vq-wf-opt__k">${esc(label)}</span>${esc(value)}</span>`;
  const reads = (o.reads || []).length
    ? `${o.reads.length} file${o.reads.length > 1 ? 's' : ''}${o.read_type ? ' · ' + o.read_type : ''}`
    : 'off';

  const overallText = overall === 'done' ? 'Completed'
                    : overall === 'partial' ? 'Completed with warnings'
                    : 'Completed with errors';

  return `
    <div class="vq-chart-card vq-wf-card" id="stats-workflow-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Pipeline Workflow</div>
          <div class="vq-chart-card__sub">
            started ${esc(wf.started_at ? new Date(wf.started_at).toLocaleString() : '—')}
            &nbsp;·&nbsp; ${ran.length} step${ran.length === 1 ? '' : 's'} run
            ${skipped ? `&nbsp;·&nbsp; ${skipped} skipped` : ''}
            ${version ? `&nbsp;·&nbsp; ViralQuest v${esc(version)}` : ''}
            &nbsp;·&nbsp; click a step for details
          </div>
        </div>
        <div class="vq-wf-head-right">
          ${_wfPill(overall, overallText)}
          <div class="vq-chart-card__big" style="margin-top:0">${_fmtDur(wf.total_seconds)}</div>
        </div>
      </div>

      <div class="vq-wf-opts">
        ${o.input ? opt('Input', o.input) : ''}
        ${o.threads != null ? opt('Threads', o.threads) : ''}
        ${opt('CAP3', o.cap3 ? 'on' : 'off', !!o.cap3)}
        ${opt('NR', o.nr_db || 'off', !!o.nr_db)}
        ${opt('BLASTn', o.blastn || 'off', !!o.blastn)}
        ${opt('Reads', reads, !!(o.reads || []).length)}
        ${o.transcriptome ? opt('Transcriptome', o.transcriptome) : ''}
        ${opt('LLM', o.llm || 'off', !!o.llm)}
        ${o.force ? opt('Force', 'on') : ''}
      </div>

      ${issues.length ? `
      <div class="vq-wf-issues">
        ${issues.map(s => `
          <button type="button" class="vq-wf-issue vq-wf-issue--${s.status}" data-wf-key="${esc(s.key)}">
            <strong>${esc(_WF_SHORT[s.key] || s.label)}</strong>
            <span>${esc(s.message || _WF_STATUS[s.status].label)}</span>
          </button>`).join('')}
      </div>` : ''}

      <div class="vq-wf-flow" id="stats-workflow-flow"></div>
      <div class="vq-wf-timebar" id="stats-workflow-timebar"></div>
      <div class="vq-wf-detail" id="stats-workflow-detail"></div>
    </div>`;
}

function _renderWorkflow(wf) {
  const wrap = document.getElementById('stats-workflow-flow');
  if (!wrap) return;
  const esc   = VQ.esc;
  const steps = wf.steps || [];
  const total = wf.total_seconds || steps.reduce((a, s) => a + (s.seconds || 0), 0);

  // ── Flowchart (SVG) ──────────────────────────────────────────────────
  const SLOT   = Math.max((wrap.clientWidth || 0) / steps.length, 86);
  const W      = SLOT * steps.length;
  const Y_MAIN = 26, Y_SKIP = 100, H = 146;
  const xOf    = i => SLOT * i + SLOT / 2;

  const svg = d3.create('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');

  // Success path: connects consecutive executed steps along the main lane;
  // each segment takes the colour of the step it leads into.
  const ran = steps.map((s, i) => ({ s, i })).filter(d => d.s.status !== 'skipped');
  ran.forEach((d, k) => {
    if (!k) return;
    const prev = ran[k - 1];
    svg.append('line')
      .attr('x1', xOf(prev.i)).attr('y1', Y_MAIN)
      .attr('x2', xOf(d.i)).attr('y2', Y_MAIN)
      .attr('stroke', _WF_STATUS[d.s.status].color)
      .attr('stroke-width', 3).attr('stroke-linecap', 'round')
      .attr('opacity', d.s.status === 'done' ? 0.55 : 0.85);
  });

  // Skipped steps hang below the path on a dashed branch (option not chosen).
  steps.forEach((s, i) => {
    if (s.status !== 'skipped') return;
    svg.append('path')
      .attr('d', `M${xOf(i)},${Y_MAIN + 4} C${xOf(i)},${Y_MAIN + 40} ${xOf(i)},${Y_SKIP - 40} ${xOf(i)},${Y_SKIP - 9}`)
      .attr('fill', 'none')
      .attr('stroke', 'var(--vq-border-dark)')
      .attr('stroke-width', 1.2).attr('stroke-dasharray', '3,3');
  });

  const nodes = svg.selectAll('g.vq-wf-node').data(steps).join('g')
    .attr('class', s => `vq-wf-node vq-wf-node--${s.status}`)
    .attr('data-wf-key', s => s.key)
    .attr('transform', (s, i) => `translate(${xOf(i)},${s.status === 'skipped' ? Y_SKIP : Y_MAIN})`)
    .style('cursor', 'pointer');

  nodes.append('circle').attr('class', 'vq-wf-node__halo')
    .attr('r', s => s.status === 'skipped' ? 13 : 17)
    .attr('fill', 'none').attr('stroke', 'var(--vq-accent)').attr('stroke-width', 2)
    .attr('opacity', 0);

  nodes.append('circle')
    .attr('r', s => s.status === 'skipped' ? 8 : 11)
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-surface)' : _WF_STATUS[s.status].color)
    .attr('stroke', s => s.status === 'skipped' ? 'var(--vq-border-dark)' : 'var(--vq-surface)')
    .attr('stroke-width', s => s.status === 'skipped' ? 1.4 : 2)
    .attr('stroke-dasharray', s => s.status === 'skipped' ? '2,2' : null);

  nodes.filter(s => s.status !== 'skipped').append('text')
    .attr('text-anchor', 'middle').attr('dy', '0.35em')
    .attr('font-size', 11).attr('font-weight', 700).attr('fill', '#fff')
    .attr('pointer-events', 'none')
    .text(s => _WF_STATUS[s.status].glyph);

  nodes.append('text')
    .attr('text-anchor', 'middle')
    .attr('y', s => s.status === 'skipped' ? 22 : 30)
    .attr('font-size', 11)
    .attr('font-weight', s => s.status === 'skipped' ? 400 : 600)
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-text-3)' : 'var(--vq-text)')
    .text(s => _WF_SHORT[s.key] || s.key);

  nodes.append('text')
    .attr('text-anchor', 'middle')
    .attr('y', s => s.status === 'skipped' ? 34 : 44)
    .attr('font-size', 10)
    .attr('font-family', 'var(--vq-font-mono)')
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-text-3)' : _WF_STATUS[s.status].color)
    .text(s => s.status === 'skipped' ? 'skipped' : _fmtDur(s.seconds));

  nodes
    .on('mousemove', (evt, s) => VQ.tooltipShow(`
      <div class="vq-tooltip__title">${esc(s.label)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Status</span><span>${_WF_STATUS[s.status].label}</span>
        ${s.seconds != null ? `<span class="vq-tooltip__key">Time</span><span>${_fmtDur(s.seconds)}</span>` : ''}
        ${s.message ? `<span class="vq-tooltip__key">Note</span><span>${esc(s.message)}</span>` : ''}
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide)
    .on('click', (evt, s) => select(s.key));

  wrap.innerHTML = '';
  wrap.appendChild(svg.node());

  // ── Time breakdown bar ───────────────────────────────────────────────
  const bar = document.getElementById('stats-workflow-timebar');
  const timed = steps.filter(s => s.seconds != null && s.seconds > 0);
  if (bar && timed.length && total > 0) {
    const longest = timed.reduce((a, s) => (s.seconds > a.seconds ? s : a));
    bar.innerHTML = `
      <div class="vq-wf-timebar__track">
        ${timed.map((s, i) => `
          <div class="vq-wf-timebar__seg" data-wf-key="${esc(s.key)}"
               style="flex:${s.seconds} 1 0;opacity:${i % 2 ? 0.55 : 0.9}"></div>`).join('')}
      </div>
      <div class="vq-wf-timebar__legend">
        <span>time per step</span>
        <span>longest: <strong>${esc(_WF_SHORT[longest.key] || longest.label)}</strong>
          · ${_fmtDur(longest.seconds)} (${_pct(longest.seconds, total)})</span>
      </div>`;
    bar.querySelectorAll('.vq-wf-timebar__seg').forEach(seg => {
      const s = steps.find(x => x.key === seg.dataset.wfKey);
      seg.addEventListener('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(s.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Time</span><span>${_fmtDur(s.seconds)}</span>
          <span class="vq-tooltip__key">of total</span><span>${_pct(s.seconds, total)}</span>
        </div>`, evt));
      seg.addEventListener('mouseleave', VQ.tooltipHide);
      seg.addEventListener('click', () => select(s.key));
    });
  }

  // ── Step details ─────────────────────────────────────────────────────
  const detail = document.getElementById('stats-workflow-detail');
  function select(key) {
    const s = steps.find(x => x.key === key);
    if (!s || !detail) return;
    svg.selectAll('g.vq-wf-node').each(function (d) {
      d3.select(this).select('.vq-wf-node__halo').attr('opacity', d.key === key ? 0.9 : 0);
    });
    const rows = Object.entries(s.details || {});
    detail.innerHTML = `
      <div class="vq-wf-detail__head">
        <div class="vq-wf-detail__title">${esc(s.label)}</div>
        ${_wfPill(s.status)}
        ${s.seconds != null
          ? `<span class="vq-wf-detail__time">${_fmtDur(s.seconds)} · ${_pct(s.seconds, total)} of total</span>`
          : ''}
      </div>
      ${s.message ? `<div class="vq-wf-detail__msg vq-wf-detail__msg--${s.status}">${esc(s.message)}</div>` : ''}
      ${rows.length ? `
        <div class="vq-wf-detail__grid">
          ${rows.map(([k, v]) => _miniRow(k,
              esc(typeof v === 'number' ? v.toLocaleString() : String(v ?? '—')))).join('')}
        </div>` : ''}`;
  }

  document.querySelectorAll('#stats-workflow-card .vq-wf-issue').forEach(btn =>
    btn.addEventListener('click', () => select(btn.dataset.wfKey)));

  // Open on the first problem, otherwise on the slowest step.
  const first = steps.find(s => s.status === 'error')
             || steps.find(s => s.status === 'partial')
             || (timed.length ? timed.reduce((a, s) => (s.seconds > a.seconds ? s : a)) : steps[0]);
  select(first.key);
}

function _pct(a, b) {
  if (!b || isNaN(a) || isNaN(b)) return '—';
  return ((a / b) * 100).toFixed(1) + '%';
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 1 — Half-donut: family/phylum distribution
// ────────────────────────────────────────────────────────────────────────

const _TAX_DONUT_PALETTE = [
  'var(--vq-tax-flaviviridae)','var(--vq-tax-parvoviridae)',
  'var(--vq-tax-phenuiviridae)','var(--vq-tax-nodaviridae)',
  'var(--vq-tax-rhabdoviridae)','var(--vq-tax-tombusviridae)',
  '#0891b2','#7c3aed','#b45309','#be185d','var(--vq-tax-other)',
];

function _familyColor(family, i) {
  const cssVar = '--vq-tax-' + (family || '').toLowerCase();
  const computed = getComputedStyle(document.documentElement).getPropertyValue(cssVar).trim();
  return computed || _TAX_DONUT_PALETTE[i % _TAX_DONUT_PALETTE.length];
}

function _renderTaxDonut(sequences, rank) {
  const wrap = document.getElementById('stats-tax-svg');
  const lgd  = document.getElementById('stats-tax-legend');
  const sub  = document.getElementById('stats-tax-sub');
  if (!wrap || !lgd) return;
  wrap.innerHTML = '';
  lgd.innerHTML  = '';

  const viral = sequences.filter(s => s.is_viral);
  const counts = {};
  viral.forEach(s => {
    const key = s.taxonomy?.[rank] || 'Unclassified';
    counts[key] = (counts[key] || 0) + 1;
  });

  const entries = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  // Collapse the long tail into "Other" so the donut stays readable.
  const MAX_SLICES = 7;
  let display = entries;
  if (entries.length > MAX_SLICES) {
    const top  = entries.slice(0, MAX_SLICES - 1);
    const rest = entries.slice(MAX_SLICES - 1);
    const otherCount = rest.reduce((a, [, n]) => a + n, 0);
    display = [...top, ['Other (' + rest.length + ')', otherCount]];
  }

  const total = entries.reduce((a, [, n]) => a + n, 0);
  if (sub) sub.textContent = `across ${total} confirmed sequence${total !== 1 ? 's' : ''}`;

  if (!total) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:20px">No viral sequences.</div>';
    return;
  }

  // Geometry
  const W   = 220;
  const H   = 130;
  const cx  = W / 2;
  const cy  = H - 14;
  const rOuter = 110;
  const rInner = 68;

  const arc = d3.arc()
    .innerRadius(rInner).outerRadius(rOuter)
    .cornerRadius(4)
    .padAngle(0.015);

  // Half-donut: π/2 from -π/2 to π/2 won't work for half-pie at bottom;
  // we want a top-opening semicircle: start at -π/2 (top), through 0, to π/2.
  // Actually we want the FLAT side at the bottom (semicircle on top).
  // That means start = -π/2, end = π/2 (sweep clockwise across the top).
  const pie = d3.pie()
    .value(d => d[1])
    .sort(null)
    .startAngle(-Math.PI / 2)
    .endAngle(Math.PI / 2);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('width',  W)
    .attr('height', H)
    .style('overflow', 'visible');

  const g = svg.append('g').attr('transform', `translate(${cx}, ${cy})`);

  const slices = pie(display);
  slices.forEach((s, i) => {
    g.append('path')
      .attr('d', arc(s))
      .attr('fill', _familyColor(s.data[0], i))
      .attr('opacity', 0.92)
      .style('cursor', 'default')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(s.data[0])}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${s.data[1]}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${((s.data[1] / total) * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Center label (big number above the flat side, inside the donut hole)
  const labelG = svg.append('g').attr('transform', `translate(${cx}, ${cy - 20})`);
  labelG.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 22)
    .attr('font-weight', 700)
    .attr('fill', 'var(--vq-primary)')
    .attr('y', 0)
    .text(entries.length);
  labelG.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 10)
    .attr('fill', 'var(--vq-text-3)')
    .attr('y', 14)
    .text(rank === 'family' ? 'families' : 'phyla');

  wrap.appendChild(svg.node());

  // Legend rows
  display.forEach(([name, n], i) => {
    const pct = (n / total) * 100;
    const row = document.createElement('div');
    row.className = 'vq-chart-legend__row';
    row.innerHTML = `
      <span class="vq-chart-legend__dot" style="background:${_familyColor(name, i)}"></span>
      <span style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap"
            title="${VQ.esc(name)}">${VQ.esc(name)}</span>
      <span class="vq-chart-legend__count">${n}</span>
      <span class="vq-chart-legend__pct">${pct.toFixed(0)}%</span>
    `;
    lgd.appendChild(row);
  });
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 2 — Horizontal bar: HMM database hits
// ────────────────────────────────────────────────────────────────────────

// ────────────────────────────────────────────────────────────────────────
//  Top domains per HMM model — vertical bars, one bank at a time
// ────────────────────────────────────────────────────────────────────────

// Same bank colours as the HMM Database Hits card.
const _HMM_DB_COLORS = {
  RVDB:   'var(--vq-accent)',
  Vfam:   'var(--vq-primary-light)',
  EggNOG: 'var(--vq-accent-dark)',
  Pfam:   'var(--vq-success)',
};
const _HMM_DB_ORDER = ['Pfam', 'RVDB', 'Vfam', 'EggNOG'];   // Pfam first: most readable names
const _TOP_DOMAINS_N = 10;

/* HMM banks that have at least one domain, in display order. */
function _domainDatabases(sequences) {
  const present = new Set();
  sequences.forEach(s => (s.orfs || []).forEach(o =>
    (o.domains || []).forEach(d => d.database && present.add(d.database))));
  return [
    ..._HMM_DB_ORDER.filter(db => present.has(db)),
    ...[...present].filter(db => !_HMM_DB_ORDER.includes(db)).sort(),
  ];
}

/* Per model of one bank: sequences carrying it (ranking) and total hits. */
function _domainCounts(sequences, database) {
  const byModel = new Map();
  sequences.forEach(s => {
    (s.orfs || []).forEach(o => (o.domains || []).forEach(d => {
      if (d.database !== database) return;
      let m = byModel.get(d.target);
      if (!m) {
        m = { target: d.target, description: d.description || '', type: d.type || '',
              seqs: new Set(), hits: 0 };
        byModel.set(d.target, m);
      }
      m.seqs.add(s.id);
      m.hits += 1;
    }));
  });
  return [...byModel.values()]
    .map(m => ({ ...m, nSeqs: m.seqs.size }))
    .sort((a, b) => b.nSeqs - a.nSeqs || b.hits - a.hits || a.target.localeCompare(b.target));
}

function _renderTopDomains(sequences, database) {
  const wrap = document.getElementById('stats-domains-svg');
  if (!wrap) return;
  // Bars are replaced on every switch, so their mouseleave never fires —
  // drop any tooltip still showing from the previous bank.
  VQ.tooltipHide();
  wrap.innerHTML = '';

  const all  = _domainCounts(sequences, database);
  const data = all.slice(0, _TOP_DOMAINS_N);
  const sub  = document.getElementById('stats-domains-sub');
  if (sub) sub.textContent =
    `${database} · top ${data.length} of ${all.length} model${all.length === 1 ? '' : 's'} · sequences per model`;

  if (!data.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No domains.</div>';
    return;
  }

  const color = _HMM_DB_COLORS[database] || 'var(--vq-accent)';
  const W     = Math.max(wrap.clientWidth || 0, 220);
  const PAD_L = 28; const PAD_R = 6;
  const PAD_T = 14; const PAD_B = 58;           // room for rotated model names
  const H     = 200;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  const xScale = d3.scaleBand(data.map(d => d.target), [0, drawW]).padding(0.28);
  const yMax   = d3.max(data, d => d.nSeqs) || 1;
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();
  const short  = t => (t.length > 13 ? t.slice(0, 12) + '…' : t);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + integer tick labels (counts of sequences)
  yScale.ticks(Math.min(4, yMax)).filter(Number.isInteger).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text')
      .attr('x', -5).attr('y', yScale(tick) + 3)
      .attr('text-anchor', 'end').attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  data.forEach(d => {
    const x  = xScale(d.target);
    const bw = xScale.bandwidth();
    const bh = Math.max(drawH - yScale(d.nSeqs), 1);
    const r  = Math.min(Math.ceil(bw / 2), 5);
    const tip = evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">${VQ.esc(d.target)}</div>
      ${d.description ? `<div style="margin-bottom:4px">${VQ.esc(d.description)}</div>` : ''}
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Sequences</span><span>${d.nSeqs.toLocaleString()}</span>
        <span class="vq-tooltip__key">Hits</span><span>${d.hits.toLocaleString()}</span>
        ${d.type ? `<span class="vq-tooltip__key">Type</span><span>${VQ.esc(d.type)}</span>` : ''}
        <span class="vq-tooltip__key">Database</span><span>${VQ.esc(database)}</span>
      </div>`, evt);

    g.append('path')
      .attr('d', _barPath(x, yScale(d.nSeqs), bw, bh, r))
      .attr('fill', color).attr('opacity', 0.85)
      .on('mousemove', tip)
      .on('mouseleave', VQ.tooltipHide);

    // Count above the bar
    g.append('text')
      .attr('x', x + bw / 2).attr('y', yScale(d.nSeqs) - 4)
      .attr('text-anchor', 'middle').attr('font-size', 10).attr('font-weight', 600)
      .attr('fill', 'var(--vq-text)').attr('font-variant-numeric', 'tabular-nums')
      .text(d.nSeqs);

    // Model name under the bar, rotated; full name + description on hover
    g.append('text')
      .attr('transform', `translate(${x + bw / 2},${drawH + 8}) rotate(-38)`)
      .attr('text-anchor', 'end').attr('font-size', 9.5)
      .attr('font-family', 'var(--vq-font-mono)')
      .attr('fill', 'var(--vq-text-2)')
      .style('cursor', 'default')
      .text(short(d.target))
      .on('mousemove', tip)
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);

  wrap.appendChild(svg.node());
}

function _renderHMMBars(hmm) {
  const wrap = document.getElementById('stats-hmm-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const eggnog = hmm.eggnog_hits ?? hmm.eggnong_hits ?? 0;
  const data = [
    { name: 'RVDB',   value: hmm.rvdb_hits || 0,   color: 'var(--vq-accent)'        },
    { name: 'Vfam',   value: hmm.vfam_hits || 0,   color: 'var(--vq-primary-light)' },
    { name: 'EggNOG', value: eggnog,               color: 'var(--vq-accent-dark)'   },
    { name: 'Pfam',   value: hmm.pfam_hits || 0,   color: 'var(--vq-success)'       },
  ];

  const W = Math.max(wrap.clientWidth || 0, 220);
  const PAD_L = 64;
  const PAD_R = 36;
  const ROW_H = 38;
  const BAR_H = 16;
  const H = data.length * ROW_H + 16;
  const drawW = W - PAD_L - PAD_R;

  const maxV = d3.max(data, d => d.value) || 1;
  const xScale = d3.scaleLinear([0, maxV], [0, drawW]).nice();

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  data.forEach((d, i) => {
    const y = 8 + i * ROW_H;

    svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + BAR_H / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 12)
      .attr('font-weight', 500)
      .attr('fill', 'var(--vq-text-2)')
      .text(d.name);

    // Track
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', drawW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', 'var(--vq-bg)');

    // Value bar
    const bw = Math.max(xScale(d.value), 4);
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', bw).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', d.color)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(d.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Hits</span><span>${d.value.toLocaleString()}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    // Value text (after the bar)
    svg.append('text')
      .attr('x', PAD_L + bw + 8).attr('y', y + BAR_H / 2 + 4)
      .attr('font-size', 12)
      .attr('font-weight', 600)
      .attr('fill', 'var(--vq-text)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(d.value.toLocaleString());
  });

  wrap.appendChild(svg.node());
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 3 — Salmon: mapping-rate gauge + reads breakdown
// ────────────────────────────────────────────────────────────────────────

function _renderSalmonChart(salmon) {
  const wrap = document.getElementById('stats-salmon-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const W    = Math.max(wrap.clientWidth || 0, 200);
  const H    = 72;
  const rate = Math.max(0, Math.min(100, salmon.mapping_rate ?? 0));

  const rateColor = rate >= 30 ? 'var(--vq-success)'
                  : rate >= 15 ? 'var(--vq-accent)'
                              : 'var(--vq-warning)';

  const PAD_L = 4;
  const PAD_R = 4;
  const BAR_Y = 46;
  const BAR_H = 12;
  const drawW = W - PAD_L - PAD_R;
  const fillW = Math.max(drawW * rate / 100, 4);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMidYMin meet')
    .style('width', '100%').style('height', 'auto');

  // Large rate % + label
  svg.append('text')
    .attr('x', PAD_L).attr('y', 22)
    .attr('font-size', 26).attr('font-weight', 700)
    .attr('letter-spacing', '-0.02em')
    .attr('fill', 'var(--vq-primary)')
    .text(rate.toFixed(1) + '%');
  svg.append('text')
    .attr('x', PAD_L).attr('y', 33)
    .attr('font-size', 11).attr('fill', 'var(--vq-text-3)')
    .text('mapping rate');

  // Track
  svg.append('rect')
    .attr('x', PAD_L).attr('y', BAR_Y)
    .attr('width', drawW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', 'var(--vq-bg)');

  // Fill — grows left to right
  svg.append('rect')
    .attr('x', PAD_L).attr('y', BAR_Y)
    .attr('width', fillW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', rateColor).attr('opacity', 0.9)
    .on('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">Mapping rate</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Reads aligned</span><span>${rate.toFixed(1)}%</span>
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // End labels
  svg.append('text')
    .attr('x', PAD_L).attr('y', H - 3)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)').text('0');
  svg.append('text')
    .attr('x', W - PAD_R).attr('y', H - 3)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--vq-text-3)').text('100%');

  wrap.appendChild(svg.node());

  // Mini-rows below the gauge
  const rows = document.getElementById('stats-salmon-rows');
  if (rows) {
    const _fmt = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();
    const mappedReads = salmon.total_reads != null && salmon.mapping_rate != null
      ? Math.round(salmon.total_reads * salmon.mapping_rate / 100)
      : null;
    rows.innerHTML = [
      mappedReads != null
        ? _miniRow('Mapped reads',         _fmt(mappedReads))
        : '',
      salmon.viral_expressed != null
        ? _miniRow('Viral seqs expressed', _fmt(salmon.viral_expressed))
        : '',
      salmon.ref_hk_count
        ? _miniRow('HK ref genes',         _fmt(salmon.ref_hk_count))
        : '',
      salmon.conserved_count
        ? _miniRow('Conserved HK detected',
                   `${_fmt(salmon.conserved_detected)} / ${_fmt(salmon.conserved_count)}`)
        : '',
      salmon.pfam_hk_count
        ? _miniRow('Pfam HK seqs',         _fmt(salmon.pfam_hk_count))
        : '',
    ].join('');
  }
}

function _renderSalmonEmpty() {
  const wrap = document.getElementById('stats-salmon-svg');
  if (wrap) wrap.innerHTML =
    '<div class="vq-empty" style="padding:40px 16px;font-size:12px">No Salmon quantification data.</div>';
}

// ── Rounded-top bar path (flat bottom to sit on x-axis) ─────────────────
function _barPath(x, y, w, h, r) {
  if (h <= 0 || w <= 0) return '';
  r = Math.min(r, w / 2, h);
  return `M${x},${y+h} L${x},${y+r} Q${x},${y} ${x+r},${y} L${x+w-r},${y} Q${x+w},${y} ${x+w},${y+r} L${x+w},${y+h} Z`;
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 4 — Pipeline funnel
// ────────────────────────────────────────────────────────────────────────

function _renderFunnel(blast) {
  const wrap = document.getElementById('stats-funnel-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const steps = [
    { label: 'Input sequences',    value: blast.total_input          ?? null, color: 'var(--vq-primary)'       },
    { label: 'RefSeq BLASTx hits', value: blast.refseq_unique_seqs   ?? 0,   color: 'var(--vq-accent)'        },
    { label: 'Viral flagged',      value: blast.total_viral_flagged  ?? null, color: 'var(--vq-primary-light)' },
    { label: 'Confirmed viral',    value: blast.total_confirmed       ?? 0,   color: 'var(--vq-success)'       },
    { label: 'BLASTn annotated',   value: blast.blastn_unique_seqs   ?? null, color: 'var(--vq-accent-dark)'   },
  ].filter(s => s.value !== null);

  if (!steps.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:20px;font-size:12px">No pipeline data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 260);
  const ROW_H = 27;
  const BAR_H = 10;
  const H     = steps.length * ROW_H + 8;
  const base  = steps[0].value || 1;

  // Column positions: label | number | bar | pct
  const DOT_X   = 8;
  const LBL_X   = 20;
  const NUM_END = Math.round(W * 0.60);
  const BAR_X   = Math.round(W * 0.62);
  const BAR_W   = Math.round(W * 0.25);
  const PCT_X   = BAR_X + BAR_W + 5;

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  steps.forEach((step, i) => {
    const mid    = 4 + i * ROW_H + ROW_H / 2;
    const pctNum = i === 0 ? 100 : (base > 0 ? step.value / base * 100 : 0);
    const pctStr = i === 0 ? '100%' : pctNum.toFixed(1) + '%';
    const fillW  = Math.max(pctNum / 100 * BAR_W, 2);

    // Dot indicator
    svg.append('circle').attr('cx', DOT_X).attr('cy', mid).attr('r', 4)
      .attr('fill', step.color);

    // Connector to next step
    if (i < steps.length - 1)
      svg.append('line')
        .attr('x1', DOT_X).attr('y1', mid + 5)
        .attr('x2', DOT_X).attr('y2', mid + ROW_H - 3)
        .attr('stroke', 'var(--vq-border-dark)')
        .attr('stroke-width', 1.2).attr('stroke-dasharray', '2,2');

    // Step label (left)
    svg.append('text')
      .attr('x', LBL_X).attr('y', mid + 4)
      .attr('font-size', 11).attr('fill', 'var(--vq-text-2)')
      .text(step.label);

    // Count — bold, colour-matched (right-aligned at NUM_END)
    svg.append('text')
      .attr('x', NUM_END).attr('y', mid + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 12).attr('font-weight', 700)
      .attr('font-family', 'var(--vq-font-mono)')
      .attr('font-variant-numeric', 'tabular-nums')
      .attr('fill', step.color)
      .text(step.value.toLocaleString());

    // Bar track
    svg.append('rect')
      .attr('x', BAR_X).attr('y', mid - BAR_H / 2)
      .attr('width', BAR_W).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', 'var(--vq-bg)');

    // Bar fill — proportional to % of input
    svg.append('rect')
      .attr('x', BAR_X).attr('y', mid - BAR_H / 2)
      .attr('width', fillW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', step.color).attr('opacity', 0.8)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(step.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${step.value.toLocaleString()}</span>
          <span class="vq-tooltip__key">of input</span><span>${pctStr}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    // Pct label after bar
    svg.append('text')
      .attr('x', PCT_X).attr('y', mid + 4)
      .attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(pctStr);
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 5 — BLASTx identity histogram (green card)
// ────────────────────────────────────────────────────────────────────────

// ── Data-driven x axis for the histograms (starts at the smallest value) ──

const _PCT_BINS = 10;   // BLAST histograms: [min, 100%] split into equal bins

/* x domain from the data; widened when every value is the same, so the
   single bar still has width. */
function _dataDomain(lo, hi, minSpan) {
  return hi - lo >= minSpan ? [lo, hi] : [Math.max(0, hi - minSpan), hi];
}

/* Percent label: integers stay integers, others keep one decimal. */
function _fmtPct(v) {
  const r = Math.round(v * 10) / 10;
  return (Number.isInteger(r) ? String(r) : r.toFixed(1)) + '%';
}

/* Tick values for [lo, hi]: both ends always labelled, plus the round ticks
   between them that are not too close to either end. */
function _edgeTicks(scale, lo, hi, n) {
  const gap  = (hi - lo) * 0.12;
  const mids = scale.ticks(n).filter(t => t > lo + gap && t < hi - gap);
  return [lo, ...mids, hi];
}

/* End labels hug the axis ends so they are never clipped at the card edge. */
function _tickAnchor(t, lo, hi) {
  return t === lo ? 'start' : t === hi ? 'end' : 'middle';
}

/* Equal-width percent bins over [lo, 100]. */
function _pctBins(values, lo) {
  const width = (100 - lo) / _PCT_BINS;
  const thresholds = d3.range(1, _PCT_BINS).map(i => lo + i * width);
  return d3.bin().domain([lo, 100]).thresholds(thresholds)(values);
}

function _renderIdentityHistogram(sequences, mode = 'identity') {
  const wrap = document.getElementById('stats-identity-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const isCoverage = mode === 'coverage';
  const values = sequences
    .filter(s => s.is_viral)
    .map(s => {
      const hit = (s.blastx_nr_hits && s.blastx_nr_hits[0]) || (s.blastx_hits && s.blastx_hits[0]);
      if (!hit) return null;
      return isCoverage ? (hit.query_coverage ?? null) : hit.pct_identity;
    })
    .filter(v => v != null);

  if (!values.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No BLASTx data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 160);
  const PAD_L = 36; const PAD_R = 6;
  const PAD_T = 6;  const PAD_B = 20;
  const H     = 110;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  // x axis from the lowest observed value up to 100%
  const [lo]   = _dataDomain(d3.min(values), 100, 1);
  const bins   = _pctBins(values, lo);
  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([lo, 100], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  const metricLabel = isCoverage ? 'coverage' : 'identity';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + tick labels
  yScale.ticks(4).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text')
      .attr('x', -5).attr('y', yScale(tick) + 3)
      .attr('text-anchor', 'end').attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  // Bars — single accent colour, rounded top / flat bottom
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    const r = Math.min(Math.ceil(bw / 2), 5);
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, r))
      .attr('fill', 'var(--vq-accent)').attr('opacity', 0.82)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${_fmtPct(bin.x0)} – ${_fmtPct(bin.x1)} ${metricLabel}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / values.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks (keep x-axis)
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
  _edgeTicks(xScale, lo, 100, 3).forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 3)
      .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 12)
      .attr('text-anchor', _tickAnchor(tick, lo, 100)).attr('font-size', 9).attr('fill', 'var(--vq-c-tick)')
      .text(_fmtPct(tick));
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 6 — Sequence length histogram, viral only (blue card)
// ────────────────────────────────────────────────────────────────────────

/* Bin width for the length histogram, scaled to the longest sequence so short
   assemblies are not squeezed into a couple of 500 bp bars. */
function _lengthBinStep(maxLen) {
  if (maxLen <= 2000) return 100;
  if (maxLen <= 4000) return 300;
  return 500;
}

function _renderLengthHistogram(sequences) {
  const wrap = document.getElementById('stats-length-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const lengths = sequences
    .filter(s => s.is_viral && s.length > 0)
    .map(s => s.length);

  if (!lengths.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No sequence data.</div>';
    return;
  }

  function _fmtLen(v) {
    return v >= 1000 ? (v / 1000).toFixed(v % 1000 === 0 ? 0 : 1) + 'k' : String(Math.round(v));
  }

  const minL = d3.min(lengths);
  const maxL = d3.max(lengths);
  const STEP = _lengthBinStep(maxL);

  const sub = document.getElementById('stats-length-sub');
  if (sub) sub.textContent =
    `${lengths.length.toLocaleString()} viral sequences · median ${_fmtLen(d3.median(lengths))} bp · ${STEP} bp bins`;

  // Same frame as the BLASTx / BLASTn histograms; height fills the card
  const W     = Math.max(wrap.clientWidth || 0, 200);
  const PAD_L = 36; const PAD_R = 6;
  const PAD_T = 6;  const PAD_B = 28;
  const H     = _fillHeight(wrap, 130);
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  // x axis from the shortest to the longest sequence; fixed-width bins
  // (width from _lengthBinStep) starting exactly at the shortest one.
  const [lo, hi]   = _dataDomain(minL, maxL, STEP);
  const thresholds = d3.range(lo + STEP, hi, STEP);
  const bins       = d3.bin().domain([lo, hi]).thresholds(thresholds)(lengths);

  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([lo, hi], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + tick labels
  yScale.ticks(4).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text')
      .attr('x', -5).attr('y', yScale(tick) + 3)
      .attr('text-anchor', 'end').attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  // Bars — single accent colour, rounded top / flat bottom
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    const r  = Math.min(Math.ceil(bw / 2), 5);
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, r))
      .attr('fill', 'var(--vq-accent)').attr('opacity', 0.82)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${_fmtLen(bin.x0)} – ${_fmtLen(bin.x1)} bp</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / lengths.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
  _edgeTicks(xScale, lo, hi, Math.min(5, bins.length)).forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 3)
      .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 12)
      .attr('text-anchor', _tickAnchor(tick, lo, hi)).attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(_fmtLen(tick));
  });

  // X axis label
  svg.append('text')
    .attr('x', PAD_L + drawW / 2).attr('y', H - 2)
    .attr('text-anchor', 'middle').attr('font-size', 9)
    .attr('fill', 'var(--vq-c-tick)').text('length (bp)');

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 6b — BLASTn identity / coverage histogram (compact, white)
// ────────────────────────────────────────────────────────────────────────

function _renderBLAstnHistogram(sequences, mode) {
  const wrap = document.getElementById('stats-blastn-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const isIdentity = mode !== 'coverage';
  const values = sequences
    .filter(s => s.is_viral && s.blastn_hits && s.blastn_hits.length > 0)
    .map(s => {
      const h = s.blastn_hits[0];
      return isIdentity ? (h.pident ?? null) : (h.qcovhsp ?? null);
    })
    .filter(v => v != null);

  if (!values.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px;font-size:12px">No BLASTn data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 200);
  const PAD_L = 40; const PAD_R = 10;
  const PAD_T = 8;  const PAD_B = 28;
  const H     = 130;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  // x axis from the lowest observed value up to 100%
  const [lo]   = _dataDomain(d3.min(values), 100, 1);
  const bins   = _pctBins(values, lo);
  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([lo, 100], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  const barColor = 'var(--vq-primary-light)';
  const metricLabel = isIdentity ? 'identity' : 'coverage';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + tick labels
  yScale.ticks(4).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.8);
    g.append('text')
      .attr('x', -8).attr('y', yScale(tick) + 3.5)
      .attr('text-anchor', 'end').attr('font-size', 9.5)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  // Bars
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, Math.min(Math.ceil(bw / 2), 5)))
      .attr('fill', barColor).attr('opacity', 0.82)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${_fmtPct(bin.x0)} – ${_fmtPct(bin.x1)} ${metricLabel}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / values.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
  _edgeTicks(xScale, lo, 100, 4).forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 4)
      .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 14)
      .attr('text-anchor', _tickAnchor(tick, lo, 100)).attr('font-size', 9.5).attr('fill', 'var(--vq-c-tick)')
      .text(_fmtPct(tick));
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 8 — NR classification: Virus vs Phage (bars only)
// ────────────────────────────────────────────────────────────────────────

const _NR_CLASS_COLORS = {
  'Virus': 'var(--vq-accent)',
  'Phage': 'var(--vq-accent-dark)',
};

function _renderNRClassification(sequences) {
  const lgd   = document.getElementById('stats-nrclass-legend');
  const sub   = document.getElementById('stats-nrclass-sub');
  const title = document.getElementById('stats-nrclass-title');
  if (!lgd) return;
  lgd.innerHTML = '';

  const viral  = sequences.filter(s => s.is_viral);
  const nrUsed = viral.some(s => s.blastx_nr_hits && s.blastx_nr_hits.length > 0);
  if (title) title.textContent = nrUsed ? 'NR Classification' : 'BLASTx Classification';

  let virus = 0, phage = 0;
  viral.forEach(s => {
    const hit = (s.blastx_nr_hits && s.blastx_nr_hits[0]) || (s.blastx_hits && s.blastx_hits[0]);
    if (!hit) return;
    if (/phage|bacteriophage/i.test(hit.subject_title)) phage++;
    else virus++;
  });

  const total   = virus + phage;
  const display = [['Virus', virus], ['Phage', phage]].filter(([, n]) => n > 0);

  const card = document.getElementById('stats-nrclass-card');
  if (display.length < 2) {
    if (card) card.style.display = 'none';
    return;
  }
  if (card) card.style.display = '';

  if (sub) sub.textContent = `${total} classified sequence${total !== 1 ? 's' : ''}`;

  const color = name => _NR_CLASS_COLORS[name] || _TAX_DONUT_PALETTE[0];
  display.forEach(([name, n]) => {
    const pct = (n / total * 100).toFixed(1);
    const col = color(name);
    const row = document.createElement('div');
    row.style.cssText = 'display:flex;align-items:center;gap:10px;padding:6px 0;border-bottom:1px dashed var(--vq-border)';
    row.innerHTML = `
      <span style="width:14px;height:14px;border-radius:50%;background:${col};flex-shrink:0"></span>
      <div style="flex:1;min-width:0">
        <div style="font-size:15px;font-weight:600;color:var(--vq-text)">${VQ.esc(name)}</div>
        <div style="height:4px;background:var(--vq-bg);border-radius:2px;margin-top:5px">
          <div style="height:100%;width:${pct}%;background:${col};border-radius:2px;opacity:.85"></div>
        </div>
      </div>
      <div style="text-align:right;flex-shrink:0">
        <div style="font-size:15px;font-weight:700;font-variant-numeric:tabular-nums;color:var(--vq-text)">${n}</div>
        <div style="font-size:12px;color:var(--vq-text-3)">${pct}%</div>
      </div>
    `;
    lgd.appendChild(row);
  });
}


// ────────────────────────────────────────────────────────────────────────
//  Heuristic score gauge — stacked classification bar + avg score track
// ────────────────────────────────────────────────────────────────────────

function _renderHeuristicGauge(heur) {
  const wrap = document.getElementById('stats-heur-gauge');
  if (!wrap) return;
  wrap.innerHTML = '';

  const total   = heur.scored || 1;
  const known   = heur.viral_known   || 0;
  const unknown = heur.viral_unknown || 0;
  const nonvir  = heur.non_viral     || 0;
  const avg     = heur.avg_score     ?? 0;

  const W    = Math.max(wrap.clientWidth || 0, 180);
  const H    = 62;
  const BAR_Y = 6;
  const BAR_H = 10;

  // Colour-code the avg-score track
  const scoreColor = avg >= 70 ? 'var(--vq-success)'
                   : avg >= 40 ? 'var(--vq-accent)'
                               : 'var(--vq-warning)';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // ── Top: stacked classification bar ──────────────────────────────────
  const segments = [
    { label: 'Viral known',   n: known,   color: 'var(--vq-success)'      },
    { label: 'Viral unknown', n: unknown, color: 'var(--vq-accent)'       },
    { label: 'Non-viral',     n: nonvir,  color: 'var(--vq-warning)'      },
  ].filter(s => s.n > 0);

  // Clip segments to a rounded rect so the stacked fill gets rounded caps
  // without square corners poking out behind an overlay stroke.
  const clipId = `vq-heur-bar-clip-${Math.random().toString(36).slice(2)}`;
  svg.append('defs').append('clipPath')
    .attr('id', clipId)
    .append('rect')
      .attr('x', 0).attr('y', BAR_Y)
      .attr('width', W).attr('height', BAR_H)
      .attr('rx', BAR_H / 2);

  const barG = svg.append('g').attr('clip-path', `url(#${clipId})`);

  let cx = 0;
  segments.forEach(seg => {
    const segW = Math.max((seg.n / total) * W, 2);
    barG.append('rect')
      .attr('x', cx).attr('y', BAR_Y)
      .attr('width', segW).attr('height', BAR_H)
      .attr('fill', seg.color).attr('opacity', 0.85)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(seg.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${seg.n}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(seg.n / total * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
    cx += segW;
  });

  // ── Bottom: avg score track ───────────────────────────────────────────
  const TRACK_Y = BAR_Y + BAR_H + 10;
  const fillW   = Math.max((avg / 100) * W, 3);

  svg.append('text')
    .attr('x', 0).attr('y', TRACK_Y - 2)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)')
    .text('avg score');

  // Track background
  svg.append('rect')
    .attr('x', 0).attr('y', TRACK_Y + 4)
    .attr('width', W).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', 'var(--vq-bg)');

  // Fill
  svg.append('rect')
    .attr('x', 0).attr('y', TRACK_Y + 4)
    .attr('width', fillW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', scoreColor).attr('opacity', 0.9)
    .on('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">Average heuristic score</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Score</span><span>${avg.toFixed(1)}</span>
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // End labels
  svg.append('text')
    .attr('x', 0).attr('y', TRACK_Y + BAR_H + 14)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)').text('0');
  svg.append('text')
    .attr('x', W).attr('y', TRACK_Y + BAR_H + 14)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--vq-text-3)').text('100');

  wrap.appendChild(svg.node());
}


window.vqInitStats = vqInitStats;
})();
