/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_viruses.js — Virome tab: what was found
   Virome KPIs, the detected-viruses table, then grouped cards
   (Taxonomy · Similarity to known viruses · Sequences & proteins).
   The Run tab (pipeline / technical) lives in section_stats.js.
   ============================================================ */

function vqInitViruses(report) {
  _initViromeTab(report);
}

// Overview helpers shared with the other overview tab (defined in export.js).
const _chip           = VQ.statChip;
const _fmtNum         = VQ.fmtNum;
const _group          = VQ.cardGroup;
const _pct            = VQ.pct;
const _redrawOnResize = VQ.redrawOnResize;

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
  const known     = viruses.filter(v => !v.noHit && v.known > 0).length;

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
            <button class="vq-toggle__btn" type="button" data-status="viral-known">Known</button>
            <button class="vq-toggle__btn" type="button" data-status="viral-unknown">Unknown</button>
            <button class="vq-toggle__btn" type="button" data-status="non-viral">Non-viral</button>
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
        ${_chip('Species', _fmtNum(viruses.filter(v => !v.noHit).length), 'accent', 'distinct best-hit species')}
        ${_chip('Families', _fmtNum(families.size), '', 'with resolved taxonomy')}
        ${_chip('Known Viruses', _fmtNum(known), 'success',
                'species with ≥1 viral-known contig')}
      </div>

      ${_group('Viruses', 'one row per species of the best BLASTx hit · classes from the heuristic score', table)}
      ${_group('Taxonomy', 'what kinds of viruses were found', taxonomy)}
      ${_group('Similarity to known viruses', 'how close each sequence is to its best reference', similarity)}
      ${_group('Sequences & proteins', 'contig sizes and protein domains',
               sequences, 1 + (domainDbs.length ? 2 : 0))}
    </div>`;

  // ── Detected viruses table (search + class filter + sorting) ──
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
//  Virome helpers — per-contig classes, per-species rows, genome types, scatter
// ────────────────────────────────────────────────────────────────────────

// Per-contig classes are the heuristic scorer's (score_heuristic._classify):
// the same classes the Sequence Viewer, the LLM prompt and the stats use.
// "Known" = >=90% identity AND >=70% query coverage to the best BLASTx hit
// (aa) or to the best *viral* BLASTn hit (nt); "non-viral" = no viral evidence,
// or a dominant non-viral BLASTn hit (likely host / contaminant).
const _KNOWN_ID  = 90;
const _KNOWN_COV = 70;
const _CLASS_ORDER = ['viral-known', 'viral-unknown', 'non-viral'];
const _CLASS_LABEL = {
  'viral-known':   'Known',
  'viral-unknown': 'Unknown',
  'non-viral':     'Non-viral',
};
// Same palette as the Viewer's badges / score meters (.vq-fill--*)
const _CLASS_COLOR = {
  'viral-known':   '#1f7a3c',
  'viral-unknown': '#9b6a09',
  'non-viral':     '#a8302b',
};
const _CLASS_CRITERION =
  `known = ≥${_KNOWN_ID}% identity and ≥${_KNOWN_COV}% coverage to the best BLASTx hit (aa) ` +
  `or to a viral BLASTn hit (nt)`;

// Viral-subject matcher — mirrors score_heuristic._VIRAL_RE.
const _VIRAL_RE = /vir(?:us|al|idae|ales|inae|oid|ion|aceae)|phage|bacteriophage/i;

/* Best BLASTx hit of a sequence: NR first, RefSeq as fallback (= best_blastx). */
function _bestBlastx(s) {
  const hits = (s.blastx_nr_hits && s.blastx_nr_hits.length) ? s.blastx_nr_hits : (s.blastx_hits || []);
  return hits.length ? hits.reduce((a, h) => (h.bit_score > a.bit_score ? h : a)) : null;
}

function _bestBlastn(s) {
  const hits = s.blastn_hits || [];
  return hits.length ? hits.reduce((a, h) => ((h.bit_score ?? 0) > (a.bit_score ?? 0) ? h : a)) : null;
}

function _bestViralBlastn(s) {
  const hits = (s.blastn_hits || []).filter(h => _VIRAL_RE.test(h.stitle || ''));
  return hits.length ? hits.reduce((a, h) => ((h.bit_score ?? 0) > (a.bit_score ?? 0) ? h : a)) : null;
}

/* Class of one contig. Reports from versions without the heuristic scorer get
   the same rule recomputed here (minus its host / contaminant check). */
function _seqClass(s) {
  const c = s.heuristic_output?.classification;
  if (_CLASS_LABEL[c]) return c;
  const bn = _bestViralBlastn(s);
  if (bn && bn.pident >= _KNOWN_ID && (bn.qcovhsp ?? 0) >= _KNOWN_COV) return 'viral-known';
  const bx = _bestBlastx(s);
  if (bx && bx.pct_identity >= _KNOWN_ID && (bx.query_coverage ?? 0) >= _KNOWN_COV) return 'viral-known';
  return 'viral-unknown';
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

const _NO_BLASTX = 'No BLASTx hit';

/* One row per best-hit species, plus one row for contigs with no BLASTx hit
   (viral evidence from HMM / BLASTn only — often the most novel candidates). */
function _virusRows(viral) {
  const groups = new Map();
  viral.forEach(s => {
    const hit = _bestBlastx(s);
    const key = hit ? _hitSpecies(hit) : _NO_BLASTX;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push({ seq: s, hit, bn: _bestBlastn(s), cls: _seqClass(s) });
  });

  return [...groups.entries()].map(([species, items]) => {
    const noHit  = species === _NO_BLASTX;
    const counts = { 'viral-known': 0, 'viral-unknown': 0, 'non-viral': 0 };
    items.forEach(i => { counts[i.cls]++; });
    const withHit = items.filter(i => i.hit);
    const best = withHit.length
      ? withHit.reduce((a, b) => (b.hit.pct_identity > a.hit.pct_identity ? b : a)) : null;
    const bnBest = items.map(i => i.bn).filter(Boolean)
      .reduce((a, b) => (!a || b.pident > a.pident ? b : a), null);
    // Row click opens a known contig when there is one, else the closest match,
    // else (no BLASTx) the highest-scoring contig.
    const score = i => i.seq.heuristic_output?.vq_score ?? 0;
    const open = items.find(i => i.cls === 'viral-known' && i === best)
      || items.filter(i => i.cls === 'viral-known').sort((a, b) => score(b) - score(a))[0]
      || best
      || items.slice().sort((a, b) => score(b) - score(a))[0];
    return {
      species, noHit, counts,
      known:    counts['viral-known'],
      family:   _mostCommon(items.map(i => i.seq.taxonomy?.family)) || '',
      genome:   _mostCommon(items.map(i => i.seq.taxonomy?.genome)) || '',
      contigs:  items.length,
      length:   items.reduce((a, i) => a + (i.seq.length || 0), 0),
      identity: best ? best.hit.pct_identity : null,
      coverage: best ? (best.hit.query_coverage ?? 0) : null,
      blastn:   bnBest ? bnBest.pident : null,
      accession: best ? best.hit.subject_id : null,
      seqId:    open.seq.id,
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

/* Contigs of a row split by class: stacked bar + "8 known · 23 unknown". */
function _classCell(counts, total) {
  const parts = _CLASS_ORDER.filter(c => counts[c]);
  const bar = parts.map(c =>
    `<span class="vq-fill--${c}" style="width:${(counts[c] / total * 100).toFixed(1)}%"></span>`).join('');
  const text = parts.map(c =>
    `<span class="vq-vt-cls vq-vt-cls--${c}">${counts[c]} ${_CLASS_LABEL[c].toLowerCase()}</span>`)
    .join('<span class="vq-vt-dot">·</span>');
  return `<span class="vq-vt-classes"><span class="vq-vt-classes__bar">${bar}</span>${text}</span>`;
}

function _renderVirusTable(rows, st) {
  const wrap = document.getElementById('stats-viruses-table');
  if (!wrap) return;
  const esc = VQ.esc;

  // Filter counts = species having at least one contig of that class.
  const counts = { all: rows.length };
  _CLASS_ORDER.forEach(c => { counts[c] = rows.filter(r => r.counts[c] > 0).length; });
  document.querySelectorAll('#stats-viruses-status [data-status]').forEach(b => {
    const k = b.dataset.status;
    b.textContent = `${k === 'all' ? 'All' : _CLASS_LABEL[k]} · ${counts[k]}`;
    b.hidden = k !== 'all' && counts[k] === 0;
  });

  let view = rows.filter(r =>
    (st.status === 'all' || r.counts[st.status] > 0) &&
    (!st.q || `${r.species} ${r.family} ${r.genome}`.toLowerCase().includes(st.q)));
  const key = st.sort;
  view = view.sort((a, b) => {
    if (a.noHit !== b.noHit) return a.noHit ? 1 : -1;          // HMM-only row last
    const x = a[key], y = b[key];
    if (typeof x === 'string' || typeof y === 'string')
      return st.dir * String(x || '').localeCompare(String(y || ''));
    return st.dir * ((x ?? -1) - (y ?? -1)) || a.species.localeCompare(b.species);
  });

  const nSpecies = rows.filter(r => !r.noHit).length;
  const noHit = rows.find(r => r.noHit);
  const sub = document.getElementById('stats-viruses-sub');
  if (sub) sub.textContent =
    `${nSpecies} species${noHit ? ` + ${noHit.contigs} contigs without a BLASTx hit` : ''} · ` +
    `classes per contig from the heuristic score (${_CLASS_CRITERION}) · click a row to open its best contig`;

  if (!rows.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px;font-size:12px">No viral sequences.</div>';
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
        ${th('known',    'Contigs by class')}
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
          <tr data-seq="${esc(r.seqId)}"${r.noHit ? ' class="vq-vt-nohit"' : ''}
              title="Open ${esc(r.seqId)} in the Sequence Viewer">
            <td class="vq-vt-species">${r.noHit
              ? `${esc(r.species)}<span class="vq-vt-note">viral evidence from HMM / BLASTn only</span>`
              : esc(r.species)}</td>
            <td>${_classCell(r.counts, r.contigs)}</td>
            <td>${r.family ? esc(r.family) : '<span class="vq-vt-na">—</span>'}</td>
            <td>${r.genome ? esc(r.genome) : '<span class="vq-vt-na">—</span>'}</td>
            <td class="vq-vt-num">${r.contigs.toLocaleString()}</td>
            <td class="vq-vt-num">${r.length.toLocaleString()}</td>
            <td>${_pctCell(r.identity)}</td>
            <td>${_pctCell(r.coverage)}</td>
            ${hasBn ? `<td>${_pctCell(r.blastn)}</td>` : ''}
            <td class="vq-td--mono">${url
              ? `<a class="vq-acc-link" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(r.accession)}</a>`
              : `<span class="vq-vt-na">${esc(r.accession || '—')}</span>`}</td>
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
      s, id, cov, cls: _seqClass(s),
      species: isX ? _hitSpecies(h) : (h.stitle || ''),
    };
  }).filter(Boolean);

  const sub = document.getElementById('stats-scatter-sub');
  if (sub) sub.textContent =
    `${pts.length} sequences · best ${isX ? 'BLASTx' : 'BLASTn'} hit · colour = heuristic class · ` +
    `shaded area = ≥${_KNOWN_ID}% id and ≥${_KNOWN_COV}% cov · click a point to open it`;

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

  // Points — rarer classes drawn last so they stay visible on top
  const order = { 'viral-unknown': 0, 'viral-known': 1, 'non-viral': 2 };
  pts.sort((a, b) => order[a.cls] - order[b.cls]).forEach(p => {
    g.append('circle')
      .attr('cx', x(p.id)).attr('cy', y(p.cov)).attr('r', 4.2)
      .attr('fill', _CLASS_COLOR[p.cls]).attr('fill-opacity', 0.72)
      .attr('stroke', 'var(--vq-surface)').attr('stroke-width', 1)
      .style('cursor', 'pointer')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(p.s.id)}</div>
        <div style="margin-bottom:4px;font-style:italic">${esc(p.species)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Identity</span><span>${p.id.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Coverage</span><span>${Number(p.cov).toFixed(1)}%</span>
          <span class="vq-tooltip__key">Family</span><span>${esc(p.s.taxonomy?.family || '—')}</span>
          <span class="vq-tooltip__key">Class</span><span>${_CLASS_LABEL[p.cls]}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide)
      .on('click', () => { VQ.tooltipHide(); VQ.jumpToViewer(p.s.id); });
  });

  // Legend
  const counts = { 'viral-known': 0, 'viral-unknown': 0, 'non-viral': 0 };
  pts.forEach(p => counts[p.cls]++);
  const lg = svg.append('g').attr('transform', `translate(${PAD_L + 4},10)`);
  let lx = 0;
  _CLASS_ORDER.filter(k => counts[k]).forEach(k => {
    const item = lg.append('g').attr('transform', `translate(${lx},0)`);
    item.append('circle').attr('r', 4).attr('cy', 0).attr('fill', _CLASS_COLOR[k]);
    const t = item.append('text').attr('x', 8).attr('y', 3.5).attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)').text(`${_CLASS_LABEL[k]} ${counts[k]}`);
    lx += 8 + (t.node().getComputedTextLength?.() || 60) + 16;
  });

  wrap.appendChild(svg.node());
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

// ── Rounded-top bar path (flat bottom to sit on x-axis) ─────────────────
function _barPath(x, y, w, h, r) {
  if (h <= 0 || w <= 0) return '';
  r = Math.min(r, w / 2, h);
  return `M${x},${y+h} L${x},${y+r} Q${x},${y} ${x+r},${y} L${x+w-r},${y} Q${x+w},${y} ${x+w},${y+r} L${x+w},${y+h} Z`;
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


window.vqInitViruses = vqInitViruses;
})();
