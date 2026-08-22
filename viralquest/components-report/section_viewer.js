/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_viewer.js — Section 3: Sequence Viewer
   Searchable / filterable list of sequence cards.
   Expanded card shows:
     • genome track (lane-packed ORFs, lane-packed domains)
     • BLASTn / BLASTx-RefSeq / BLASTx-NR hits (subtabs)
     • LLM analysis
     • FASTA preview
   Export per-card: PNG / SVG / PDF / FASTA
   Batch export: composite PNG/SVG (proportions kept) + combined FASTA
   ============================================================ */

const _VW = {
  sequences: [],
  filtered:  [],
  selected:  new Set(),
  domainMap: {},
  colourIdx: 0,
  // Multi-select taxonomy filter state (checkbox dropdowns)
  tax: { family: new Set(), phylum: new Set(), genus: new Set() },
  // Low-complexity highlighting (genome-map bands + FASTA highlight) is per
  // sequence and on by default, so this holds only the ids the user switched
  // OFF. Deliberately NOT reset by vqInitViewer so the choices survive a
  // re-init (the multi-sample report re-mounts the viewer on every sample
  // filter change). Exports read the live SVG, so they follow these flags.
  lowCxOff: new Set(),
};

/* Should this sequence's low-complexity regions be highlighted? */
function _showLowCx(seq) {
  return !_VW.lowCxOff.has(seq.id);
}


/* Build a checkbox-dropdown filter field (multi-select).
   Lives inside the retractable filter panel; toggling is handled by the
   global .vq-menu open/close logic, with stopPropagation inside the panel
   so checking several boxes doesn't dismiss the dropdown. */
function _msField(id, label, allLabel, values) {
  const esc = VQ.esc;
  const items = values.length
    ? values.map(v =>
        `<label class="vq-ms__item">
           <input type="checkbox" value="${esc(v)}">${esc(v)}
         </label>`).join('')
    : '<div class="vq-ms__empty">None available</div>';
  return `
    <div class="vq-filter-field">
      <span class="vq-filter-field__label">${esc(label)}</span>
      <div class="vq-menu vq-ms" id="${id}">
        <button class="vq-ms__toggle" data-menu-toggle type="button">
          <span class="vq-ms__value" data-ms-label data-all-label="${esc(allLabel)}">${esc(allLabel)}</span>
          <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2.5" aria-hidden="true">
            <polyline points="6 9 12 15 18 9"/>
          </svg>
        </button>
        <div class="vq-menu__panel vq-ms__panel" role="menu">${items}</div>
      </div>
    </div>`;
}

const _DOM_COLOURS = [
  'var(--vq-dom-1)','var(--vq-dom-2)','var(--vq-dom-3)','var(--vq-dom-4)',
  'var(--vq-dom-5)','var(--vq-dom-6)','var(--vq-dom-7)','var(--vq-dom-8)',
];

// ── Init ────────────────────────────────────────────────────────────────────

function vqInitViewer(sequences, mountId) {
  const el = document.getElementById(mountId || 'section-viewer');
  if (!el) return;
  const esc = VQ.esc;

  _VW.sequences = sequences;
  _VW.filtered  = [...sequences];
  _VW.tax = { family: new Set(), phylum: new Set(), genus: new Set() };

  const phyla    = [...new Set(sequences.map(s => s.taxonomy?.phylum).filter(Boolean))].sort();
  const families = [...new Set(sequences.map(s => s.taxonomy?.family).filter(Boolean))].sort();
  const genera   = [...new Set(sequences.map(s => s.taxonomy?.genus).filter(Boolean))].sort();
  const hasHeur  = sequences.some(s => s.heuristic_output != null);
  const hasLlm   = sequences.some(s => s.llm_output != null);

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">
          Sequence Viewer
          <span class="vq-count-label" id="viewer-count">
            ${sequences.length} sequence${sequences.length !== 1 ? 's' : ''}
          </span>
        </div>
      </div>
      <div class="vq-section-actions">
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-expand-all" type="button">
          Expand all
        </button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-collapse-all" type="button">
          Collapse all
        </button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-select-all" type="button">
          Select all
        </button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-deselect" type="button">
          Deselect
        </button>
        <div class="vq-menu" id="viewer-batch-menu">
          <button class="vq-btn vq-btn--primary vq-btn--sm" data-menu-toggle id="viewer-batch-btn"
                  type="button" disabled>
            <span id="viewer-batch-label">Export selected</span>
            <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="2.5">
              <polyline points="6 9 12 15 18 9"/>
            </svg>
          </button>
          <div class="vq-menu__panel" role="menu">
            <button class="vq-menu__item" data-batch="png" role="menuitem" type="button">
              Combined PNG (one image)
            </button>
            <button class="vq-menu__item" data-batch="svg" role="menuitem" type="button">
              Combined SVG (one file)
            </button>
            <button class="vq-menu__item" data-batch="png-each" role="menuitem" type="button">
              PNG — one per sequence
            </button>
            <button class="vq-menu__item" data-batch="svg-each" role="menuitem" type="button">
              SVG — one per sequence
            </button>
            <button class="vq-menu__item" data-batch="fasta" role="menuitem" type="button">
              FASTA (selected sequences)
            </button>
          </div>
        </div>
      </div>
    </div>

    <div class="vq-toolbar">
      <div class="vq-search-wrap">
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2" aria-hidden="true">
          <circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/>
        </svg>
        <input class="vq-input" id="viewer-search" type="search"
               placeholder="Search ID, species, family, genus, phylum…"
               aria-label="Search sequences">
      </div>
      <select class="vq-select" id="viewer-sort" aria-label="Sort sequences">
        <option value="score-desc">Score: high → low</option>
        <option value="score-asc">Score: low → high</option>
        <option value="len-desc">Length: long → short</option>
        <option value="len-asc">Length: short → long</option>
        <option value="id-asc">Sequence ID (A→Z)</option>
      </select>
      <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-filter-toggle"
              type="button" aria-expanded="false" aria-controls="viewer-filter-panel">
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2" aria-hidden="true">
          <line x1="4" y1="6" x2="20" y2="6"/><circle cx="9" cy="6" r="2"/>
          <line x1="4" y1="12" x2="20" y2="12"/><circle cx="15" cy="12" r="2"/>
          <line x1="4" y1="18" x2="20" y2="18"/><circle cx="9" cy="18" r="2"/>
        </svg>
        Filters
        <span class="vq-filter-badge" id="viewer-filter-badge" hidden>0</span>
      </button>
    </div>

    <div class="vq-filter-panel" id="viewer-filter-panel" hidden>
      <div class="vq-filter-panel__grid">

        <div class="vq-filter-field">
          <span class="vq-filter-field__label">Length (nt)</span>
          <div class="vq-range">
            <input class="vq-input vq-input--sm" type="number" id="viewer-len-min"
                   min="0" step="1" placeholder="min" aria-label="Min sequence length">
            <span class="vq-filter-sep">–</span>
            <input class="vq-input vq-input--sm" type="number" id="viewer-len-max"
                   min="0" step="1" placeholder="max" aria-label="Max sequence length">
          </div>
        </div>

        <div class="vq-filter-field vq-filter-field--wide">
          <span class="vq-filter-field__label">BLAST hits</span>
          <div class="vq-blast-filter">
            <span class="vq-filter-sub">Identity&nbsp;%</span>
            <div class="vq-range">
              <input class="vq-input vq-input--sm" type="number" id="viewer-ident-min"
                     min="0" max="100" step="0.1" placeholder="min" aria-label="Min identity %">
              <span class="vq-filter-sep">–</span>
              <input class="vq-input vq-input--sm" type="number" id="viewer-ident-max"
                     min="0" max="100" step="0.1" placeholder="max" aria-label="Max identity %">
            </div>
            <span class="vq-filter-sub">Coverage&nbsp;%</span>
            <div class="vq-range">
              <input class="vq-input vq-input--sm" type="number" id="viewer-cov-min"
                     min="0" max="100" step="0.1" placeholder="min" aria-label="Min coverage %">
              <span class="vq-filter-sep">–</span>
              <input class="vq-input vq-input--sm" type="number" id="viewer-cov-max"
                     min="0" max="100" step="0.1" placeholder="max" aria-label="Max coverage %">
            </div>
            <label class="vq-filter-check">
              <input type="checkbox" id="viewer-filter-blastn"> BLASTn
            </label>
            <label class="vq-filter-check">
              <input type="checkbox" id="viewer-filter-blastx"> BLASTx NR
            </label>
          </div>
        </div>

        ${_msField('viewer-ms-family', 'Family', 'All families', families)}
        ${_msField('viewer-ms-phylum', 'Phylum', 'All phyla',    phyla)}
        ${_msField('viewer-ms-genus',  'Genus',  'All genera',   genera)}

        ${hasHeur ? `
        <div class="vq-filter-field">
          <span class="vq-filter-field__label">Heuristic score</span>
          <div class="vq-range">
            <input class="vq-input vq-input--sm" type="number" id="viewer-heur-min"
                   min="0" max="100" step="1" placeholder="min" aria-label="Min heuristic score">
            <span class="vq-filter-sep">–</span>
            <input class="vq-input vq-input--sm" type="number" id="viewer-heur-max"
                   min="0" max="100" step="1" placeholder="max" aria-label="Max heuristic score">
          </div>
        </div>` : ''}

        ${hasLlm ? `
        <div class="vq-filter-field">
          <span class="vq-filter-field__label">LLM score</span>
          <div class="vq-range">
            <input class="vq-input vq-input--sm" type="number" id="viewer-llm-min"
                   min="0" max="100" step="1" placeholder="min" aria-label="Min LLM score">
            <span class="vq-filter-sep">–</span>
            <input class="vq-input vq-input--sm" type="number" id="viewer-llm-max"
                   min="0" max="100" step="1" placeholder="max" aria-label="Max LLM score">
          </div>
        </div>` : ''}

      </div>
      <div class="vq-filter-panel__foot">
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-filter-clear" type="button">
          Clear filters
        </button>
      </div>
    </div>

    <div id="viewer-list"></div>
  `;

  ['viewer-search','viewer-sort',
   'viewer-ident-min','viewer-ident-max','viewer-cov-min','viewer-cov-max',
   'viewer-len-min','viewer-len-max',
   'viewer-heur-min','viewer-heur-max','viewer-llm-min','viewer-llm-max']
    .forEach(id => document.getElementById(id)?.addEventListener('input', _applyFilters));
  ['viewer-filter-blastn','viewer-filter-blastx']
    .forEach(id => document.getElementById(id)?.addEventListener('change', _applyFilters));

  _wireMultiSelect('viewer-ms-family', _VW.tax.family);
  _wireMultiSelect('viewer-ms-phylum', _VW.tax.phylum);
  _wireMultiSelect('viewer-ms-genus',  _VW.tax.genus);

  // Retractable filter panel
  const filterToggle = document.getElementById('viewer-filter-toggle');
  const filterPanel  = document.getElementById('viewer-filter-panel');
  filterToggle?.addEventListener('click', () => {
    const willOpen = filterPanel.hidden;
    filterPanel.hidden = !willOpen;
    filterToggle.setAttribute('aria-expanded', String(willOpen));
    filterToggle.classList.toggle('active', willOpen);
  });
  document.getElementById('viewer-filter-clear')?.addEventListener('click', () => {
    filterPanel.querySelectorAll('input[type="number"]').forEach(i => { i.value = ''; });
    filterPanel.querySelectorAll('input[type="checkbox"]').forEach(c => { c.checked = false; });
    Object.values(_VW.tax).forEach(set => set.clear());
    filterPanel.querySelectorAll('[data-ms-label]').forEach(l => { l.textContent = l.dataset.allLabel; });
    filterPanel.querySelectorAll('.vq-ms').forEach(m => m.classList.remove('vq-ms--active'));
    _applyFilters();
  });

  document.getElementById('viewer-select-all')?.addEventListener('click', () => {
    _VW.filtered.forEach(s => _VW.selected.add(s.id));
    _renderList(); _updateBatchBtn();
  });
  document.getElementById('viewer-deselect')?.addEventListener('click', () => {
    _VW.selected.clear();
    _renderList(); _updateBatchBtn();
  });
  document.getElementById('viewer-expand-all')?.addEventListener('click', () => {
    document.querySelectorAll('.vq-seq-card').forEach(c => _openCard(c));
  });
  document.getElementById('viewer-collapse-all')?.addEventListener('click', () => {
    document.querySelectorAll('.vq-seq-card.open').forEach(c => {
      c.classList.remove('open');
      c.querySelector('.vq-seq-card__head')?.setAttribute('aria-expanded', 'false');
    });
  });

  const batchMenu = document.getElementById('viewer-batch-menu');
  const batchBtn  = document.getElementById('viewer-batch-btn');
  batchBtn?.addEventListener('click', e => {
    if (batchBtn.disabled) return;
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => {
      if (m !== batchMenu) m.classList.remove('open');
    });
    batchMenu.classList.toggle('open');
  });
  batchMenu?.querySelectorAll('[data-batch]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      batchMenu.classList.remove('open');
      _runBatch(item.dataset.batch);
    });
  });

  _sortFiltered();
  _renderList();
}

// ── Scoring helpers ─────────────────────────────────────────────────────────

/* Heuristic-primary: the rule-based score leads; the LLM score is the
   fallback only when no heuristic score is present. */
function _primaryScore(seq) {
  const h = seq.heuristic_output, l = seq.llm_output;
  if (h) return { score: h.vq_score, cls: h.classification, src: 'heur' };
  if (l) return { score: l.vq_score, cls: l.classification, src: 'llm'  };
  return null;
}

/* Compact header chip: number + a meter whose WIDTH encodes magnitude (0–100)
   and whose COLOUR encodes classification (reuses the badge palette). */
function _scoreChip(srcLabel, score, cls) {
  const esc = VQ.esc;
  const w   = Math.max(0, Math.min(100, score ?? 0));
  return `
    <span class="vq-score-chip" title="${esc(srcLabel)} VQ score ${score} · ${esc(cls)}">
      <span class="vq-score-chip__src">${esc(srcLabel)}</span>
      <span class="vq-score-chip__val">${score}</span>
      <span class="vq-score-chip__meter">
        <span class="vq-fill--${esc(cls)}" style="width:${w}%"></span>
      </span>
    </span>`;
}

function _scoreChips(seq) {
  const esc = VQ.esc;
  const h = seq.heuristic_output, l = seq.llm_output;
  if (!h && !l) return '';
  const chips = [];
  if (h) chips.push(_scoreChip('H',  h.vq_score, h.classification));
  if (l) chips.push(_scoreChip('AI', l.vq_score, l.classification));
  const disagree = (h && l && h.classification !== l.classification)
    ? `<span class="vq-score-disagree" title="Heuristic and AI disagree: ${esc(h.classification)} vs ${esc(l.classification)}">⚠</span>`
    : '';
  return `<span class="vq-score-chips">${chips.join('')}${disagree}</span>`;
}

/* One heuristic component sub-score bar (width = magnitude; null = absent). */
function _compBar(label, val, weight) {
  if (val == null) {
    return `<div class="vq-comp"><span class="vq-comp__lbl">${label}</span>
      <span class="vq-comp__na">not present</span></div>`;
  }
  const w  = Math.max(0, Math.min(100, val));
  const wt = weight != null ? ` title="weight ${(weight * 100).toFixed(0)}%"` : '';
  return `
    <div class="vq-comp"${wt}>
      <span class="vq-comp__lbl">${label}</span>
      <span class="vq-comp__track"><span style="width:${w}%"></span></span>
      <span class="vq-comp__val">${Math.round(val)}</span>
    </div>`;
}

function _heurColumn(h, safe) {
  const esc = VQ.esc;
  const c   = h.components || null;
  const wt  = c?.weights || {};
  const comps = c ? `
    ${_compBar('BLASTn', c.blastn, wt.blastn)}
    ${_compBar('BLASTx', c.blastx, wt.blastx)}
    ${_compBar('HMM',    c.hmm,    wt.hmm)}` : '';
  const penalty = (c && c.penalty != null && c.penalty < 1)
    ? `<div class="vq-score-col__penalty">Score penalty ×${c.penalty.toFixed(2)} — see analysis</div>`
    : '';
  return `
    <div class="vq-score-col">
      <div class="vq-score-col__head">
        <span class="vq-score-col__title">Heuristic · rule-based</span>
        <span id="heurinfo-${safe}" role="img" tabindex="0" aria-label="How the heuristic score is computed"
              style="display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;
                     border-radius:50%;border:1px solid var(--vq-border-dark);color:var(--vq-text-3);
                     font-size:10px;font-weight:700;cursor:help;flex:0 0 auto;margin-right:auto">?</span>
        <span class="vq-badge vq-badge--${esc(h.classification)}">${esc(h.classification)}</span>
      </div>
      <div class="vq-score-col__num">${h.vq_score}<small> / 100</small></div>
      ${comps}
      ${penalty}
      ${h.analysis ? `<div class="vq-score-col__analysis">${esc(h.analysis)}</div>` : ''}
    </div>`;
}

function _llmColumn(l) {
  const esc = VQ.esc;
  return `
    <div class="vq-score-col">
      <div class="vq-score-col__head">
        <span class="vq-score-col__title">AI · LLM analysis</span>
        <span class="vq-badge vq-badge--${esc(l.classification)}">${esc(l.classification)}</span>
      </div>
      <div class="vq-score-col__num">${l.vq_score}<small> / 100</small></div>
      ${l.analysis ? `<div class="vq-score-col__analysis">${esc(l.analysis)}</div>` : ''}
    </div>`;
}

function _scorePanel(seq, safe) {
  const h = seq.heuristic_output, l = seq.llm_output;
  if (!h && !l) return '';
  const cols = [];
  if (h) cols.push(_heurColumn(h, safe));   // heuristic-primary: rule-based leads
  if (l) cols.push(_llmColumn(l));
  return `
    <div class="vq-body-label">VQ scores</div>
    <div class="vq-score-panel">${cols.join('')}</div>`;
}

// ── Filtering ──────────────────────────────────────────────────────────────

function _hitPassesBlastFilter(hit, identMin, identMax, covMin, covMax) {
  const pident = hit.pct_identity ?? hit.pident;
  const qcov   = hit.query_coverage ?? hit.qcovhsp;
  if (identMin !== null && pident != null && pident < identMin) return false;
  if (identMax !== null && pident != null && pident > identMax) return false;
  if (covMin   !== null && qcov   != null && qcov   < covMin)   return false;
  if (covMax   !== null && qcov   != null && qcov   > covMax)   return false;
  return true;
}

/* Read a numeric filter input; '' / absent → null (i.e. unbounded). */
function _numVal(id) {
  const v = document.getElementById(id)?.value;
  return v != null && v !== '' ? parseFloat(v) : null;
}

/* Wire a checkbox-dropdown filter to its backing Set. */
function _wireMultiSelect(rootId, set) {
  const root = document.getElementById(rootId);
  if (!root) return;
  const toggle = root.querySelector('[data-menu-toggle]');
  const label  = root.querySelector('[data-ms-label]');
  toggle?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => { if (m !== root) m.classList.remove('open'); });
    root.classList.toggle('open');
  });
  // Keep the dropdown open while ticking several boxes (the global document
  // click handler closes any open .vq-menu).
  root.querySelector('.vq-ms__panel')?.addEventListener('click', e => e.stopPropagation());
  root.querySelectorAll('input[type="checkbox"]').forEach(cb => {
    cb.addEventListener('change', () => {
      if (cb.checked) set.add(cb.value); else set.delete(cb.value);
      const n = set.size;
      label.textContent = n === 0 ? label.dataset.allLabel
                        : n === 1 ? [...set][0]
                        : `${n} selected`;
      root.classList.toggle('vq-ms--active', n > 0);
      _applyFilters();
    });
  });
}

/* Recompute the count badge on the Filters toggle. With the panel collapsed
   this is the only cue that filters are silently narrowing the results. */
function _updateFilterBadge() {
  const badge = document.getElementById('viewer-filter-badge');
  if (!badge) return;
  const has = id => { const v = document.getElementById(id)?.value; return v != null && v !== ''; };
  const chk = id => document.getElementById(id)?.checked ?? false;
  let n = 0;
  if (has('viewer-len-min')  || has('viewer-len-max'))  n++;
  // Mirror _applyFilters: BLAST only filters when a source is checked AND a range is set.
  const blastRange = has('viewer-ident-min') || has('viewer-ident-max') ||
                     has('viewer-cov-min')   || has('viewer-cov-max');
  if ((chk('viewer-filter-blastn') || chk('viewer-filter-blastx')) && blastRange) n++;
  if (_VW.tax.family.size) n++;
  if (_VW.tax.phylum.size) n++;
  if (_VW.tax.genus.size)  n++;
  if (has('viewer-heur-min') || has('viewer-heur-max')) n++;
  if (has('viewer-llm-min')  || has('viewer-llm-max'))  n++;
  badge.textContent = String(n);
  badge.hidden = n === 0;
}

function _applyFilters() {
  const q       = document.getElementById('viewer-search')?.value.toLowerCase() || '';
  const famSet  = _VW.tax.family;
  const phySet  = _VW.tax.phylum;
  const genSet  = _VW.tax.genus;
  const heurMin = _numVal('viewer-heur-min');
  const heurMax = _numVal('viewer-heur-max');
  const llmMin  = _numVal('viewer-llm-min');
  const llmMax  = _numVal('viewer-llm-max');

  const identMinRaw = document.getElementById('viewer-ident-min')?.value;
  const identMaxRaw = document.getElementById('viewer-ident-max')?.value;
  const covMinRaw   = document.getElementById('viewer-cov-min')?.value;
  const covMaxRaw   = document.getElementById('viewer-cov-max')?.value;
  const identMin = identMinRaw !== '' ? parseFloat(identMinRaw) : null;
  const identMax = identMaxRaw !== '' ? parseFloat(identMaxRaw) : null;
  const covMin   = covMinRaw   !== '' ? parseFloat(covMinRaw)   : null;
  const covMax   = covMaxRaw   !== '' ? parseFloat(covMaxRaw)   : null;
  const filterBlastn = document.getElementById('viewer-filter-blastn')?.checked ?? false;
  const filterBlastx = document.getElementById('viewer-filter-blastx')?.checked ?? false;
  const hasBlastFilter = (filterBlastn || filterBlastx) &&
    (identMin !== null || identMax !== null || covMin !== null || covMax !== null);

  const lenMinRaw = document.getElementById('viewer-len-min')?.value;
  const lenMaxRaw = document.getElementById('viewer-len-max')?.value;
  const lenMin = lenMinRaw !== '' ? parseInt(lenMinRaw, 10) : null;
  const lenMax = lenMaxRaw !== '' ? parseInt(lenMaxRaw, 10) : null;

  _VW.filtered = _VW.sequences.filter(s => {
    const t = s.taxonomy || {};
    if (q) {
      const hay = [s.id, t.species, t.family, t.genus, t.phylum, t.order]
        .filter(Boolean).join(' ').toLowerCase();
      if (!hay.includes(q)) return false;
    }
    if (phySet.size && !phySet.has(t.phylum)) return false;
    if (famSet.size && !famSet.has(t.family)) return false;
    if (genSet.size && !genSet.has(t.genus))  return false;
    // Score ranges: a bound is set but the score is absent → exclude.
    if (heurMin !== null || heurMax !== null) {
      const hs = s.heuristic_output?.vq_score;
      if (hs == null) return false;
      if (heurMin !== null && hs < heurMin) return false;
      if (heurMax !== null && hs > heurMax) return false;
    }
    if (llmMin !== null || llmMax !== null) {
      const ls = s.llm_output?.vq_score;
      if (ls == null) return false;
      if (llmMin !== null && ls < llmMin) return false;
      if (llmMax !== null && ls > llmMax) return false;
    }
    if (lenMin !== null && (s.length ?? 0) < lenMin) return false;
    if (lenMax !== null && (s.length ?? 0) > lenMax) return false;
    if (hasBlastFilter) {
      if (filterBlastn) {
        const hits = s.blastn_hits || [];
        if (!hits.some(h => _hitPassesBlastFilter(h, identMin, identMax, covMin, covMax))) return false;
      }
      if (filterBlastx) {
        const hits = s.blastx_nr_hits || [];
        if (!hits.some(h => _hitPassesBlastFilter(h, identMin, identMax, covMin, covMax))) return false;
      }
    }
    return true;
  });

  const countEl = document.getElementById('viewer-count');
  if (countEl) countEl.textContent =
    `${_VW.filtered.length} of ${_VW.sequences.length} sequence${_VW.sequences.length !== 1 ? 's' : ''}`;

  _updateFilterBadge();
  _sortFiltered();
  _renderList();
}

// ── Sorting ──────────────────────────────────────────────────────────────────

function _sortFiltered() {
  const mode = document.getElementById('viewer-sort')?.value || 'score-desc';
  const score = s => _primaryScore(s)?.score ?? -1;
  const len   = s => s.length ?? 0;
  // Stable, deterministic ordering: ties on the primary key fall back to id.
  const byId  = (a, b) => a.id.localeCompare(b.id);
  const cmp = {
    'score-desc': (a, b) => score(b) - score(a) || byId(a, b),
    'score-asc':  (a, b) => score(a) - score(b) || byId(a, b),
    'len-desc':   (a, b) => len(b) - len(a)     || byId(a, b),
    'len-asc':    (a, b) => len(a) - len(b)     || byId(a, b),
    'id-asc':     byId,
  }[mode] || ((a, b) => score(b) - score(a) || byId(a, b));
  _VW.filtered.sort(cmp);
}

// ── List rendering ─────────────────────────────────────────────────────────

function _renderList() {
  const list = document.getElementById('viewer-list');
  if (!list) return;
  list.innerHTML = '';

  if (!_VW.filtered.length) {
    list.innerHTML = '<div class="vq-empty">No sequences match the current filters.</div>';
    return;
  }

  const frag = document.createDocumentFragment();
  _VW.filtered.forEach(seq => frag.appendChild(_seqCard(seq)));
  list.appendChild(frag);
}

/* Re-draw one card's low-complexity highlighting in place, without collapsing
   it. The genome SVG is only rebuilt when it already exists (a card that was
   never opened builds it lazily and picks up the current flag then).
   Because the export helpers read the live SVG out of the DOM, redrawing here
   is what makes this card's PNG/SVG/PDF exports honour its checkbox. */
function _refreshLowCx(card, seq) {
  if (!card || !seq) return;

  const wrap = card.querySelector('.vq-genome-wrap');
  if (wrap && wrap.querySelector('svg')) {
    wrap.innerHTML = '';
    wrap.appendChild(_genomeSVG(seq, wrap.clientWidth));
  }
  const fasta = card.querySelector('.vq-fasta');
  if (fasta) fasta.innerHTML = _fastaHighlightedHTML(seq);
}

// ── Sequence card ──────────────────────────────────────────────────────────

function _seqCard(seq) {
  const esc  = VQ.esc;
  const safe = VQ.safeId(seq.id);
  const card = document.createElement('div');
  card.className = 'vq-seq-card';
  card.id        = 'seq-card-' + safe;
  card.dataset.seqId = seq.id;
  card.style.marginBottom = '6px';

  const tax = seq.taxonomy;

  const orfCount     = (seq.orfs || []).length;
  const blastnCount  = (seq.blastn_hits || []).length;
  const refseqCount  = (seq.blastx_hits || []).length;
  const nrCount      = (seq.blastx_nr_hits || []).length;

  // Per-sequence low-complexity toggle — only offered when this sequence has
  // dust-masked regions of its own.
  const lowCxCount = ((seq.seq_quality || {}).low_complexity_regions || []).length;
  const hasLowCx   = lowCxCount > 0;

  // FASTA preview (header + sequence, seq-quality regions colour-highlighted).
  const hasFasta  = !!(seq.sequence || seq.sequence_nt);
  const fastaHtml = hasFasta ? _fastaHighlightedHTML(seq) : '';
  const topHitContent = _topHitContent(seq);   // middle block: best BLASTx hit + taxonomy
  // Column widths adapt: with no BLASTx hit (HMM-only) the middle block is absent,
  // so the FASTA pane grows and the sequence-quality pane shrinks.
  const hasTopHit = !!topHitContent;
  // Dot-plot column hugs its content (no grow) so there's no empty space to the
  // right of the legend; the freed width goes to the top-hit / FASTA panes.
  const dotCol    = 'flex:0 1 auto;min-width:0';
  const fastaCol  = hasTopHit ? 'flex:1.1 1 260px;min-width:250px' : 'flex:2 1 380px;min-width:320px';

  card.innerHTML = `
    <div class="vq-seq-card__head" role="button" tabindex="0" aria-expanded="false">
      <input type="checkbox" class="vq-checkbox" data-seq-checkbox
             ${_VW.selected.has(seq.id) ? 'checked' : ''}
             aria-label="Select ${esc(seq.id)}">
      <span class="vq-seq-card__id">${esc(seq.id)}</span>
      <div class="vq-seq-card__meta">
        ${tax?.family ? `
          <span class="vq-seq-card__meta-item">
            <span class="vq-seq-card__meta-label">Family</span>${esc(tax.family)}
          </span>` : ''}
        ${tax?.genus ? `
          <span class="vq-seq-card__meta-item">
            <span class="vq-seq-card__meta-label">Genus</span>${esc(tax.genus)}
          </span>` : ''}
        <span class="vq-seq-card__meta-item">
          ${(seq.length || 0).toLocaleString()} nt
        </span>
        <span class="vq-seq-card__meta-item">
          GC ${seq.gc_content != null ? seq.gc_content.toFixed(1) : '—'}%
        </span>
        <span class="vq-seq-card__meta-item">
          ${orfCount} ORF${orfCount !== 1 ? 's' : ''}
        </span>
        ${seq.cluster_id ? `<span class="vq-badge vq-badge--cluster">${esc(seq.cluster_id)}</span>` : ''}
        ${_scoreChips(seq)}
      </div>
      <div class="vq-section-actions" style="margin-left:auto;">
        ${_seqExportMenu(safe)}
      </div>
      <svg class="vq-seq-card__chevron" width="14" height="14" viewBox="0 0 24 24"
           fill="none" stroke="currentColor" stroke-width="2.5" aria-hidden="true">
        <polyline points="9 18 15 12 9 6"/>
      </svg>
    </div>

    <div class="vq-seq-card__body" id="body-${safe}">

      ${_scorePanel(seq, safe)}

      <div class="vq-body-label" style="display:flex;align-items:center;gap:7px">
        <span>Genome map</span>
        <span id="mapinfo-${safe}" role="img" tabindex="0" aria-label="What the genome map shows"
              style="display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;
                     border-radius:50%;border:1px solid var(--vq-border-dark);color:var(--vq-text-3);
                     font-size:10px;font-weight:700;cursor:help">?</span>
        ${hasLowCx ? `
        <label class="vq-filter-check vq-filter-check--inline" id="lowcx-wrap-${safe}"
               title="Amber shading over the ${lowCxCount} dustmasker-flagged low-complexity region${lowCxCount !== 1 ? 's' : ''} of this sequence, in the genome map and the FASTA preview. This sequence's exports follow the setting.">
          <input type="checkbox" id="lowcx-${safe}" ${_showLowCx(seq) ? 'checked' : ''}>
          Highlight low-complexity regions
        </label>` : ''}
      </div>
      <div class="vq-genome-wrap" id="genome-wrap-${safe}"></div>

      ${(seq.seq_quality || hasFasta || topHitContent) ? `
      <div style="display:flex;gap:14px;flex-wrap:wrap;align-items:flex-start;margin-top:8px">

        ${seq.seq_quality ? `
        <div style="${dotCol};display:flex;flex-direction:column">
          <div class="vq-body-label" style="display:flex;align-items:center;gap:7px;min-height:24px;margin:0 0 6px">
            <span>Sequence quality</span>
            <span id="siginfo-${safe}" role="img" tabindex="0" aria-label="What each signal means"
                  style="display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;
                         border-radius:50%;border:1px solid var(--vq-border-dark);color:var(--vq-text-3);
                         font-size:10px;font-weight:700;cursor:help">?</span>
          </div>
          <div class="vq-qual-block" style="height:360px;display:flex;flex-direction:row;align-items:center;
                                            justify-content:flex-start;gap:10px">
            <div class="vq-dotplot-wrap" id="dotplot-wrap-${safe}"
                 style="flex:0 1 auto;min-width:0;height:100%;aspect-ratio:1;max-width:260px;
                        display:flex;align-items:center;justify-content:center"></div>
            <div class="vq-seqqual-meta"
                 style="flex:0 0 auto;display:flex;flex-direction:column;gap:11px;font-size:11.5px;color:var(--vq-text-2)">
              ${_repeatSummary(seq.seq_quality)}
            </div>
          </div>
        </div>` : ''}

        ${topHitContent ? `
        <div style="flex:1 1 240px;min-width:230px;display:flex;flex-direction:column">
          <div class="vq-body-label" style="display:flex;align-items:center;gap:7px;min-height:24px;margin:0 0 6px">
            <span>Top hit &amp; taxonomy</span>
            <span id="taxinfo-${safe}" role="img" tabindex="0" aria-label="Taxonomy source"
                  style="display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;
                         border-radius:50%;border:1px solid var(--vq-border-dark);color:var(--vq-text-3);
                         font-size:10px;font-weight:700;cursor:help">?</span>
          </div>
          <div class="vq-qual-block" style="height:360px;overflow:auto;display:flex;flex-direction:column;
                                            justify-content:flex-start;gap:2px;font-size:11.5px;line-height:1.75">
            ${topHitContent}
          </div>
        </div>` : ''}

        ${hasFasta ? `
        <div style="${fastaCol};display:flex;flex-direction:column">
          <div class="vq-body-label" style="display:flex;align-items:center;justify-content:space-between;
                                            gap:10px;min-height:24px;margin:0 0 6px">
            <span style="display:inline-flex;align-items:center;gap:7px">
              <span>FASTA preview</span>
              <span id="fastainfo-${safe}" role="img" tabindex="0" aria-label="What the highlight colours mean"
                    style="display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;
                           border-radius:50%;border:1px solid var(--vq-border-dark);color:var(--vq-text-3);
                           font-size:10px;font-weight:700;cursor:help">?</span>
            </span>
            <button type="button" class="vq-btn vq-btn--ghost vq-btn--sm" id="fasta-copy-${safe}">Copy FASTA</button>
          </div>
          <div class="vq-fasta" style="height:360px;max-height:360px;margin:0">${fastaHtml}</div>
        </div>` : ''}

      </div>` : ''}

      <div class="vq-body-label">BLAST hits</div>
      <div class="vq-panel">
        <div class="vq-subtabs" role="tablist" aria-label="BLAST hit sources">
          <button class="vq-subtab active" type="button" data-blast="blastn" role="tab">
            BLASTn<span class="vq-subtab__count">${blastnCount}</span>
          </button>
          <button class="vq-subtab" type="button" data-blast="refseq" role="tab">
            BLASTx · RefSeq<span class="vq-subtab__count">${refseqCount}</span>
          </button>
          <button class="vq-subtab" type="button" data-blast="nr" role="tab">
            BLASTx · NR<span class="vq-subtab__count">${nrCount}</span>
          </button>
        </div>
        <div class="vq-panel__body" id="blast-body-${safe}"></div>
      </div>

    </div>`;

  // Info (?) hovers — bind a tooltip to a badge by id. `src` may be a function
  // so a tooltip whose content depends on viewer state (the low-complexity
  // toggle) is built at hover time rather than frozen at card creation.
  const _bindInfo = (id, src) => {
    const el = card.querySelector('#' + id + '-' + safe);
    if (!el) return;
    const html = () => (typeof src === 'function' ? src() : src);
    el.addEventListener('mousemove', e => VQ.tooltipShow(html(), e));
    el.addEventListener('mouseleave', VQ.tooltipHide);
    el.addEventListener('focus', () => {
      const r = el.getBoundingClientRect();
      VQ.tooltipShow(html(), { clientX: r.right, clientY: r.bottom });
    });
    el.addEventListener('blur', VQ.tooltipHide);
  };
  _bindInfo('mapinfo',   () => _mapInfoHTML(seq));
  _bindInfo('fastainfo', () => _fastaInfoHTML(seq));
  _bindInfo('heurinfo',  _heurInfoHTML());
  _bindInfo('siginfo', _signalsInfoHTML());
  _bindInfo('taxinfo', `
    <div class="vq-tooltip__title">Top hit &amp; taxonomy</div>
    <div style="max-width:240px">Best BLASTx hit (NR preferred, else RefSeq).
    Taxonomic lineage resolved from <strong>ICTV</strong> and <strong>NCBI</strong> taxonomy.</div>`);

  // Per-sequence low-complexity highlight toggle
  card.querySelector('#lowcx-' + safe)?.addEventListener('change', e => {
    if (e.target.checked) _VW.lowCxOff.delete(seq.id);
    else                  _VW.lowCxOff.add(seq.id);
    _refreshLowCx(card, seq);
  });

  // Copy FASTA (header + sequence) to clipboard
  card.querySelector('#fasta-copy-' + safe)?.addEventListener('click', function () {
    navigator.clipboard?.writeText(VQ.buildFasta([seq])).then(() => {
      this.textContent = 'Copied!';
      setTimeout(() => { this.textContent = 'Copy FASTA'; }, 1500);
    }).catch(() => {
      this.textContent = 'Copy failed';
      setTimeout(() => { this.textContent = 'Copy FASTA'; }, 1500);
    });
  });

  // Checkbox
  card.querySelector('[data-seq-checkbox]')?.addEventListener('change', e => {
    if (e.target.checked) _VW.selected.add(seq.id);
    else                  _VW.selected.delete(seq.id);
    _updateBatchBtn();
  });

  // Export menu wiring
  _wireExportMenu(card, 'seq-menu-' + safe, () => card.querySelector('.vq-genome-wrap svg'), seq);

  // BLAST subtabs (lazy-render hit tables)
  const subtabs   = card.querySelectorAll('.vq-subtab');
  const blastBody = card.querySelector('#blast-body-' + safe);
  const renderBlast = (kind) => {
    let hits = [];
    if      (kind === 'blastn') hits = seq.blastn_hits        || [];
    else if (kind === 'refseq') hits = seq.blastx_hits || [];
    else                        hits = seq.blastx_nr_hits     || [];
    blastBody.innerHTML = _blastTable(hits, kind);
  };
  subtabs.forEach(t => t.addEventListener('click', () => {
    subtabs.forEach(s => s.classList.toggle('active', s === t));
    renderBlast(t.dataset.blast);
  }));

  // Expand / collapse
  const head = card.querySelector('.vq-seq-card__head');
  head.addEventListener('click', e => {
    if (e.target.closest('button') || e.target.closest('input') ||
        e.target.closest('.vq-menu')) return;
    _toggleCard(card, seq);
  });
  head.addEventListener('keydown', e => {
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); _toggleCard(card, seq); }
  });

  return card;
}

function _openCard(card) {
  if (card.classList.contains('open')) return;
  const id  = card.dataset.seqId;
  const seq = _VW.sequences.find(s => s.id === id);
  if (seq) _toggleCard(card, seq, /*forceOpen*/ true);
}

function _toggleCard(card, seq, forceOpen = false) {
  const isOpen = forceOpen ? true : !card.classList.contains('open');
  card.classList.toggle('open', isOpen);
  card.querySelector('.vq-seq-card__head')?.setAttribute('aria-expanded', String(isOpen));
  if (!isOpen) return;

  const safe = VQ.safeId(seq.id);
  const wrap = card.querySelector('#genome-wrap-' + safe);
  if (wrap && !wrap.querySelector('svg')) wrap.appendChild(_genomeSVG(seq, wrap.clientWidth));

  // Native-SVG self-similarity dot plot (rendered lazily on first open).
  if (seq.seq_quality) {
    const dpWrap = card.querySelector('#dotplot-wrap-' + safe);
    if (dpWrap && !dpWrap.querySelector('svg')) {
      const node = _dotPlotSVG(seq.seq_quality);
      if (node) dpWrap.appendChild(node);
      else dpWrap.innerHTML = '<div class="vq-empty" style="padding:8px 0;font-size:11px">No dot plot.</div>';
    }
  }

  // Render the default BLASTn table
  const blastBody = card.querySelector('#blast-body-' + safe);
  if (blastBody && !blastBody.innerHTML) {
    blastBody.innerHTML = _blastTable(seq.blastn_hits || [], 'blastn');
  }
}

// ── BLAST hits table ────────────────────────────────────────────────────────

function _blastTable(hits, kind) {
  if (!hits.length) {
    return `<div class="vq-empty" style="padding:24px 16px">No ${kind} hits.</div>`;
  }
  const esc = VQ.esc;
  const rows = hits.map(h => {
    const accession = h.accession || h.subject_id || '—';
    const title     = h.stitle || h.subject_title || '—';
    const pident    = h.pct_identity ?? h.pident;
    const qcov      = h.query_coverage ?? h.qcovhsp;
    const evalue    = h.e_value ?? h.evalue;
    const bitscore  = h.bit_score ?? h.bitscore;
    const qrange    = (h.query_start != null && h.query_end != null)
      ? `${h.query_start}–${h.query_end}`
      : '—';
    return `
    <tr>
      <td class="vq-td--mono">${esc(accession)}</td>
      <td>${esc(title)}</td>
      <td class="vq-td--num">${pident != null ? pident.toFixed(1) + '%' : '—'}</td>
      <td class="vq-td--num">${qcov != null ? (+qcov).toFixed(1) + '%' : '—'}</td>
      <td class="vq-td--num">${evalue != null ? evalue.toExponential(1) : '—'}</td>
      <td class="vq-td--num">${bitscore ?? '—'}</td>
      <td class="vq-td--num">${qrange}</td>
    </tr>`;
  }).join('');
  return `
    <table class="vq-table">
      <thead>
        <tr>
          <th>Accession</th><th>Title</th><th>Identity</th>
          <th>Coverage</th><th>E-value</th><th>Bit score</th><th>Query range</th>
        </tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>`;
}

// ── Per-card export menu ───────────────────────────────────────────────────

function _seqExportMenu(safeId) {
  return `
    <div class="vq-menu" id="seq-menu-${safeId}">
      <button class="vq-btn vq-btn--sm vq-btn--ghost" data-menu-toggle type="button"
              aria-label="Export options">
        Export
        <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2.5">
          <polyline points="6 9 12 15 18 9"/>
        </svg>
      </button>
      <div class="vq-menu__panel" role="menu">
        <button class="vq-menu__item" data-fmt="png" role="menuitem" type="button">PNG image</button>
        <button class="vq-menu__item" data-fmt="svg" role="menuitem" type="button">SVG vector</button>
        <button class="vq-menu__item" data-fmt="pdf" role="menuitem" type="button">PDF (print)</button>
        <button class="vq-menu__item" data-fmt="fasta" role="menuitem" type="button">FASTA sequence</button>
      </div>
    </div>`;
}

function _wireExportMenu(card, menuId, getSvg, seq) {
  const menu = card.querySelector('#' + menuId);
  if (!menu) return;
  const toggle = menu.querySelector('[data-menu-toggle]');
  toggle?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => {
      if (m !== menu) m.classList.remove('open');
    });
    menu.classList.toggle('open');
  });
  menu.querySelectorAll('[data-fmt]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      menu.classList.remove('open');
      const fmt = item.dataset.fmt;
      if (fmt === 'fasta') {
        VQ.downloadText(VQ.buildFasta([seq]), `${seq.id}.fasta`);
        return;
      }
      // Make sure the card is open & the SVG rendered
      _openCard(card);
      const svg = getSvg();
      if (!svg) return;
      if (fmt === 'png') VQ.exportPNG(svg, seq.id + '.png');
      if (fmt === 'svg') VQ.exportSVG(svg, seq.id + '.svg');
      if (fmt === 'pdf') VQ.exportPDF(svg, seq.id);
    });
  });
}

// ── Batch ──────────────────────────────────────────────────────────────────

function _updateBatchBtn() {
  const n     = _VW.selected.size;
  const btn   = document.getElementById('viewer-batch-btn');
  const label = document.getElementById('viewer-batch-label');
  if (!btn) return;
  btn.disabled = n === 0;
  if (label) label.textContent = n ? `Export ${n} selected` : 'Export selected';
}

function _ensureSvgForSeq(seqId) {
  const safe = VQ.safeId(seqId);
  const card = document.getElementById('seq-card-' + safe);
  if (!card) return null;
  _openCard(card);
  return card.querySelector('.vq-genome-wrap svg');
}

function _selectedSeqs() {
  return [..._VW.selected].map(id => _VW.sequences.find(s => s.id === id)).filter(Boolean);
}

async function _runBatch(kind) {
  const seqs = _selectedSeqs();
  if (!seqs.length) return;

  if (kind === 'fasta') {
    VQ.downloadText(VQ.buildFasta(seqs), `viralquest_selected_${seqs.length}.fasta`);
    return;
  }

  // Make sure every selected sequence has its SVG built
  const items = seqs.map(s => ({
    svgEl: _ensureSvgForSeq(s.id),
    label: `${s.id}  ·  ${s.taxonomy?.species || ''}  ·  ${(s.length||0).toLocaleString()} nt`,
    filename: s.id,
  })).filter(it => it.svgEl);

  if (!items.length) return;

  if (kind === 'png') {
    VQ.exportCompositePNG(items, `viralquest_${items.length}_sequences.png`);
  } else if (kind === 'svg') {
    VQ.exportCompositeSVG(items, `viralquest_${items.length}_sequences.svg`);
  } else if (kind === 'png-each') {
    await VQ.exportBatchPNG(items.map(it => ({ svgEl: it.svgEl, filename: it.filename + '.png' })), 350);
  } else if (kind === 'svg-each') {
    for (const it of items) { VQ.exportSVG(it.svgEl, it.filename + '.svg'); await _sleep(180); }
  }
}

const _sleep = ms => new Promise(r => setTimeout(r, ms));

// ── Genome SVG ────────────────────────────────────────────────────────────────

function _genomeSVG(seq, containerWidth) {
  const orfs   = seq.orfs || [];
  const seqLen = seq.length || 1;

  const PAD_L = 50;   // wider left gutter for frame labels
  const PAD_R = 24;

  // Long genomes get more horizontal pixels so the bars stay readable and
  // the wrap scrolls horizontally (CSS: overflow-x:auto on .vq-genome-wrap).
  const SCROLL_THRESHOLD = 20000;
  const containerW = Math.max(containerWidth || 0, 720);
  const fluidMode  = seqLen <= SCROLL_THRESHOLD;
  const W = fluidMode
    ? containerW
    : Math.max(containerW, Math.round(seqLen / 10));   // ~10 nt per pixel

  const drawW = W - PAD_L - PAD_R;
  const SCALE = d3.scaleLinear([0, seqLen], [0, drawW]);

  // ── Read-coverage band (only when --reads produced a profile) ───────────
  const cov        = (seq.coverage && (seq.coverage.bins || []).length) ? seq.coverage : null;
  // Strand-split coverage needs room for a two-sided (sense up / antisense down) plot.
  const covStrand  = cov && (cov.bins_fwd || []).length === (cov.bins || []).length
                          && (cov.bins_rev || []).length === (cov.bins || []).length;
  const COV_AREA_H = covStrand ? 56 : 36;         // height of the coverage area
  // Room reserved = legend block row + bars + the caption//gap row below.
  const COV_H      = cov ? COV_LEG_H + COV_AREA_H + 18 : 0;

  // ── Frame-based lane assignment ─────────────────────────────────────────
  const FRAME_INDEX = { '+1': 0, '+2': 1, '+3': 2, '-1': 3, '-2': 4, '-3': 5 };
  const FRAME_LABEL = ['+1','+2','+3','-1','-2','-3'];

  function frameLane(orf) {
    const f = orf.frame;
    // frame may be a number (2, -2) or already a string ('+2', '-2')
    const key = typeof f === 'number'
      ? (f > 0 ? '+' : '') + f
      : String(f ?? '').replace(/\s+/g, '');
    return FRAME_INDEX[key] ?? (orf.strand === '-' ? 3 : 0);
  }
  const nLanes = 6;

  // Pack domains into sub-lanes per ORF (deduplicated: best score per target)
  const orfsWithFrame = orfs.map(o => ({
    ...o,
    frameLane:    frameLane(o),
    _domainLanes: _assignDomainLanes(_bestDomainPerTarget(o.domains || [])),
  }));
  orfsWithFrame.forEach(o => {
    o._nDomLanes = Math.max(0, ...o._domainLanes.map(d => d.lane + 1));
  });

  const AXIS_Y = 28;
  const ORF_H  = 16;
  const DOM_H  = 7;
  const DOM_G  = 1.5;
  const LANE_G = 12;
  const STRAND_GAP = 14;   // gap between + and - frame groups

  // Each lane: enough vertical room for the ORF bar + its domain stack.
  const laneDomMax = new Array(nLanes).fill(0);
  orfsWithFrame.forEach(o => {
    if (o._nDomLanes > laneDomMax[o.frameLane]) laneDomMax[o.frameLane] = o._nDomLanes;
  });
  const laneHeights = laneDomMax.map(n => ORF_H + n * (DOM_H + DOM_G) + LANE_G);
  const laneYs      = [];
  let cy = AXIS_Y + 14 + COV_H;
  laneHeights.forEach((h, i) => {
    laneYs.push(cy);
    cy += h;
    if (i === 2) cy += STRAND_GAP;   // separator between + and - strand frames
  });

  const totalH = cy + 36;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg' + (fluidMode ? ' vq-genome-svg--fluid' : ''))
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('preserveAspectRatio', 'xMinYMin meet');
  if (fluidMode) {
    svg.style('width', '100%').style('height', 'auto');
  } else {
    svg.attr('width', W).attr('height', totalH).style('min-width', W + 'px');
  }

  // ── Axis ────────────────────────────────────────────────────────────────
  // Very long genomes scroll horizontally, so a fixed ~14 labels would leave
  // the reader thousands of nt from the nearest coordinate. Past the threshold
  // the tick count follows the drawn width (a label roughly every LONG_TICK_PX)
  // and switches to thousands separators so the bigger numbers stay legible.
  const LONG_SEQ_NT  = 50000;
  const LONG_TICK_PX = 120;   // d3 snaps to nice steps → a label every ~100 px
  const isLong = seqLen > LONG_SEQ_NT;
  const nTicks = isLong
    ? Math.max(8, Math.round(drawW / LONG_TICK_PX))
    : Math.min(14, Math.max(4, Math.round(W / 120)));

  const gridH = totalH - AXIS_Y - 26;
  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, ${AXIS_Y})`)
    .call(d3.axisTop(SCALE).ticks(nTicks).tickSize(0)
            .tickFormat(isLong ? d3.format(',') : null));
  axisG.selectAll('.tick line')
    .attr('y1', 0).attr('y2', gridH)
    .attr('stroke', 'var(--vq-border)').attr('stroke-dasharray', '3,3');
  axisG.select('.domain').attr('stroke', 'var(--vq-border-dark)');
  axisG.selectAll('.tick text')
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)');

  // ── Read-coverage track ─────────────────────────────────────────────────
  if (cov) _drawCoverageTrack(svg, cov, SCALE, seqLen, PAD_L, drawW, AXIS_Y, COV_AREA_H);

  // ── Low-complexity bands (dustmasker) — behind the ORF lanes ─────────────
  // Suppressed when this card's "Highlight low-complexity regions" checkbox is
  // off; exports take the SVG as drawn, so the bands are omitted there too.
  const sqBands = _showLowCx(seq)
    && seq.seq_quality && (seq.seq_quality.low_complexity_regions || []).length
    ? seq.seq_quality.low_complexity_regions : null;
  if (sqBands) _drawLowComplexityBands(svg, sqBands, SCALE, PAD_L, AXIS_Y, AXIS_Y + gridH);

  // ── Frame labels (left gutter) + lane backgrounds ───────────────────────
  laneYs.forEach((y, i) => {
    // Subtle striping so the lanes are visually distinct
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y - 2)
      .attr('width', drawW)
      .attr('height', ORF_H + laneDomMax[i] * (DOM_H + DOM_G) + 4)
      .attr('fill', i % 2 ? 'var(--vq-surface-2)' : 'transparent')
      .attr('opacity', 0.5);

    svg.append('text')
      .attr('x', PAD_L - 8)
      .attr('y', y + ORF_H / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 10)
      .attr('font-family', 'var(--vq-font-mono)')
      .attr('fill', i < 3 ? 'var(--vq-accent-dark)' : 'var(--vq-warning)')
      .text(FRAME_LABEL[i]);
  });

  // Backbone for + strand and - strand (between frame groups)
  const plusBackboneY  = laneYs[2] + ORF_H + laneDomMax[2] * (DOM_H + DOM_G) + 6;
  svg.append('line')
    .attr('x1', PAD_L).attr('x2', PAD_L + drawW)
    .attr('y1', plusBackboneY).attr('y2', plusBackboneY)
    .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1.5)
    .attr('stroke-dasharray', '2,3');

  // ── ORFs ────────────────────────────────────────────────────────────────
  orfsWithFrame.forEach(orf => {
    const x1  = PAD_L + SCALE(orf.start_position);
    const x2  = PAD_L + SCALE(orf.stop_position);
    const bw  = Math.max(x2 - x1, 4);
    const y   = laneYs[orf.frameLane];
    const col = orf.strand === '-' ? 'var(--vq-warning)' : 'var(--vq-accent)';

    const g = svg.append('g');

    const aa        = orf.aa_sequence || '';
    const copyHint  = aa ? `<div style="margin-top:5px;padding-top:5px;
        border-top:1px solid rgba(255,255,255,.12);font-size:10px;
        color:rgba(255,255,255,.55)">Click to copy AA sequence</div>` : '';

    g.append('polygon')
      .attr('points', _orfArrowPoints(x1, x1 + bw, y, ORF_H, orf.strand))
      .attr('fill', col).attr('opacity', .88)
      .style('cursor', aa ? 'pointer' : 'default')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(orf.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Strand</span><span>${VQ.esc(orf.strand)}</span>
          <span class="vq-tooltip__key">Frame</span><span>${VQ.esc(orf.frame)}</span>
          <span class="vq-tooltip__key">Length</span><span>${orf.length_aa} aa</span>
          <span class="vq-tooltip__key">Position</span>
          <span>${orf.start_position}–${orf.stop_position}</span>
          <span class="vq-tooltip__key">Type</span><span>${VQ.esc(orf.orf_type)}</span>
          <span class="vq-tooltip__key">Domains</span><span>${(orf._domainLanes||[]).length}</span>
        </div>
        ${copyHint}`, evt))
      .on('mouseleave', VQ.tooltipHide)
      .on('click', evt => {
        if (!aa) return;
        navigator.clipboard.writeText(aa).then(() => {
          VQ.tooltipShow(`<div class="vq-tooltip__title">✓ AA sequence copied</div>
            <div style="font-size:10px;color:rgba(255,255,255,.6);margin-top:3px">${orf.length_aa} aa</div>`, evt);
          setTimeout(VQ.tooltipHide, 1500);
        }).catch(() => {
          VQ.tooltipShow(`<div class="vq-tooltip__title">⚠ Copy not available</div>`, evt);
          setTimeout(VQ.tooltipHide, 1500);
        });
      });

    // Show label only when the full name fits in the body of the arrow (excluding tip)
    const tip        = Math.min(ORF_H * 0.75, bw * 0.25);
    const labelW     = bw - tip;
    const estimatedW = orf.name.length * 6.2 + 10;
    if (labelW > estimatedW) {
      const labelX = orf.strand === '+'
        ? x1 + (bw - tip) / 2
        : x1 + tip + (bw - tip) / 2;
      g.append('text')
        .attr('x', labelX).attr('y', y + ORF_H / 2 + 3.5)
        .attr('text-anchor', 'middle')
        .attr('font-size', 9.5)
        .attr('fill', 'rgba(255,255,255,.92)')
        .attr('pointer-events', 'none')
        .text(orf.name);
    }

    // HMM domains — lane-packed below the ORF bar within the same frame lane
    (orf._domainLanes || []).forEach(d => {
      const domStart = orf.start_position + d.dom.start * 3;
      const domEnd   = orf.start_position + d.dom.stop  * 3;
      const dx1      = PAD_L + SCALE(Math.min(domStart, domEnd));
      const dx2      = PAD_L + SCALE(Math.max(domStart, domEnd));
      const dw       = Math.max(dx2 - dx1, 3);
      const domCol   = _domainColor(d.dom.target);
      const dy       = y + ORF_H + d.lane * (DOM_H + DOM_G) + 1;

      g.append('rect')
        .attr('x', dx1).attr('width', dw)
        .attr('y', dy).attr('height', DOM_H)
        .attr('rx', 2)
        .attr('fill', domCol).attr('opacity', .9)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(d.dom.target)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Database</span><span>${VQ.esc(d.dom.database)}</span>
            <span class="vq-tooltip__key">Score</span><span>${d.dom.score}</span>
            <span class="vq-tooltip__key">E-value</span><span>${d.dom.e_value?.toExponential(2) ?? '—'}</span>
            <span class="vq-tooltip__key">aa range</span><span>${d.dom.start}–${d.dom.stop}</span>
          </div>
          ${d.dom.description ? `<div style="margin-top:5px;padding-top:5px;
              border-top:1px solid rgba(255,255,255,.12);font-size:10px;
              color:rgba(255,255,255,.72)">${VQ.esc(d.dom.description)}</div>` : ''}`,
          evt))
        .on('mouseleave', VQ.tooltipHide);
    });
  });

  // ── Legend ──────────────────────────────────────────────────────────────
  const legendG = svg.append('g').attr('transform', `translate(${PAD_L}, ${totalH - 18})`);
  legendG.append('polygon')
    .attr('points', _orfArrowPoints(0, 18, 0, 9, '+'))
    .attr('fill', 'var(--vq-accent)');
  legendG.append('text').attr('x', 24).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('(+) ORF');
  legendG.append('polygon')
    .attr('points', _orfArrowPoints(70, 88, 0, 9, '-'))
    .attr('fill', 'var(--vq-warning)');
  legendG.append('text').attr('x', 94).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('(−) ORF');
  legendG.append('rect').attr('x', 140).attr('width', 12).attr('height', 5).attr('rx', 1)
    .attr('fill', 'var(--vq-dom-1)').attr('y', 2);
  legendG.append('text').attr('x', 158).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('HMM domain');
  // Coverage is NOT listed here — it has its own legend block at the top of the
  // coverage band (_drawCovLegend), next to the data it describes.
  const legendX = 240;
  if (!fluidMode) {
    legendG.append('text').attr('x', legendX).attr('y', 8.5)
      .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)')
      .text(`↔ ${seqLen.toLocaleString()} nt · scroll horizontally`);
  }

  return svg.node();
}

// ── Read-coverage track ─────────────────────────────────────────────────────

// Colour a coverage bin by its mean read base quality (Phred). Green = high
// confidence, amber = borderline, red = low quality (possible sequencing error).
function _qualColor(q) {
  if (q >= 30) return 'var(--vq-success)';
  if (q >= 20) return 'var(--vq-warning)';
  if (q > 0)   return 'var(--vq-danger)';
  return 'var(--vq-text-3)';
}

// Height of the legend row reserved above the coverage bars. Shared by
// _genomeSVG (which budgets the vertical space) and _drawCoverageTrack (which
// draws into it) — the two must agree or the bars drift out of their band.
const COV_LEG_H = 17;

// SVG has no text metrics before layout, so the legend advances by an estimate
// of each label's width at 9.5px sans-serif.
const _legendTextW = s => s.length * 5.35 + 2;

/* Legend block for the coverage band — drawn in its own row above the bars, so
   the swatches sit next to the data they describe instead of in the shared
   bottom legend. Each swatch carries its own "Q…" label: a single free-floating
   "Q" reads as if it labelled the first colour. Returns the block width. */
function _drawCovLegend(g, x, y, qual) {
  const PADX = 7, SW = 10, SW_GAP = 4, ITEM_GAP = 11, H = 14;

  const items = [{ label: 'Read coverage', lead: true }];
  if (qual) {
    items.push({ sw: 'var(--vq-success)', label: 'Q≥30'   });
    items.push({ sw: 'var(--vq-warning)', label: 'Q20–29' });
    items.push({ sw: 'var(--vq-danger)',  label: 'Q<20'   });
  } else {
    items.push({ sw: 'var(--vq-accent)', label: 'depth' });
  }

  // Lay out left→right, recording each item's x before drawing anything.
  let w = PADX;
  const xs = items.map(it => {
    const at = w;
    w += (it.sw ? SW + SW_GAP : 0) + _legendTextW(it.label) + ITEM_GAP;
    return at;
  });
  w = w - ITEM_GAP + PADX;

  const blk = g.append('g').attr('class', 'vq-cov-legend')
    .attr('transform', `translate(${x}, ${y})`);
  blk.append('rect')
    .attr('width', w).attr('height', H).attr('rx', 4)
    .attr('fill', 'var(--vq-surface-2)').attr('stroke', 'var(--vq-border)');

  items.forEach((it, i) => {
    let cx = xs[i];
    if (it.sw) {
      blk.append('rect')
        .attr('x', cx).attr('y', H / 2 - 2.5).attr('width', SW).attr('height', 5)
        .attr('rx', 1).attr('fill', it.sw).attr('opacity', it.label === 'depth' ? 0.35 : 0.9);
      cx += SW + SW_GAP;
    }
    blk.append('text').attr('x', cx).attr('y', H / 2 + 3.2)
      .style('font-size', '9.5px')
      .style('font-weight', it.lead ? 600 : 400)
      .style('fill', it.lead ? 'var(--vq-text-2)' : 'var(--vq-text-3)')
      .text(it.label);
  });
  return w;
}

function _drawCoverageTrack(svg, cov, SCALE, seqLen, PAD_L, drawW, AXIS_Y, areaH) {
  const legY   = AXIS_Y + 12;             // legend block row (reserved: COV_LEG_H)
  const top    = legY + COV_LEG_H;        // bars start below it
  const bottom = top + areaH;
  const bins   = cov.bins;
  const nb     = bins.length;
  const maxD   = cov.max_depth > 0 ? cov.max_depth : 1;

  // Per-base quality (from the same BAM) aligned 1:1 with the depth bins.
  const qual   = (cov.quality_bins && cov.quality_bins.length === nb) ? cov.quality_bins : null;
  // Strand-split depth (sense = forward, antisense = reverse), same binning.
  const fwd    = (cov.bins_fwd && cov.bins_fwd.length === nb) ? cov.bins_fwd : null;
  const rev    = (cov.bins_rev && cov.bins_rev.length === nb) ? cov.bins_rev : null;
  const strand = fwd && rev;

  const binNt  = i => ((i + 0.5) / nb) * seqLen;       // bin centre in nt
  const xOf    = i => PAD_L + SCALE(binNt(i));
  const g      = svg.append('g').attr('class', 'vq-cov-track');
  const barCol = i => qual ? _qualColor(qual[i]) : 'var(--vq-accent)';
  const bw     = Math.max(0.6, drawW / nb);

  if (strand) {
    // Two-sided track: sense reads grow up, antisense grow down from a central
    // baseline. Both bars are coloured by base quality (same column colour).
    const mid  = top + areaH / 2;
    const maxS = Math.max(1, d3.max(fwd) || 0, d3.max(rev) || 0);
    const yUp  = d3.scaleLinear([0, maxS], [mid, top]);      // forward → up
    const yDn  = d3.scaleLinear([0, maxS], [mid, bottom]);   // reverse → down

    g.selectAll('rect.vq-cov-fwd').data(fwd).join('rect')
      .attr('x', (d, i) => xOf(i) - bw / 2).attr('y', d => yUp(d))
      .attr('width', bw).attr('height', d => mid - yUp(d))
      .attr('fill', (d, i) => barCol(i)).attr('opacity', 0.85);
    g.selectAll('rect.vq-cov-rev').data(rev).join('rect')
      .attr('x', (d, i) => xOf(i) - bw / 2).attr('y', mid)
      .attr('width', bw).attr('height', d => yDn(d) - mid)
      .attr('fill', (d, i) => barCol(i)).attr('opacity', 0.5);

    g.append('line').attr('x1', PAD_L).attr('x2', PAD_L + drawW)
      .attr('y1', mid).attr('y2', mid)
      .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1);

    // Gutter: strand markers + scale.
    g.append('text').attr('x', PAD_L - 8).attr('y', top + 7).attr('text-anchor', 'end')
      .attr('font-size', 8).attr('font-weight', 700).attr('fill', 'var(--vq-text-2)').text('＋');
    g.append('text').attr('x', PAD_L - 8).attr('y', mid + 3).attr('text-anchor', 'end')
      .attr('font-size', 9).attr('font-family', 'var(--vq-font-mono)').attr('font-weight', 700)
      .attr('fill', 'var(--vq-text-2)').text('Cov');
    g.append('text').attr('x', PAD_L - 8).attr('y', bottom - 1).attr('text-anchor', 'end')
      .attr('font-size', 8).attr('font-weight', 700).attr('fill', 'var(--vq-text-2)').text('－');
  } else {
    const yScale = d3.scaleLinear([0, maxD], [bottom, top]);
    if (qual) {
      g.selectAll('rect.vq-cov-bar').data(bins).join('rect')
        .attr('class', 'vq-cov-bar')
        .attr('x', (d, i) => xOf(i) - bw / 2).attr('y', d => yScale(d))
        .attr('width', bw).attr('height', d => bottom - yScale(d))
        .attr('fill', (d, i) => _qualColor(qual[i])).attr('opacity', 0.75);
    } else {
      const area = d3.area().x((d, i) => xOf(i)).y0(bottom).y1(d => yScale(d)).curve(d3.curveMonotoneX);
      g.append('path').datum(bins).attr('d', area)
        .attr('fill', 'var(--vq-accent)').attr('opacity', 0.30);
    }
    const line = d3.line().x((d, i) => xOf(i)).y(d => yScale(d)).curve(d3.curveMonotoneX);
    g.append('path').datum(bins).attr('d', line).attr('fill', 'none')
      .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1).attr('opacity', 0.85);
    g.append('line').attr('x1', PAD_L).attr('x2', PAD_L + drawW)
      .attr('y1', bottom).attr('y2', bottom)
      .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1);
    g.append('text').attr('x', PAD_L - 8).attr('y', top + areaH / 2 - 3).attr('text-anchor', 'end')
      .attr('font-size', 9).attr('font-family', 'var(--vq-font-mono)').attr('font-weight', 700)
      .attr('fill', 'var(--vq-text-2)').text('Cov');
    g.append('text').attr('x', PAD_L - 8).attr('y', top + areaH / 2 + 8).attr('text-anchor', 'end')
      .attr('font-size', 8).attr('fill', 'var(--vq-text-3)').text(`${Math.round(maxD)}×`);
  }

  // Legend block (left) and summary caption (right) share the reserved row.
  _drawCovLegend(g, PAD_L, legY, qual);
  const qCap = qual ? ` · Q̄ ${(cov.mean_quality ?? 0).toFixed(0)}` : '';
  const cap  = strand
    ? `sense ${(cov.mean_depth_fwd ?? 0).toFixed(1)}× · anti ${(cov.mean_depth_rev ?? 0).toFixed(1)}× · CV ${cov.cv.toFixed(2)}${qCap}`
    : `mean ${cov.mean_depth.toFixed(1)}× · breadth ${(cov.breadth_1x * 100).toFixed(0)}% · CV ${cov.cv.toFixed(2)}${qCap}`;
  g.append('text')
    .attr('x', PAD_L + drawW).attr('y', legY + 10)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--vq-text-3)').text(cap);

  // Hover overlay → depth tooltip at the pointer position.
  g.append('rect')
    .attr('x', PAD_L).attr('y', top - 2)
    .attr('width', drawW).attr('height', areaH + 4)
    .attr('fill', 'transparent').style('cursor', 'crosshair')
    .on('mousemove', evt => {
      const [mx] = d3.pointer(evt);
      const nt   = Math.max(0, Math.min(seqLen, SCALE.invert(mx - PAD_L)));
      const bi   = Math.max(0, Math.min(nb - 1, Math.floor((nt / seqLen) * nb)));
      const qRow = qual
        ? `<span class="vq-tooltip__key">Quality</span><span>Q${(qual[bi] ?? 0).toFixed(0)}</span>`
        : '';
      const strandRows = strand
        ? `<span class="vq-tooltip__key">Sense (＋)</span><span>${fwd[bi].toFixed(1)}×</span>
           <span class="vq-tooltip__key">Antisense (－)</span><span>${rev[bi].toFixed(1)}×</span>`
        : '';
      VQ.tooltipShow(`
        <div class="vq-tooltip__title">Read coverage</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Position</span><span>~${Math.round(nt).toLocaleString()} nt</span>
          <span class="vq-tooltip__key">Depth</span><span>${bins[bi].toFixed(1)}×</span>
          ${strandRows}
          ${qRow}
          <span class="vq-tooltip__key">Mean</span><span>${cov.mean_depth.toFixed(1)}×</span>
          <span class="vq-tooltip__key">Max</span><span>${cov.max_depth.toFixed(0)}×</span>
        </div>`, evt);
    })
    .on('mouseleave', VQ.tooltipHide);
}

// ── Domain helpers ─────────────────────────────────────────────────────────

function _orfArrowPoints(x1, x2, y, h, strand) {
  const tip = Math.min(h * 0.75, (x2 - x1) * 0.25);
  const mid  = y + h / 2;
  if (strand === '+') {
    const ax = x2;
    const bx = x2 - tip;
    return `${x1},${y} ${bx},${y} ${ax},${mid} ${bx},${y + h} ${x1},${y + h}`;
  } else {
    const ax = x1;
    const bx = x1 + tip;
    return `${ax},${mid} ${bx},${y} ${x2},${y} ${x2},${y + h} ${bx},${y + h}`;
  }
}

function _bestDomainPerTarget(domains) {
  const best = {};
  domains.forEach(d => {
    const t = d.target || '';
    if (!best[t] || d.score > best[t].score) best[t] = d;
  });
  return Object.values(best);
}

function _assignDomainLanes(domains) {
  const sorted = [...domains].sort((a, b) => a.start - b.start);
  const lanes  = [];
  return sorted.map(d => {
    for (let lane = 0; lane < 20; lane++) {
      if (!lanes[lane]) lanes[lane] = [];
      const overlaps = lanes[lane].some(([s, e]) =>
        d.start < e && d.stop > s);
      if (!overlaps) {
        lanes[lane].push([d.start, d.stop]);
        return { dom: d, lane };
      }
    }
    return { dom: d, lane: 0 };
  });
}

function _domainColor(target) {
  if (_VW.domainMap[target]) return _VW.domainMap[target];
  const col = _DOM_COLOURS[_VW.colourIdx % _DOM_COLOURS.length];
  _VW.colourIdx++;
  _VW.domainMap[target] = col;
  return col;
}


// ── Sequence-quality overlays ────────────────────────────────────────────────

function _drawLowComplexityBands(svg, regions, SCALE, PAD_L, yTop, yBottom) {
  // Purely visual marker — pointer-events disabled so it never intercepts the
  // coverage track's hover (depth/quality tooltip) underneath it. The region
  // detail is available in the FASTA preview highlight and the dot plot.
  const g = svg.append('g').attr('class', 'vq-lowcx-bands').style('pointer-events', 'none');
  regions.forEach(([a, b]) => {
    const x1 = PAD_L + SCALE(a);
    const x2 = PAD_L + SCALE(b + 1);
    g.append('rect')
      .attr('x', x1).attr('y', yTop)
      .attr('width', Math.max(1, x2 - x1)).attr('height', Math.max(0, yBottom - yTop))
      .attr('fill', 'var(--vq-warning)').attr('opacity', 0.10);
  });
}

// Shared square size (px) for the dot plot and the FASTA-preview pane, so the
// two blocks are visually symmetric side by side.
const _SEQQUAL_SIZE = 300;

// Native-SVG self-similarity dot plot. Main diagonal = identity; off-diagonal
// segments (slope +1) = direct repeats; anti-diagonal (slope −1) = inverted.
function _dotPlotSVG(sq, size = _SEQQUAL_SIZE) {
  const dp = sq.dotplot;
  if (!dp || !(dp.segments || []).length) return null;

  const L   = dp.length || 1;
  const PAD = 34;
  const W   = size, H = size;
  const s   = d3.scaleLinear([0, L], [PAD, W - 8]);
  const sy  = d3.scaleLinear([0, L], [H - PAD, 8]);   // y grows upward

  const svg = d3.create('svg')
    .attr('class', 'vq-dotplot-svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMidYMid meet')
    .style('width', '100%').style('height', '100%').style('max-width', '100%');

  // Plot frame — white fill so the plot stands out against the block background.
  svg.append('rect')
    .attr('x', PAD).attr('y', 8).attr('width', W - PAD - 8).attr('height', H - PAD - 8)
    .attr('fill', 'var(--vq-surface)').attr('stroke', 'var(--vq-border)').attr('stroke-width', 1);

  // Axis labels (bold).
  svg.append('text').attr('x', (W + PAD) / 2).attr('y', H - 8)
    .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('font-weight', 700)
    .attr('fill', 'var(--vq-text-2)').text('query (nt)');
  svg.append('text').attr('x', 10).attr('y', (H - PAD) / 2)
    .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('font-weight', 700)
    .attr('fill', 'var(--vq-text-2)')
    .attr('transform', `rotate(-90 10 ${(H - PAD) / 2})`).text('subject (nt)');

  // Classify every segment with the shared classifier, then draw background
  // layers first (diagonal, low-complexity) and the structural repeats on top,
  // so what is drawn matches exactly what the legend counts.
  const low  = sq.low_complexity_regions || [];
  const g    = svg.append('g');
  const line = (seg, stroke, width, opacity, dash, cap) => {
    const [x1, y1, x2, y2] = seg;
    const l = g.append('line')
      .attr('x1', s(x1)).attr('y1', sy(y1)).attr('x2', s(x2)).attr('y2', sy(y2))
      .attr('stroke', stroke).attr('stroke-width', width).attr('opacity', opacity);
    if (dash) l.attr('stroke-dasharray', dash);
    if (cap)  l.attr('stroke-linecap', cap);
  };
  const repeats = [];
  dp.segments.forEach(seg => {
    const c = _segClass(seg, low);
    if      (c === 'diag')  line(seg, 'var(--vq-text-3)', 1,   0.45, '4,3');
    else if (c === 'lowcx') line(seg, 'var(--vq-warning)', 1.4, 0.5);
    else repeats.push(seg);
  });
  repeats.forEach(seg => {
    const inverted = _segClass(seg, low) === 'inverted';
    line(seg, inverted ? 'var(--vq-danger)' : 'var(--vq-accent)', 2.6, 0.95, null, 'round');
  });
  return svg.node();
}

// Inner content for the middle "Top hit & taxonomy" block: best BLASTx hit
// (NR preferred, else RefSeq) + taxonomy lineage. Returns '' when there is
// neither a hit nor taxonomy to show.
function _topHitContent(seq) {
  const esc  = VQ.esc;
  const t    = seq.taxonomy || {};
  const pool = (seq.blastx_nr_hits && seq.blastx_nr_hits.length)
    ? seq.blastx_nr_hits : (seq.blastx_hits || []);
  const best = pool.length ? pool.reduce((a, b) => (b.bit_score > a.bit_score ? b : a)) : null;
  const hasTax = t.family || t.genus || t.species || t.scientific_name;
  if (!best && !hasTax) return '';

  const kv = (k, v) => (v || v === 0) ? `
    <div style="display:flex;justify-content:space-between;gap:10px">
      <span style="color:var(--vq-text-3)">${esc(k)}</span>
      <span style="text-align:right;word-break:break-word">${esc(String(v))}</span>
    </div>` : '';

  const hit = best ? `
    <div style="font-weight:600;line-height:1.35;margin-bottom:4px;word-break:break-word">
      ${esc(best.subject_title || best.species || best.subject_id || '—')}
    </div>
    ${kv('Accession', best.subject_id)}
    ${kv('Identity',  best.pct_identity != null ? best.pct_identity.toFixed(1) + '%' : null)}
    ${kv('Coverage',  best.query_coverage ? best.query_coverage.toFixed(0) + '%' : null)}
    ${kv('E-value',   best.e_value != null ? best.e_value.toExponential(1) : null)}
    ${kv('Bit score', best.bit_score != null ? Math.round(best.bit_score) : null)}`
    : `<div style="color:var(--vq-text-3)">No BLASTx hit</div>`;

  const tax = hasTax ? `
    <div style="border-top:1px solid var(--vq-border);margin-top:8px;padding-top:8px">
      ${kv('Kingdom',   t.kingdom)}
      ${kv('Phylum',    t.phylum)}
      ${kv('Class',     t['class'])}
      ${kv('Order',     t.order)}
      ${kv('Family',    t.family)}
      ${kv('Subfamily', t.subfamily)}
      ${kv('Genus',     t.genus)}
      ${kv('Species',   t.species || t.scientific_name)}
      ${t.genome ? `
      <div style="border-top:1px dashed var(--vq-border);margin-top:6px;padding-top:6px">
        ${kv('Genome', t.genome)}
      </div>` : ''}
    </div>` : '';

  return hit + tax;
}

// One classifier shared by the dot plot, the legend count and the FASTA
// highlight, so all three agree. A segment is:
//   diag     — the identity diagonal
//   lowcx    — a self-match inside a dust-masked low-complexity region
//              (counted under low-complexity %, not as a structural repeat)
//   direct   — off-diagonal, slope +1 (duplicated segment)
//   inverted — anti-diagonal, slope −1 (reverse-complement segment)
function _inLowCx(pos, regions) {
  for (let k = 0; k < regions.length; k++)
    if (pos >= regions[k][0] && pos <= regions[k][1]) return true;
  return false;
}
function _segClass(seg, low) {
  const [x1, y1, x2, y2] = seg;
  if (Math.abs(y1 - x1) < 1e-6 && Math.abs(y2 - x2) < 1e-6) return 'diag';
  if (_inLowCx((x1 + x2) / 2, low) || _inLowCx((y1 + y2) / 2, low)) return 'lowcx';
  return (y2 - y1) * (x2 - x1) < 0 ? 'inverted' : 'direct';
}

// Structural repeats (from the drawn dot-plot segments) — the single source of
// truth for the legend count and the FASTA highlight.
function _dotplotRepeats(sq) {
  const low = sq.low_complexity_regions || [];
  const direct = [], inverted = [];
  (sq.dotplot && sq.dotplot.segments || []).forEach(seg => {
    const c = _segClass(seg, low);
    if (c === 'direct') direct.push(seg);
    else if (c === 'inverted') inverted.push(seg);
  });
  return { direct, inverted };
}

// Compact textual summary of repeats (from the dot plot) + k-mer repetitiveness.
function _repeatSummary(sq) {
  const rep      = _dotplotRepeats(sq);
  const direct   = rep.direct.length;
  const inverted = rep.inverted.length;
  const score    = ((sq.kmer_repeat_score || 0) * 100).toFixed(1);
  const lowcx   = ((sq.low_complexity_frac || 0) * 100).toFixed(1);
  const dot = c => `<span style="display:inline-block;width:8px;height:8px;border-radius:2px;`
    + `background:${c};margin-right:7px;vertical-align:middle;flex:0 0 auto"></span>`;
  // Full-name legend items, stacked in a column to the right of the dot plot.
  const rows = [
    `${dot('var(--vq-accent)')}${direct} direct repeat${direct !== 1 ? 's' : ''}`,
    `${dot('var(--vq-danger)')}${inverted} inverted repeat${inverted !== 1 ? 's' : ''}`,
    `${dot('var(--vq-warning)')}low-complexity ${lowcx}%`,
    `k${sq.kmer_size} repeat score ${score}%`,
  ];
  return rows.map(r => `<div style="display:flex;align-items:center;white-space:nowrap">${r}</div>`).join('');
}

// Colour codes matching the legend / dot plot: 1 low-complexity (amber),
// 2 direct repeat (blue), 3 inverted repeat (red). Higher priority wins overlaps.
const _HL_STYLE = {
  1: 'background:rgba(224,146,28,.32)',   // --vq-warning
  2: 'background:rgba(47,134,214,.30)',   // --vq-accent
  3: 'background:rgba(210,74,61,.34)',    // --vq-danger
};
const _HL_TITLE = { 1: 'low complexity', 2: 'direct repeat', 3: 'inverted repeat' };

// FASTA header + sequence, with seq-quality regions wrapped in coloured <span>s.
function _fastaHighlightedHTML(seq) {
  const esc = VQ.esc;
  const raw = VQ.buildFasta([seq]).replace(/\n+$/, '');
  const nl  = raw.indexOf('\n');
  const header = nl >= 0 ? raw.slice(0, nl) : raw;
  const body   = nl >= 0 ? raw.slice(nl + 1) : '';
  const sq = seq.seq_quality;
  if (!sq || !body) return esc(raw);

  const n    = body.length;
  const code = new Uint8Array(n);
  const mark = (a, b, v) => {
    for (let i = Math.max(0, a | 0); i <= Math.min(n - 1, b | 0); i++)
      if (v > code[i]) code[i] = v;
  };
  // Low-complexity spans follow this card's checkbox; repeats are always shown.
  if (_showLowCx(seq))
    (sq.low_complexity_regions || []).forEach(([a, b]) => mark(a, b, 1));
  // Repeats from the same dot-plot segments the legend counts — both copies of
  // each segment (query x-range and subject y-range) are marked.
  const rep = _dotplotRepeats(sq);
  const markSeg = (segs, v) => segs.forEach(([x1, y1, x2, y2]) => {
    mark(Math.min(x1, x2), Math.max(x1, x2), v);
    mark(Math.min(y1, y2), Math.max(y1, y2), v);
  });
  markSeg(rep.direct, 2);
  markSeg(rep.inverted, 3);

  let html = esc(header) + '\n';
  let i = 0;
  while (i < n) {
    const c = code[i];
    let j = i + 1;
    while (j < n && code[j] === c) j++;
    const chunk = esc(body.slice(i, j));
    html += c === 0
      ? chunk
      : `<span style="${_HL_STYLE[c]};border-radius:2px" title="${_HL_TITLE[c]}">${chunk}</span>`;
    i = j;
  }
  return html;
}

// ── Info (?) tooltip bodies ────────────────────────────────────────────────

// Shared row builders: `_tipRow(colour, text)` — pass '' as colour for a
// bullet-less row that still aligns with the swatched ones.
const _tipDot = c => `<span style="display:inline-block;width:9px;height:9px;border-radius:2px;`
  + `background:${c};margin-right:7px;flex:0 0 auto;margin-top:3px"></span>`;
const _tipRow = (c, t) =>
  `<div style="display:flex;align-items:flex-start;margin:3px 0;max-width:250px">`
  + `${c ? _tipDot(c) : '<span style="width:16px;flex:0 0 auto"></span>'}<span>${t}</span></div>`;
const _tipNote = t =>
  `<div style="margin-top:5px;color:var(--vq-text-3);font-size:11px;max-width:250px">${t}</div>`;
const _tipLead = t =>
  `<div style="margin-bottom:4px;color:var(--vq-text-3);font-size:11px;max-width:250px">${t}</div>`;

// Genome map: what each track/glyph is. Coverage rows are gated on the profile
// actually present for this sequence, so the tooltip never promises a track the
// card did not draw.
function _mapInfoHTML(seq) {
  const cov    = (seq.coverage && (seq.coverage.bins || []).length) ? seq.coverage : null;
  const nb     = cov ? cov.bins.length : 0;
  const qual   = cov && (cov.quality_bins || []).length === nb;
  const strand = cov && (cov.bins_fwd || []).length === nb
                     && (cov.bins_rev || []).length === nb;
  const lowcx  = ((seq.seq_quality || {}).low_complexity_regions || []).length;
  return `
    <div class="vq-tooltip__title">Genome map</div>
    ${_tipLead('ORFs, HMM domains and read support along the contig, in nucleotide coordinates.')}
    ${_tipRow('var(--vq-accent)',  'ORF (+) — arrow points 5′→3′ on the forward strand; ORFs are laned by reading frame.')}
    ${_tipRow('var(--vq-warning)', 'ORF (−) — arrow points 3′→5′ on the reverse strand.')}
    ${_tipRow('var(--vq-dom-1)',   'HMM domain — profile hit (RVDB / Vfam / Pfam / EggNOG) drawn inside its ORF; each target gets its own colour.')}
    ${cov ? _tipRow(qual ? 'var(--vq-success)' : 'var(--vq-accent)',
        qual ? 'Read coverage — bar height is depth, bar colour is mean base quality (green Q≥30, amber Q20–29, red Q&lt;20).'
             : 'Read coverage — depth per bin, from the reads aligned back to the contig.') : ''}
    ${strand ? _tipRow('', 'Sense reads grow up from the centre line, antisense grow down (＋ / － in the gutter).') : ''}
    ${lowcx && _showLowCx(seq) ? _tipRow('var(--vq-warning)', 'Amber band — low-complexity region; shading only, it never blocks the hover below it.') : ''}
    ${lowcx && !_showLowCx(seq) ? _tipRow('', 'Low-complexity shading is off for this sequence — re-enable it with the checkbox next to the Genome map label.') : ''}
    ${_tipNote('Hover any ORF, domain or the coverage band for details.')}`;
}

// FASTA preview: the highlight palette (same regions the dot plot shows).
function _fastaInfoHTML(seq) {
  return `
    <div class="vq-tooltip__title">FASTA preview</div>
    ${_tipLead('Sequence-quality regions are highlighted in place, in the same colours as the dot plot.')}
    ${_tipRow('var(--vq-accent)',  'Direct repeat — both copies of a segment duplicated in the same orientation.')}
    ${_tipRow('var(--vq-danger)',  'Inverted repeat — a segment and the reverse complement it matches.')}
    ${_showLowCx(seq)
        ? _tipRow('var(--vq-warning)', 'Low complexity — dustmasker-flagged low-information region.')
        : _tipRow('', 'Low-complexity highlighting is off for this sequence — re-enable it with the checkbox next to the Genome map label.')}
    ${_tipNote('Where regions overlap, the strongest signal wins. Hover a highlight for its label.')}`;
}

// Heuristic score: the rules, mirroring score_heuristic.py's CALIBRATION block.
function _heurInfoHTML() {
  return `
    <div class="vq-tooltip__title">Heuristic score · basic rules</div>
    ${_tipLead('Deterministic 0–100 score computed from the evidence itself — no LLM involved.')}
    ${_tipRow('', 'Three components: <strong>BLASTn 30%</strong>, <strong>BLASTx 30%</strong>, <strong>HMM 40%</strong>. Each BLAST component mixes identity (60%) and query coverage (40%); HMM mixes best bit score with hit multiplicity.')}
    ${_tipRow('', 'A component that is <em>absent</em> is dropped and its weight redistributed over the rest — absent is not the same as zero.')}
    ${_tipRow('', 'Low-coverage penalty ×0.5 inside a BLAST component when the hit covers &lt;30% of the query.')}
    ${_tipRow('', 'Single-evidence penalty ×0.7 when only one component is present (no corroboration).')}
    ${_tipRow('', 'False-positive penalty ×0.2 when the dominant BLASTn hit is non-viral at ≥90% identity and ≥80% coverage — this one also forces the <strong>non-viral</strong> class.')}
    ${_tipNote('Class: viral-known needs ≥90% identity and ≥70% coverage in a BLAST channel; otherwise viral-unknown.')}`;
}

// Shared tooltip body for the sequence-quality legend info (?) badge.
function _signalsInfoHTML() {
  return `
    <div class="vq-tooltip__title">Sequence-quality signals</div>
    ${_tipRow('var(--vq-accent)',  'Direct repeat — a segment duplicated elsewhere in the same orientation.')}
    ${_tipRow('var(--vq-danger)',  'Inverted repeat — a segment matching the reverse complement of another.')}
    ${_tipRow('var(--vq-warning)', 'Low complexity — dustmasker-flagged low-information region (AT-rich, homopolymer, …).')}
    ${_tipRow('', 'Repeat score — fraction of k-mers (length k) that occur more than once in the sequence.')}
    ${_tipNote('The same colours mark these regions in the FASTA sequence.')}`;
}

window.vqInitViewer = vqInitViewer;
})();
