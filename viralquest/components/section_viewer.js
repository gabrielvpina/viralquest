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
};

const _DOM_COLOURS = [
  'var(--vq-dom-1)','var(--vq-dom-2)','var(--vq-dom-3)','var(--vq-dom-4)',
  'var(--vq-dom-5)','var(--vq-dom-6)','var(--vq-dom-7)','var(--vq-dom-8)',
];

// ── Init ────────────────────────────────────────────────────────────────────

function vqInitViewer(sequences) {
  const el = document.getElementById('section-viewer');
  if (!el) return;
  const esc = VQ.esc;

  _VW.sequences = sequences;
  _VW.filtered  = [...sequences];

  const phyla    = [...new Set(sequences.map(s => s.taxonomy?.phylum).filter(Boolean))].sort();
  const families = [...new Set(sequences.map(s => s.taxonomy?.family).filter(Boolean))].sort();
  const genera   = [...new Set(sequences.map(s => s.taxonomy?.genus).filter(Boolean))].sort();

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
               placeholder="Search sequence ID or species…" aria-label="Search sequences">
      </div>
      <div class="vq-toolbar__filters">
        <select class="vq-select" id="viewer-phylum" aria-label="Filter by phylum">
          <option value="">All phyla</option>
          ${phyla.map(p => `<option value="${esc(p)}">${esc(p)}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-family" aria-label="Filter by family">
          <option value="">All families</option>
          ${families.map(f => `<option value="${esc(f)}">${esc(f)}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-genus" aria-label="Filter by genus">
          <option value="">All genera</option>
          ${genera.map(g => `<option value="${esc(g)}">${esc(g)}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-score" aria-label="Filter by VQ score">
          <option value="">Any VQ score</option>
          <option value="80">VQ ≥ 80</option>
          <option value="60">VQ ≥ 60</option>
          <option value="40">VQ ≥ 40</option>
        </select>
        <select class="vq-select" id="viewer-classification" aria-label="Filter by classification">
          <option value="">Any classification</option>
          <option value="viral-known">Viral known</option>
          <option value="viral-unknown">Viral unknown</option>
          <option value="non-viral">Non-viral</option>
        </select>
      </div>
    </div>

    <div id="viewer-list"></div>
  `;

  ['viewer-search','viewer-phylum','viewer-family','viewer-genus','viewer-score','viewer-classification']
    .forEach(id => document.getElementById(id)?.addEventListener('input', _applyFilters));

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

  _renderList();
}

// ── Filtering ──────────────────────────────────────────────────────────────

function _applyFilters() {
  const q       = document.getElementById('viewer-search')?.value.toLowerCase() || '';
  const phylum  = document.getElementById('viewer-phylum')?.value || '';
  const family  = document.getElementById('viewer-family')?.value || '';
  const genus   = document.getElementById('viewer-genus')?.value  || '';
  const cls     = document.getElementById('viewer-classification')?.value || '';
  const minScore = parseInt(document.getElementById('viewer-score')?.value || '0', 10);

  _VW.filtered = _VW.sequences.filter(s => {
    if (q && !s.id.toLowerCase().includes(q) &&
        !(s.taxonomy?.species || '').toLowerCase().includes(q)) return false;
    if (phylum && s.taxonomy?.phylum !== phylum) return false;
    if (family && s.taxonomy?.family !== family) return false;
    if (genus  && s.taxonomy?.genus  !== genus)  return false;
    if (cls    && s.llm_output?.classification !== cls) return false;
    if (minScore && (s.llm_output?.vq_score ?? 0) < minScore) return false;
    return true;
  });

  const countEl = document.getElementById('viewer-count');
  if (countEl) countEl.textContent =
    `${_VW.filtered.length} of ${_VW.sequences.length} sequence${_VW.sequences.length !== 1 ? 's' : ''}`;

  _renderList();
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

// ── Sequence card ──────────────────────────────────────────────────────────

function _seqCard(seq) {
  const esc  = VQ.esc;
  const safe = VQ.safeId(seq.id);
  const card = document.createElement('div');
  card.className = 'vq-seq-card';
  card.id        = 'seq-card-' + safe;
  card.dataset.seqId = seq.id;
  card.style.marginBottom = '6px';

  const llm = seq.llm_output;
  const tax = seq.taxonomy;
  const cls = llm?.classification || 'non-viral';

  const orfCount     = (seq.orfs || []).length;
  const blastnCount  = (seq.blastn_hits || []).length;
  const refseqCount  = (seq.blastx_hits || []).length;
  const nrCount      = (seq.blastx_nr_hits || []).length;

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
        ${llm ? `<span class="vq-badge vq-badge--${esc(cls)}">${esc(cls)} · ${llm.vq_score}</span>` : ''}
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

      <div class="vq-body-label">Genome map</div>
      <div class="vq-genome-wrap" id="genome-wrap-${safe}"></div>

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

      ${llm?.analysis ? `
        <div class="vq-body-label">LLM analysis</div>
        <div style="font-size:12.5px;color:var(--vq-text-2);line-height:1.55;
                    background:var(--vq-surface-2);border:1px solid var(--vq-border);
                    border-radius:var(--vq-radius-sm);padding:10px 12px">
          <div style="font-size:11px;color:var(--vq-text-3);margin-bottom:4px">
            ${esc(llm.classification)} · score ${llm.vq_score}
          </div>
          ${esc(llm.analysis)}
        </div>` : ''}

      ${seq.sequence ? `
        <div class="vq-body-label">FASTA preview</div>
        <div class="vq-fasta">${esc(seq.sequence)}</div>` : ''}

    </div>`;

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
    if      (kind === 'blastn') hits = seq.blastn_hits    || [];
    else if (kind === 'refseq') hits = seq.blastx_hits    || [];
    else                        hits = seq.blastx_nr_hits || [];
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
  const esc     = VQ.esc;
  const isBlastx = kind === 'refseq' || kind === 'nr';

  const rows = hits.map(h => {
    if (isBlastx) {
      const pct = h.pct_identity  != null ? h.pct_identity.toFixed(1)   : '—';
      const cov = h.query_coverage != null ? h.query_coverage.toFixed(1) : '—';
      const ev  = h.e_value        != null ? h.e_value.toExponential(1)  : '—';
      const bs  = h.bit_score      != null ? h.bit_score.toFixed(0)      : '—';
      const rng = (h.query_start != null && h.query_end != null)
                    ? `${h.query_start}–${h.query_end}` : '—';
      return `
        <tr>
          <td class="vq-td--mono">${esc(h.subject_id   ?? '—')}</td>
          <td>${esc(h.subject_title ?? '—')}</td>
          <td class="vq-td--num">${pct}%</td>
          <td class="vq-td--num">${cov}%</td>
          <td class="vq-td--num">${ev}</td>
          <td class="vq-td--num">${bs}</td>
          <td class="vq-td--num">${rng}</td>
        </tr>`;
    } else {
      const pct = h.pident  != null ? h.pident.toFixed(1)       : '—';
      const ev  = h.evalue  != null ? h.evalue.toExponential(1) : '—';
      const bs  = h.bit_score != null ? h.bit_score.toFixed(0)  : '—';
      return `
        <tr>
          <td class="vq-td--mono">—</td>
          <td>${esc(h.stitle ?? '—')}</td>
          <td class="vq-td--num">${pct}%</td>
          <td class="vq-td--num">${h.qcovhsp ?? '—'}%</td>
          <td class="vq-td--num">${ev}</td>
          <td class="vq-td--num">${bs}</td>
          <td class="vq-td--num">—</td>
        </tr>`;
    }
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

  // ── Frame-based lane assignment ─────────────────────────────────────────
  // orf.frame is serialised as an integer (1,2,3,-1,-2,-3) from Python.
  const FRAME_LABEL = ['+1','+2','+3','-1','-2','-3'];

  function frameLane(orf) {
    const n = Number(orf.frame);
    if (Number.isFinite(n) && n !== 0) {
      return n > 0 ? n - 1 : 3 + (-n - 1);
    }
    return orf.strand === '-' ? 3 : 0;
  }
  const nLanes = 6;

  // For each ORF keep only the highest-scoring domain per database.
  function _bestPerDb(domains) {
    const best = {};
    for (const d of domains) {
      const db = d.database || '';
      if (!best[db] || d.score > best[db].score) best[db] = d;
    }
    return Object.values(best);
  }

  // Pack domains into sub-lanes per ORF
  const orfsWithFrame = orfs.map(o => ({
    ...o,
    frameLane:    frameLane(o),
    _domainLanes: _assignDomainLanes(_bestPerDb(o.domains || [])),
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
  let cy = AXIS_Y + 14;
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
  const gridH = totalH - AXIS_Y - 26;
  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, ${AXIS_Y})`)
    .call(d3.axisTop(SCALE).ticks(Math.min(14, Math.max(4, Math.round(W / 120)))).tickSize(0));
  axisG.selectAll('.tick line')
    .attr('y1', 0).attr('y2', gridH)
    .attr('stroke', 'var(--vq-border)').attr('stroke-dasharray', '3,3');
  axisG.select('.domain').attr('stroke', 'var(--vq-border-dark)');
  axisG.selectAll('.tick text')
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)');

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

    // Build tooltip with AA sequence preview if available
    const aa = orf.aa_sequence || '';
    const aaPreview = aa.length > 60 ? aa.slice(0, 60) + '…' : aa;
    const aaRow = aa
      ? `<div style="margin-top:6px;padding-top:6px;border-top:1px solid rgba(255,255,255,.12)">
           <div class="vq-tooltip__key" style="margin-bottom:3px">Sequence (${aa.length} aa)</div>
           <div style="font-family:var(--vq-font-mono);font-size:10px;word-break:break-all;
                       max-width:240px;color:rgba(255,255,255,.92)">${VQ.esc(aaPreview)}</div>
         </div>`
      : '';

    g.append('rect')
      .attr('x', x1).attr('width', bw)
      .attr('y', y).attr('height', ORF_H)
      .attr('rx', 3)
      .attr('fill', col).attr('opacity', .88)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(orf.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Strand</span><span>${VQ.esc(orf.strand)}</span>
          <span class="vq-tooltip__key">Frame</span><span>${VQ.esc(orf.frame)}</span>
          <span class="vq-tooltip__key">Length</span><span>${orf.length_aa} aa</span>
          <span class="vq-tooltip__key">Position</span>
          <span>${orf.start_position}–${orf.stop_position}</span>
          <span class="vq-tooltip__key">Type</span><span>${VQ.esc(orf.orf_type)}</span>
          <span class="vq-tooltip__key">Domains</span><span>${(orf.domains||[]).length}</span>
        </div>
        ${aaRow}`, evt))
      .on('mouseleave', VQ.tooltipHide);

    if (bw > 38) {
      g.append('text')
        .attr('x', x1 + bw / 2).attr('y', y + ORF_H / 2 + 3.5)
        .attr('text-anchor', 'middle')
        .attr('font-size', 9.5)
        .attr('fill', 'rgba(255,255,255,.92)')
        .attr('pointer-events', 'none')
        .text((orf.strand === '+' ? '▶ ' : '◀ ') + orf.name);
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
            <span class="vq-tooltip__key">In</span><span>${VQ.esc(orf.name)}</span>
            <span class="vq-tooltip__key">DB</span><span>${VQ.esc(d.dom.database)}</span>
            <span class="vq-tooltip__key">Description</span><span>${VQ.esc(d.dom.description)}</span>
            <span class="vq-tooltip__key">Score</span><span>${d.dom.score}</span>
            <span class="vq-tooltip__key">E-value</span><span>${d.dom.e_value?.toExponential(2) ?? '—'}</span>
            <span class="vq-tooltip__key">aa range</span><span>${d.dom.start}–${d.dom.stop}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    });
  });

  // ── Legend ──────────────────────────────────────────────────────────────
  const legendG = svg.append('g').attr('transform', `translate(${PAD_L}, ${totalH - 18})`);
  legendG.append('rect').attr('width', 12).attr('height', 9).attr('rx', 2)
    .attr('fill', 'var(--vq-accent)');
  legendG.append('text').attr('x', 18).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('(+) ORF');
  legendG.append('rect').attr('x', 70).attr('width', 12).attr('height', 9).attr('rx', 2)
    .attr('fill', 'var(--vq-warning)');
  legendG.append('text').attr('x', 88).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('(−) ORF');
  legendG.append('rect').attr('x', 140).attr('width', 12).attr('height', 5).attr('rx', 1)
    .attr('fill', 'var(--vq-dom-1)').attr('y', 2);
  legendG.append('text').attr('x', 158).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text('HMM domain');
  if (!fluidMode) {
    legendG.append('text').attr('x', 230).attr('y', 8.5)
      .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)')
      .text(`↔ ${seqLen.toLocaleString()} nt · scroll horizontally`);
  }

  return svg.node();
}

// ── Domain lane assignment (within an ORF) ─────────────────────────────────

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


window.vqInitViewer = vqInitViewer;
})();
