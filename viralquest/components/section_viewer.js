/* ============================================================
   section_viewer.js — Section 3: All Viruses Viewer
   Searchable / filterable list of sequence cards.
   Each card expands to a genome track SVG showing all ORFs
   (vertical lane packing for overlaps) plus HMM domains.
   ============================================================ */

'use strict';

// ── Module state ─────────────────────────────────────────────────────────────
const _VW = {
  sequences:  [],
  filtered:   [],
  selected:   new Set(),
  domainMap:  {},          // domain target → colour
  colourIdx:  0,
};

const _DOM_COLOURS = [
  'var(--vq-dom-1)','var(--vq-dom-2)','var(--vq-dom-3)','var(--vq-dom-4)',
  'var(--vq-dom-5)','var(--vq-dom-6)','var(--vq-dom-7)','var(--vq-dom-8)',
];

// ── Init ─────────────────────────────────────────────────────────────────────

function vqInitViewer(sequences) {
  const el = document.getElementById('section-viewer');
  if (!el) return;
  _VW.sequences = sequences;
  _VW.filtered  = [...sequences];

  // Collect taxonomy options
  const phyla    = [...new Set(sequences.map(s => s.taxonomy?.phylum).filter(Boolean))].sort();
  const families = [...new Set(sequences.map(s => s.taxonomy?.family).filter(Boolean))].sort();
  const genera   = [...new Set(sequences.map(s => s.taxonomy?.genus).filter(Boolean))].sort();

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">All Viruses Viewer</div>
        <div class="vq-section-sub" id="viewer-count">
          ${sequences.length} sequence${sequences.length !== 1 ? 's' : ''}
        </div>
      </div>
      <div style="display:flex;gap:var(--vq-space-2);flex-wrap:wrap">
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-select-all">Select all</button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-deselect">Deselect</button>
        <button class="vq-btn vq-btn--primary vq-btn--sm" id="viewer-export-batch" disabled>
          Export selected (PNG)
        </button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="viewer-export-svg-batch" disabled>
          Export selected (SVG)
        </button>
      </div>
    </div>

    <div class="vq-toolbar">
      <div class="vq-search-wrap">
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2">
          <circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/>
        </svg>
        <input class="vq-input" id="viewer-search" placeholder="Search sequence ID or species…" type="search">
      </div>
      <div class="vq-toolbar__filters">
        <select class="vq-select" id="viewer-phylum">
          <option value="">All phyla</option>
          ${phyla.map(p => `<option value="${p}">${p}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-family">
          <option value="">All families</option>
          ${families.map(f => `<option value="${f}">${f}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-genus">
          <option value="">All genera</option>
          ${genera.map(g => `<option value="${g}">${g}</option>`).join('')}
        </select>
        <select class="vq-select" id="viewer-score">
          <option value="">Any VQ score</option>
          <option value="80">VQ ≥ 80</option>
          <option value="60">VQ ≥ 60</option>
          <option value="40">VQ ≥ 40</option>
        </select>
      </div>
    </div>

    <div id="viewer-list"></div>
  `;

  // Wire up controls
  ['viewer-search','viewer-phylum','viewer-family','viewer-genus','viewer-score'].forEach(id => {
    document.getElementById(id)?.addEventListener('input', _applyFilters);
  });

  document.getElementById('viewer-select-all')?.addEventListener('click', () => {
    _VW.filtered.forEach(s => _VW.selected.add(s.id));
    _renderList();
    _updateBatchBtn();
  });

  document.getElementById('viewer-deselect')?.addEventListener('click', () => {
    _VW.selected.clear();
    _renderList();
    _updateBatchBtn();
  });

  document.getElementById('viewer-export-batch')?.addEventListener('click', _batchExportPNG);
  document.getElementById('viewer-export-svg-batch')?.addEventListener('click', _batchExportSVG);

  _renderList();
}

// ── Filtering ────────────────────────────────────────────────────────────────

function _applyFilters() {
  const q      = document.getElementById('viewer-search')?.value.toLowerCase() || '';
  const phylum = document.getElementById('viewer-phylum')?.value || '';
  const family = document.getElementById('viewer-family')?.value || '';
  const genus  = document.getElementById('viewer-genus')?.value  || '';
  const minScore = parseInt(document.getElementById('viewer-score')?.value || '0', 10);

  _VW.filtered = _VW.sequences.filter(s => {
    if (q && !s.id.toLowerCase().includes(q) &&
        !(s.taxonomy?.species || '').toLowerCase().includes(q)) return false;
    if (phylum && s.taxonomy?.phylum !== phylum) return false;
    if (family && s.taxonomy?.family !== family) return false;
    if (genus  && s.taxonomy?.genus  !== genus)  return false;
    if (minScore && (s.llm_output?.vq_score ?? 0) < minScore) return false;
    return true;
  });

  const countEl = document.getElementById('viewer-count');
  if (countEl) countEl.textContent =
    `${_VW.filtered.length} of ${_VW.sequences.length} sequence${_VW.sequences.length !== 1 ? 's' : ''}`;

  _renderList();
}

// ── List rendering ───────────────────────────────────────────────────────────

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

// ── Sequence card ────────────────────────────────────────────────────────────

function _seqCard(seq) {
  const card = document.createElement('div');
  card.className  = 'vq-seq-card vq-card';
  card.id         = 'seq-card-' + seq.id;
  card.style.marginBottom = 'var(--vq-space-3)';

  const llm    = seq.llm_output;
  const tax    = seq.taxonomy;
  const cls    = llm?.classification || 'non-viral';
  const badgeCls = `vq-badge--${cls.replace('-','–').replace(/\s/g,'_') || 'non-viral'}`;

  card.innerHTML = `
    <div class="vq-seq-card__head">
      <input type="checkbox" class="vq-checkbox" ${_VW.selected.has(seq.id) ? 'checked' : ''}
             onchange="vqViewerToggle('${seq.id}', this.checked)">
      <span class="vq-seq-card__id">${seq.id}</span>
      <div class="vq-seq-card__meta">
        ${tax?.family ? `<span>📂 ${tax.family}</span>` : ''}
        ${tax?.genus  ? `<span>🔬 ${tax.genus}</span>`  : ''}
        <span>${(seq.length || 0).toLocaleString()} nt</span>
        <span>GC ${seq.gc_content?.toFixed(1) ?? '—'}%</span>
        <span>${(seq.orfs || []).length} ORF${(seq.orfs || []).length !== 1 ? 's' : ''}</span>
        ${seq.cluster_id ? `<span class="vq-badge vq-badge--cluster">${seq.cluster_id}</span>` : ''}
        ${llm ? `<span class="vq-badge vq-badge--${cls.replace('-','–')}">${cls} · ${llm.vq_score}</span>` : ''}
      </div>
      <div style="display:flex;gap:var(--vq-space-2);margin-left:auto;align-items:center">
        <button class="vq-btn vq-btn--sm vq-btn--ghost" onclick="vqExportSeqPNG('${seq.id}')">PNG</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost" onclick="vqExportSeqSVG('${seq.id}')">SVG</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost" onclick="vqExportSeqPDF('${seq.id}')">PDF</button>
      </div>
      <svg class="vq-seq-card__chevron" width="14" height="14" viewBox="0 0 24 24"
           fill="none" stroke="currentColor" stroke-width="2.5">
        <polyline points="9 18 15 12 9 6"/>
      </svg>
    </div>
    <div class="vq-seq-card__body" id="body-${seq.id}">
      <div class="vq-genome-wrap" id="genome-wrap-${seq.id}"></div>
      ${llm?.analysis ? `
        <div style="margin-top:var(--vq-space-4);padding-top:var(--vq-space-4);
                    border-top:1px solid var(--vq-border);">
          <div style="font-size:var(--vq-text-xs);font-weight:600;color:var(--vq-text-3);
                      text-transform:uppercase;letter-spacing:.06em;margin-bottom:4px">
            LLM Analysis
          </div>
          <div style="font-size:var(--vq-text-sm);color:var(--vq-text-2);line-height:1.6">
            ${llm.analysis}
          </div>
        </div>` : ''}
    </div>`;

  // Toggle open/close
  const head = card.querySelector('.vq-seq-card__head');
  head.addEventListener('click', e => {
    if (e.target.closest('button') || e.target.closest('input')) return;
    const open = !card.classList.contains('open');
    card.classList.toggle('open', open);
    if (open) {
      const wrap = document.getElementById('genome-wrap-' + seq.id);
      if (wrap && !wrap.querySelector('svg')) {
        wrap.appendChild(_genomeSVG(seq));
      }
    }
  });

  return card;
}

// ── Public helpers (called from HTML onclick) ────────────────────────────────

function vqViewerToggle(seqId, checked) {
  if (checked) _VW.selected.add(seqId);
  else         _VW.selected.delete(seqId);
  _updateBatchBtn();
}

function vqExportSeqPNG(seqId) {
  const svg = _getSVG(seqId);
  if (svg) vqExportPNG(svg, `${seqId}.png`);
}

function vqExportSeqSVG(seqId) {
  const svg = _getSVG(seqId);
  if (svg) vqExportSVG(svg, `${seqId}.svg`);
}

function vqExportSeqPDF(seqId) {
  const svg = _getSVG(seqId);
  if (svg) vqExportPDF(svg, `ViralQuest — ${seqId}`);
}

function _getSVG(seqId) {
  // Ensure card is open and SVG rendered
  const card = document.getElementById('seq-card-' + seqId);
  if (!card) return null;
  if (!card.classList.contains('open')) {
    card.classList.add('open');
    const wrap = document.getElementById('genome-wrap-' + seqId);
    const seq  = _VW.sequences.find(s => s.id === seqId);
    if (wrap && seq && !wrap.querySelector('svg')) {
      wrap.appendChild(_genomeSVG(seq));
    }
  }
  return document.getElementById('genome-wrap-' + seqId)?.querySelector('svg');
}

function _updateBatchBtn() {
  const n = _VW.selected.size;
  ['viewer-export-batch','viewer-export-svg-batch'].forEach(id => {
    const btn = document.getElementById(id);
    if (btn) {
      btn.disabled = n === 0;
      if (id === 'viewer-export-batch')
        btn.textContent = n ? `Export ${n} selected (PNG)` : 'Export selected (PNG)';
      else
        btn.textContent = n ? `Export ${n} selected (SVG)` : 'Export selected (SVG)';
    }
  });
}

async function _batchExportPNG() {
  const items = [..._VW.selected].map(id => ({
    svgEl:    _getSVG(id) || document.createElement('svg'),
    filename: `${id}.png`,
  })).filter(x => x.svgEl);
  await vqExportBatchPNG(items, 400);
}

async function _batchExportSVG() {
  for (const id of _VW.selected) {
    const svg = _getSVG(id);
    if (svg) { vqExportSVG(svg, `${id}.svg`); await _sleep(200); }
  }
}

const _sleep = ms => new Promise(r => setTimeout(r, ms));

// ── Genome SVG ────────────────────────────────────────────────────────────────

function _genomeSVG(seq) {
  const orfs   = seq.orfs || [];
  const seqLen = seq.length || 1;

  const PAD_L  = 20;
  const PAD_R  = 20;
  const W      = Math.max(900, seqLen / 8);  // scale down very long sequences
  const drawW  = W - PAD_L - PAD_R;
  const SCALE  = d3.scaleLinear([0, seqLen], [0, drawW]);

  // Assign ORFs to vertical lanes
  const orfsWithLanes = _assignLanes(orfs);
  const nLanes        = Math.max(...orfsWithLanes.map(o => o.lane), 0) + 1;

  const AXIS_Y  = 30;
  const LANE_H  = 22;
  const LANE_G  = 8;
  const ORF_H   = 16;
  const DOM_H   = 8;
  const totalH  = AXIS_Y + nLanes * (LANE_H + LANE_G) + 50;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('width',   W)
    .attr('height',  totalH);

  // --- axis ---
  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, ${AXIS_Y})`)
    .call(d3.axisTop(SCALE).ticks(8).tickSize(-(totalH - AXIS_Y - 20)));

  axisG.select('.domain').attr('stroke', 'var(--vq-border-dark)');
  axisG.selectAll('.tick line').attr('stroke', 'var(--vq-border)').attr('stroke-dasharray', '3,3');
  axisG.selectAll('.tick text').style('font-size', '9px').style('fill', 'var(--vq-text-3)');

  // --- backbone ---
  svg.append('line')
    .attr('x1', PAD_L).attr('x2', PAD_L + drawW)
    .attr('y1', AXIS_Y + 2).attr('y2', AXIS_Y + 2)
    .attr('stroke', 'var(--vq-border-dark)')
    .attr('stroke-width', 2);

  // --- ORFs ---
  orfsWithLanes.forEach(orf => {
    const x1  = PAD_L + SCALE(orf.start_position);
    const x2  = PAD_L + SCALE(orf.stop_position);
    const bw  = Math.max(x2 - x1, 4);
    const y   = AXIS_Y + 10 + orf.lane * (LANE_H + LANE_G);
    const col = orf.strand === '-' ? 'var(--vq-warning)' : 'var(--vq-accent)';

    const g = svg.append('g');

    // ORF box
    g.append('rect')
      .attr('x', x1).attr('width', bw)
      .attr('y', y).attr('height', ORF_H)
      .attr('rx', 3)
      .attr('fill', col)
      .attr('opacity', .8)
      .on('mousemove', evt => {
        vqTooltipShow(`
          <div class="vq-tooltip__title">${orf.name}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Strand</span><span>${orf.strand}</span>
            <span class="vq-tooltip__key">Frame</span><span>${orf.frame}</span>
            <span class="vq-tooltip__key">Length</span><span>${orf.length_aa} aa</span>
            <span class="vq-tooltip__key">Position</span>
            <span>${orf.start_position}–${orf.stop_position}</span>
            <span class="vq-tooltip__key">Type</span><span>${orf.orf_type}</span>
            <span class="vq-tooltip__key">Domains</span><span>${(orf.domains||[]).length}</span>
          </div>`, evt);
      })
      .on('mouseleave', vqTooltipHide);

    // Strand arrow indicator
    const arrowX = orf.strand === '+' ? x1 + bw - 8 : x1 + 8;
    g.append('text')
      .attr('x', x1 + bw / 2).attr('y', y + ORF_H / 2 + 4)
      .attr('text-anchor', 'middle')
      .attr('font-size', 8)
      .attr('fill', 'rgba(255,255,255,.85)')
      .attr('pointer-events', 'none')
      .text(bw > 40 ? (orf.strand === '+' ? '▶' : '◀') + ' ' + orf.name : '');

    // --- HMM domains (overlaid below ORF) ---
    (orf.domains || []).forEach(dom => {
      const domStart = orf.start_position + dom.start * 3;
      const domEnd   = orf.start_position + dom.stop  * 3;
      const dx1      = PAD_L + SCALE(Math.min(domStart, domEnd));
      const dx2      = PAD_L + SCALE(Math.max(domStart, domEnd));
      const dw       = Math.max(dx2 - dx1, 3);
      const domCol   = _domainColor(dom.target);

      g.append('rect')
        .attr('x', dx1).attr('width', dw)
        .attr('y', y + ORF_H)
        .attr('height', DOM_H)
        .attr('rx', 2)
        .attr('fill', domCol)
        .attr('opacity', .9)
        .on('mousemove', evt => {
          vqTooltipShow(`
            <div class="vq-tooltip__title">${dom.target}</div>
            <div class="vq-tooltip__row">
              <span class="vq-tooltip__key">DB</span><span>${dom.database}</span>
              <span class="vq-tooltip__key">Description</span><span>${dom.description}</span>
              <span class="vq-tooltip__key">Score</span><span>${dom.score}</span>
              <span class="vq-tooltip__key">E-value</span><span>${dom.e_value?.toExponential(2)}</span>
            </div>`, evt);
        })
        .on('mouseleave', vqTooltipHide);
    });
  });

  // --- Legend ---
  const legendG = svg.append('g').attr('transform', `translate(${PAD_L}, ${totalH - 20})`);
  legendG.append('rect').attr('width', 12).attr('height', 10).attr('rx', 2)
    .attr('fill', 'var(--vq-accent)');
  legendG.append('text').attr('x', 16).attr('y', 9)
    .style('font-size', '9px').style('fill', 'var(--vq-text-3)').text('(+) ORF');
  legendG.append('rect').attr('x', 60).attr('width', 12).attr('height', 10).attr('rx', 2)
    .attr('fill', 'var(--vq-warning)');
  legendG.append('text').attr('x', 76).attr('y', 9)
    .style('font-size', '9px').style('fill', 'var(--vq-text-3)').text('(−) ORF');
  legendG.append('text').attr('x', 120).attr('y', 9)
    .style('font-size', '9px').style('fill', 'var(--vq-text-3)').text('■ HMM domain');

  return svg.node();
}

// ── Lane assignment (greedy interval packing) ────────────────────────────────

function _assignLanes(orfs) {
  const sorted = [...orfs].sort((a, b) => a.start_position - b.start_position);
  const lanes  = [];   // each lane: array of [start, stop]

  return sorted.map(orf => {
    let lane = 0;
    while (true) {
      if (!lanes[lane]) lanes[lane] = [];
      const overlaps = lanes[lane].some(([s, e]) =>
        orf.start_position < e && orf.stop_position > s
      );
      if (!overlaps) {
        lanes[lane].push([orf.start_position, orf.stop_position]);
        return { ...orf, lane };
      }
      lane++;
    }
  });
}

// ── Domain colour (deterministic per target name) ────────────────────────────

function _domainColor(target) {
  if (_VW.domainMap[target]) return _VW.domainMap[target];
  const col = _DOM_COLOURS[_VW.colourIdx % _DOM_COLOURS.length];
  _VW.colourIdx++;
  _VW.domainMap[target] = col;
  return col;
}
