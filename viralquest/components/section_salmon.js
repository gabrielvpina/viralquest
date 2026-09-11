/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_salmon.js — Section 5: Salmon Quantification

   Two pathways:
     reference  — viral + bundled HK + ref-HK + transcriptome
     de_novo    — assembled contigs with viral/HK-matched prefixes

   Left card:   viral sequences — boxplot if clustered, dot if singleton.
   Right card:  housekeeping genes (reference) OR HK-matched contigs (de novo).
   Ref-HK card: reference housekeeping genes (reference pathway + --hk-genes only).
   Toggle:      TPM ↔ NumReads.
   Toggle:      overlay per-kingdom housekeeping maxima in the main plot.
   ============================================================ */

const _SQ = {
  data:      null,
  clusters:  null,
  metric:    'tpm',
  pathway:   'reference',
  showHK:    false,
  groupMode: 'cluster',
  minVal:    0,
  maxVal:    null,
  resizeT:   null,
};

// Display helpers for kingdom names coming from conserved_quant entries.
// Keys are the uppercase strings stored in the JSON (VQ_CONS_{KINGDOM}_*).
const _KG_LABEL = {
  MAMMALS:    'Mammals',
  ARTHROPODS: 'Arthropods',
  PLANTS:     'Plants',
  FISH:       'Fish',
  FUNGI:      'Fungi',
  BACTERIA:   'Bacteria',
  NEMATODES:  'Nematodes',
  HUMAN:      'Human',
  MOUSE:      'Mouse',
};

// Short abbreviations used on the viral-panel overlay (≤4 chars).
const _KG_ABBR = {
  MAMMALS:    'MAM',
  ARTHROPODS: 'ART',
  PLANTS:     'PLT',
  FISH:       'FSH',
  FUNGI:      'FNG',
  BACTERIA:   'BAC',
  NEMATODES:  'NEM',
  HUMAN:      'HSA',   // Homo sapiens
  MOUSE:      'MMU',   // Mus musculus
};

function _kgLabel(k) { return _KG_LABEL[k] || (k.charAt(0) + k.slice(1).toLowerCase()); }
function _kgAbbr(k)  { return _KG_ABBR[k]  || k.slice(0, 3).toUpperCase(); }

function vqInitSalmon(salmonQuant, clusters) {
  const el = document.getElementById('section-salmon');
  if (!el) return;

  _SQ.data     = salmonQuant;
  _SQ.clusters = clusters || [];
  _SQ.pathway  = salmonQuant.pathway || 'reference';

  const isRef      = _SQ.pathway === 'reference';
  const hkTitle    = isRef ? 'Housekeeping Genes' : 'HK-Matched Contigs';
  const refHkData  = salmonQuant.ref_hk_quant || [];
  const showRefHk  = isRef && refHkData.length > 0;

  const pathwayColor = isRef
    ? 'background:#dbeafe;color:#1d4ed8'
    : 'background:#dcfce7;color:#166534';
  const pathwayLabel = isRef ? 'Reference' : 'De Novo';
  const pathwayBadge = `<span style="display:inline-block;padding:2px 8px;border-radius:10px;`
    + `font-size:11px;font-weight:600;vertical-align:middle;${pathwayColor};margin-left:8px;">`
    + `${pathwayLabel}</span>`;

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Salmon Quantification ${pathwayBadge}</div>
        <div class="vq-section-sub">
          Mapping rate: <strong>${(salmonQuant.mapping_rate ?? 0).toFixed(1)}%</strong>
          &nbsp;·&nbsp; ${(salmonQuant.total_reads ?? 0).toLocaleString()} total reads
        </div>
      </div>
      <div class="vq-section-actions">
        <div class="vq-toggle" role="tablist" aria-label="Group mode">
          <button class="vq-toggle__btn active"
                  id="sq-btn-cluster" data-group="cluster" type="button"
                  role="tab" aria-selected="true">By Cluster</button>
          <button class="vq-toggle__btn"
                  id="sq-btn-individual" data-group="individual" type="button"
                  role="tab" aria-selected="false">By Sequence</button>
        </div>
        <div class="vq-toggle" role="tablist" aria-label="Metric">
          <button class="vq-toggle__btn ${_SQ.metric==='tpm'?'active':''}"
                  id="sq-btn-tpm" data-metric="tpm" type="button"
                  role="tab" aria-selected="${_SQ.metric==='tpm'}">TPM</button>
          <button class="vq-toggle__btn ${_SQ.metric==='num_reads'?'active':''}"
                  id="sq-btn-reads" data-metric="num_reads" type="button"
                  role="tab" aria-selected="${_SQ.metric==='num_reads'}">Reads</button>
        </div>
        <label class="vq-checkbox-label" style="display:inline-flex;align-items:center;gap:6px;
               font-size:12px;color:var(--vq-text-2);cursor:pointer">
          <input type="checkbox" class="vq-checkbox" id="sq-show-hk">
          Overlay housekeeping
        </label>
        <div class="vq-menu" id="sq-export-menu">
          <button class="vq-btn vq-btn--sm vq-btn--ghost" data-menu-toggle type="button">
            Export
            <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="2.5">
              <polyline points="6 9 12 15 18 9"/>
            </svg>
          </button>
          <div class="vq-menu__panel" role="menu">
            <button class="vq-menu__item" data-fmt="png" role="menuitem" type="button">PNG image</button>
            <button class="vq-menu__item" data-fmt="svg" role="menuitem" type="button">SVG vector</button>
          </div>
        </div>
      </div>
    </div>

    <div class="vq-toolbar__filter-row" style="margin-bottom:var(--vq-space-3)">
      <div class="vq-filter-group">
        <span class="vq-filter-label" id="sq-min-label">Min TPM</span>
        <input class="vq-input vq-input--sm" type="number" id="sq-min-val"
               min="0" step="0.01" placeholder="0"
               aria-label="Minimum value threshold">
        <span class="vq-filter-sub">— hide rows below this value</span>
      </div>
      <div class="vq-filter-group">
        <span class="vq-filter-label" id="sq-max-label">Max TPM</span>
        <input class="vq-input vq-input--sm" type="number" id="sq-max-val"
               min="0" step="0.01" placeholder="—"
               aria-label="Maximum value threshold">
        <span class="vq-filter-sub">— hide rows above this value &amp; cap scale</span>
      </div>
    </div>

    <div style="display:grid;grid-template-columns:1fr 320px;gap:var(--vq-space-3);
                align-items:stretch" id="sq-grid">
      <div class="vq-card" style="display:flex;flex-direction:column">
        <div class="vq-card__header">
          <div class="vq-card__title">Viral Sequences</div>
          <span class="vq-panel__count" id="sq-viral-count">0</span>
        </div>
        <div class="vq-card__body" style="padding:12px 14px;flex:1">
          <div id="sq-viral-wrap" style="width:100%"></div>
        </div>
      </div>
      <div class="vq-card" style="display:flex;flex-direction:column">
        <div class="vq-card__header">
          <div class="vq-card__title">${hkTitle}</div>
          <span class="vq-panel__count" id="sq-cons-count">0</span>
        </div>
        <div class="vq-card__body" id="sq-cons-wrap" style="padding:12px 14px;flex:1"></div>
      </div>
    </div>

    ${showRefHk ? `
    <div class="vq-card" style="margin-top:var(--vq-space-3)">
      <div class="vq-card__header">
        <div class="vq-card__title">Reference Housekeeping Genes</div>
        <span class="vq-panel__count" id="sq-refhk-count">0</span>
      </div>
      <div class="vq-card__body" id="sq-refhk-wrap" style="padding:12px 14px"></div>
    </div>` : ''}

    ${_SQ.pathway === 'reference' ? _renderHostHits(salmonQuant.host_viral_hits || []) : ''}

    ${(salmonQuant.pfam_hk_quant || []).length ? `
    <div class="vq-card" style="margin-top:var(--vq-space-3)">
      <div class="vq-card__header">
        <div class="vq-card__title">Housekeeping Domains (Pfam)</div>
        <span class="vq-panel__count" id="sq-pfamhk-count">0</span>
      </div>
      <div class="vq-card__body" style="padding:12px 14px">
        <div id="sq-pfamhk-wrap" style="width:100%"></div>
      </div>
    </div>` : ''}
  `;

  // Min-value filter
  document.getElementById('sq-min-val')?.addEventListener('input', e => {
    const v = parseFloat(e.target.value);
    _SQ.minVal = isNaN(v) || v < 0 ? 0 : v;
    _drawViralPanel();
    _drawPfamHkPanel();
  });

  // Max-value filter
  document.getElementById('sq-max-val')?.addEventListener('input', e => {
    const v = parseFloat(e.target.value);
    _SQ.maxVal = isNaN(v) || v <= 0 ? null : v;
    _drawViralPanel();
    _drawPfamHkPanel();
  });

  // Wire group-mode toggle
  document.querySelectorAll('[data-group]').forEach(btn => {
    btn.addEventListener('click', () => _setGroupMode(btn.dataset.group));
  });

  // Wire metric toggle
  document.querySelectorAll('[data-metric]').forEach(btn => {
    btn.addEventListener('click', () => _setMetric(btn.dataset.metric));
  });

  // Overlay toggle
  document.getElementById('sq-show-hk')?.addEventListener('change', e => {
    _SQ.showHK = e.target.checked;
    _drawViralPanel();
  });

  // Export menu
  const exportMenu = document.getElementById('sq-export-menu');
  const exportTog  = exportMenu?.querySelector('[data-menu-toggle]');
  exportTog?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => {
      if (m !== exportMenu) m.classList.remove('open');
    });
    exportMenu.classList.toggle('open');
  });
  exportMenu?.querySelectorAll('[data-fmt]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      exportMenu.classList.remove('open');
      const svg = document.querySelector('#sq-viral-wrap svg');
      if (!svg) return;
      if (item.dataset.fmt === 'png') VQ.exportPNG(svg, 'salmon_viral.png');
      else                            VQ.exportSVG(svg, 'salmon_viral.svg');
    });
  });

  // Responsive grid
  const grid  = document.getElementById('sq-grid');
  const apply = () => {
    if (window.innerWidth < 980) grid.style.gridTemplateColumns = '1fr';
    else                         grid.style.gridTemplateColumns = '1fr 320px';
  };
  apply();
  window.addEventListener('resize', () => {
    apply();
    clearTimeout(_SQ.resizeT);
    _SQ.resizeT = setTimeout(_draw, 150);
  });

  // ResizeObserver — re-draw once the wraps actually have real widths
  let lastW = 0;
  const ro = new ResizeObserver(() => {
    const w = document.getElementById('sq-viral-wrap')?.clientWidth || 0;
    if (!w || w === lastW) return;
    lastW = w;
    clearTimeout(_SQ.resizeT);
    _SQ.resizeT = setTimeout(_draw, 60);
  });
  const viralWrap = document.getElementById('sq-viral-wrap');
  if (viralWrap) ro.observe(viralWrap);

  _draw();
}

function _setGroupMode(mode) {
  _SQ.groupMode = mode;
  document.querySelectorAll('[data-group]').forEach(btn => {
    const active = btn.dataset.group === mode;
    btn.classList.toggle('active', active);
    btn.setAttribute('aria-selected', String(active));
  });
  _drawViralPanel();
  _drawPfamHkPanel();
}

function _setMetric(metric) {
  _SQ.metric = metric;
  document.querySelectorAll('[data-metric]').forEach(btn => {
    const active = btn.dataset.metric === metric;
    btn.classList.toggle('active', active);
    btn.setAttribute('aria-selected', String(active));
  });
  const lbl  = document.getElementById('sq-min-label');
  const inp  = document.getElementById('sq-min-val');
  const lbl2 = document.getElementById('sq-max-label');
  const inp2 = document.getElementById('sq-max-val');
  if (lbl)  lbl.textContent  = metric === 'tpm' ? 'Min TPM' : 'Min Reads';
  if (lbl2) lbl2.textContent = metric === 'tpm' ? 'Max TPM' : 'Max Reads';
  if (inp)  inp.step  = metric === 'tpm' ? '0.01' : '1';
  if (inp2) inp2.step = metric === 'tpm' ? '0.01' : '1';
  _draw();
}

function _draw() {
  _drawViralPanel();
  _drawConsPanel();
  if (_SQ.pathway === 'reference') _drawRefHkPanel();
  _drawPfamHkPanel();
}

// ── Housekeeping summary per kingdom ───────────────────────────────────────
//
// The conserved panel is dominated by zeros: every kingdom ships thousands of
// genes and only the handful belonging to the actual host ever collects reads
// (a typical arthropod library detects ~10 of ~1700 ARTHROPODS genes). A plain
// median over every gene is therefore 0.0 for *every* kingdom, which collapsed
// both this panel and the viral-panel overlay onto x=0.
//
// So: rank kingdoms by `sum`, never by a measure of position. A lone stray gene
// at high TPM (one FISH transcript at ~31 TPM) would otherwise outrank the real
// host (ARTHROPODS, ~0.6 median but ~2500 TPM total).
//
// The viral-panel overlay marks each kingdom's `max` — the single most
// expressed housekeeping gene. That is the line a viral contig has to clear to
// be "more transcribed than the host's housekeeping", which is the comparison
// the overlay exists to make; a median sits far below it and understates the
// bar. `median` is kept for the tooltip, computed over detected genes only.

function _hkSummary() {
  const cons   = (_SQ.data.conserved_quant || []);
  const metric = _SQ.metric;
  const byKingdom = {};
  cons.forEach(c => {
    const k = c.kingdom || 'Unknown';
    (byKingdom[k] ??= []).push(c[metric] ?? 0);
  });
  return Object.keys(byKingdom).map(k => {
    const values = byKingdom[k];
    const nz     = values.filter(v => v > 0).sort(d3.ascending);
    return {
      kingdom:  k,
      values,
      detected: nz.length,
      total:    values.length,
      sum:      d3.sum(values),
      max:      d3.max(values) ?? 0,
      median:   nz.length ? (d3.quantile(nz, 0.5) ?? 0) : 0,
    };
  }).sort((a, b) => b.sum - a.sum);
}

// Kingdoms worth drawing a reference line for — one with nothing detected has
// no expression level to mark, and would just stack another line on x=0.
const _hkOverlayRows = () => _hkSummary().filter(h => h.detected > 0);

// Kingdom totals span 0 → thousands; keep small values readable without
// letting the big ones run past the label gutter.
function _fmtVal(v) {
  if (v == null || isNaN(v)) return '—';
  if (v === 0)    return '0';
  if (v >= 1000)  return v.toFixed(0);
  if (v >= 10)    return v.toFixed(1);
  return v.toFixed(2);
}

// ── Viral panel ────────────────────────────────────────────────────────────

function _drawViralPanel() {
  const wrap = document.getElementById('sq-viral-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const viral  = (_SQ.data.viral_quant || []);
  const metric = _SQ.metric;
  const label  = metric === 'tpm' ? 'TPM' : 'Reads';

  // Apply min/max-value filters
  const minVal    = _SQ.minVal || 0;
  const maxFilter = _SQ.maxVal;
  const viralFiltered = viral.filter(v => {
    const val = v[metric] ?? 0;
    return val >= minVal && (maxFilter == null || val <= maxFilter);
  });

  const countEl = document.getElementById('sq-viral-count');
  if (countEl) {
    const filtered = minVal > 0 || maxFilter != null;
    countEl.textContent = filtered
      ? `${viralFiltered.length} of ${viral.length} sequence${viral.length !== 1 ? 's' : ''}`
      : `${viral.length} sequence${viral.length !== 1 ? 's' : ''}`;
  }

  if (!viralFiltered.length) {
    wrap.innerHTML = viral.length
      ? `<div class="vq-empty">All ${viral.length} sequence${viral.length !== 1 ? 's' : ''} hidden by the active ${label} filter.</div>`
      : '<div class="vq-empty">No viral quantification data.</div>';
    return;
  }

  // Build groups depending on view mode
  const groups = {};

  if (_SQ.groupMode === 'individual') {
    // Every viral entry is its own row
    viralFiltered.forEach(v => {
      const key = v.name;
      groups[key] = { label: key, values: [v[metric] ?? 0], seqIds: [v.name] };
    });
  } else {
    // Group by cluster_id when available, fall back to individual
    const seqCluster = {};
    (_SQ.clusters || []).forEach(c => {
      c.members.forEach(m => { seqCluster[m.seq_id] = c.cluster_id; });
    });
    viralFiltered.forEach(v => {
      const cid = seqCluster[v.name];
      const key = cid != null ? 'cluster_' + cid : v.name;
      (groups[key] ??= { label: key, values: [], seqIds: [] }).values.push(v[metric] ?? 0);
      groups[key].seqIds.push(v.name);
    });
  }

  // Sort by median descending
  const keys = Object.keys(groups).sort((a, b) => {
    const ma = d3.quantile(groups[a].values.slice().sort(d3.ascending), 0.5) ?? 0;
    const mb = d3.quantile(groups[b].values.slice().sort(d3.ascending), 0.5) ?? 0;
    return mb - ma;
  });

  // Geometry
  const containerW = wrap.clientWidth || 700;
  const W     = Math.max(containerW, 560);
  const PAD_L = 200;
  const PAD_R = 30;
  const PAD_T = 26;
  const ROW_H = 30;
  const drawW = W - PAD_L - PAD_R;

  const HK = _SQ.showHK ? _hkOverlayRows() : [];

  const allVals = viralFiltered.map(v => v[metric] ?? 0).slice();
  HK.forEach(h => allVals.push(h.max));
  const dataMax  = d3.max(allVals) || 1;
  const scaleMax = maxFilter != null ? maxFilter : dataMax;

  const xScale = d3.scaleLinear([0, scaleMax], [0, drawW]);

  const overlayH = HK.length ? 24 : 0;
  const PAD_B   = 60 + overlayH;
  const H       = PAD_T + keys.length * ROW_H + PAD_B;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // X axis
  const gridH = H - PAD_T - PAD_B;
  const axG = svg.append('g')
    .attr('transform', `translate(${PAD_L},${PAD_T})`)
    .call(d3.axisTop(xScale).ticks(Math.max(4, Math.round(W / 140))).tickSize(0));
  axG.selectAll('.tick line')
    .attr('y1', 0).attr('y2', gridH)
    .attr('stroke', '#dde3ec').attr('stroke-dasharray', '3,3');
  axG.select('.domain').remove();
  axG.selectAll('.tick text')
    .style('font-size', '10px').style('fill', '#94a3b8');

  // X axis label
  svg.append('text')
    .attr('x', PAD_L + drawW / 2)
    .attr('y', H - 12 - overlayH)
    .attr('text-anchor', 'middle')
    .attr('font-size', 11)
    .attr('fill', 'var(--vq-text-2)')
    .text(label);

  // Rows
  keys.forEach((key, i) => {
    const grp  = groups[key];
    const vals = grp.values.slice().sort(d3.ascending);
    const y    = PAD_T + i * ROW_H + ROW_H / 2;
    const isSingleton = vals.length === 1;
    const gEl  = svg.append('g').attr('transform', `translate(${PAD_L},0)`);

    const fullLabel = grp.label;
    const truncated = fullLabel.length > 28 ? fullLabel.slice(0, 27) + '…' : fullLabel;
    const lbl = svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + 4)
      .attr('text-anchor', 'end').attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)').attr('font-family', 'var(--vq-font-mono)')
      .text(truncated);
    if (truncated !== fullLabel) {
      lbl.style('cursor', 'help')
        .on('mousemove', evt => VQ.tooltipShow(
          `<div class="vq-tooltip__title">${VQ.esc(fullLabel)}</div>
           <div class="vq-tooltip__row">
             <span class="vq-tooltip__key">Members</span><span>${vals.length}</span>
           </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    }

    if (isSingleton) {
      const cx = xScale(vals[0]);
      gEl.append('circle')
        .attr('cx', cx).attr('cy', y).attr('r', 5)
        .attr('fill', 'var(--vq-accent)')
        .attr('stroke', '#fff').attr('stroke-width', 1)
        .attr('cursor', 'pointer')
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(grp.seqIds[0])}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">${label}</span><span>${vals[0].toFixed(2)}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide)
        .on('click', () => VQ.jumpToViewer(grp.seqIds[0]));
    } else {
      const q1  = d3.quantile(vals, 0.25);
      const med = d3.quantile(vals, 0.5);
      const q3  = d3.quantile(vals, 0.75);
      const iqr = q3 - q1;
      const lo  = Math.max(d3.min(vals), q1 - 1.5 * iqr);
      const hi  = Math.min(d3.max(vals), q3 + 1.5 * iqr);
      const bh  = 14;

      gEl.append('line')
        .attr('x1', xScale(lo)).attr('x2', xScale(hi))
        .attr('y1', y).attr('y2', y)
        .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5);
      [lo, hi].forEach(x => {
        gEl.append('line')
          .attr('x1', xScale(x)).attr('x2', xScale(x))
          .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
          .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5);
      });
      gEl.append('rect')
        .attr('x', xScale(q1)).attr('width', Math.max(xScale(q3) - xScale(q1), 2))
        .attr('y', y - bh / 2).attr('height', bh)
        .attr('rx', 2)
        .attr('fill', 'var(--vq-accent)').attr('opacity', 0.25)
        .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(grp.label)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Median ${label}</span><span>${med.toFixed(2)}</span>
            <span class="vq-tooltip__key">Q1–Q3</span><span>${q1.toFixed(2)}–${q3.toFixed(2)}</span>
            <span class="vq-tooltip__key">Members</span><span>${vals.length}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
      gEl.append('line')
        .attr('x1', xScale(med)).attr('x2', xScale(med))
        .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
        .attr('stroke', 'var(--vq-primary)').attr('stroke-width', 2);

      vals.filter(v => v < lo || v > hi).forEach(v => {
        gEl.append('circle')
          .attr('cx', xScale(v)).attr('cy', y).attr('r', 3)
          .attr('fill', 'var(--vq-danger)').attr('opacity', 0.7);
      });
    }
  });

  // Housekeeping overlay: vertical reference lines + bottom legend chips
  if (HK.length) {
    const overlayY = H - PAD_B + 36;
    HK.forEach((h, i) => {
      const x = PAD_L + xScale(h.max);
      svg.append('line')
        .attr('x1', x).attr('x2', x)
        .attr('y1', PAD_T).attr('y2', H - PAD_B)
        .attr('stroke', 'var(--vq-success)')
        .attr('stroke-dasharray', '4,3')
        .attr('stroke-width', 1.2)
        .attr('opacity', 0.7)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(_kgLabel(h.kingdom))} housekeeping</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Max ${label}</span><span>${h.max.toFixed(2)}</span>
            <span class="vq-tooltip__key">Median ${label}</span><span>${h.median.toFixed(2)}</span>
            <span class="vq-tooltip__key">Total ${label}</span><span>${h.sum.toFixed(2)}</span>
            <span class="vq-tooltip__key">Detected</span><span>${h.detected} / ${h.total}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);

      svg.append('text')
        .attr('x', x).attr('y', PAD_T - 8)
        .attr('text-anchor', 'middle')
        .attr('font-size', 9.5)
        .attr('fill', 'var(--vq-success)')
        .text(_kgAbbr(h.kingdom));
    });

    svg.append('text')
      .attr('x', PAD_L).attr('y', overlayY)
      .attr('font-size', 11).attr('fill', 'var(--vq-text-2)')
      .text('Housekeeping (max ' + label + '):');
    let lx = PAD_L + 200;
    HK.forEach(h => {
      svg.append('line')
        .attr('x1', lx).attr('x2', lx + 18)
        .attr('y1', overlayY - 4).attr('y2', overlayY - 4)
        .attr('stroke', 'var(--vq-success)')
        .attr('stroke-dasharray', '4,3').attr('stroke-width', 1.5);
      svg.append('text')
        .attr('x', lx + 22).attr('y', overlayY)
        .attr('font-size', 11)
        .attr('fill', 'var(--vq-text-3)')
        .text(`${_kgLabel(h.kingdom)} (${h.max.toFixed(1)})`);
      lx += 22 + (h.kingdom.length + 8) * 6.2;
    });
  }

  wrap.appendChild(svg.node());
}

// ── Housekeeping / HK-matched contigs panel ────────────────────────────────

function _drawConsPanel() {
  const wrap = document.getElementById('sq-cons-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const cons   = (_SQ.data.conserved_quant || []);
  const metric = _SQ.metric;
  const label  = metric === 'tpm' ? 'TPM' : 'Reads';

  const kingdoms  = _hkSummary();
  const detected  = d3.sum(kingdoms, k => k.detected);

  const countEl = document.getElementById('sq-cons-count');
  if (countEl) countEl.textContent = cons.length
    ? `${detected} / ${cons.length} detected`
    : '0 genes';

  if (!cons.length) {
    const msg = _SQ.pathway === 'de_novo'
      ? 'No HK-matched contigs found.'
      : 'No housekeeping data.';
    wrap.innerHTML = `<div class="vq-empty" style="font-size:12px">${msg}</div>`;
    return;
  }

  const maxVal = d3.max(kingdoms, x => x.sum) || 1;
  const W      = Math.max(wrap.clientWidth || 280, 250);
  const PAD_L  = 96;
  const PAD_R  = 56;
  const ROW_H  = 34;
  const BAR_H  = 12;
  const H      = kingdoms.length * ROW_H + 36;

  // Symlog, not linear: the host kingdom outweighs the incidental ones by three
  // orders of magnitude, and a linear axis flattens every non-host kingdom to
  // nothing. Symlog (not log) because kingdoms with zero detected genes are the
  // common case and must still land on a real zero.
  const xScale = d3.scaleSymlog([0, maxVal], [0, W - PAD_L - PAD_R]);

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  svg.append('text')
    .attr('x', PAD_L).attr('y', 14)
    .attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
    .text(`Total ${label} per kingdom · log`);

  kingdoms.forEach((row, i) => {
    const y    = 26 + i * ROW_H;
    const isOn = row.detected > 0;

    svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + BAR_H / 2 + 3.5)
      .attr('text-anchor', 'end')
      .attr('font-size', 11)
      .attr('fill', isOn ? 'var(--vq-text-2)' : 'var(--vq-text-3)')
      .text(_kgLabel(row.kingdom));

    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', W - PAD_L - PAD_R).attr('height', BAR_H)
      .attr('rx', 2)
      .attr('fill', 'var(--vq-bg)');

    // A kingdom with nothing detected draws no bar at all — an undetected
    // kingdom is the signal here, so it must not look like a small measurement.
    const bw = isOn ? Math.max(xScale(row.sum), 2) : 0;
    if (bw) {
      svg.append('rect')
        .attr('x', PAD_L).attr('y', y)
        .attr('width', bw).attr('height', BAR_H)
        .attr('rx', 2)
        .attr('fill', 'var(--vq-success)').attr('opacity', 0.85);
    }

    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', W - PAD_L - PAD_R).attr('height', BAR_H)
      .attr('fill', 'transparent')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(_kgLabel(row.kingdom))}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Total ${label}</span><span>${_fmtVal(row.sum)}</span>
          <span class="vq-tooltip__key">Detected</span><span>${row.detected} / ${row.total}</span>
          <span class="vq-tooltip__key">Median ${label}</span><span>${isOn ? _fmtVal(row.median) : '—'}</span>
          <span class="vq-tooltip__key">Max ${label}</span><span>${_fmtVal(row.max)}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    svg.append('text')
      .attr('x', W - PAD_R + 6).attr('y', y + BAR_H / 2 + 3.5)
      .attr('font-size', 10)
      .attr('fill', isOn ? 'var(--vq-text-2)' : 'var(--vq-text-3)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(_fmtVal(row.sum));

    svg.append('text')
      .attr('x', PAD_L).attr('y', y + BAR_H + 11)
      .attr('font-size', 9)
      .attr('fill', 'var(--vq-text-3)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(`${row.detected}/${row.total} detected`);
  });

  wrap.appendChild(svg.node());
}

// ── Reference HK genes panel (reference pathway + --hk-genes only) ─────────

function _drawRefHkPanel() {
  const wrap = document.getElementById('sq-refhk-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const refHk  = (_SQ.data.ref_hk_quant || []);
  const metric = _SQ.metric;
  const label  = metric === 'tpm' ? 'TPM' : 'Reads';

  const countEl = document.getElementById('sq-refhk-count');
  if (countEl) countEl.textContent = `${refHk.length} gene${refHk.length !== 1 ? 's' : ''}`;

  if (!refHk.length) {
    wrap.innerHTML = '<div class="vq-empty" style="font-size:12px">No reference HK data.</div>';
    return;
  }

  const sorted = [...refHk].sort((a, b) => (b[metric] ?? 0) - (a[metric] ?? 0));
  const maxVal = (sorted[0]?.[metric] ?? 0) || 1;

  const W     = Math.max(wrap.clientWidth || 600, 400);
  const PAD_L = 180;
  const PAD_R = 60;
  const ROW_H = 24;
  const BAR_H = 10;
  const H     = sorted.length * ROW_H + 36;

  const xScale = d3.scaleLinear([0, maxVal], [0, W - PAD_L - PAD_R]);

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  svg.append('text')
    .attr('x', PAD_L).attr('y', 14)
    .attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
    .text(`${label} per reference HK gene`);

  sorted.forEach((entry, i) => {
    const y     = 26 + i * ROW_H;
    const val   = entry[metric] ?? 0;
    const gName = entry.name.replace(/^VQ_REFHK_/, '');
    const truncated = gName.length > 26 ? gName.slice(0, 25) + '…' : gName;

    const lbl = svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + BAR_H / 2 + 3.5)
      .attr('text-anchor', 'end')
      .attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)')
      .attr('font-family', 'var(--vq-font-mono)')
      .text(truncated);
    if (truncated !== gName) {
      lbl.style('cursor', 'help')
        .on('mousemove', evt => VQ.tooltipShow(
          `<div class="vq-tooltip__title">${VQ.esc(gName)}</div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    }

    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', W - PAD_L - PAD_R).attr('height', BAR_H)
      .attr('rx', 2).attr('fill', 'var(--vq-bg)');

    const bw = Math.max(xScale(val), val > 0 ? 2 : 0);
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', bw).attr('height', BAR_H)
      .attr('rx', 2)
      .attr('fill', 'var(--vq-primary)').attr('opacity', 0.7)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(gName)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">${label}</span><span>${val.toFixed(2)}</span>
          <span class="vq-tooltip__key">Reads</span><span>${(entry.num_reads ?? 0).toFixed(0)}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    svg.append('text')
      .attr('x', PAD_L + bw + 6).attr('y', y + BAR_H / 2 + 3.5)
      .attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(val.toFixed(1));
  });

  wrap.appendChild(svg.node());
}

// ── Pfam HK panel ──────────────────────────────────────────────────────────

function _stripHtml(html) {
  return (html || '').replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim();
}

function _drawPfamHkPanel() {
  const wrap = document.getElementById('sq-pfamhk-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const raw    = (_SQ.data.pfam_hk_quant || []);
  const metric = _SQ.metric;
  const label  = metric === 'tpm' ? 'TPM' : 'Reads';

  const minVal    = _SQ.minVal || 0;
  const maxFilter = _SQ.maxVal;
  const filtered  = raw.filter(e => {
    const val = e[metric] ?? 0;
    return val >= minVal && (maxFilter == null || val <= maxFilter);
  });

  const countEl = document.getElementById('sq-pfamhk-count');
  if (countEl) {
    const isFiltered = minVal > 0 || maxFilter != null;
    countEl.textContent = isFiltered
      ? `${filtered.length} of ${raw.length} domain${raw.length !== 1 ? 's' : ''}`
      : `${raw.length} domain${raw.length !== 1 ? 's' : ''}`;
  }

  if (!filtered.length) return;

  // Build groups: cluster mode → by pfam_target; individual → each seq alone
  const groups = {};
  if (_SQ.groupMode === 'individual') {
    filtered.forEach(e => {
      groups[e.seq_id] = {
        label:       e.seq_id,
        values:      [e[metric] ?? 0],
        pfam_target: e.pfam_target,
        pfam_desc:   e.pfam_desc,
        pfam_details:e.pfam_details,
        seqIds:      [e.seq_id],
      };
    });
  } else {
    filtered.forEach(e => {
      const key = e.pfam_target;
      if (!groups[key]) {
        groups[key] = {
          label:       key,
          values:      [],
          pfam_target: e.pfam_target,
          pfam_desc:   e.pfam_desc,
          pfam_details:e.pfam_details,
          seqIds:      [],
        };
      }
      groups[key].values.push(e[metric] ?? 0);
      groups[key].seqIds.push(e.seq_id);
    });
  }

  const keys = Object.keys(groups).sort((a, b) => {
    const ma = d3.quantile(groups[a].values.slice().sort(d3.ascending), 0.5) ?? 0;
    const mb = d3.quantile(groups[b].values.slice().sort(d3.ascending), 0.5) ?? 0;
    return mb - ma;
  });

  const containerW = wrap.clientWidth || 700;
  const W     = Math.max(containerW, 560);
  const PAD_L = 200;
  const PAD_R = 30;
  const PAD_T = 26;
  const ROW_H = 30;
  const PAD_B = 48;
  const drawW = W - PAD_L - PAD_R;

  const allVals  = filtered.map(e => e[metric] ?? 0);
  const dataMax  = d3.max(allVals) || 1;
  const scaleMax = maxFilter != null ? maxFilter : dataMax;
  const xScale   = d3.scaleLinear([0, scaleMax], [0, drawW]);
  const H       = PAD_T + keys.length * ROW_H + PAD_B;

  const COLOR = 'var(--vq-success)';

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // X axis
  const gridH = H - PAD_T - PAD_B;
  const axG = svg.append('g')
    .attr('transform', `translate(${PAD_L},${PAD_T})`)
    .call(d3.axisTop(xScale).ticks(Math.max(4, Math.round(W / 140))).tickSize(0));
  axG.selectAll('.tick line')
    .attr('y1', 0).attr('y2', gridH)
    .attr('stroke', '#dde3ec').attr('stroke-dasharray', '3,3');
  axG.select('.domain').remove();
  axG.selectAll('.tick text')
    .style('font-size', '10px').style('fill', '#94a3b8');

  svg.append('text')
    .attr('x', PAD_L + drawW / 2).attr('y', H - 12)
    .attr('text-anchor', 'middle').attr('font-size', 11)
    .attr('fill', 'var(--vq-text-2)').text(label);

  keys.forEach((key, i) => {
    const grp  = groups[key];
    const vals = grp.values.slice().sort(d3.ascending);
    const y    = PAD_T + i * ROW_H + ROW_H / 2;
    const gEl  = svg.append('g').attr('transform', `translate(${PAD_L},0)`);

    // Label
    const truncated = grp.label.length > 28 ? grp.label.slice(0, 27) + '…' : grp.label;
    const detailText = _stripHtml(grp.pfam_details);
    svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + 4)
      .attr('text-anchor', 'end').attr('font-size', 11)
      .attr('fill', 'var(--vq-success)').attr('font-family', 'var(--vq-font-mono)')
      .attr('cursor', 'help')
      .text(truncated)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(grp.pfam_target)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Description</span>
          <span>${VQ.esc(grp.pfam_desc)}</span>
        </div>
        ${detailText ? `<div style="margin-top:6px;padding-top:6px;border-top:1px solid rgba(255,255,255,.15);
                                    font-size:10px;color:rgba(255,255,255,.8);max-width:260px;line-height:1.45">
          ${VQ.esc(detailText.slice(0, 240))}${detailText.length > 240 ? '…' : ''}
        </div>` : ''}`, evt))
      .on('mouseleave', VQ.tooltipHide);

    if (vals.length === 1) {
      gEl.append('circle')
        .attr('cx', xScale(vals[0])).attr('cy', y).attr('r', 5)
        .attr('fill', COLOR).attr('stroke', '#fff').attr('stroke-width', 1)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(grp.seqIds[0])}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">${label}</span>
            <span>${vals[0].toFixed(2)}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    } else {
      const q1  = d3.quantile(vals, 0.25);
      const med = d3.quantile(vals, 0.5);
      const q3  = d3.quantile(vals, 0.75);
      const iqr = q3 - q1;
      const lo  = Math.max(d3.min(vals), q1 - 1.5 * iqr);
      const hi  = Math.min(d3.max(vals), q3 + 1.5 * iqr);
      const bh  = 14;

      gEl.append('line')
        .attr('x1', xScale(lo)).attr('x2', xScale(hi))
        .attr('y1', y).attr('y2', y)
        .attr('stroke', COLOR).attr('stroke-width', 1.5);
      [lo, hi].forEach(x => {
        gEl.append('line')
          .attr('x1', xScale(x)).attr('x2', xScale(x))
          .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
          .attr('stroke', COLOR).attr('stroke-width', 1.5);
      });
      gEl.append('rect')
        .attr('x', xScale(q1)).attr('width', Math.max(xScale(q3) - xScale(q1), 2))
        .attr('y', y - bh / 2).attr('height', bh).attr('rx', 2)
        .attr('fill', COLOR).attr('opacity', 0.25)
        .attr('stroke', COLOR).attr('stroke-width', 1.5)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(grp.label)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Median ${label}</span><span>${med.toFixed(2)}</span>
            <span class="vq-tooltip__key">Q1–Q3</span><span>${q1.toFixed(2)}–${q3.toFixed(2)}</span>
            <span class="vq-tooltip__key">Sequences</span><span>${vals.length}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
      gEl.append('line')
        .attr('x1', xScale(med)).attr('x2', xScale(med))
        .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
        .attr('stroke', 'var(--vq-primary)').attr('stroke-width', 2);

      vals.filter(v => v < lo || v > hi).forEach(v => {
        gEl.append('circle')
          .attr('cx', xScale(v)).attr('cy', y).attr('r', 3)
          .attr('fill', 'var(--vq-danger)').attr('opacity', 0.7);
      });
    }
  });

  wrap.appendChild(svg.node());
}

// ── Host–viral hits table ──────────────────────────────────────────────────

function _renderHostHits(hits) {
  if (!hits.length) return '';
  const esc = VQ.esc;

  const rows = hits.map(h => `
    <tr>
      <td class="vq-td--mono">${esc(h.transcript?.name ?? '—')}</td>
      <td class="vq-td--mono">${esc(h.blast_hits?.[0]?.viral_seq_id ?? '—')}</td>
      <td class="vq-td--num">${h.blast_hits?.[0]?.pident?.toFixed(1) ?? '—'}%</td>
      <td class="vq-td--num">${h.blast_hits?.[0]?.qcovhsp ?? '—'}%</td>
      <td class="vq-td--num">${h.transcript?.tpm?.toFixed(2) ?? '—'}</td>
    </tr>`).join('');

  return `
    <div class="vq-card" style="margin-top:var(--vq-space-3)">
      <div class="vq-card__header">
        <div class="vq-card__title">Host–Viral Transcriptome Hits</div>
        <span class="vq-panel__count">${hits.length}</span>
      </div>
      <div class="vq-card__body" style="overflow-x:auto;padding:0">
        <table class="vq-table">
          <thead>
            <tr>
              ${['Host Transcript','Viral Match','Identity','Coverage','TPM']
                .map(h => `<th>${h}</th>`).join('')}
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </div>`;
}


window.vqInitSalmon = vqInitSalmon;
})();
