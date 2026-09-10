/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_quant.js — Section 4: RNA Quantification
   Built on the cross-sample clusters.  Pick clusters (checkboxes):
     · one cluster  → vertical dot plot of each member's TPM
     · ≥2 clusters  → one boxplot per cluster (members overlaid)
   Members from samples without Salmon (TPM absent) are omitted and the
   count is annotated under each cluster.
   ============================================================ */

let _clusters   = [];
let _selected   = new Set();
let _logScale   = false;
let _hmSamples  = [];   // sample summaries carrying quantification, for the heatmap

function _tpmMembers(c) {
  return c.members.filter(m => m.tpm != null && !isNaN(m.tpm));
}
function _quantifiable(c) {
  return _tpmMembers(c).length > 0;
}

function vqInitQuant(clusters, samples) {
  const el = document.getElementById('section-quant');
  if (!el) return;

  _clusters = (clusters || []).filter(_quantifiable);

  // Heatmap columns: every sample that produced quantification. Cluster members
  // carry the only per-sample TPM the boxplot needs, but the heatmap also puts
  // housekeeping genes on the axis, so it needs the sample summaries too.
  const quantSamples = new Set(_clusters.flatMap(c => _tpmMembers(c).map(m => m.sample)));
  _hmSamples = (samples || []).filter(s => s.salmon || quantSamples.has(s.sample));
  _hkGenes   = _collectHkGenes(_hmSamples);

  // Hide the whole section when no cluster carries any TPM (no Salmon anywhere).
  if (!_clusters.length) {
    el.hidden = true;
    document.getElementById('tab-quant')?.style.setProperty('display', 'none');
    return;
  }

  const esc = VQ.esc;
  // Pre-select up to the first 3 clusters so the plot isn't empty on open.
  _selected = new Set(_clusters.slice(0, Math.min(3, _clusters.length)).map(c => c.gid));

  // Heatmap opens populated but not crowded: up to 8 clusters, and the 5 most
  // expressed housekeeping genes as a baseline to read the viral rows against.
  _hmClusters = new Set(_clusters.slice(0, Math.min(8, _clusters.length)).map(c => c.gid));
  _hmGenes = new Set(
    _hkGenes.slice()
      .sort((a, b) => (d3.mean(Object.values(b.bySample)) || 0) - (d3.mean(Object.values(a.bySample)) || 0))
      .slice(0, 5).map(g => g.name)
  );

  const rows = _clusters.map(c => {
    const n = _tpmMembers(c).length;
    const omitted = c.members.length - n;
    return `
      <label class="vr-sample" style="display:flex">
        <input type="checkbox" value="${esc(c.gid)}" ${_selected.has(c.gid) ? 'checked' : ''}>
        <span class="vq-badge vq-badge--cluster">${esc(c.gid)}</span>
        <span class="clu-species">${esc(c.species)}</span>
        <span class="vr-sample__n">${n} TPM${omitted ? ` · ${omitted} no salmon` : ''}</span>
      </label>`;
  }).join('');

  const hkAvailable = _hkGenes.length > 0;

  const hmCluRows = _clusters.map(c => `
      <label class="vr-sample" style="display:flex">
        <input type="checkbox" value="${esc(c.gid)}" ${_hmClusters.has(c.gid) ? 'checked' : ''}>
        <span class="vq-badge vq-badge--cluster">${esc(c.gid)}</span>
        <span class="clu-species">${esc(c.species)}</span>
        <span class="vr-sample__n">${_tpmMembers(c).length} TPM</span>
      </label>`).join('');

  const hmHkRows = _hkGenes.map(g => {
    const mean = d3.mean(Object.values(g.bySample)) || 0;
    return `
      <label class="vr-sample" style="display:flex">
        <input type="checkbox" value="${esc(g.name)}" ${_hmGenes.has(g.name) ? 'checked' : ''}>
        <span class="clu-species" style="font-family:var(--vq-font-mono)">${esc(g.name)}</span>
        <span class="vr-sample__n">${_fmtTpm(mean)} mean</span>
      </label>`;
  }).join('');

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">RNA Quantification</div>
        <div class="vq-section-sub">Salmon TPM per cluster member · select clusters to compare</div>
      </div>
      <div class="vq-section-actions">
        <div class="vq-toggle" role="tablist" aria-label="Y scale">
          <button class="vq-toggle__btn active" type="button" data-scale="lin">Linear</button>
          <button class="vq-toggle__btn"        type="button" data-scale="log">Log</button>
        </div>
      </div>
    </div>

    <div class="vq-stats-page">
      <div class="vq-card" style="margin-bottom:var(--vq-space-4)">
        <div class="vq-card__header">
          <div class="vq-card__title">Clusters</div>
          <button class="vq-btn vq-btn--sm vq-btn--ghost" id="qt-clear" type="button">Clear selection</button>
        </div>
        <div class="vq-card__body">
          <div class="qt-cluster-list" id="qt-list">${rows}</div>
        </div>
      </div>

      <div class="vq-chart-card">
        <div class="vq-chart-card__head">
          <div>
            <div class="vq-chart-card__title">TPM Distribution</div>
            <div class="vq-chart-card__sub" id="qt-sub"></div>
          </div>
          <div class="vq-section-actions">${_exportBtn('qt-export')}</div>
        </div>
        <div class="vq-chart-card__body"><div id="qt-plot" style="width:100%"></div></div>
      </div>

      <div class="vq-section-header" style="margin-top:var(--vq-space-6)">
        <div>
          <div class="vq-section-title">Cluster Heatmap</div>
          <div class="vq-section-sub">
            TPM per cluster across samples${hkAvailable ? ', against the curated housekeeping panel' : ''}
            &nbsp;·&nbsp; pick the rows, the scale and the palette, then export the figure
          </div>
        </div>
      </div>

      <div class="hm-picker">
        <div class="vq-card">
          <div class="vq-card__header">
            <div class="vq-card__title">Clusters</div>
            <div class="hm-picker__acts">
              <button class="vq-btn vq-btn--sm vq-btn--ghost" data-hm-all="clu" type="button">All</button>
              <button class="vq-btn vq-btn--sm vq-btn--ghost" data-hm-none="clu" type="button">None</button>
            </div>
          </div>
          <div class="vq-card__body">
            <div class="qt-cluster-list hm-list" id="hm-clu-list">${hmCluRows}</div>
          </div>
        </div>

        ${hkAvailable ? `
        <div class="vq-card">
          <div class="vq-card__header">
            <div class="vq-card__title">Housekeeping genes</div>
            <div class="hm-picker__acts">
              <button class="vq-btn vq-btn--sm vq-btn--ghost" data-hm-all="hk" type="button">All</button>
              <button class="vq-btn vq-btn--sm vq-btn--ghost" data-hm-none="hk" type="button">None</button>
            </div>
          </div>
          <div class="vq-card__body">
            <div class="qt-cluster-list hm-list" id="hm-hk-list">${hmHkRows}</div>
          </div>
        </div>` : ''}
      </div>

      <div class="vq-chart-card">
        <div class="vq-chart-card__head">
          <div>
            <div class="vq-chart-card__title">Heatmap</div>
            <div class="vq-chart-card__sub" id="hm-sub"></div>
          </div>
          <div class="vq-section-actions">
            <button class="vq-btn vq-btn--sm vq-btn--ghost" id="hm-export-png" type="button">Export PNG</button>
            <button class="vq-btn vq-btn--sm vq-btn--ghost" id="hm-export-svg" type="button">Export SVG</button>
          </div>
        </div>

        <div class="hm-controls">
          <div class="vq-toggle" role="tablist" aria-label="Heatmap scale">
            <button class="vq-toggle__btn" type="button" data-hm-scale="lin">Linear</button>
            <button class="vq-toggle__btn active" type="button" data-hm-scale="log">Log</button>
            <button class="vq-toggle__btn" type="button" data-hm-scale="z">Z-score</button>
          </div>
          <label class="hm-ctl">Palette
            <select class="vq-input vq-input--sm" id="hm-palette">
              ${Object.entries(_PALETTES).map(([k, p]) =>
                `<option value="${k}">${esc(p.label)}</option>`).join('')}
            </select>
          </label>
          <label class="hm-ctl">Sort rows
            <select class="vq-input vq-input--sm" id="hm-sort">
              <option value="mean">Mean TPM</option>
              <option value="name">Name</option>
              <option value="type">Clusters first</option>
            </select>
          </label>
          <label class="hm-ctl hm-ctl--check">
            <input type="checkbox" class="vq-checkbox" id="hm-labels" checked> Cell values
          </label>
        </div>

        <div class="vq-chart-card__body"><div id="hm-plot" style="width:100%"></div></div>
      </div>
    </div>
  `;

  document.getElementById('qt-list').addEventListener('change', e => {
    if (e.target.matches('input[type="checkbox"]')) {
      if (e.target.checked) _selected.add(e.target.value);
      else _selected.delete(e.target.value);
      _renderPlot();
    }
  });
  el.querySelectorAll('[data-scale]').forEach(btn => btn.addEventListener('click', () => {
    el.querySelectorAll('[data-scale]').forEach(b => b.classList.toggle('active', b === btn));
    _logScale = btn.dataset.scale === 'log';
    _renderPlot();
  }));
  document.getElementById('qt-export')?.addEventListener('click', () => {
    const svg = document.querySelector('#qt-plot svg');
    if (svg) VQ.exportPNG(svg, 'quantification.png');
  });
  document.getElementById('qt-clear')?.addEventListener('click', () => {
    _selected.clear();
    document.querySelectorAll('#qt-list input[type="checkbox"]').forEach(cb => { cb.checked = false; });
    _renderPlot();
  });

  // ── Heatmap controls ────────────────────────────────────────────────────
  const bindList = (id, set) => document.getElementById(id)?.addEventListener('change', e => {
    if (!e.target.matches('input[type="checkbox"]')) return;
    if (e.target.checked) set.add(e.target.value); else set.delete(e.target.value);
    _renderHeatmap();
  });
  bindList('hm-clu-list', _hmClusters);
  bindList('hm-hk-list',  _hmGenes);

  el.querySelectorAll('[data-hm-all],[data-hm-none]').forEach(btn => {
    btn.addEventListener('click', () => {
      const on   = btn.hasAttribute('data-hm-all');
      const kind = btn.getAttribute(on ? 'data-hm-all' : 'data-hm-none');
      const listId = kind === 'clu' ? 'hm-clu-list' : 'hm-hk-list';
      const set    = kind === 'clu' ? _hmClusters : _hmGenes;
      set.clear();
      document.querySelectorAll(`#${listId} input[type="checkbox"]`).forEach(cb => {
        cb.checked = on;
        if (on) set.add(cb.value);
      });
      _renderHeatmap();
    });
  });

  el.querySelectorAll('[data-hm-scale]').forEach(btn => btn.addEventListener('click', () => {
    el.querySelectorAll('[data-hm-scale]').forEach(b => b.classList.toggle('active', b === btn));
    _hmScale = btn.dataset.hmScale;
    _renderHeatmap();
  }));
  document.getElementById('hm-palette')?.addEventListener('change', e => {
    _hmPalette = e.target.value; _renderHeatmap();
  });
  document.getElementById('hm-sort')?.addEventListener('change', e => {
    _hmSort = e.target.value; _renderHeatmap();
  });
  document.getElementById('hm-labels')?.addEventListener('change', e => {
    _hmLabels = e.target.checked; _renderHeatmap();
  });
  document.getElementById('hm-export-png')?.addEventListener('click', () => {
    const svg = document.querySelector('#hm-plot svg');
    if (svg) VQ.exportPNG(svg, 'cluster_heatmap.png');
  });
  document.getElementById('hm-export-svg')?.addEventListener('click', () => {
    const svg = document.querySelector('#hm-plot svg');
    if (svg) VQ.exportSVG(svg, 'cluster_heatmap.svg');
  });

  // Re-render when the card gains or changes width (e.g. when this tab first
  // opens, or the window resizes) so the plot always fills the card.
  let _lastW = 0;
  const plotHost = document.getElementById('qt-plot');
  new ResizeObserver(() => {
    const w = Math.round(plotHost.clientWidth);
    if (w && w !== _lastW) { _lastW = w; _renderPlot(); }
  }).observe(plotHost);

  let _lastHmW = 0;
  const hmHost = document.getElementById('hm-plot');
  new ResizeObserver(() => {
    const w = Math.round(hmHost.clientWidth);
    if (w && w !== _lastHmW) { _lastHmW = w; _renderHeatmap(); }
  }).observe(hmHost);

  _renderPlot();
  _renderHeatmap();
}

// ── Plot ────────────────────────────────────────────────────────────────────

function _renderPlot() {
  const host = document.getElementById('qt-plot');
  const sub  = document.getElementById('qt-sub');
  if (!host) return;
  host.innerHTML = '';

  const chosen = _clusters.filter(c => _selected.has(c.gid));
  if (!chosen.length) {
    host.innerHTML = `<div class="vq-empty">Select one or more clusters to plot TPM.</div>`;
    if (sub) sub.textContent = '';
    return;
  }

  const showBox = chosen.length >= 2;
  const cols = chosen.map(c => {
    const members = _tpmMembers(c).slice().sort((a, b) => a.tpm - b.tpm);
    return { gid: c.gid, species: c.species, members,
             omitted: c.members.length - members.length };
  });
  if (sub) {
    const totalOmit = cols.reduce((a, c) => a + c.omitted, 0);
    sub.textContent = `${chosen.length} cluster${chosen.length > 1 ? 's' : ''} · `
      + (showBox ? 'boxplot' : 'dot plot')
      + ` · scale: ${_logScale ? 'log₁₀(TPM + 1)' : 'linear TPM'}`
      + (totalOmit ? ` · ${totalOmit} member(s) omitted (no Salmon)` : '');
  }

  const tx = _logScale ? (v => Math.log10(v + 1)) : (v => v);
  const allVals = cols.flatMap(c => c.members.map(m => tx(m.tpm)));
  const maxV = d3.max(allVals) || 1;

  const padL = 54, padR = 18, padT = 14, padB = 70;
  const H = 320;

  // Fill 100% of the card: the viewBox width tracks the container's pixel
  // width (height fixed), so width:100% renders at exactly that size — no
  // distortion, no letter-boxing, and it spans the whole card.
  const W = Math.max(320, Math.round(host.clientWidth) || 900);
  const plotW = W - padL - padR;
  const colW  = showBox ? plotW / cols.length : plotW;

  const svg = d3.select(host).append('svg')
    .attr('viewBox', `0 0 ${W} ${H}`).attr('preserveAspectRatio', 'xMidYMin meet')
    .style('width', '100%').style('height', 'auto').style('display', 'block');

  const y = d3.scaleLinear().domain([0, maxV]).nice().range([H - padB, padT]);

  // y axis + grid
  y.ticks(5).forEach(t => {
    // Inline stroke (not the .ov-grid class) so the gridlines survive PNG/SVG
    // export, where external stylesheet rules are not carried over.
    svg.append('line').attr('x1', padL).attr('x2', W - padR)
      .attr('y1', y(t)).attr('y2', y(t))
      .attr('stroke', 'var(--vq-border)').attr('stroke-width', 1)
      .attr('shape-rendering', 'crispEdges');
    const lbl = _logScale ? (Math.pow(10, t) - 1) : t;
    svg.append('text').attr('x', padL - 8).attr('y', y(t))
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(_fmtTpm(lbl));
  });
  svg.append('text').attr('x', 14).attr('y', H / 2)
    .attr('text-anchor', 'middle').attr('class', 'ov-axis-label')
    .attr('transform', `rotate(-90 14 ${H / 2})`)
    .text(_logScale ? 'log₁₀(TPM + 1)' : 'TPM');

  if (showBox) {
    cols.forEach((c, i) => {
      const cx = padL + i * colW + colW / 2;
      const vals = c.members.map(m => tx(m.tpm));
      const boxW = Math.min(colW * 0.5, 110);
      if (vals.length >= 1) _drawBox(svg, cx, boxW, vals, y);
      c.members.forEach(m => {
        const jitter = vals.length > 1 ? (_hash(m.gid) - 0.5) * Math.min(colW * 0.4, 80) : 0;
        _point(svg, cx + jitter, y(tx(m.tpm)), m);
      });
      svg.append('text').attr('x', cx).attr('y', H - padB + 16)
        .attr('text-anchor', 'middle').attr('font-size', 10)
        .attr('font-family', 'var(--vq-font-mono)').attr('fill', 'var(--vq-text-2)')
        .text(c.gid);
      svg.append('text').attr('x', cx).attr('y', H - padB + 30)
        .attr('text-anchor', 'middle').attr('class', 'ov-axis-label')
        .text(`n=${c.members.length}${c.omitted ? ` (+${c.omitted})` : ''}`);
    });
  } else {
    // Single cluster: one x position per member so points never overlap.
    const c = cols[0];
    const xb = d3.scalePoint().domain(c.members.map((_, j) => j))
      .range([padL, W - padR]).padding(0.7);
    c.members.forEach((m, j) => {
      const px = xb(j);
      _point(svg, px, y(tx(m.tpm)), m);
      svg.append('text').attr('x', px).attr('y', H - padB + 14)
        .attr('text-anchor', 'end').attr('class', 'ov-axis-label')
        .attr('transform', `rotate(-40 ${px} ${H - padB + 14})`)
        .text(_trunc(m.sample, 12));
    });
    svg.append('text').attr('x', (padL + W - padR) / 2).attr('y', H - 8)
      .attr('text-anchor', 'middle').attr('font-size', 10)
      .attr('font-family', 'var(--vq-font-mono)').attr('fill', 'var(--vq-text-2)')
      .text(`${c.gid} · n=${c.members.length}${c.omitted ? ` (+${c.omitted})` : ''}`);
  }
}

function _drawBox(svg, cx, boxW, vals, y) {
  const s = vals.slice().sort((a, b) => a - b);
  const q1 = d3.quantile(s, 0.25), med = d3.quantile(s, 0.5), q3 = d3.quantile(s, 0.75);
  const lo = s[0], hi = s[s.length - 1];
  const x0 = cx - boxW / 2, x1 = cx + boxW / 2;

  svg.append('line').attr('x1', cx).attr('x2', cx).attr('y1', y(lo)).attr('y2', y(hi))
    .attr('stroke', 'var(--vq-text-3)').attr('stroke-width', 1);
  svg.append('rect').attr('x', x0).attr('y', y(q3)).attr('width', boxW)
    .attr('height', Math.max(1, y(q1) - y(q3)))
    .attr('fill', 'var(--vq-accent-subtle, #e8eef7)')
    .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.4).attr('rx', 2);
  svg.append('line').attr('x1', x0).attr('x2', x1).attr('y1', y(med)).attr('y2', y(med))
    .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 2);
}

// ── Heatmap ─────────────────────────────────────────────────────────────────
/* Complements the boxplot: the boxplot answers "how does one cluster spread
   across its members", the heatmap answers "which virus is where, and how does
   that compare to the host's own housekeeping baseline".

   Rows are chosen by the user — viral clusters and, when the reference pathway
   ran with a curated housekeeping panel, individual HK genes on the same axis.
   Putting them in one matrix is the point: a viral row brighter than the HK
   rows means the virus out-transcribes the host's own reference genes. */

const _PALETTES = {
  viridis: { label: 'Viridis', fn: v => d3.interpolateViridis(v),        diverging: false },
  magma:   { label: 'Magma',   fn: v => d3.interpolateMagma(v),          diverging: false },
  blues:   { label: 'Blues',   fn: v => d3.interpolateBlues(0.12 + v * 0.88), diverging: false },
};

let _hmClusters = new Set();
let _hmGenes    = new Set();
let _hkGenes    = [];        // [{ name, bySample: {sample: tpm} }]
let _hmScale    = 'log';     // 'lin' | 'log' | 'z'
let _hmPalette  = 'viridis';
let _hmLabels   = true;
let _hmSort     = 'mean';    // 'mean' | 'name' | 'type'

/* Union of the curated HK panel across samples, with each sample's TPM. */
function _collectHkGenes(samples) {
  const byName = new Map();
  samples.forEach(s => {
    ((s.salmon && s.salmon.ref_hk && s.salmon.ref_hk.genes) || []).forEach(g => {
      if (!byName.has(g.name)) byName.set(g.name, { name: g.name, bySample: {} });
      byName.get(g.name).bySample[s.sample] = g.tpm;
    });
  });
  return [...byName.values()];
}

/* One row per selected cluster / gene, keyed by sample.

   A cluster's members in one sample are SUMMED, not averaged: Salmon splits
   reads between near-identical contigs of the same virus, so the fragments'
   TPMs are partial counts of one thing. Averaging would under-report a virus
   purely because its assembly fragmented. */
function _hmRows() {
  const rows = [];

  _clusters.filter(c => _hmClusters.has(c.gid)).forEach(c => {
    const bySample = {};
    _tpmMembers(c).forEach(m => {
      bySample[m.sample] = (bySample[m.sample] || 0) + m.tpm;
    });
    rows.push({
      key: c.gid, kind: 'cluster',
      label: c.species || c.gid, sub: c.gid,
      n: _tpmMembers(c).length, bySample,
    });
  });

  _hkGenes.filter(g => _hmGenes.has(g.name)).forEach(g => {
    rows.push({
      key: 'hk:' + g.name, kind: 'hk',
      label: g.name, sub: 'housekeeping',
      n: null, bySample: g.bySample,
    });
  });

  const meanOf = r => {
    const v = _hmSamples.map(s => r.bySample[s.sample]).filter(x => x != null);
    return v.length ? d3.mean(v) : -1;
  };
  const cmp = {
    mean: (a, b) => meanOf(b) - meanOf(a),
    name: (a, b) => a.label.localeCompare(b.label),
    type: (a, b) => (a.kind === b.kind ? meanOf(b) - meanOf(a) : (a.kind === 'cluster' ? -1 : 1)),
  }[_hmSort];
  return rows.sort(cmp);
}

function _renderHeatmap() {
  const host = document.getElementById('hm-plot');
  const sub  = document.getElementById('hm-sub');
  if (!host) return;
  host.innerHTML = '';

  const rows = _hmRows();
  if (!rows.length) {
    host.innerHTML = `<div class="vq-empty">Select clusters or housekeeping genes to build the heatmap.</div>`;
    if (sub) sub.textContent = '';
    return;
  }
  const cols = _hmSamples;

  // Per-row transform. z-score is per row, so rows with very different absolute
  // levels can still be compared by shape; it needs ≥2 samples to mean anything.
  const canZ = cols.length >= 2;
  const mode = (_hmScale === 'z' && !canZ) ? 'log' : _hmScale;
  const raw  = (r, s) => r.bySample[s.sample];

  let valueOf, domain, interp, legendTicks;
  if (mode === 'z') {
    const stats = new Map(rows.map(r => {
      const v = cols.map(s => raw(r, s)).filter(x => x != null);
      const mu = d3.mean(v) ?? 0;
      const sd = v.length > 1 ? (d3.deviation(v) || 0) : 0;
      return [r.key, { mu, sd }];
    }));
    valueOf = (r, s) => {
      const v = raw(r, s); if (v == null) return null;
      const { mu, sd } = stats.get(r.key);
      return sd ? (v - mu) / sd : 0;
    };
    const vals   = rows.flatMap(r => cols.map(s => valueOf(r, s))).filter(v => v != null);
    const maxAbs = Math.max(1e-6, d3.max(vals.map(Math.abs)) || 1);
    domain = [-maxAbs, maxAbs];
    interp = v => d3.interpolateRdBu(1 - v);          // red = high, blue = low
    legendTicks = [-maxAbs, 0, maxAbs].map(v => [v, v.toFixed(1)]);
  } else {
    valueOf = (r, s) => {
      const v = raw(r, s); if (v == null) return null;
      return mode === 'log' ? Math.log10(v + 1) : v;
    };
    const vals = rows.flatMap(r => cols.map(s => valueOf(r, s))).filter(v => v != null);
    const maxV = d3.max(vals) || 1;
    domain = [0, maxV];
    interp = _PALETTES[_hmPalette].fn;
    legendTicks = [0, maxV / 2, maxV].map(v =>
      [v, _fmtTpm(mode === 'log' ? Math.pow(10, v) - 1 : v)]);
  }

  const colour = d3.scaleSequential(domain, interp);

  if (sub) {
    const nClu = rows.filter(r => r.kind === 'cluster').length;
    const nHk  = rows.length - nClu;
    sub.textContent =
      `${rows.length} row${rows.length > 1 ? 's' : ''} (${nClu} cluster${nClu === 1 ? '' : 's'}`
      + (nHk ? `, ${nHk} housekeeping` : '') + `) × ${cols.length} sample${cols.length > 1 ? 's' : ''}`
      + ` · ${mode === 'z' ? 'z-score per row' : mode === 'log' ? 'log₁₀(TPM + 1)' : 'linear TPM'}`
      + (_hmScale === 'z' && !canZ ? ' · z-score needs ≥2 samples, showing log' : '');
  }

  // Geometry. Labels sit in a left gutter; sample names are rotated above.
  const W     = Math.max(360, Math.round(host.clientWidth) || 900);
  const padL  = Math.min(230, Math.max(120, Math.round(W * 0.24)));
  const padR  = 16;
  const headH = 74;
  const rowH  = rows.length > 26 ? 18 : rows.length > 14 ? 22 : 26;
  const legH  = 46;
  const gridW = W - padL - padR;
  const cellW = gridW / cols.length;
  const H     = headH + rows.length * rowH + legH;

  const svg = d3.select(host).append('svg')
    .attr('viewBox', `0 0 ${W} ${H}`).attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto').style('display', 'block');

  // Exports carry no external stylesheet, so every fill/font is inline here.
  const FONT = 'var(--vq-font)';
  const txt = (x, y, s, opts = {}) => svg.append('text')
    .attr('x', x).attr('y', y)
    .attr('font-family', opts.mono ? 'var(--vq-font-mono)' : FONT)
    .attr('font-size', opts.size || 10)
    .attr('font-weight', opts.weight || 400)
    .attr('fill', opts.fill || 'var(--vq-text-2)')
    .attr('text-anchor', opts.anchor || 'start')
    .attr('dominant-baseline', opts.baseline || 'central')
    .text(s);

  // Column headers (rotated so long sample names stay readable)
  cols.forEach((s, j) => {
    const cx = padL + j * cellW + cellW / 2;
    svg.append('text')
      .attr('x', cx).attr('y', headH - 8)
      .attr('font-family', FONT).attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)').attr('text-anchor', 'start')
      .attr('transform', `rotate(-42 ${cx} ${headH - 8})`)
      .text(_trunc(s.sample, 18));
  });

  rows.forEach((r, i) => {
    const y = headH + i * rowH;

    // Row label: virus name (or gene), with the cluster id / kind underneath
    // when there is vertical room for it.
    txt(padL - 8, y + rowH / 2 - (rowH >= 22 ? 5 : 0), _trunc(r.label, Math.floor(padL / 6.4)),
        { anchor: 'end', weight: 600, size: rowH >= 22 ? 11 : 10,
          fill: r.kind === 'hk' ? 'var(--vq-text-2)' : 'var(--vq-text)' })
      .append('title').text(r.label);
    if (rowH >= 22) {
      txt(padL - 8, y + rowH / 2 + 7,
          r.kind === 'cluster' ? `${r.sub} · n=${r.n}` : r.sub,
          { anchor: 'end', size: 9, fill: 'var(--vq-text-3)', mono: r.kind === 'cluster' });
    }

    cols.forEach((s, j) => {
      const x  = padL + j * cellW;
      const v  = valueOf(r, s);
      const rv = raw(r, s);

      if (v == null) {
        // Absent ≠ zero: no member of this cluster in this sample (or Salmon
        // never ran there). Left unpainted so it cannot read as "measured 0".
        svg.append('rect').attr('x', x).attr('y', y)
          .attr('width', cellW - 1).attr('height', rowH - 1)
          .attr('fill', 'var(--vq-surface-2)')
          .attr('stroke', 'var(--vq-border)').attr('stroke-width', 1)
          .attr('stroke-dasharray', '2,2')
          .on('mousemove', e => VQ.tooltipShow(
            `<b>${VQ.esc(r.label)}</b><br>${VQ.esc(s.sample)}<br>not quantified`, e))
          .on('mouseleave', () => VQ.tooltipHide());
        return;
      }

      const fill = colour(v);
      svg.append('rect').attr('x', x).attr('y', y)
        .attr('width', cellW - 1).attr('height', rowH - 1)
        .attr('fill', fill)
        .on('mousemove', e => VQ.tooltipShow(
          `<b>${VQ.esc(r.label)}</b><br>${VQ.esc(s.sample)}<br>`
          + `TPM ${rv.toFixed(2)}${mode === 'z' ? ` · z ${v.toFixed(2)}` : ''}`, e))
        .on('mouseleave', () => VQ.tooltipHide());

      if (_hmLabels && cellW >= 46 && rowH >= 18) {
        // Flip the label to white on dark cells so it stays legible in both
        // the light and dark ends of every palette.
        const lab = d3.lab(fill);
        txt(x + cellW / 2, y + rowH / 2,
            mode === 'z' ? v.toFixed(1) : _fmtTpm(rv),
            { anchor: 'middle', size: 9.5, fill: lab.l < 62 ? '#ffffff' : '#111827' });
      }
    });
  });

  // Colour legend
  const legY = headH + rows.length * rowH + 16;
  const legW = Math.min(220, gridW);
  const gradId = 'hm-grad-' + _hmPalette + '-' + mode;
  const grad = svg.append('defs').append('linearGradient')
    .attr('id', gradId).attr('x1', '0%').attr('x2', '100%');
  d3.range(0, 1.0001, 0.05).forEach(t => {
    grad.append('stop').attr('offset', `${t * 100}%`)
      .attr('stop-color', colour(domain[0] + t * (domain[1] - domain[0])));
  });
  svg.append('rect').attr('x', padL).attr('y', legY)
    .attr('width', legW).attr('height', 10).attr('rx', 2)
    .attr('fill', `url(#${gradId})`)
    .attr('stroke', 'var(--vq-border)').attr('stroke-width', 1);
  legendTicks.forEach(([v, label]) => {
    const t = (v - domain[0]) / (domain[1] - domain[0] || 1);
    txt(padL + t * legW, legY + 22, label,
        { anchor: t === 0 ? 'start' : t >= 1 ? 'end' : 'middle', size: 9, fill: 'var(--vq-text-3)' });
  });
  txt(padL + legW + 12, legY + 5,
      mode === 'z' ? 'z-score (row)' : mode === 'log' ? 'TPM (log scale)' : 'TPM',
      { size: 10, fill: 'var(--vq-text-3)' });
}

// ── helpers ─────────────────────────────────────────────────────────────────

function _point(svg, cx, cy, m) {
  svg.append('circle').attr('cx', cx).attr('cy', cy).attr('r', 4)
    .attr('fill', 'var(--vq-accent)').attr('opacity', 0.8)
    .attr('stroke', '#fff').attr('stroke-width', 1)
    .on('mousemove', e => VQ.tooltipShow(
      `<b>${VQ.esc(m.sample)}</b> · ${VQ.esc(m.seq_id)}<br>TPM ${m.tpm.toFixed(2)}`, e))
    .on('mouseleave', () => VQ.tooltipHide());
}

function _trunc(s, n) {
  s = String(s ?? '');
  return s.length > n ? s.slice(0, n - 1) + '…' : s;
}

function _exportBtn(id) {
  return `<button class="vq-btn vq-btn--sm vq-btn--ghost" id="${id}" type="button">Export PNG</button>`;
}

// Deterministic 0..1 jitter from a string so points don't jump on re-render.
function _hash(str) {
  let h = 0;
  for (let i = 0; i < str.length; i++) h = (h * 31 + str.charCodeAt(i)) | 0;
  return (Math.abs(h) % 1000) / 1000;
}

function _fmtTpm(v) {
  if (v >= 1000) return (v / 1000).toFixed(1) + 'k';
  if (v >= 10)   return v.toFixed(0);
  return v.toFixed(1);
}

window.vqInitQuant = vqInitQuant;
})();
