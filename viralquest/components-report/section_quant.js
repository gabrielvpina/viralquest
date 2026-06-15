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

let _clusters = [];
let _selected = new Set();
let _logScale = false;

function _tpmMembers(c) {
  return c.members.filter(m => m.tpm != null && !isNaN(m.tpm));
}
function _quantifiable(c) {
  return _tpmMembers(c).length > 0;
}

function vqInitQuant(clusters) {
  const el = document.getElementById('section-quant');
  if (!el) return;

  _clusters = (clusters || []).filter(_quantifiable);

  // Hide the whole section when no cluster carries any TPM (no Salmon anywhere).
  if (!_clusters.length) {
    el.hidden = true;
    document.getElementById('tab-quant')?.style.setProperty('display', 'none');
    return;
  }

  const esc = VQ.esc;
  // Pre-select up to the first 3 clusters so the plot isn't empty on open.
  _selected = new Set(_clusters.slice(0, Math.min(3, _clusters.length)).map(c => c.gid));

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

  // Re-render when the card gains or changes width (e.g. when this tab first
  // opens, or the window resizes) so the plot always fills the card.
  let _lastW = 0;
  const plotHost = document.getElementById('qt-plot');
  new ResizeObserver(() => {
    const w = Math.round(plotHost.clientWidth);
    if (w && w !== _lastW) { _lastW = w; _renderPlot(); }
  }).observe(plotHost);

  _renderPlot();
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
