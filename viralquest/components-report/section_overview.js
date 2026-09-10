/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_overview.js — Section 1: multi-sample Overview
   Cross-sample comparison of every loaded ViralQuest result.
   Every chart lives in a dynamic card that adapts from a few
   samples to many, and degrades gracefully when data is absent.
   ============================================================ */

const STEP_LABELS = [
  ['nr',     'NR'],
  ['blastn', 'BLASTn'],
  ['salmon', 'Salmon'],
  ['llm',    'LLM'],
  ['cap3',   'CAP3'],
];

const PALETTE = [
  '#2563eb', '#16a34a', '#db2777', '#d97706', '#7c3aed',
  '#0891b2', '#dc2626', '#65a30d', '#c026d3', '#0d9488',
  '#ea580c', '#4f46e5', '#059669', '#e11d48', '#9333ea',
];
const colorFor = i => PALETTE[i % PALETTE.length];

// Single accent colour for all quantitative-per-sample charts; only the
// family-distribution chart distinguishes samples/categories by colour.
const QUANT = 'var(--vq-accent)';

function vqInitOverview(report) {
  const el = document.getElementById('section-overview');
  if (!el) return;

  const esc      = VQ.esc;
  const samples  = report.samples || [];
  const meta     = report.meta || {};

  if (!samples.length) {
    el.innerHTML = `<div class="vq-empty">No samples loaded.</div>`;
    return;
  }

  // ── Aggregates ──────────────────────────────────────────────────────────
  const totalConfirmed = samples.reduce((a, s) => a + (s.n_confirmed || 0), 0);
  const familyUnion = new Set();
  samples.forEach(s => Object.keys(s.families || {}).forEach(f => familyUnion.add(f)));
  const totalClusters = (report.clusters || []).length;
  const anyLLM  = samples.some(s => (s.llm_scores || []).length);
  const anyHeur = samples.some(s => (s.heuristic_scores || []).length);
  const salmonSamples = samples.filter(s => s.salmon);
  const anySalmon = salmonSamples.length > 0;
  const anyConserved = salmonSamples.some(s => Object.keys(s.salmon.conserved || {}).length);

  const fmt = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();

  // ── Scaffold ────────────────────────────────────────────────────────────
  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Overview</div>
        <div class="vq-section-sub">
          ${samples.length} sample${samples.length > 1 ? 's' : ''} compared
          &nbsp;·&nbsp; ${fmt(totalConfirmed)} confirmed viral sequences
          &nbsp;·&nbsp; ${familyUnion.size} viral famil${familyUnion.size === 1 ? 'y' : 'ies'}
        </div>
      </div>
    </div>

    <div class="vq-stats-page">
      <div class="vq-stats-row vq-stats-row--top">
        ${_chip('Samples',            fmt(samples.length), 'accent', 'ViralQuest runs')}
        ${_chip('Confirmed Viral',    fmt(totalConfirmed), 'success', 'across all samples')}
        ${_chip('Viral Families',     fmt(familyUnion.size), '', 'distinct, union')}
        ${_chip('Cross-sample Clusters', fmt(totalClusters), 'accent',
                report.has_clusters ? 'shared between samples' : 'none detected')}
      </div>

      <!-- Pipeline uniformity matrix — always full width -->
      <div class="vq-chart-card" id="ov-steps-card" style="min-height:auto">
        <div class="vq-chart-card__head">
          <div>
            <div class="vq-chart-card__title">Pipeline Steps per Sample</div>
            <div class="vq-chart-card__sub">which optional stages ran — analysis uniformity</div>
          </div>
        </div>
        <div class="vq-chart-card__body" id="ov-steps-body"></div>
      </div>

      ${anyConserved ? `
      <!-- Host attribution — full width, one column per sample -->
      <div class="vq-chart-card" id="ov-kingdoms-card" style="min-height:auto">
        <div class="vq-chart-card__head">
          <div>
            <div class="vq-chart-card__title">Conserved Housekeeping by Kingdom</div>
            <div class="vq-chart-card__sub">
              total TPM and detected genes per kingdom — host provenance and contamination
            </div>
          </div>
        </div>
        <div class="vq-chart-card__body" id="ov-kingdoms-body"></div>
      </div>` : ''}

      <div class="vq-masonry">
        ${anySalmon ? _card('ov-maprate', 'Salmon Mapping Rate per Sample', '% of reads mapped') : ''}
        ${_card('ov-seqs',    'Confirmed Sequences per Sample', 'final viral contigs')}
        ${_card('ov-fams',    'Viral Family Diversity per Sample', 'distinct families')}
        ${_card('ov-famdist', 'Family Distribution per Sample', 'composition, all samples')}
        ${_card('ov-lengths', 'Sequence Length per Sample', 'min · median · max (nt)')}
        ${anyHeur ? _card('ov-heur', 'Heuristic Score per Sample', 'VQ score distribution') : ''}
        ${anyLLM  ? _card('ov-llm',  'LLM Score per Sample', 'VQ score distribution') : ''}
        ${report.has_clusters ? _card('ov-clusters', 'Cross-sample Clusters per Sample', 'shared sequence groups') : ''}
      </div>
    </div>
  `;

  // ── Render charts ───────────────────────────────────────────────────────
  _renderStepsMatrix(document.getElementById('ov-steps-body'), samples);
  if (anyConserved) _renderKingdomMatrix(_body('ov-kingdoms'), salmonSamples);
  if (anySalmon)
    _barChart(_body('ov-maprate'),
              salmonSamples.map(s => ({ label: s.sample, value: s.salmon.mapping_rate ?? 0 })),
              { unit: '%' });
  _barChart(_body('ov-seqs'),    samples.map(s => ({ label: s.sample, value: s.n_confirmed || 0 })));
  _lineDotChart(_body('ov-fams'), samples.map(s => ({ label: s.sample, value: Object.keys(s.families || {}).length })));
  _stackedFamilies(_body('ov-famdist'), samples);
  _lengthRanges(_body('ov-lengths'), samples);
  if (anyHeur) _scoreBoxes(_body('ov-heur'), samples, 'heuristic_scores');
  if (anyLLM)  _scoreBoxes(_body('ov-llm'),  samples, 'llm_scores');
  if (report.has_clusters)
    _barChart(_body('ov-clusters'), samples.map(s => ({ label: s.sample, value: s.n_clusters || 0 })));
}

// ── Card scaffolding ────────────────────────────────────────────────────────

function _chip(label, value, mod, sub) {
  const esc = VQ.esc;
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${esc(label)}</div>
      <div class="vq-stat-chip__value" title="${esc(value)}">${esc(value)}</div>
      ${sub ? `<div class="vq-stat-chip__sub" title="${esc(sub)}">${esc(sub)}</div>` : ''}
    </div>`;
}

function _card(id, title, sub, cls) {
  const esc = VQ.esc;
  return `
    <div class="vq-chart-card${cls ? ' ' + cls : ''}" id="${id}-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">${esc(title)}</div>
          <div class="vq-chart-card__sub">${esc(sub)}</div>
        </div>
      </div>
      <div class="vq-chart-card__body" id="${id}-body" style="width:100%"></div>
    </div>`;
}

const _body = id => document.getElementById(id + '-body');

// ── Pipeline steps matrix ───────────────────────────────────────────────────

function _renderStepsMatrix(host, samples) {
  if (!host) return;
  const esc = VQ.esc;
  const yes = `<span class="ov-step ov-step--yes" title="ran">✓</span>`;
  const no  = `<span class="ov-step ov-step--no"  title="not run">·</span>`;

  const rows = samples.map(s => `
    <tr>
      <td class="ov-matrix__sample" title="${esc(s.sample)}">${esc(s.sample)}</td>
      ${STEP_LABELS.map(([k]) => `<td>${(s.steps || {})[k] ? yes : no}</td>`).join('')}
    </tr>`).join('');

  host.innerHTML = `
    <div class="ov-matrix-wrap">
      <table class="ov-matrix">
        <thead>
          <tr>
            <th class="ov-matrix__sample">Sample</th>
            ${STEP_LABELS.map(([, lbl]) => `<th>${esc(lbl)}</th>`).join('')}
          </tr>
        </thead>
        <tbody>${rows}</tbody>
      </table>
    </div>`;
}

// ── Generic SVG sizing ──────────────────────────────────────────────────────

function _svg(host, width, height) {
  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('viewBox', `0 0 ${width} ${height}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%')
    .style('height', 'auto')
    .style('display', 'block');
  return svg;
}

// Row height scales with sample count so few/many both look right.
const _rowH = n => Math.max(16, Math.min(34, Math.round(360 / Math.max(n, 1))));

// ── Horizontal bar chart ────────────────────────────────────────────────────

function _barChart(host, data, opts) {
  if (!host) return;
  if (!data.length) { host.innerHTML = `<div class="vq-empty">No data.</div>`; return; }
  const unit = (opts && opts.unit) || '';

  const W = 440, padL = 110, padR = 48, padT = 8;
  const rh = _rowH(data.length), gap = 6;
  const H  = padT + data.length * (rh + gap);
  const max = d3.max(data, d => d.value) || 1;
  const x = d3.scaleLinear().domain([0, max]).range([padL, W - padR]);

  const svg = _svg(host, W, H);
  data.forEach((d, i) => {
    const y = padT + i * (rh + gap);
    const g = svg.append('g');
    g.append('text').attr('x', padL - 8).attr('y', y + rh / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(_trunc(d.label, 16));
    g.append('rect').attr('x', padL).attr('y', y)
      .attr('width', Math.max(1, x(d.value) - padL)).attr('height', rh)
      .attr('rx', 3).attr('fill', QUANT);
    g.append('text').attr('x', x(d.value) + 6).attr('y', y + rh / 2)
      .attr('dominant-baseline', 'central').attr('class', 'ov-value-label')
      .text(d.value.toLocaleString() + unit);
  });
}

// ── Connected dots / line chart (per sample) ────────────────────────────────

function _lineDotChart(host, data) {
  if (!host) return;
  if (!data.length) { host.innerHTML = `<div class="vq-empty">No data.</div>`; return; }

  const W = 440, padL = 40, padR = 24, padT = 16, padB = 64;
  const H = 240;
  const max = d3.max(data, d => d.value) || 1;
  const x = d3.scalePoint().domain(data.map((_, i) => i)).range([padL, W - padR]).padding(0.5);
  const y = d3.scaleLinear().domain([0, max]).nice().range([H - padB, padT]);

  const svg = _svg(host, W, H);

  // y gridlines
  y.ticks(4).forEach(t => {
    svg.append('line').attr('x1', padL).attr('x2', W - padR)
      .attr('y1', y(t)).attr('y2', y(t)).attr('class', 'ov-grid');
    svg.append('text').attr('x', padL - 6).attr('y', y(t))
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(t);
  });

  const line = d3.line().x((_, i) => x(i)).y(d => y(d.value));
  svg.append('path').datum(data).attr('fill', 'none')
    .attr('stroke', QUANT).attr('stroke-width', 2).attr('d', line);

  data.forEach((d, i) => {
    svg.append('circle').attr('cx', x(i)).attr('cy', y(d.value)).attr('r', 4)
      .attr('fill', QUANT).attr('stroke', '#fff').attr('stroke-width', 1.5)
      .on('mousemove', e => VQ.tooltipShow(`<b>${VQ.esc(d.label)}</b><br>${d.value}`, e))
      .on('mouseleave', () => VQ.tooltipHide());
    svg.append('text').attr('x', x(i)).attr('y', H - padB + 12)
      .attr('text-anchor', 'end').attr('class', 'ov-axis-label')
      .attr('transform', `rotate(-40 ${x(i)} ${H - padB + 12})`)
      .text(_trunc(d.label, 12));
  });
}

// ── Stacked family distribution (one row per sample) ────────────────────────

function _stackedFamilies(host, samples) {
  if (!host) return;
  const esc = VQ.esc;

  // Union of families ranked by total abundance → stable colour mapping.
  const totals = {};
  samples.forEach(s => Object.entries(s.families || {}).forEach(([f, n]) => {
    totals[f] = (totals[f] || 0) + n;
  }));
  const families = Object.keys(totals).sort((a, b) => totals[b] - totals[a]);
  if (!families.length) { host.innerHTML = `<div class="vq-empty">No family data.</div>`; return; }
  const cIdx = Object.fromEntries(families.map((f, i) => [f, i]));

  const W = 440, padL = 110, padR = 12, padT = 6;
  const rh = _rowH(samples.length), gap = 6;
  const H  = padT + samples.length * (rh + gap);
  const svg = _svg(host, W, H);

  samples.forEach((s, i) => {
    const y = padT + i * (rh + gap);
    const total = Object.values(s.families || {}).reduce((a, n) => a + n, 0) || 1;
    const x = d3.scaleLinear().domain([0, total]).range([padL, W - padR]);

    svg.append('text').attr('x', padL - 8).attr('y', y + rh / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(_trunc(s.sample, 16));

    // Rounded outer corners (matching the other overview bars) via a per-row
    // clip; inner segment boundaries stay crisp.
    const clipId = `ov-fam-clip-${i}`;
    svg.append('clipPath').attr('id', clipId)
      .append('rect').attr('x', padL).attr('y', y)
      .attr('width', (W - padR) - padL).attr('height', rh).attr('rx', 3).attr('ry', 3);
    const rowG = svg.append('g').attr('clip-path', `url(#${clipId})`);

    let acc = 0;
    // Draw segments in the global family order for visual consistency.
    families.forEach(f => {
      const n = (s.families || {})[f];
      if (!n) return;
      const x0 = x(acc), x1 = x(acc + n);
      rowG.append('rect').attr('x', x0).attr('y', y)
        .attr('width', Math.max(0.5, x1 - x0)).attr('height', rh)
        .attr('fill', colorFor(cIdx[f]))
        .on('mousemove', e => VQ.tooltipShow(`<b>${esc(f)}</b><br>${esc(s.sample)}: ${n}`, e))
        .on('mouseleave', () => VQ.tooltipHide());
      acc += n;
    });
  });

  // Legend (top families; rest folded into "other" note)
  const legendMax = 12;
  const shown = families.slice(0, legendMax);
  const legend = shown.map(f =>
    `<span class="ov-legend__item"><span class="ov-legend__dot" style="background:${colorFor(cIdx[f])}"></span>${esc(f)}</span>`
  ).join('');
  const extra = families.length > legendMax ? `<span class="ov-legend__item">+${families.length - legendMax} more</span>` : '';
  const wrap = document.createElement('div');
  wrap.className = 'ov-legend';
  wrap.innerHTML = legend + extra;
  host.appendChild(wrap);
}

// ── Length ranges per sample (min · median · max) ───────────────────────────

function _lengthRanges(host, samples) {
  if (!host) return;
  const rows = samples.map(s => {
    const L = (s.lengths || []).slice().sort((a, b) => a - b);
    if (!L.length) return { label: s.sample, min: null };
    return {
      label: s.sample,
      min: L[0], max: L[L.length - 1],
      med: L[Math.floor(L.length / 2)],
    };
  });
  const valid = rows.filter(r => r.min != null);
  if (!valid.length) { host.innerHTML = `<div class="vq-empty">No length data.</div>`; return; }

  const W = 440, padL = 110, padR = 56, padT = 8;
  const rh = _rowH(rows.length), gap = 6;
  const H  = padT + rows.length * (rh + gap);
  const max = d3.max(valid, r => r.max) || 1;
  const x = d3.scaleLinear().domain([0, max]).nice().range([padL, W - padR]);
  const svg = _svg(host, W, H);

  rows.forEach((r, i) => {
    const y = padT + i * (rh + gap), cy = y + rh / 2;
    svg.append('text').attr('x', padL - 8).attr('y', cy)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(_trunc(r.label, 16));
    if (r.min == null) return;
    svg.append('line').attr('x1', x(r.min)).attr('x2', x(r.max))
      .attr('y1', cy).attr('y2', cy).attr('stroke', QUANT).attr('stroke-width', 2);
    [['min', r.min], ['max', r.max]].forEach(([, v]) =>
      svg.append('circle').attr('cx', x(v)).attr('cy', cy).attr('r', 2.5).attr('fill', QUANT));
    svg.append('circle').attr('cx', x(r.med)).attr('cy', cy).attr('r', 4)
      .attr('fill', '#fff').attr('stroke', QUANT).attr('stroke-width', 2)
      .on('mousemove', e => VQ.tooltipShow(
        `<b>${VQ.esc(r.label)}</b><br>min ${r.min} · median ${r.med} · max ${r.max} nt`, e))
      .on('mouseleave', () => VQ.tooltipHide());
    svg.append('text').attr('x', x(r.max) + 6).attr('y', cy)
      .attr('dominant-baseline', 'central').attr('class', 'ov-value-label')
      .text(r.max.toLocaleString());
  });
}

// ── Score "boxes" per sample (min–max range + mean dot) ─────────────────────

function _scoreBoxes(host, samples, key) {
  if (!host) return;
  const rows = samples.map(s => {
    const v = (s[key] || []).slice().sort((a, b) => a - b);
    if (!v.length) return { label: s.sample, n: 0 };
    return {
      label: s.sample, n: v.length,
      min: v[0], max: v[v.length - 1],
      mean: v.reduce((a, b) => a + b, 0) / v.length,
    };
  });
  const valid = rows.filter(r => r.n);
  if (!valid.length) { host.innerHTML = `<div class="vq-empty">No score data.</div>`; return; }

  const W = 440, padL = 110, padR = 48, padT = 8;
  const rh = _rowH(rows.length), gap = 6;
  const H  = padT + rows.length * (rh + gap);
  const x = d3.scaleLinear().domain([0, 100]).range([padL, W - padR]);
  const svg = _svg(host, W, H);

  [0, 25, 50, 75, 100].forEach(t =>
    svg.append('line').attr('x1', x(t)).attr('x2', x(t))
      .attr('y1', padT).attr('y2', H - 2).attr('class', 'ov-grid'));

  rows.forEach((r, i) => {
    const y = padT + i * (rh + gap), cy = y + rh / 2;
    svg.append('text').attr('x', padL - 8).attr('y', cy)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .attr('class', 'ov-axis-label').text(_trunc(r.label, 16));
    if (!r.n) return;
    svg.append('line').attr('x1', x(r.min)).attr('x2', x(r.max))
      .attr('y1', cy).attr('y2', cy).attr('stroke', QUANT).attr('stroke-width', 2);
    svg.append('circle').attr('cx', x(r.mean)).attr('cy', cy).attr('r', 4)
      .attr('fill', QUANT).attr('stroke', '#fff').attr('stroke-width', 1.5)
      .on('mousemove', e => VQ.tooltipShow(
        `<b>${VQ.esc(r.label)}</b><br>n=${r.n} · mean ${r.mean.toFixed(1)}<br>range ${r.min.toFixed(0)}–${r.max.toFixed(0)}`, e))
      .on('mouseleave', () => VQ.tooltipHide());
  });
}

// ── Conserved housekeeping: kingdom × sample ────────────────────────────────

/* Kingdom totals span three orders of magnitude between the host kingdom and
   the incidental ones, so cell shading is symlog — a linear ramp paints every
   non-host kingdom the same blank white and throws away the contamination
   signal.  A kingdom with nothing detected is left deliberately unshaded: an
   undetected kingdom is the finding, and must not read as a small measurement. */
function _renderKingdomMatrix(host, samples) {
  if (!host) return;
  const esc = VQ.esc;

  const kingdoms = [...new Set(
    samples.flatMap(s => Object.keys(s.salmon.conserved || {}))
  )];
  if (!kingdoms.length) { host.innerHTML = `<div class="vq-empty">No conserved housekeeping data.</div>`; return; }

  const cell = (s, k) => (s.salmon.conserved || {})[k] || null;
  const totalFor = k => d3.sum(samples, s => (cell(s, k) || {}).tpm_sum || 0);
  kingdoms.sort((a, b) => totalFor(b) - totalFor(a) || a.localeCompare(b));

  const maxTpm = d3.max(samples.flatMap(s => kingdoms.map(k => (cell(s, k) || {}).tpm_sum || 0))) || 1;
  const shade  = d3.scaleSymlog().domain([0, maxTpm]).range([0, 0.85]);

  const rows = kingdoms.map(k => {
    const cells = samples.map(s => {
      const c = cell(s, k);
      if (!c) return `<td class="ov-kg-cell" title="${esc(k)} — not quantified in ${esc(s.sample)}">·</td>`;
      const on = c.detected > 0;
      const bg = on ? `background:color-mix(in srgb, var(--vq-success) ${(shade(c.tpm_sum) * 100).toFixed(1)}%, transparent)` : '';
      const tip = `${k} · ${s.sample}\nTotal TPM: ${_fmtTpm(c.tpm_sum)}\nDetected: ${c.detected} / ${c.n} genes`;
      return `<td class="ov-kg-cell${on ? '' : ' ov-kg-cell--off'}" style="${bg}" title="${esc(tip)}">
                <span class="ov-kg-cell__v">${on ? _fmtTpm(c.tpm_sum) : '0'}</span>
                <span class="ov-kg-cell__n">${c.detected}/${c.n}</span>
              </td>`;
    }).join('');
    return `<tr><td class="ov-matrix__sample" title="${esc(k)}">${esc(_kgLabel(k))}</td>${cells}</tr>`;
  }).join('');

  host.innerHTML = `
    <div class="ov-matrix-wrap">
      <table class="ov-matrix">
        <thead>
          <tr>
            <th class="ov-matrix__sample">Kingdom</th>
            ${samples.map(s => `<th title="${esc(s.sample)}">${esc(_trunc(s.sample, 14))}</th>`).join('')}
          </tr>
        </thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
    <div class="ov-legend">
      <span class="ov-legend__item">cell: total TPM · detected/total genes</span>
      <span class="ov-legend__item">shading is log-scaled; unshaded = nothing detected</span>
    </div>`;
}

/* Kingdom keys arrive upper-case from the bundled FASTA prefixes. */
function _kgLabel(k) {
  return String(k || '').charAt(0) + String(k || '').slice(1).toLowerCase();
}

function _fmtTpm(v) {
  if (v == null || isNaN(v)) return '—';
  if (v === 0)   return '0';
  if (v >= 1000) return v.toFixed(0);
  if (v >= 10)   return v.toFixed(1);
  return v.toFixed(2);
}

// ── utils ───────────────────────────────────────────────────────────────────

function _trunc(s, n) {
  s = String(s ?? '');
  return s.length > n ? s.slice(0, n - 1) + '…' : s;
}

window.vqInitOverview = vqInitOverview;
})();
