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

      <!-- Pipeline workflow per sample — always full width -->
      ${_workflowCardHtml(samples)}

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
  _renderWorkflowMatrix(samples);
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

// ── Pipeline workflow per sample ────────────────────────────────────────────
// Same data as the per-sample "Pipeline Workflow" card (pipeline_stats.workflow):
// one lane per sample, one column per step, success path between executed
// steps, per-step time, run options and a click-through detail panel.
// Samples from runs older than v3.0.2 carry no workflow record; their lane
// falls back to the detected optional stages (NR, BLASTn, Salmon, LLM).

const _WF_ORDER = [
  'parse', 'orfs', 'refseq', 'hmm', 'nr', 'blastn', 'pfam', 'taxonomy',
  'clusters', 'salmon', 'seq_quality', 'coverage', 'heuristic', 'llm', 'export',
];
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
const _WF_OVERALL_TEXT = { done: 'Completed', partial: 'Warnings', error: 'Errors', legacy: 'No record' };

function _fmtDur(sec) {
  if (sec == null || isNaN(sec)) return '—';
  if (sec < 1)    return (sec * 1000).toFixed(0) + ' ms';
  if (sec < 60)   return sec.toFixed(1) + ' s';
  const r = Math.round(sec), m = Math.floor(r / 60), s = r % 60;
  if (r < 3600) return `${m}m ${String(s).padStart(2, '0')}s`;
  return `${Math.floor(m / 60)}h ${String(m % 60).padStart(2, '0')}m`;
}

const _wfStatus = s => _WF_STATUS[s] || _WF_STATUS.done;

function _wfPill(status, text) {
  const cls = _WF_STATUS[status] ? status : 'skipped';
  return `<span class="vq-wf-pill vq-wf-pill--${cls}">${VQ.esc(text ?? _wfStatus(status).label)}</span>`;
}

function _wfMiniRow(label, value, wide = false) {
  return `
    <div class="ov-wf-kv${wide ? ' ov-wf-kv--wide' : ''}">
      <span class="ov-wf-kv__k">${VQ.esc(label)}</span>
      <span class="ov-wf-kv__v">${value}</span>
    </div>`;
}

/** One lane per sample: {sample, legacy, steps: {key: step}, total, overall, options, ...}. */
function _wfLanes(samples) {
  return samples.map(s => {
    const wf = s.workflow;
    if (wf && (wf.steps || []).length) {
      const steps = {};
      wf.steps.forEach(st => { steps[st.key] = st; });
      const list  = wf.steps;
      const total = wf.total_seconds || list.reduce((a, st) => a + (st.seconds || 0), 0);
      const overall = list.some(st => st.status === 'error')   ? 'error'
                    : list.some(st => st.status === 'partial') ? 'partial' : 'done';
      return { sample: s.sample, legacy: false, steps, list, total, overall,
               options: wf.options || {}, started: wf.started_at, version: wf.version };
    }
    // Legacy run: only the optional stages detected from the results are known.
    const flags = s.steps || {};
    const steps = {};
    ['nr', 'blastn', 'salmon', 'llm'].forEach(k => {
      if (!(k in flags)) return;
      steps[k] = { key: k, label: _WF_SHORT[k], status: flags[k] ? 'done' : 'skipped',
                   seconds: null, details: {},
                   message: 'detected from results — no workflow record (run before v3.0.2)' };
    });
    return { sample: s.sample, legacy: true, steps, list: Object.values(steps), total: null,
             overall: 'legacy', options: { cap3: !!flags.cap3 }, started: null, version: null };
  });
}

function _wfColumns(lanes) {
  const present = new Set();
  lanes.forEach(l => Object.keys(l.steps).forEach(k => present.add(k)));
  const cols = _WF_ORDER.filter(k => present.has(k));
  present.forEach(k => { if (!cols.includes(k)) cols.push(k); });
  return cols;
}

/** Option chips for one sample — same set as the per-sample workflow card. */
function _wfOptionChips(o, legacy) {
  const esc = VQ.esc;
  const opt = (label, value, on = true) =>
    `<span class="vq-wf-opt${on ? '' : ' vq-wf-opt--off'}">
       <span class="vq-wf-opt__k">${esc(label)}</span>${esc(value)}</span>`;
  if (legacy) return opt('CAP3', o.cap3 ? 'on' : 'off', !!o.cap3);
  const reads = (o.reads || []).length
    ? `${o.reads.length} file${o.reads.length > 1 ? 's' : ''}${o.read_type ? ' · ' + o.read_type : ''}`
    : 'off';
  return [
    o.input ? opt('Input', o.input) : '',
    o.threads != null ? opt('Threads', o.threads) : '',
    opt('CAP3', o.cap3 ? 'on' : 'off', !!o.cap3),
    opt('NR', o.nr_db || 'off', !!o.nr_db),
    opt('BLASTn', o.blastn || 'off', !!o.blastn),
    opt('Reads', reads, !!(o.reads || []).length),
    o.transcriptome ? opt('Transcriptome', o.transcriptome) : '',
    opt('LLM', o.llm || 'off', !!o.llm),
    o.force ? opt('Force', 'on') : '',
  ].join('');
}

/** Cross-sample option uniformity: one chip per option, "mixed" when samples differ. */
function _wfSharedOptions(lanes) {
  const esc = VQ.esc;
  const recs = lanes.filter(l => !l.legacy);
  if (!recs.length) return '';
  const fields = [
    ['Threads', o => o.threads != null ? String(o.threads) : '—'],
    ['CAP3',    o => o.cap3 ? 'on' : 'off'],
    ['NR',      o => o.nr_db || 'off'],
    ['BLASTn',  o => o.blastn || 'off'],
    ['Reads',   o => (o.reads || []).length ? (o.read_type || 'on') : 'off'],
    ['LLM',     o => o.llm || 'off'],
  ];
  return fields.map(([label, get]) => {
    const byVal = new Map();
    recs.forEach(l => {
      const v = get(l.options);
      if (!byVal.has(v)) byVal.set(v, []);
      byVal.get(v).push(l.sample);
    });
    if (byVal.size === 1) {
      const v = [...byVal.keys()][0];
      const off = v === 'off' || v === '—';
      return `<span class="vq-wf-opt${off ? ' vq-wf-opt--off' : ''}">
                <span class="vq-wf-opt__k">${esc(label)}</span>${esc(v)}</span>`;
    }
    const tip = [...byVal.entries()]
      .map(([v, ss]) => `${v}: ${ss.length} sample${ss.length > 1 ? 's' : ''}`).join(' · ');
    return `<span class="vq-wf-opt vq-wf-opt--mixed" title="${esc(tip)}">
              <span class="vq-wf-opt__k">${esc(label)}</span>mixed (${byVal.size})</span>`;
  }).join('');
}

function _workflowCardHtml(samples) {
  const esc    = VQ.esc;
  const lanes  = _wfLanes(samples);
  const nRec   = lanes.filter(l => !l.legacy).length;
  const counts = { done: 0, partial: 0, error: 0, legacy: 0 };
  lanes.forEach(l => { counts[l.overall]++; });
  const issues = [];
  lanes.forEach((l, i) => l.list.forEach(st => {
    if (st.status === 'error' || st.status === 'partial') issues.push({ i, l, st });
  }));
  const shared = _wfSharedOptions(lanes);

  return `
    <div class="vq-chart-card vq-wf-card" id="ov-steps-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Pipeline Workflow per Sample</div>
          <div class="vq-chart-card__sub">
            ${nRec} of ${lanes.length} sample${lanes.length === 1 ? '' : 's'} with a workflow record
            &nbsp;·&nbsp; step status, time and options per run
            &nbsp;·&nbsp; click a step or a sample for details
          </div>
        </div>
        <div class="vq-wf-head-right ov-wf-head-pills">
          ${counts.done    ? _wfPill('done',    `${counts.done} completed`) : ''}
          ${counts.partial ? _wfPill('partial', `${counts.partial} with warnings`) : ''}
          ${counts.error   ? _wfPill('error',   `${counts.error} with errors`) : ''}
          ${counts.legacy  ? _wfPill('skipped', `${counts.legacy} without record`) : ''}
        </div>
      </div>

      ${shared ? `
      <div class="vq-wf-opts">
        <span class="ov-wf-opts__label">Run options</span>${shared}
      </div>` : ''}

      ${issues.length ? `
      <div class="vq-wf-issues">
        ${issues.map(({ i, l, st }) => `
          <button type="button" class="vq-wf-issue vq-wf-issue--${st.status}"
                  data-wf-lane="${i}" data-wf-key="${esc(st.key)}">
            <strong>${esc(l.sample)} · ${esc(_WF_SHORT[st.key] || st.label)}</strong>
            <span>${esc(st.message || _wfStatus(st.status).label)}</span>
          </button>`).join('')}
      </div>` : ''}

      <div class="ov-wf-matrix-wrap" id="ov-wf-matrix"></div>
      <div class="vq-wf-detail" id="ov-wf-detail"></div>
    </div>`;
}

function _renderWorkflowMatrix(samples) {
  const host = document.getElementById('ov-wf-matrix');
  if (!host) return;
  const esc   = VQ.esc;
  const lanes = _wfLanes(samples);
  const cols  = _wfColumns(lanes);
  const maxTotal = Math.max(0, ...lanes.map(l => l.total || 0));

  const ran = st => st && st.status !== 'skipped';

  const laneCells = (l, li) => {
    const idx = cols.map((k, c) => (ran(l.steps[k]) ? c : -1)).filter(c => c >= 0);
    const nextRan = c => idx.find(x => x >= c);
    return cols.map((k, c) => {
      const st = l.steps[k];
      // Success path: a segment between consecutive executed steps takes the
      // colour of the step it leads into (as in the per-sample flowchart).
      const seg = (from, to) => {
        const n = nextRan(to);
        if (n == null || !idx.some(x => x <= from)) return '';
        const s2 = l.steps[cols[n]].status;
        return `style="--ov-wf-c:${_wfStatus(s2).color};opacity:${s2 === 'done' ? 0.55 : 0.85}"`;
      };
      const left  = idx.some(x => x < c) && idx.some(x => x >= c) ? seg(c - 1, c) : '';
      const right = idx.some(x => x <= c) && idx.some(x => x > c) ? seg(c, c + 1) : '';
      const lines = `${left ? `<span class="ov-wf-line ov-wf-line--l" ${left}></span>` : ''}
                     ${right ? `<span class="ov-wf-line ov-wf-line--r" ${right}></span>` : ''}`;
      if (!st) return `<div class="ov-wf-cell ov-wf-cell--none">${lines}<span class="ov-wf-dot"></span></div>`;
      const time = st.status === 'skipped' ? 'skipped' : (l.legacy ? 'ran' : _fmtDur(st.seconds));
      return `
        <button type="button" class="ov-wf-cell ov-wf-cell--${st.status}${l.legacy ? ' ov-wf-cell--legacy' : ''}"
                data-wf-lane="${li}" data-wf-key="${esc(k)}">
          ${lines}
          <span class="ov-wf-dot" style="--ov-wf-c:${_wfStatus(st.status).color}">${_wfStatus(st.status).glyph}</span>
          <span class="ov-wf-time">${esc(time)}</span>
        </button>`;
    }).join('');
  };

  host.innerHTML = `
    <div class="ov-wf-matrix" style="grid-template-columns:minmax(120px,max-content) repeat(${cols.length},minmax(62px,1fr)) minmax(96px,max-content)">
      <div class="ov-wf-hd ov-wf-hd--sample">Sample</div>
      ${cols.map(k => `<div class="ov-wf-hd">${esc(_WF_SHORT[k] || k)}</div>`).join('')}
      <div class="ov-wf-hd ov-wf-hd--total">Total</div>
      ${lanes.map((l, li) => `
        <button type="button" class="ov-wf-sample" data-wf-lane="${li}" title="${esc(l.sample)}">
          <span class="ov-wf-sample__name">${esc(l.sample)}</span>
          ${_wfPill(l.overall === 'legacy' ? 'skipped' : l.overall, _WF_OVERALL_TEXT[l.overall])}
        </button>
        ${laneCells(l, li)}
        <div class="ov-wf-total" data-wf-lane="${li}">
          <span>${l.total != null ? _fmtDur(l.total) : '—'}</span>
          ${l.total != null && maxTotal > 0
            ? `<span class="ov-wf-total__bar"><span style="width:${(100 * l.total / maxTotal).toFixed(1)}%"></span></span>`
            : ''}
        </div>`).join('')}
    </div>`;

  // ── Hover tooltip ────────────────────────────────────────────────────────
  host.querySelectorAll('button.ov-wf-cell').forEach(btn => {
    const l  = lanes[+btn.dataset.wfLane];
    const st = l.steps[btn.dataset.wfKey];
    btn.addEventListener('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">${esc(l.sample)} · ${esc(st.label)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Status</span><span>${_wfStatus(st.status).label}</span>
        ${st.seconds != null ? `<span class="vq-tooltip__key">Time</span><span>${_fmtDur(st.seconds)}</span>` : ''}
        ${st.message ? `<span class="vq-tooltip__key">Note</span><span>${esc(st.message)}</span>` : ''}
      </div>`, evt));
    btn.addEventListener('mouseleave', VQ.tooltipHide);
  });

  // ── Detail panel: sample run (options, time per step) + selected step ────
  const detail = document.getElementById('ov-wf-detail');
  const card   = document.getElementById('ov-steps-card');

  function select(li, key) {
    const l = lanes[li];
    if (!l || !detail) return;
    card.querySelectorAll('.ov-wf-cell, .ov-wf-sample').forEach(el =>
      el.classList.toggle('is-active',
        +el.dataset.wfLane === li && (el.classList.contains('ov-wf-sample') ? !key : el.dataset.wfKey === key)));
    card.querySelectorAll('.ov-wf-sample').forEach(el =>
      el.classList.toggle('is-lane', +el.dataset.wfLane === li));

    const st = key ? l.steps[key] : null;
    const timed = l.list.filter(x => x.seconds != null && x.seconds > 0);
    const longest = timed.length ? timed.reduce((a, x) => (x.seconds > a.seconds ? x : a)) : null;
    const runN  = l.list.filter(ran).length;

    // Same step across the other samples (time spread).
    const across = key
      ? lanes.map(x => x.steps[key]).filter(x => x && x.seconds != null && x.status !== 'skipped')
             .map(x => x.seconds).sort((a, b) => a - b)
      : [];
    const med = across.length ? across[Math.floor((across.length - 1) / 2)] : null;

    detail.innerHTML = `
      <div class="vq-wf-detail__head">
        <div class="vq-wf-detail__title">${esc(l.sample)}</div>
        ${_wfPill(l.overall === 'legacy' ? 'skipped' : l.overall,
                  l.legacy ? 'No workflow record' : _WF_OVERALL_TEXT[l.overall])}
        <span class="vq-wf-detail__time">
          ${l.legacy ? 'run before v3.0.2' : `
            ${l.started ? 'started ' + esc(new Date(l.started).toLocaleString()) + ' · ' : ''}
            ${runN} step${runN === 1 ? '' : 's'} run · ${_fmtDur(l.total)}
            ${l.version ? ' · ViralQuest v' + esc(l.version) : ''}`}
        </span>
      </div>
      <div class="vq-wf-opts">${_wfOptionChips(l.options, l.legacy)}</div>
      ${timed.length && l.total > 0 ? `
      <div class="vq-wf-timebar">
        <div class="vq-wf-timebar__track">
          ${timed.map((x, i) => `
            <div class="vq-wf-timebar__seg${x.key === key ? ' is-active' : ''}" data-wf-key="${esc(x.key)}"
                 style="flex:${x.seconds} 1 0;opacity:${x.key === key ? 1 : (i % 2 ? 0.55 : 0.9)}"></div>`).join('')}
        </div>
        <div class="vq-wf-timebar__legend">
          <span>time per step</span>
          <span>longest: <strong>${esc(_WF_SHORT[longest.key] || longest.label)}</strong>
            · ${_fmtDur(longest.seconds)} (${VQ.pct(longest.seconds, l.total)})</span>
        </div>
      </div>` : ''}
      ${st ? `
      <div class="ov-wf-step">
        <div class="vq-wf-detail__head">
          <div class="vq-wf-detail__title">${esc(st.label)}</div>
          ${_wfPill(st.status)}
          ${st.seconds != null
            ? `<span class="vq-wf-detail__time">${_fmtDur(st.seconds)}${l.total ? ' · ' + VQ.pct(st.seconds, l.total) + ' of total' : ''}</span>`
            : ''}
        </div>
        ${st.message ? `<div class="vq-wf-detail__msg vq-wf-detail__msg--${st.status}">${esc(st.message)}</div>` : ''}
        <div class="vq-wf-detail__grid">
          ${Object.entries(st.details || {}).map(([k, v]) => _wfMiniRow(k,
              esc(typeof v === 'number' ? v.toLocaleString() : String(v ?? '—')))).join('')}
          ${across.length > 1 ? _wfMiniRow(`Time across ${across.length} samples`,
              `min ${_fmtDur(across[0])} · median ${_fmtDur(med)} · max ${_fmtDur(across[across.length - 1])}`, true) : ''}
        </div>
      </div>` : ''}`;

    detail.querySelectorAll('.vq-wf-timebar__seg').forEach(seg => {
      const x = l.steps[seg.dataset.wfKey];
      seg.addEventListener('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(x.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Time</span><span>${_fmtDur(x.seconds)}</span>
          <span class="vq-tooltip__key">of total</span><span>${VQ.pct(x.seconds, l.total)}</span>
        </div>`, evt));
      seg.addEventListener('mouseleave', VQ.tooltipHide);
      seg.addEventListener('click', () => { VQ.tooltipHide(); select(li, x.key); });
    });
  }

  card.querySelectorAll('.ov-wf-cell[data-wf-key], .vq-wf-issue').forEach(el =>
    el.addEventListener('click', () => select(+el.dataset.wfLane, el.dataset.wfKey)));
  card.querySelectorAll('.ov-wf-sample').forEach(el =>
    el.addEventListener('click', () => select(+el.dataset.wfLane, null)));

  // Open on the first problem, otherwise on the slowest sample's run.
  let li0 = lanes.findIndex(l => l.overall === 'error');
  if (li0 < 0) li0 = lanes.findIndex(l => l.overall === 'partial');
  if (li0 >= 0) {
    const st = lanes[li0].list.find(x => x.status === 'error')
            || lanes[li0].list.find(x => x.status === 'partial');
    select(li0, st.key);
  } else {
    const slowest = lanes.reduce((a, l, i) => ((l.total || 0) > (lanes[a].total || 0) ? i : a), 0);
    select(slowest, null);
  }
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
