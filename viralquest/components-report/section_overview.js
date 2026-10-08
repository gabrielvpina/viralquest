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

function vqInitOverview(report) {
  const el = document.getElementById('section-overview');
  if (!el) return;

  const esc      = VQ.esc;
  const samples  = report.samples || [];

  if (!samples.length) {
    el.innerHTML = `<div class="vq-empty">No samples loaded.</div>`;
    return;
  }

  // ── Aggregates ──────────────────────────────────────────────────────────
  const NOV = window.vqNovelty;
  const seqs = report.sequences || [];
  const bySample = new Map(samples.map(s => [s.sample, []]));
  seqs.forEach(q => { if (bySample.has(q.sample)) bySample.get(q.sample).push(q); });

  const totalConfirmed = samples.reduce((a, s) => a + (s.n_confirmed || 0), 0);
  const familyUnion = new Set();
  samples.forEach(s => Object.keys(s.families || {}).forEach(f => familyUnion.add(f)));
  const coreFamilies = [...familyUnion].filter(f =>
    f !== 'Unclassified' && samples.every(s => (s.families || {})[f]));

  // Species: best BLASTx hit (NR first, RefSeq fallback), as in report_data.
  const speciesSamples = new Map();
  seqs.forEach(q => {
    const hits = (q.blastx_nr_hits || []).length ? q.blastx_nr_hits : (q.blastx_hits || []);
    const sp = hits.length ? hits[0].species : null;
    if (!sp) return;
    if (!speciesSamples.has(sp)) speciesSamples.set(sp, new Set());
    speciesSamples.get(sp).add(q.sample);
  });
  const sharedSpecies = [...speciesSamples.values()].filter(ss => ss.size > 1).length;

  const novelOf = list => NOV ? list.filter(q => NOV.novel.includes(NOV.tier(q))).length : null;
  samples.forEach(s => { s._novel = novelOf(bySample.get(s.sample) || []); });
  const totalNovel = NOV ? samples.reduce((a, s) => a + (s._novel || 0), 0) : null;

  const confirmedPer = samples.map(s => s.n_confirmed || 0).sort((a, b) => a - b);
  const medConfirmed = confirmedPer[Math.floor((confirmedPer.length - 1) / 2)];

  const anyLLM  = samples.some(s => (s.llm_scores || []).length);
  const anyHeur = samples.some(s => (s.heuristic_scores || []).length);
  const salmonSamples = samples.filter(s => s.salmon);
  const anySalmon = salmonSamples.length > 0;
  const anyConserved = salmonSamples.some(s => Object.keys(s.salmon.conserved || {}).length);
  const nWorkflow = samples.filter(s => s.workflow).length;

  const fmt = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();

  const group = (title, hint, body) => `
    <section class="vq-group">
      <div class="vq-group__head">
        <h3 class="vq-group__title">${esc(title)}</h3>
        ${hint ? `<span class="vq-group__hint">${esc(hint)}</span>` : ''}
      </div>
      <div class="ov-stack">${body}</div>
    </section>`;

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
        ${_chip('Samples', fmt(samples.length), 'accent',
                [`${nWorkflow} with workflow record`, anySalmon ? `${salmonSamples.length} with Salmon` : '']
                  .filter(Boolean).join(' · '))}
        ${_chip('Confirmed Viral', fmt(totalConfirmed), 'success',
                `median ${fmt(medConfirmed)} per sample (${fmt(confirmedPer[0])}–${fmt(confirmedPer[confirmedPer.length - 1])})`)}
        ${_chip('Viral Species', fmt(speciesSamples.size), '',
                `best BLASTx hit · ${fmt(sharedSpecies)} shared by ≥2 samples`)}
        ${totalNovel != null ? _chip('Putative Novel', fmt(totalNovel), 'accent',
                `${VQ.pct(totalNovel, totalConfirmed)} of confirmed · novel species + divergent`) : ''}
        ${_chip('Viral Families', fmt(familyUnion.size), '',
                `${fmt(coreFamilies.length)} found in every sample`)}
      </div>

      ${group('Run', 'pipeline steps, status and options per sample', _workflowCardHtml(samples))}

      ${group('Detection & composition', 'from input contigs to putative novel viruses, and what they are', `
        ${_funnelCardHtml(NOV)}
        ${_bubbleCardHtml(samples)}`)}

      ${anySalmon ? group('Read support', 'Salmon quantification — how much of each library is viral', `
        ${_salmonCardHtml(salmonSamples)}
        ${anyConserved ? `
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
`) : ''}

      ${group('Sequence properties', 'length and score distributions per sample', `
        ${_card('ov-lengths', 'Sequence Length per Sample',
                 'density of confirmed contig lengths (log scale) · solid = median · dashed = mean')}
        ${anyHeur || anyLLM ? `
        <div class="vq-grid${anyHeur && anyLLM ? ' vq-grid--2' : ' ov-grid--1'}">
          ${anyHeur ? _card('ov-heur', 'Heuristic Score per Sample',
                             'VQ score · box = Q1–Q3, whiskers 1.5 × IQR, ◆ mean, dots = every contig') : ''}
          ${anyLLM  ? _card('ov-llm',  'LLM Score per Sample',
                             'VQ score · box = Q1–Q3, whiskers 1.5 × IQR, ◆ mean, dots = every contig') : ''}
        </div>` : ''}`)}
    </div>
  `;

  // ── Render charts ───────────────────────────────────────────────────────
  _renderWorkflowMatrix(samples);
  _renderFunnel(samples, NOV);
  _renderBubbleMatrix(samples, bySample, NOV);
  if (anySalmon) _renderSalmon(salmonSamples);
  if (anyConserved) _renderKingdomMatrix(_body('ov-kingdoms'), salmonSamples);
  _renderLengthDensity(samples);
  if (anyHeur) _renderScoreBoxes('ov-heur', samples, 'heuristic_scores');
  if (anyLLM)  _renderScoreBoxes('ov-llm',  samples, 'llm_scores');
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
  'parse', 'orfs', 'refseq', 'hmm', 'nr', 'blastn', 'pfam', 'synteny', 'taxonomy',
  'clusters', 'salmon', 'seq_quality', 'coverage', 'heuristic', 'llm', 'export',
];
const _WF_SHORT = {
  parse: 'Parse FASTA', orfs: 'ORFs', refseq: 'RefSeq', hmm: 'HMM filter',
  nr: 'Diamond NR', blastn: 'BLASTn', pfam: 'Pfam', synteny: 'Synteny', taxonomy: 'Taxonomy',
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

// ── Detection funnel per sample ─────────────────────────────────────────────
// One horizontal funnel per sample: input contigs → viral-flagged (RefSeq +
// HMM) → confirmed → putative novel.  Bar height is log-scaled against the
// largest input so samples stay comparable; the band between two stages
// carries the retention rate.

function _funnelStages(NOV) {
  return [
    { key: 'n_input',     label: 'Input contigs',  color: 'var(--vq-text-3)' },
    { key: 'n_viral',     label: 'Viral-flagged',  color: 'var(--vq-primary-light)' },
    { key: 'n_confirmed', label: 'Confirmed',      color: 'var(--vq-success)' },
    ...(NOV ? [{ key: '_novel', label: 'Putative novel', color: NOV.color['novel-species'] }] : []),
  ];
}

function _funnelCardHtml(NOV) {
  return `
    <div class="vq-chart-card" id="ov-funnel-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Detection Funnel per Sample</div>
          <div class="vq-chart-card__sub">
            input contigs → viral-flagged (RefSeq + HMM) → confirmed${NOV ? ' → putative novel' : ''}
            &nbsp;·&nbsp; bar height on a log scale, band label = share kept from the previous stage
          </div>
        </div>
      </div>
      <div class="vq-chart-card__body ov-funnel" id="ov-funnel-body"></div>
    </div>`;
}

function _renderFunnel(samples, NOV) {
  const host = _body('ov-funnel');
  if (!host) return;
  const draw = () => _drawFunnel(host, samples, NOV);
  draw();
  VQ.redrawOnResize(host, draw);
}

function _drawFunnel(host, samples, NOV) {
  const esc    = VQ.esc;
  const stages = _funnelStages(NOV);
  const val    = (s, st) => { const v = s[st.key]; return v == null || isNaN(v) ? null : +v; };
  const maxN   = d3.max(samples, s => d3.max(stages, st => val(s, st) || 0)) || 1;

  const W = Math.max(host.clientWidth || 0, 640);
  const padL = 150, padR = 16, headH = 28, rowH = 48, barW = 8;
  const H = headH + samples.length * rowH + 4;
  const lastW = 110;                               // last stage only carries its count
  const colW = (W - padL - padR - lastW) / Math.max(1, stages.length - 1);
  const xBar = i => padL + i * colW + 6;
  const hOf  = n => n == null ? 0 : Math.max(3, (Math.log10(n + 1) / Math.log10(maxN + 1)) * (rowH - 10));

  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');

  stages.forEach((st, i) => {
    svg.append('text').attr('class', 'ov-funnel__stage')
      .attr('x', xBar(i)).attr('y', 14).text(st.label);
    svg.append('rect').attr('x', xBar(i)).attr('y', 19).attr('width', 18).attr('height', 3)
      .attr('rx', 1.5).attr('fill', st.color);
  });

  samples.forEach((s, r) => {
    const mid = headH + r * rowH + rowH / 2;
    const g = svg.append('g').attr('class', 'ov-funnel__row');
    g.append('rect').attr('class', 'ov-funnel__hover')
      .attr('x', 0).attr('y', mid - rowH / 2).attr('width', W).attr('height', rowH);
    if (r) g.append('line').attr('class', 'ov-grid')
      .attr('x1', 0).attr('x2', W).attr('y1', mid - rowH / 2).attr('y2', mid - rowH / 2);

    g.append('text').attr('class', 'ov-funnel__sample')
      .attr('x', 8).attr('y', mid).attr('dominant-baseline', 'central')
      .text(_trunc(s.sample, 20));

    const vals = stages.map(st => val(s, st));
    // Bands between consecutive stages (drawn first, under the bars).
    stages.forEach((st, i) => {
      if (!i || vals[i] == null || vals[i - 1] == null) return;
      const h0 = hOf(vals[i - 1]), h1 = hOf(vals[i]);
      const x0 = xBar(i - 1) + barW, x1 = xBar(i);
      g.append('path')
        .attr('d', `M${x0},${mid - h0 / 2} L${x1},${mid - h1 / 2} L${x1},${mid + h1 / 2} L${x0},${mid + h0 / 2} Z`)
        .attr('fill', st.color).attr('opacity', 0.13);
      g.append('text').attr('class', 'ov-funnel__rate')
        .attr('x', x1 - 8).attr('y', mid).attr('text-anchor', 'end').attr('dominant-baseline', 'central')
        .text(VQ.pct(vals[i], vals[i - 1]));
    });
    stages.forEach((st, i) => {
      const x = xBar(i), h = hOf(vals[i]);
      if (vals[i] != null) g.append('rect')
        .attr('x', x).attr('y', mid - h / 2).attr('width', barW).attr('height', h)
        .attr('rx', 2).attr('fill', st.color);
      g.append('text').attr('class', 'ov-funnel__count')
        .attr('x', x + barW + 6).attr('y', mid).attr('dominant-baseline', 'central')
        .text(vals[i] == null ? '—' : vals[i].toLocaleString());
    });

    g.on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(s.sample)}</div>
        <div class="vq-tooltip__row">
          ${stages.map((st, i) => `
            <span class="vq-tooltip__key">${esc(st.label)}</span>
            <span>${vals[i] == null ? '—' : vals[i].toLocaleString()}${
              i && vals[i] != null && vals[i - 1] != null ? ` · ${VQ.pct(vals[i], vals[i - 1])}` : ''}</span>`).join('')}
        </div>`, evt))
     .on('mouseleave', VQ.tooltipHide);
  });
}

// ── Family × sample bubble matrix ───────────────────────────────────────────
// Rows = viral families (top N by prevalence, then abundance), columns =
// samples.  Bubble area ∝ the chosen measure; colour = the family's dominant
// novelty tier in that sample.  The right margin shows prevalence (in how many
// samples the family occurs) — core vs sample-specific families at a glance.

const _BUBBLE_TOP = 15;

function _bubbleCardHtml(samples) {
  const anyTpm = samples.some(s => s.salmon);
  return `
    <div class="vq-chart-card" id="ov-bubbles-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Viral Families across Samples</div>
          <div class="vq-chart-card__sub">
            bubble area ∝ value · colour = dominant novelty tier · right: prevalence across samples
          </div>
        </div>
        <div class="vq-toggle" id="ov-bubbles-mode" role="tablist" aria-label="Bubble size">
          <button class="vq-toggle__btn active" type="button" data-mode="n">Sequences</button>
          <button class="vq-toggle__btn" type="button" data-mode="pct">% of sample</button>
          ${anyTpm ? `<button class="vq-toggle__btn" type="button" data-mode="tpm">TPM</button>` : ''}
        </div>
      </div>
      <div class="vq-chart-card__body ov-bubbles" id="ov-bubbles-body"></div>
      <div class="ov-legend" id="ov-bubbles-legend"></div>
    </div>`;
}

function _bubbleData(samples, bySample, NOV) {
  const cells = new Map();                       // family → sample → cell
  const totals = new Map(samples.map(s => [s.sample, 0]));
  samples.forEach(s => (bySample.get(s.sample) || []).forEach(q => {
    const fam = (q.taxonomy || {}).family || 'Unclassified';
    if (!cells.has(fam)) cells.set(fam, new Map());
    const row = cells.get(fam);
    if (!row.has(s.sample)) row.set(s.sample, { n: 0, tpm: 0, nov: {} });
    const c = row.get(s.sample);
    c.n += 1;
    c.tpm += q.tpm || 0;
    if (NOV) { const t = NOV.tier(q); c.nov[t] = (c.nov[t] || 0) + 1; }
    totals.set(s.sample, totals.get(s.sample) + 1);
  }));

  const famStats = [...cells.entries()].map(([fam, row]) => ({
    fam, row,
    prev:  row.size,
    total: [...row.values()].reduce((a, c) => a + c.n, 0),
  }));
  const ranked = famStats.filter(f => f.fam !== 'Unclassified')
    .sort((a, b) => b.prev - a.prev || b.total - a.total || a.fam.localeCompare(b.fam));
  const rows = ranked.slice(0, _BUBBLE_TOP);

  // Fold the tail into "Other families"; keep "Unclassified" as its own last row.
  const merge = (label, list) => {
    const row = new Map();
    list.forEach(f => f.row.forEach((c, smp) => {
      if (!row.has(smp)) row.set(smp, { n: 0, tpm: 0, nov: {} });
      const m = row.get(smp);
      m.n += c.n; m.tpm += c.tpm;
      Object.entries(c.nov).forEach(([t, k]) => { m.nov[t] = (m.nov[t] || 0) + k; });
    }));
    return { fam: label, row, prev: row.size, total: list.reduce((a, f) => a + f.total, 0), folded: list.length };
  };
  const tail = ranked.slice(_BUBBLE_TOP);
  if (tail.length) rows.push(merge(`Other families (${tail.length})`, tail));
  const uncl = famStats.find(f => f.fam === 'Unclassified');
  if (uncl) rows.push(uncl);
  return { rows, totals };
}

function _renderBubbleMatrix(samples, bySample, NOV) {
  const host = _body('ov-bubbles');
  if (!host) return;
  const data = _bubbleData(samples, bySample, NOV);
  let mode = 'n';
  const draw = () => _drawBubbles(host, samples, data, NOV, mode);

  document.querySelectorAll('#ov-bubbles-mode .vq-toggle__btn').forEach(btn =>
    btn.addEventListener('click', () => {
      mode = btn.dataset.mode;
      document.querySelectorAll('#ov-bubbles-mode .vq-toggle__btn')
        .forEach(b => b.classList.toggle('active', b === btn));
      VQ.tooltipHide();
      draw();
    }));

  // Legend: novelty tiers present + size note.
  const legend = document.getElementById('ov-bubbles-legend');
  if (legend && NOV) {
    const present = new Set();
    data.rows.forEach(r => r.row.forEach(c => Object.keys(c.nov).forEach(t => present.add(t))));
    legend.innerHTML = NOV.order.filter(t => present.has(t)).map(t =>
      `<span class="ov-legend__item" title="${VQ.esc(NOV.long[t])}">
         <span class="ov-legend__dot" style="background:${NOV.color[t]}"></span>${VQ.esc(NOV.label[t])}</span>`).join('')
      + `<span class="ov-legend__item">dominant tier of the family in that sample</span>`;
  }

  if (!data.rows.length) { host.innerHTML = `<div class="vq-empty">No family data.</div>`; return; }
  draw();
  VQ.redrawOnResize(host, draw);
}

function _drawBubbles(host, samples, data, NOV, mode) {
  const esc = VQ.esc;
  const { rows, totals } = data;
  const value = (c, smp) => !c ? 0
    : mode === 'pct' ? (totals.get(smp) ? 100 * c.n / totals.get(smp) : 0)
    : mode === 'tpm' ? c.tpm
    : c.n;
  const fmtV = v => mode === 'pct' ? v.toFixed(1) + '%'
    : mode === 'tpm' ? _fmtTpm(v) : v.toLocaleString();
  const dominant = c => {
    const e = Object.entries(c.nov);
    return e.length ? e.reduce((a, b) => (b[1] > a[1] ? b : a))[0] : null;
  };

  const avail = host.clientWidth || 900;
  const padL = 190, padR = 150, rowH = 30;
  const longest = d3.max(samples, s => s.sample.length) || 6;
  const minCol = 34;
  const colW = Math.max(minCol, Math.min(150, (avail - padL - padR) / samples.length));
  const rotate = colW < longest * 6.6 + 8;
  const headH = rotate ? Math.min(120, longest * 5.2 + 22) : 30;
  const W = Math.max(avail, padL + padR + colW * samples.length);
  const gridW = colW * samples.length;
  const H = headH + rows.length * rowH + 6;

  const maxV = d3.max(rows, r => d3.max(samples, s => value(r.row.get(s.sample), s.sample))) || 1;
  const rMax = Math.min(colW, rowH) / 2 - 1.5;
  const rOf  = d3.scaleSqrt().domain([0, maxV]).range([0, rMax]);
  const xOf  = i => padL + i * colW + colW / 2;
  const yOf  = j => headH + j * rowH + rowH / 2;

  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');

  // Column headers (rotated when samples are many / names long).
  samples.forEach((s, i) => {
    const t = svg.append('text').attr('class', 'ov-bubbles__col').text(_trunc(s.sample, 22));
    if (rotate) t.attr('transform', `translate(${xOf(i) + 3},${headH - 8}) rotate(-40)`);
    else t.attr('x', xOf(i)).attr('y', headH - 10).attr('text-anchor', 'middle');
    t.append('title').text(s.sample);
  });
  svg.append('text').attr('class', 'ov-bubbles__col ov-bubbles__col--meta')
    .attr('x', padL + gridW + 14).attr('y', headH - 10).text('prevalence · total');

  rows.forEach((r, j) => {
    const y = yOf(j);
    const g = svg.append('g').attr('class', 'ov-bubbles__row');
    g.append('rect').attr('class', 'ov-funnel__hover')
      .attr('x', 0).attr('y', y - rowH / 2).attr('width', W).attr('height', rowH);
    g.append('line').attr('class', 'ov-bubbles__guide')
      .attr('x1', padL).attr('x2', padL + gridW).attr('y1', y).attr('y2', y);
    g.append('text').attr('class', 'ov-bubbles__fam' + (r.folded || r.fam === 'Unclassified' ? ' is-muted' : ''))
      .attr('x', padL - 12).attr('y', y).attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(_trunc(r.fam, 28)).append('title').text(r.fam);

    samples.forEach((s, i) => {
      const c = r.row.get(s.sample);
      const v = value(c, s.sample);
      if (!c || v <= 0) {
        g.append('circle').attr('cx', xOf(i)).attr('cy', y).attr('r', 1.6).attr('class', 'ov-bubbles__empty');
        return;
      }
      const tier = NOV ? dominant(c) : null;
      g.append('circle')
        .attr('class', 'ov-bubbles__dot')
        .attr('cx', xOf(i)).attr('cy', y).attr('r', Math.max(2.5, rOf(v)))
        .attr('fill', tier ? NOV.color[tier] : 'var(--vq-accent)')
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(r.fam)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Sample</span><span>${esc(s.sample)}</span>
            <span class="vq-tooltip__key">Sequences</span><span>${c.n.toLocaleString()}</span>
            <span class="vq-tooltip__key">of sample</span><span>${VQ.pct(c.n, totals.get(s.sample))}</span>
            ${s.salmon ? `<span class="vq-tooltip__key">TPM</span><span>${_fmtTpm(c.tpm)}</span>` : ''}
            ${NOV ? NOV.order.filter(t => c.nov[t]).map(t =>
              `<span class="vq-tooltip__key">${esc(NOV.label[t])}</span><span>${c.nov[t]}</span>`).join('') : ''}
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
      if (rOf(v) >= 11) g.append('text').attr('class', 'ov-bubbles__val')
        .attr('x', xOf(i)).attr('y', y).attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
        .text(mode === 'pct' ? Math.round(v) : (mode === 'tpm' ? '' : v));
    });

    // Prevalence bar + total in the selected measure.
    const px = padL + gridW + 14, pw = 60;
    g.append('rect').attr('x', px).attr('y', y - 3).attr('width', pw).attr('height', 6)
      .attr('rx', 3).attr('class', 'ov-bubbles__track');
    g.append('rect').attr('x', px).attr('y', y - 3).attr('width', pw * r.prev / samples.length)
      .attr('height', 6).attr('rx', 3)
      .attr('fill', r.prev === samples.length ? 'var(--vq-success)' : 'var(--vq-accent)');
    const tot = mode === 'pct' ? null
      : samples.reduce((a, s) => a + value(r.row.get(s.sample), s.sample), 0);
    g.append('text').attr('class', 'ov-bubbles__prev')
      .attr('x', px + pw + 8).attr('y', y).attr('dominant-baseline', 'central')
      .text(`${r.prev}/${samples.length}${tot != null ? ' · ' + fmtV(tot) : ''}`);
  });
}

// ── Salmon read support per sample ──────────────────────────────────────────
// Where each library's reads went: viral contigs / other indexed targets
// (housekeeping, host) / unmapped, plus the viral load in reads per million.
// Viral reads are usually a tiny share, so the load is shown on a log scale.

function _salmonCardHtml(salmonSamples) {
  return `
    <div class="vq-chart-card" id="ov-salmon-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Read Mapping &amp; Viral Load</div>
          <div class="vq-chart-card__sub">
            Salmon · ${salmonSamples.length} sample${salmonSamples.length === 1 ? '' : 's'}
            &nbsp;·&nbsp; read fate per library and viral reads per million (log scale)
          </div>
        </div>
      </div>
      <div class="ov-salmon" id="ov-salmon-body"></div>
      <div class="ov-legend">
        <span class="ov-legend__item"><span class="ov-legend__dot" style="background:var(--vq-accent)"></span>viral contigs</span>
        <span class="ov-legend__item"><span class="ov-legend__dot" style="background:var(--vq-primary-light);opacity:.45"></span>other targets (housekeeping · host)</span>
        <span class="ov-legend__item"><span class="ov-legend__dot" style="background:var(--vq-c-track);border:1px solid var(--vq-border)"></span>unmapped</span>
      </div>
    </div>`;
}

function _renderSalmon(salmonSamples) {
  const host = _body('ov-salmon');
  if (!host) return;
  const esc = VQ.esc;
  const rows = salmonSamples.map(s => {
    const q      = s.salmon;
    const rate   = q.mapping_rate ?? null;                 // % of input reads
    const mapped = q.total_reads ?? null;                  // reads assigned to any target
    const input  = mapped != null && rate ? mapped / (rate / 100) : null;
    const viral  = q.viral_reads ?? null;
    const rpm    = viral != null && input ? viral / input * 1e6 : null;
    return { s, rate, mapped, input, viral, rpm };
  });
  const LOG_MAX = 6;                                       // 10^6 RPM = every read viral
  const fmtN = n => n == null ? '—' : Math.round(n).toLocaleString();
  const fmtRpm = v => v == null ? '—' : v >= 100 ? Math.round(v).toLocaleString() : v.toFixed(1);

  host.innerHTML = `
    <div class="ov-salmon__grid">
      <div class="ov-salmon__hd">Sample</div>
      <div class="ov-salmon__hd ov-salmon__hd--num">Reads</div>
      <div class="ov-salmon__hd">Read fate</div>
      <div class="ov-salmon__hd ov-salmon__hd--num">Mapped</div>
      <div class="ov-salmon__hd ov-salmon__hd--num">Viral reads</div>
      <div class="ov-salmon__hd">Viral load (reads per million)</div>
      ${rows.map((r, i) => {
        const pv = r.viral != null && r.input ? 100 * r.viral / r.input : 0;
        const pm = r.rate ?? 0;
        const loadW = r.rpm ? Math.max(1, 100 * Math.min(LOG_MAX, Math.log10(r.rpm + 1)) / LOG_MAX) : 0;
        return `
        <div class="ov-salmon__sample" title="${esc(r.s.sample)}"><span>${esc(r.s.sample)}</span></div>
        <div class="ov-salmon__num">${r.input != null ? fmtN(r.input) : '—'}</div>
        <div class="ov-salmon__bar" data-i="${i}"><span class="ov-salmon__track">
          <span class="ov-salmon__seg ov-salmon__seg--viral" style="width:${pv > 0 ? Math.max(pv, 0.8) : 0}%"></span>
          <span class="ov-salmon__seg ov-salmon__seg--other" style="width:${Math.max(0, pm - Math.max(pv, pv > 0 ? 0.8 : 0))}%"></span>
        </span></div>
        <div class="ov-salmon__num"><strong>${r.rate != null ? r.rate.toFixed(1) + '%' : '—'}</strong></div>
        <div class="ov-salmon__num">${fmtN(r.viral)}
          <span class="ov-salmon__muted">${r.input ? VQ.pct(r.viral, r.input) : ''}</span></div>
        <div class="ov-salmon__load" data-i="${i}">
          <span class="ov-salmon__loadtrack">
            ${[1, 2, 3, 4, 5].map(k => `<i style="left:${100 * k / LOG_MAX}%"></i>`).join('')}
            <span class="ov-salmon__loadfill" style="width:${loadW}%"></span>
          </span>
          <span class="ov-salmon__loadval">${fmtRpm(r.rpm)}</span>
        </div>`;
      }).join('')}
      <div></div><div></div><div></div><div></div><div></div>
      <div class="ov-salmon__axis">
        ${['1', '10', '100', '1k', '10k', '100k', '1M'].map((t, k) =>
          `<span style="left:${100 * k / LOG_MAX}%">${t}</span>`).join('')}
      </div>
    </div>`;

  host.querySelectorAll('[data-i]').forEach(el => {
    const r = rows[+el.dataset.i];
    el.addEventListener('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">${esc(r.s.sample)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Input reads</span><span>${fmtN(r.input)} (estimated)</span>
        <span class="vq-tooltip__key">Mapped</span><span>${fmtN(r.mapped)} · ${r.rate != null ? r.rate.toFixed(1) + '%' : '—'}</span>
        <span class="vq-tooltip__key">Viral reads</span><span>${fmtN(r.viral)} · ${r.input ? VQ.pct(r.viral, r.input) : '—'}</span>
        <span class="vq-tooltip__key">of mapped</span><span>${r.mapped ? VQ.pct(r.viral, r.mapped) : '—'}</span>
        <span class="vq-tooltip__key">Viral load</span><span>${fmtRpm(r.rpm)} RPM</span>
        <span class="vq-tooltip__key">Viral TPM</span><span>${_fmtTpm(r.s.salmon.viral_tpm_sum)}</span>
      </div>`, evt));
    el.addEventListener('mouseleave', VQ.tooltipHide);
  });
}

// ── Shared helpers for the per-sample distribution cards ────────────────────
// Both cards draw at the card's real pixel width (no viewBox scaling), so their
// text matches the other cards whatever the column width.

function _quantiles(vals) {
  const v = vals.slice().sort((a, b) => a - b);
  return {
    v, n: v.length,
    min: v[0], max: v[v.length - 1],
    q1: d3.quantileSorted(v, 0.25), med: d3.quantileSorted(v, 0.5), q3: d3.quantileSorted(v, 0.75),
    mean: d3.mean(v),
  };
}

/** Pixel-sized SVG that fills its host; redrawn by VQ.redrawOnResize. */
function _pxSvg(host, minW, H) {
  const W = Math.max(host.clientWidth || 0, minW);
  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');
  return { svg, W };
}

function _fmtLen(v) {
  if (v == null || isNaN(v)) return '—';
  return Math.round(v).toLocaleString();
}

function _lenTick(v) {
  return v >= 1e6 ? (v / 1e6) + 'M' : v >= 1000 ? (v / 1000) + 'k' : String(v);
}

// ── Sequence length: one log-scale density ridge per sample ─────────────────

function _renderLengthDensity(samples) {
  const host = _body('ov-lengths');
  if (!host) return;
  const rows = samples.map(s => {
    const L = (s.lengths || []).filter(x => x > 0);
    if (!L.length) return { s, q: null };
    const q = _quantiles(L);
    // N50: length at which half the assembled bases are in contigs this long or longer.
    let acc = 0, n50 = q.max;
    const half = d3.sum(q.v) / 2;
    for (let i = q.v.length - 1; i >= 0; i--) { acc += q.v[i]; if (acc >= half) { n50 = q.v[i]; break; } }
    return { s, q, n50 };
  });
  if (!rows.some(r => r.q)) { host.innerHTML = `<div class="vq-empty">No length data.</div>`; return; }
  const draw = () => _drawLengthDensity(host, rows);
  draw();
  VQ.redrawOnResize(host, draw);
}

function _drawLengthDensity(host, rows) {
  const esc  = VQ.esc;
  const padL = 140, statW = 76, padR = statW * 3 + 8, rowH = 36, headH = 22, axisH = 24;
  const H = headH + rows.length * rowH + axisH;
  const { svg, W } = _pxSvg(host, 620, H);

  const all  = rows.filter(r => r.q);
  const lo   = d3.min(all, r => r.q.min), hi = d3.max(all, r => r.q.max);
  const dLo  = lo / 1.15, dHi = hi * 1.15;        // data range + a little air
  const x    = d3.scaleLog().domain([dLo, dHi]).range([padL, W - padR - 12]);
  const ticks = [];
  for (let p = Math.floor(Math.log10(dLo)); p <= Math.ceil(Math.log10(dHi)); p++)
    [1, 2, 5].forEach(m => { const t = m * Math.pow(10, p); if (t >= dLo && t <= dHi) ticks.push(t); });

  // Grid + axis
  ticks.forEach(t => {
    const major = Math.log10(t) % 1 === 0;
    svg.append('line').attr('class', major ? 'ov-grid' : 'ov-grid ov-grid--minor')
      .attr('x1', x(t)).attr('x2', x(t)).attr('y1', headH).attr('y2', H - axisH + 2);
    svg.append('text').attr('class', 'ov-axis-label')
      .attr('x', x(t)).attr('y', H - axisH + 15).attr('text-anchor', 'middle').text(_lenTick(t));
  });

  // Stat column headers
  const statX = k => W - padR + 8 + statW * (k + 1) - 6;
  ['Median', 'Mean', 'N50'].forEach((h, k) =>
    svg.append('text').attr('class', 'ov-dist__hd').attr('x', statX(k)).attr('y', 13)
      .attr('text-anchor', 'end').text(h));

  const grid = d3.range(0, 121).map(i =>
    Math.log10(dLo) + (Math.log10(dHi) - Math.log10(dLo)) * i / 120);

  rows.forEach((r, j) => {
    const top = headH + j * rowH, base = top + rowH - 5, amp = rowH - 9;
    const g = svg.append('g').attr('class', 'ov-funnel__row');
    g.append('rect').attr('class', 'ov-funnel__hover')
      .attr('x', 0).attr('y', top).attr('width', W).attr('height', rowH);
    g.append('line').attr('class', 'ov-dist__base')
      .attr('x1', padL).attr('x2', W - padR - 12).attr('y1', base).attr('y2', base);
    g.append('text').attr('class', 'ov-funnel__sample')
      .attr('x', 8).attr('y', top + rowH / 2).attr('dominant-baseline', 'central')
      .text(_trunc(r.s.sample, 20)).append('title').text(r.s.sample);

    if (!r.q) {
      g.append('text').attr('class', 'ov-axis-label').attr('x', padL + 6).attr('y', top + rowH / 2)
        .attr('dominant-baseline', 'central').text('no sequences');
      return;
    }
    const q = r.q;
    const logs = q.v.map(Math.log10);
    // Gaussian KDE on log10(length).  Silverman's robust rule (min of SD and
    // IQR/1.34) at 0.9× — the plain rule oversmooths the long right tail of
    // contig lengths into a flat hump; the floor keeps tiny samples readable.
    const sd  = d3.deviation(logs) || 0.1;
    const iqr = (d3.quantile(logs.slice().sort((a, b) => a - b), 0.75)
               - d3.quantile(logs.slice().sort((a, b) => a - b), 0.25)) / 1.34;
    const bw  = Math.max(0.035, 0.9 * Math.min(sd, iqr || sd) * Math.pow(logs.length, -0.2));
    const dens = grid.map(gx => [gx, d3.sum(logs, l => Math.exp(-0.5 * ((gx - l) / bw) ** 2))]);
    const dMax = d3.max(dens, d => d[1]) || 1;
    const pts = dens.map(([gx, d]) => [x(Math.pow(10, gx)), base - (d / dMax) * amp]);

    const clipId = `ov-len-clip-${j}`;
    g.append('clipPath').attr('id', clipId).append('rect')
      .attr('x', padL).attr('y', top).attr('width', W - padR - 12 - padL).attr('height', rowH);
    const area = d3.area().x(d => d[0]).y0(base).y1(d => d[1]).curve(d3.curveBasis);
    const line = d3.line().x(d => d[0]).y(d => d[1]).curve(d3.curveBasis);
    const cg = g.append('g').attr('clip-path', `url(#${clipId})`);
    cg.append('path').attr('class', 'ov-dist__area').attr('d', area(pts));
    cg.append('path').attr('class', 'ov-dist__line').attr('d', line(pts));

    // Rug: every contig as a tick under the curve.
    if (q.n <= 400) q.v.forEach(v => cg.append('line').attr('class', 'ov-dist__rug')
      .attr('x1', x(v)).attr('x2', x(v)).attr('y1', base).attr('y2', base + 3));

    const yAt = v => {                       // curve height at a value (for marker lines)
      const px = x(v);
      const k = d3.bisector(d => d[0]).left(pts, px);
      const p = pts[Math.min(k, pts.length - 1)];
      return p ? p[1] : base - amp;
    };
    g.append('line').attr('class', 'ov-dist__mean')
      .attr('x1', x(q.mean)).attr('x2', x(q.mean)).attr('y1', base).attr('y2', yAt(q.mean));
    g.append('line').attr('class', 'ov-dist__median')
      .attr('x1', x(q.med)).attr('x2', x(q.med)).attr('y1', base).attr('y2', yAt(q.med));

    [_fmtLen(q.med), _fmtLen(q.mean), _fmtLen(r.n50)].forEach((t, k) =>
      g.append('text').attr('class', 'ov-dist__stat' + (k ? '' : ' ov-dist__stat--key'))
        .attr('x', statX(k)).attr('y', top + rowH / 2).attr('text-anchor', 'end')
        .attr('dominant-baseline', 'central').text(t));

    g.on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(r.s.sample)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${q.n.toLocaleString()}</span>
          <span class="vq-tooltip__key">Min</span><span>${_fmtLen(q.min)} nt</span>
          <span class="vq-tooltip__key">Q1</span><span>${_fmtLen(q.q1)} nt</span>
          <span class="vq-tooltip__key">Median</span><span>${_fmtLen(q.med)} nt</span>
          <span class="vq-tooltip__key">Mean</span><span>${_fmtLen(q.mean)} nt</span>
          <span class="vq-tooltip__key">Q3</span><span>${_fmtLen(q.q3)} nt</span>
          <span class="vq-tooltip__key">Max</span><span>${_fmtLen(q.max)} nt</span>
          <span class="vq-tooltip__key">N50</span><span>${_fmtLen(r.n50)} nt</span>
        </div>`, evt))
     .on('mouseleave', VQ.tooltipHide);
  });
}

// ── Score distribution: boxplot + every value per sample ────────────────────

function _renderScoreBoxes(id, samples, key) {
  const host = _body(id);
  if (!host) return;
  const rows = samples.map(s => {
    const vals = (s[key] || []).filter(v => v != null && !isNaN(v));
    if (!vals.length) return { s, q: null };
    const q = _quantiles(vals);
    const iqr = q.q3 - q.q1;
    const loF = q.q1 - 1.5 * iqr, hiF = q.q3 + 1.5 * iqr;
    q.wLo = d3.min(q.v.filter(v => v >= loF));
    q.wHi = d3.max(q.v.filter(v => v <= hiF));
    return { s, q };
  });
  if (!rows.some(r => r.q)) { host.innerHTML = `<div class="vq-empty">No score data.</div>`; return; }
  const draw = () => _drawScoreBoxes(host, rows);
  draw();
  VQ.redrawOnResize(host, draw);
}

function _drawScoreBoxes(host, rows) {
  const esc  = VQ.esc;
  const padL = 130, padR = 64, rowH = 34, headH = 8, axisH = 24, boxH = 16;
  const H = headH + rows.length * rowH + axisH;
  const { svg, W } = _pxSvg(host, 380, H);
  const x = d3.scaleLinear().domain([0, 100]).range([padL, W - padR]);

  // Score bands (low / mid / high) behind the grid, as in the Run tab gauge.
  [[0, 40, 'ov-box__band--lo'], [40, 70, 'ov-box__band--mid'], [70, 100, 'ov-box__band--hi']].forEach(([a, b, c]) =>
    svg.append('rect').attr('class', 'ov-box__band ' + c)
      .attr('x', x(a)).attr('width', x(b) - x(a)).attr('y', headH).attr('height', rows.length * rowH));
  [0, 20, 40, 60, 80, 100].forEach(t => {
    svg.append('line').attr('class', 'ov-grid')
      .attr('x1', x(t)).attr('x2', x(t)).attr('y1', headH).attr('y2', H - axisH + 2);
    svg.append('text').attr('class', 'ov-axis-label').attr('x', x(t)).attr('y', H - axisH + 15)
      .attr('text-anchor', 'middle').text(t);
  });
  svg.append('text').attr('class', 'ov-dist__hd').attr('x', W - 8).attr('y', H - axisH + 15)
    .attr('text-anchor', 'end').text('median');

  rows.forEach((r, j) => {
    const top = headH + j * rowH, cy = top + rowH / 2;
    const g = svg.append('g').attr('class', 'ov-funnel__row');
    g.append('rect').attr('class', 'ov-funnel__hover')
      .attr('x', 0).attr('y', top).attr('width', W).attr('height', rowH);
    g.append('text').attr('class', 'ov-funnel__sample')
      .attr('x', 8).attr('y', cy).attr('dominant-baseline', 'central')
      .text(_trunc(r.s.sample, 18)).append('title').text(r.s.sample);
    if (!r.q) {
      g.append('text').attr('class', 'ov-axis-label').attr('x', padL + 6).attr('y', cy)
        .attr('dominant-baseline', 'central').text('no scores');
      return;
    }
    const q = r.q;

    // Whiskers (1.5 × IQR) with caps
    g.append('line').attr('class', 'ov-box__whisker')
      .attr('x1', x(q.wLo)).attr('x2', x(q.q1)).attr('y1', cy).attr('y2', cy);
    g.append('line').attr('class', 'ov-box__whisker')
      .attr('x1', x(q.q3)).attr('x2', x(q.wHi)).attr('y1', cy).attr('y2', cy);
    [q.wLo, q.wHi].forEach(v => g.append('line').attr('class', 'ov-box__whisker')
      .attr('x1', x(v)).attr('x2', x(v)).attr('y1', cy - 5).attr('y2', cy + 5));

    // Box
    g.append('rect').attr('class', 'ov-box__box')
      .attr('x', x(q.q1)).attr('width', Math.max(2, x(q.q3) - x(q.q1)))
      .attr('y', cy - boxH / 2).attr('height', boxH).attr('rx', 4);

    // Every value, jittered vertically (deterministic, so redraws don't shimmer);
    // dots shrink and fade as n grows so the box stays readable.
    const ptR = q.n > 150 ? 1.6 : q.n > 50 ? 1.9 : 2.2;
    const ptA = q.n > 150 ? 0.28 : q.n > 50 ? 0.38 : 0.5;
    q.v.forEach((v, i) => {
      const jit = ((Math.sin(i * 12.9898 + v * 78.233) * 43758.5453) % 1 + 1) % 1 - 0.5;
      const out = v < q.wLo || v > q.wHi;
      g.append('circle').attr('class', 'ov-box__pt' + (out ? ' ov-box__pt--out' : ''))
        .attr('cx', x(v)).attr('cy', cy + jit * (boxH + 4)).attr('r', out ? 2.6 : ptR)
        .style('fill-opacity', out ? null : ptA);
    });

    // Median bar + mean diamond on top
    g.append('line').attr('class', 'ov-box__median')
      .attr('x1', x(q.med)).attr('x2', x(q.med)).attr('y1', cy - boxH / 2 - 2).attr('y2', cy + boxH / 2 + 2);
    g.append('path').attr('class', 'ov-box__mean')
      .attr('d', d3.symbol().type(d3.symbolDiamond).size(42)())
      .attr('transform', `translate(${x(q.mean)},${cy})`);

    g.append('text').attr('class', 'ov-dist__stat ov-dist__stat--key')
      .attr('x', W - 8).attr('y', cy - 5).attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(q.med.toFixed(1));
    g.append('text').attr('class', 'ov-axis-label')
      .attr('x', W - 8).attr('y', cy + 8).attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(`n=${q.n}`);

    g.on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(r.s.sample)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Scored</span><span>${q.n.toLocaleString()}</span>
          <span class="vq-tooltip__key">Min</span><span>${q.min.toFixed(1)}</span>
          <span class="vq-tooltip__key">Q1</span><span>${q.q1.toFixed(1)}</span>
          <span class="vq-tooltip__key">Median</span><span>${q.med.toFixed(1)}</span>
          <span class="vq-tooltip__key">Mean</span><span>${q.mean.toFixed(1)}</span>
          <span class="vq-tooltip__key">Q3</span><span>${q.q3.toFixed(1)}</span>
          <span class="vq-tooltip__key">Max</span><span>${q.max.toFixed(1)}</span>
        </div>`, evt))
     .on('mouseleave', VQ.tooltipHide);
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
