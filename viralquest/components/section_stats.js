/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_stats.js — Run tab: how the pipeline ran
   Workflow card, run KPIs, then grouped cards
   (Detection & scoring · Data & assembly).
   The Virome tab (results) lives in section_viruses.js.
   ============================================================ */

function vqInitStats(report) {
  _initRunTab(report);
}

// Overview helpers shared with the other overview tab (defined in export.js).
const _chip           = VQ.statChip;
const _fmtNum         = VQ.fmtNum;
const _group          = VQ.cardGroup;
const _pct            = VQ.pct;
const _redrawOnResize = VQ.redrawOnResize;

// ════════════════════════════════════════════════════════════════════════
//  RUN — how the pipeline ran: workflow, filtering, databases, scoring, QC
// ════════════════════════════════════════════════════════════════════════

function _initRunTab(report) {
  const el = document.getElementById('section-run');
  if (!el) return;

  const esc     = VQ.esc;
  const meta    = report.meta             || {};
  const sum     = report.summary          || {};
  const blast   = report.blast_stats      || {};
  const hmm     = report.hmm_stats        || {};
  const llm     = report.llm_stats        || {};
  const heur    = report.heuristic_stats  || {};
  const salmon  = report.salmon_stats     || {};
  const seqQual = report.seq_quality_stats || {};
  const seqs    = report.sequences        || [];
  const wf      = report.workflow || (report.pipeline_stats || {}).workflow || null;
  const hasWf   = !!(wf && (wf.steps || []).length);

  const totalSeqs      = sum.total_sequences ?? seqs.length;
  const confirmedViral = sum.confirmed_viral ?? seqs.filter(s => s.is_viral).length;
  const totalOrfs      = sum.total_orfs      ?? seqs.reduce((a, s) => a + (s.orfs || []).length, 0);
  const totalClusters  = sum.total_clusters  ?? (report.clusters || []).length;
  const hmmTotal = (hmm.rvdb_hits || 0) + (hmm.vfam_hits || 0)
                 + (hmm.eggnog_hits ?? hmm.eggnong_hits ?? 0) + (hmm.pfam_hits || 0);

  const detection = `
    <div class="vq-chart-card" id="stats-funnel-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Detection Pipeline</div>
          <div class="vq-chart-card__sub">sequences at each filtering step</div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="padding:6px 14px 10px">
        <div id="stats-funnel-svg" style="width:100%"></div>
      </div>
    </div>
    <div class="vq-chart-card" id="stats-hmm-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">HMM Database Hits</div>
          <div class="vq-chart-card__sub">threshold-passing hits per database</div>
        </div>
        <div class="vq-chart-card__big" id="stats-hmm-total">${_fmtNum(hmmTotal)}</div>
      </div>
      <div class="vq-chart-card__body">
        <div id="stats-hmm-svg" style="width:100%"></div>
      </div>
    </div>`;

  const scoring = `
    ${heur.present ? `
    <div class="vq-chart-card" id="stats-heur-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Heuristic Scoring</div>
          <div class="vq-chart-card__sub">rule-based · deterministic</div>
        </div>
        <div class="vq-chart-card__big">${heur.avg_score != null ? Number(heur.avg_score).toFixed(1) : '—'}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start;gap:0">
        <div id="stats-heur-gauge" style="width:100%;margin-bottom:10px"></div>
        ${_miniRow('Scored',        _fmtNum(heur.scored))}
        ${_miniRow('Viral known',   _fmtNum(heur.viral_known))}
        ${_miniRow('Viral unknown', _fmtNum(heur.viral_unknown))}
        ${heur.non_viral ? _miniRow('Non-viral', _fmtNum(heur.non_viral)) : ''}
      </div>
    </div>` : ''}
    ${llm.present ? `
    <div class="vq-chart-card" id="stats-llm-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">LLM Scoring</div>
          <div class="vq-chart-card__sub">${esc(llm.model || '—')}${llm.mode ? ' · ' + esc(llm.mode) + '-token mode' : ''}</div>
        </div>
        <div class="vq-chart-card__big">${llm.avg_score != null ? Number(llm.avg_score).toFixed(1) : '—'}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('Sent to LLM',   _fmtNum(llm.scored))}
        ${_miniRow('Viral known',   _fmtNum(llm.viral_known))}
        ${_miniRow('Viral unknown', _fmtNum(llm.viral_unknown))}
        ${llm.api_error   ? _miniRow('API errors',        _fmtNum(llm.api_error))   : ''}
        ${llm.parse_error ? _miniRow('Invalid responses', _fmtNum(llm.parse_error)) : ''}
      </div>
    </div>` : ''}`;

  const dataQc = `
    ${seqQual.present ? `
    <div class="vq-chart-card" id="stats-seqqual-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Sequence Quality</div>
          <div class="vq-chart-card__sub">structural signals across sequences</div>
        </div>
        <div class="vq-chart-card__big">${_fmtNum(seqQual.seqs_analyzed)}</div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('With repeats',        _fmtNum(seqQual.with_repeats))}
        ${_miniRow('With low complexity', _fmtNum(seqQual.with_low_complexity))}
        ${_miniRow('Total repeats',       _fmtNum(seqQual.total_repeats))}
        ${_miniRow(`Repeat score (k${seqQual.kmer_size ?? '?'})`,
                   seqQual.mean_repeat_score != null
                     ? (seqQual.mean_repeat_score * 100).toFixed(1) + '%'
                     + (seqQual.max_repeat_score != null
                         ? ' · max ' + (seqQual.max_repeat_score * 100).toFixed(1) + '%' : '')
                     : '—')}
      </div>
    </div>` : ''}
    ${sum.cap3_used ? `
    <div class="vq-chart-card" id="stats-cap3-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">CAP3 Assembly</div>
          <div class="vq-chart-card__sub">input re-assembled before analysis</div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
        ${_miniRow('Contigs',  _fmtNum(sum.cap3_contigs))}
        ${_miniRow('Singlets', _fmtNum(sum.cap3_singlets))}
      </div>
    </div>` : ''}
    ${salmon.present ? `
    <div class="vq-chart-card" id="stats-salmon-card">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Read Mapping · Salmon</div>
          <div class="vq-chart-card__sub">
            ${_fmtNum(salmon.total_reads)} reads
            ${salmon.pathway ? ' · ' + esc(salmon.pathway) + ' pathway' : ''}
          </div>
        </div>
      </div>
      <div class="vq-chart-card__body" style="flex-direction:column;justify-content:flex-start;gap:0">
        <div id="stats-salmon-svg" style="width:100%"></div>
        <div id="stats-salmon-rows" style="width:100%;padding:0 4px 4px"></div>
      </div>
    </div>` : ''}`;

  el.innerHTML = `
    ${hasWf ? '' : `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Run</div>
        <div class="vq-section-sub">
          ${esc(meta.input_file?.name || 'Sample report')}
          &nbsp;·&nbsp; run ${esc(meta.timestamp ? new Date(meta.timestamp).toLocaleString() : '—')}
          &nbsp;·&nbsp; ViralQuest v${esc(meta.viralquest_version || '?')}
        </div>
      </div>
    </div>`}

    <div class="vq-stats-page">
      ${hasWf ? _workflowCardHtml(wf, meta.viralquest_version) : ''}

      <div class="vq-stats-row vq-stats-row--top">
        ${_chip('Input Sequences', _fmtNum(totalSeqs), 'accent',
                meta.input_file?.name ? esc(meta.input_file.name) : '')}
        ${_chip('ORFs Predicted', _fmtNum(totalOrfs), '', 'across confirmed sequences')}
        ${_chip('Confirmed Viral', _fmtNum(confirmedViral), 'success',
                _pct(confirmedViral, totalSeqs) + ' of input')}
        ${hasWf
          ? _chip('Runtime', _fmtDur(wf.total_seconds), '',
                  `${(wf.steps || []).filter(x => x.status !== 'skipped').length} steps run`)
          : _chip('Clusters', _fmtNum(totalClusters), 'accent', 'species clusters')}
      </div>

      ${_group('Detection & scoring', 'how sequences were filtered, confirmed and scored',
               detection + scoring, 2 + !!heur.present + !!llm.present)}
      ${_group('Data & assembly', 'input processing and read-level quality',
               dataQc, !!seqQual.present + !!sum.cap3_used + !!salmon.present)}
    </div>`;

  if (hasWf) _renderWorkflow(wf);
  const draw = () => {
    _renderFunnel(blast);
    _renderHMMBars(hmm);
    if (heur.present)   _renderHeuristicGauge(heur);
    if (salmon.present) _renderSalmonChart(salmon);
  };
  draw();
  _redrawOnResize(el, () => { draw(); if (hasWf) _renderWorkflow(wf); });
}

// ────────────────────────────────────────────────────────────────────────
//  Mini-row helper (label · value line inside a card)
// ────────────────────────────────────────────────────────────────────────

function _miniRow(label, value) {
  return `
    <div style="display:flex;justify-content:space-between;align-items:baseline;
                padding:6px 0;border-bottom:1px dashed var(--vq-border);font-size:12px;
                gap:12px;min-width:0">
      <span style="color:var(--vq-text-3);flex-shrink:0">${VQ.esc(label)}</span>
      <span style="color:var(--vq-text);font-weight:600;font-variant-numeric:tabular-nums;
                   white-space:nowrap;overflow:hidden;text-overflow:ellipsis">${value}</span>
    </div>`;
}

// ────────────────────────────────────────────────────────────────────────
//  Pipeline workflow card — flowchart of every step (run / skipped),
//  success path between completed steps, per-step time and details.
// ────────────────────────────────────────────────────────────────────────

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

function _fmtDur(sec) {
  if (sec == null || isNaN(sec)) return '—';
  if (sec < 1)    return (sec * 1000).toFixed(0) + ' ms';
  if (sec < 60)   return sec.toFixed(1) + ' s';
  const r = Math.round(sec), m = Math.floor(r / 60), s = r % 60;
  if (r < 3600) return `${m}m ${String(s).padStart(2, '0')}s`;
  return `${Math.floor(m / 60)}h ${String(m % 60).padStart(2, '0')}m`;
}

function _wfOverall(steps) {
  if (steps.some(s => s.status === 'error'))   return 'error';
  if (steps.some(s => s.status === 'partial')) return 'partial';
  return 'done';
}

function _wfPill(status, text) {
  return `<span class="vq-wf-pill vq-wf-pill--${status}">${VQ.esc(text ?? _WF_STATUS[status].label)}</span>`;
}

function _workflowCardHtml(wf, version) {
  const esc     = VQ.esc;
  const steps   = wf.steps || [];
  const ran     = steps.filter(s => s.status !== 'skipped');
  const skipped = steps.length - ran.length;
  const overall = _wfOverall(steps);
  const issues  = steps.filter(s => s.status === 'error' || s.status === 'partial');
  const o       = wf.options || {};

  const opt = (label, value, on = true) =>
    `<span class="vq-wf-opt${on ? '' : ' vq-wf-opt--off'}">
       <span class="vq-wf-opt__k">${esc(label)}</span>${esc(value)}</span>`;
  const reads = (o.reads || []).length
    ? `${o.reads.length} file${o.reads.length > 1 ? 's' : ''}${o.read_type ? ' · ' + o.read_type : ''}`
    : 'off';

  const overallText = overall === 'done' ? 'Completed'
                    : overall === 'partial' ? 'Completed with warnings'
                    : 'Completed with errors';

  return `
    <div class="vq-chart-card vq-wf-card" id="stats-workflow-card" style="min-height:auto">
      <div class="vq-chart-card__head">
        <div>
          <div class="vq-chart-card__title">Pipeline Workflow</div>
          <div class="vq-chart-card__sub">
            started ${esc(wf.started_at ? new Date(wf.started_at).toLocaleString() : '—')}
            &nbsp;·&nbsp; ${ran.length} step${ran.length === 1 ? '' : 's'} run
            ${skipped ? `&nbsp;·&nbsp; ${skipped} skipped` : ''}
            ${version ? `&nbsp;·&nbsp; ViralQuest v${esc(version)}` : ''}
            &nbsp;·&nbsp; click a step for details
          </div>
        </div>
        <div class="vq-wf-head-right">
          ${_wfPill(overall, overallText)}
        </div>
      </div>

      <div class="vq-wf-opts">
        ${o.input ? opt('Input', o.input) : ''}
        ${o.threads != null ? opt('Threads', o.threads) : ''}
        ${opt('CAP3', o.cap3 ? 'on' : 'off', !!o.cap3)}
        ${opt('NR', o.nr_db || 'off', !!o.nr_db)}
        ${opt('BLASTn', o.blastn || 'off', !!o.blastn)}
        ${opt('Reads', reads, !!(o.reads || []).length)}
        ${o.transcriptome ? opt('Transcriptome', o.transcriptome) : ''}
        ${opt('LLM', o.llm || 'off', !!o.llm)}
        ${o.force ? opt('Force', 'on') : ''}
      </div>

      ${issues.length ? `
      <div class="vq-wf-issues">
        ${issues.map(s => `
          <button type="button" class="vq-wf-issue vq-wf-issue--${s.status}" data-wf-key="${esc(s.key)}">
            <strong>${esc(_WF_SHORT[s.key] || s.label)}</strong>
            <span>${esc(s.message || _WF_STATUS[s.status].label)}</span>
          </button>`).join('')}
      </div>` : ''}

      <div class="vq-wf-flow" id="stats-workflow-flow"></div>
      <div class="vq-wf-timebar" id="stats-workflow-timebar"></div>
      <div class="vq-wf-detail" id="stats-workflow-detail"></div>
    </div>`;
}

function _renderWorkflow(wf) {
  const wrap = document.getElementById('stats-workflow-flow');
  if (!wrap) return;
  const esc   = VQ.esc;
  const steps = wf.steps || [];
  const total = wf.total_seconds || steps.reduce((a, s) => a + (s.seconds || 0), 0);

  // ── Flowchart (SVG) ──────────────────────────────────────────────────
  const SLOT   = Math.max((wrap.clientWidth || 0) / steps.length, 86);
  const W      = SLOT * steps.length;
  const Y_MAIN = 26, Y_SKIP = 100, H = 146;
  const xOf    = i => SLOT * i + SLOT / 2;

  const svg = d3.create('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');

  // Success path: connects consecutive executed steps along the main lane;
  // each segment takes the colour of the step it leads into.
  const ran = steps.map((s, i) => ({ s, i })).filter(d => d.s.status !== 'skipped');
  ran.forEach((d, k) => {
    if (!k) return;
    const prev = ran[k - 1];
    svg.append('line')
      .attr('x1', xOf(prev.i)).attr('y1', Y_MAIN)
      .attr('x2', xOf(d.i)).attr('y2', Y_MAIN)
      .attr('stroke', _WF_STATUS[d.s.status].color)
      .attr('stroke-width', 3).attr('stroke-linecap', 'round')
      .attr('opacity', d.s.status === 'done' ? 0.55 : 0.85);
  });

  // Skipped steps hang below the path on a dashed branch (option not chosen).
  steps.forEach((s, i) => {
    if (s.status !== 'skipped') return;
    svg.append('path')
      .attr('d', `M${xOf(i)},${Y_MAIN + 4} C${xOf(i)},${Y_MAIN + 40} ${xOf(i)},${Y_SKIP - 40} ${xOf(i)},${Y_SKIP - 9}`)
      .attr('fill', 'none')
      .attr('stroke', 'var(--vq-border-dark)')
      .attr('stroke-width', 1.2).attr('stroke-dasharray', '3,3');
  });

  const nodes = svg.selectAll('g.vq-wf-node').data(steps).join('g')
    .attr('class', s => `vq-wf-node vq-wf-node--${s.status}`)
    .attr('data-wf-key', s => s.key)
    .attr('transform', (s, i) => `translate(${xOf(i)},${s.status === 'skipped' ? Y_SKIP : Y_MAIN})`)
    .style('cursor', 'pointer');

  nodes.append('circle').attr('class', 'vq-wf-node__halo')
    .attr('r', s => s.status === 'skipped' ? 13 : 17)
    .attr('fill', 'none').attr('stroke', 'var(--vq-accent)').attr('stroke-width', 2)
    .attr('opacity', 0);

  nodes.append('circle')
    .attr('r', s => s.status === 'skipped' ? 8 : 11)
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-surface)' : _WF_STATUS[s.status].color)
    .attr('stroke', s => s.status === 'skipped' ? 'var(--vq-border-dark)' : 'var(--vq-surface)')
    .attr('stroke-width', s => s.status === 'skipped' ? 1.4 : 2)
    .attr('stroke-dasharray', s => s.status === 'skipped' ? '2,2' : null);

  nodes.filter(s => s.status !== 'skipped').append('text')
    .attr('text-anchor', 'middle').attr('dy', '0.35em')
    .attr('font-size', 11).attr('font-weight', 700).attr('fill', '#fff')
    .attr('pointer-events', 'none')
    .text(s => _WF_STATUS[s.status].glyph);

  nodes.append('text')
    .attr('text-anchor', 'middle')
    .attr('y', s => s.status === 'skipped' ? 22 : 30)
    .attr('font-size', 11)
    .attr('font-weight', s => s.status === 'skipped' ? 400 : 600)
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-text-3)' : 'var(--vq-text)')
    .text(s => _WF_SHORT[s.key] || s.key);

  nodes.append('text')
    .attr('text-anchor', 'middle')
    .attr('y', s => s.status === 'skipped' ? 34 : 44)
    .attr('font-size', 10)
    .attr('font-family', 'var(--vq-font-mono)')
    .attr('fill', s => s.status === 'skipped' ? 'var(--vq-text-3)' : _WF_STATUS[s.status].color)
    .text(s => s.status === 'skipped' ? 'skipped' : _fmtDur(s.seconds));

  nodes
    .on('mousemove', (evt, s) => VQ.tooltipShow(`
      <div class="vq-tooltip__title">${esc(s.label)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Status</span><span>${_WF_STATUS[s.status].label}</span>
        ${s.seconds != null ? `<span class="vq-tooltip__key">Time</span><span>${_fmtDur(s.seconds)}</span>` : ''}
        ${s.message ? `<span class="vq-tooltip__key">Note</span><span>${esc(s.message)}</span>` : ''}
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide)
    .on('click', (evt, s) => select(s.key));

  wrap.innerHTML = '';
  wrap.appendChild(svg.node());

  // ── Time breakdown bar ───────────────────────────────────────────────
  const bar = document.getElementById('stats-workflow-timebar');
  const timed = steps.filter(s => s.seconds != null && s.seconds > 0);
  if (bar && timed.length && total > 0) {
    const longest = timed.reduce((a, s) => (s.seconds > a.seconds ? s : a));
    bar.innerHTML = `
      <div class="vq-wf-timebar__track">
        ${timed.map((s, i) => `
          <div class="vq-wf-timebar__seg" data-wf-key="${esc(s.key)}"
               style="flex:${s.seconds} 1 0;opacity:${i % 2 ? 0.55 : 0.9}"></div>`).join('')}
      </div>
      <div class="vq-wf-timebar__legend">
        <span>time per step</span>
        <span>longest: <strong>${esc(_WF_SHORT[longest.key] || longest.label)}</strong>
          · ${_fmtDur(longest.seconds)} (${_pct(longest.seconds, total)})</span>
      </div>`;
    bar.querySelectorAll('.vq-wf-timebar__seg').forEach(seg => {
      const s = steps.find(x => x.key === seg.dataset.wfKey);
      seg.addEventListener('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(s.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Time</span><span>${_fmtDur(s.seconds)}</span>
          <span class="vq-tooltip__key">of total</span><span>${_pct(s.seconds, total)}</span>
        </div>`, evt));
      seg.addEventListener('mouseleave', VQ.tooltipHide);
      seg.addEventListener('click', () => select(s.key));
    });
  }

  // ── Step details ─────────────────────────────────────────────────────
  const detail = document.getElementById('stats-workflow-detail');
  function select(key) {
    const s = steps.find(x => x.key === key);
    if (!s || !detail) return;
    svg.selectAll('g.vq-wf-node').each(function (d) {
      d3.select(this).select('.vq-wf-node__halo').attr('opacity', d.key === key ? 0.9 : 0);
    });
    const rows = Object.entries(s.details || {});
    detail.innerHTML = `
      <div class="vq-wf-detail__head">
        <div class="vq-wf-detail__title">${esc(s.label)}</div>
        ${_wfPill(s.status)}
        ${s.seconds != null
          ? `<span class="vq-wf-detail__time">${_fmtDur(s.seconds)} · ${_pct(s.seconds, total)} of total</span>`
          : ''}
      </div>
      ${s.message ? `<div class="vq-wf-detail__msg vq-wf-detail__msg--${s.status}">${esc(s.message)}</div>` : ''}
      ${rows.length ? `
        <div class="vq-wf-detail__grid">
          ${rows.map(([k, v]) => _miniRow(k,
              esc(typeof v === 'number' ? v.toLocaleString() : String(v ?? '—')))).join('')}
        </div>` : ''}`;
  }

  document.querySelectorAll('#stats-workflow-card .vq-wf-issue').forEach(btn =>
    btn.addEventListener('click', () => select(btn.dataset.wfKey)));

  // Open on the first problem, otherwise on the slowest step.
  const first = steps.find(s => s.status === 'error')
             || steps.find(s => s.status === 'partial')
             || (timed.length ? timed.reduce((a, s) => (s.seconds > a.seconds ? s : a)) : steps[0]);
  select(first.key);
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 2 — Horizontal bar: HMM database hits
// ────────────────────────────────────────────────────────────────────────

function _renderHMMBars(hmm) {
  const wrap = document.getElementById('stats-hmm-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const eggnog = hmm.eggnog_hits ?? hmm.eggnong_hits ?? 0;
  const data = [
    { name: 'RVDB',   value: hmm.rvdb_hits || 0,   color: 'var(--vq-accent)'        },
    { name: 'Vfam',   value: hmm.vfam_hits || 0,   color: 'var(--vq-primary-light)' },
    { name: 'EggNOG', value: eggnog,               color: 'var(--vq-accent-dark)'   },
    { name: 'Pfam',   value: hmm.pfam_hits || 0,   color: 'var(--vq-success)'       },
  ];

  const W = Math.max(wrap.clientWidth || 0, 220);
  const PAD_L = 64;
  const PAD_R = 36;
  const ROW_H = 38;
  const BAR_H = 16;
  const H = data.length * ROW_H + 16;
  const drawW = W - PAD_L - PAD_R;

  const maxV = d3.max(data, d => d.value) || 1;
  const xScale = d3.scaleLinear([0, maxV], [0, drawW]).nice();

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  data.forEach((d, i) => {
    const y = 8 + i * ROW_H;

    svg.append('text')
      .attr('x', PAD_L - 8).attr('y', y + BAR_H / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 12)
      .attr('font-weight', 500)
      .attr('fill', 'var(--vq-text-2)')
      .text(d.name);

    // Track
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', drawW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', 'var(--vq-bg)');

    // Value bar
    const bw = Math.max(xScale(d.value), 4);
    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', bw).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', d.color)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(d.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Hits</span><span>${d.value.toLocaleString()}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    // Value text (after the bar)
    svg.append('text')
      .attr('x', PAD_L + bw + 8).attr('y', y + BAR_H / 2 + 4)
      .attr('font-size', 12)
      .attr('font-weight', 600)
      .attr('fill', 'var(--vq-text)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(d.value.toLocaleString());
  });

  wrap.appendChild(svg.node());
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 3 — Salmon: mapping-rate gauge + reads breakdown
// ────────────────────────────────────────────────────────────────────────

function _renderSalmonChart(salmon) {
  const wrap = document.getElementById('stats-salmon-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const W    = Math.max(wrap.clientWidth || 0, 200);
  const H    = 72;
  const rate = Math.max(0, Math.min(100, salmon.mapping_rate ?? 0));

  const rateColor = rate >= 30 ? 'var(--vq-success)'
                  : rate >= 15 ? 'var(--vq-accent)'
                              : 'var(--vq-warning)';

  const PAD_L = 4;
  const PAD_R = 4;
  const BAR_Y = 46;
  const BAR_H = 12;
  const drawW = W - PAD_L - PAD_R;
  const fillW = Math.max(drawW * rate / 100, 4);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMidYMin meet')
    .style('width', '100%').style('height', 'auto');

  // Large rate % + label
  svg.append('text')
    .attr('x', PAD_L).attr('y', 22)
    .attr('font-size', 26).attr('font-weight', 700)
    .attr('letter-spacing', '-0.02em')
    .attr('fill', 'var(--vq-primary)')
    .text(rate.toFixed(1) + '%');
  svg.append('text')
    .attr('x', PAD_L).attr('y', 33)
    .attr('font-size', 11).attr('fill', 'var(--vq-text-3)')
    .text('mapping rate');

  // Track
  svg.append('rect')
    .attr('x', PAD_L).attr('y', BAR_Y)
    .attr('width', drawW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', 'var(--vq-bg)');

  // Fill — grows left to right
  svg.append('rect')
    .attr('x', PAD_L).attr('y', BAR_Y)
    .attr('width', fillW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', rateColor).attr('opacity', 0.9)
    .on('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">Mapping rate</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Reads aligned</span><span>${rate.toFixed(1)}%</span>
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // End labels
  svg.append('text')
    .attr('x', PAD_L).attr('y', H - 3)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)').text('0');
  svg.append('text')
    .attr('x', W - PAD_R).attr('y', H - 3)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--vq-text-3)').text('100%');

  wrap.appendChild(svg.node());

  // Mini-rows below the gauge
  const rows = document.getElementById('stats-salmon-rows');
  if (rows) {
    const _fmt = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();
    const mappedReads = salmon.total_reads != null && salmon.mapping_rate != null
      ? Math.round(salmon.total_reads * salmon.mapping_rate / 100)
      : null;
    rows.innerHTML = [
      mappedReads != null
        ? _miniRow('Mapped reads',         _fmt(mappedReads))
        : '',
      salmon.viral_expressed != null
        ? _miniRow('Viral seqs expressed', _fmt(salmon.viral_expressed))
        : '',
      salmon.ref_hk_count
        ? _miniRow('HK ref genes',         _fmt(salmon.ref_hk_count))
        : '',
      salmon.conserved_count
        ? _miniRow('Conserved HK detected',
                   `${_fmt(salmon.conserved_detected)} / ${_fmt(salmon.conserved_count)}`)
        : '',
      salmon.pfam_hk_count
        ? _miniRow('Pfam HK seqs',         _fmt(salmon.pfam_hk_count))
        : '',
    ].join('');
  }
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 4 — Pipeline funnel
// ────────────────────────────────────────────────────────────────────────

function _renderFunnel(blast) {
  const wrap = document.getElementById('stats-funnel-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const steps = [
    { label: 'Input sequences',    value: blast.total_input          ?? null, color: 'var(--vq-primary)'       },
    { label: 'RefSeq BLASTx hits', value: blast.refseq_unique_seqs   ?? 0,   color: 'var(--vq-accent)'        },
    { label: 'Viral flagged',      value: blast.total_viral_flagged  ?? null, color: 'var(--vq-primary-light)' },
    { label: 'Confirmed viral',    value: blast.total_confirmed       ?? 0,   color: 'var(--vq-success)'       },
    { label: 'BLASTn annotated',   value: blast.blastn_unique_seqs   ?? null, color: 'var(--vq-accent-dark)'   },
  ].filter(s => s.value !== null);

  if (!steps.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:20px;font-size:12px">No pipeline data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 260);
  const ROW_H = 27;
  const BAR_H = 10;
  const H     = steps.length * ROW_H + 8;
  const base  = steps[0].value || 1;

  // Column positions: label | number | bar | pct
  const DOT_X   = 8;
  const LBL_X   = 20;
  const NUM_END = Math.round(W * 0.60);
  const BAR_X   = Math.round(W * 0.62);
  const BAR_W   = Math.round(W * 0.25);
  const PCT_X   = BAR_X + BAR_W + 5;

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  steps.forEach((step, i) => {
    const mid    = 4 + i * ROW_H + ROW_H / 2;
    const pctNum = i === 0 ? 100 : (base > 0 ? step.value / base * 100 : 0);
    const pctStr = i === 0 ? '100%' : pctNum.toFixed(1) + '%';
    const fillW  = Math.max(pctNum / 100 * BAR_W, 2);

    // Dot indicator
    svg.append('circle').attr('cx', DOT_X).attr('cy', mid).attr('r', 4)
      .attr('fill', step.color);

    // Connector to next step
    if (i < steps.length - 1)
      svg.append('line')
        .attr('x1', DOT_X).attr('y1', mid + 5)
        .attr('x2', DOT_X).attr('y2', mid + ROW_H - 3)
        .attr('stroke', 'var(--vq-border-dark)')
        .attr('stroke-width', 1.2).attr('stroke-dasharray', '2,2');

    // Step label (left)
    svg.append('text')
      .attr('x', LBL_X).attr('y', mid + 4)
      .attr('font-size', 11).attr('fill', 'var(--vq-text-2)')
      .text(step.label);

    // Count — bold, colour-matched (right-aligned at NUM_END)
    svg.append('text')
      .attr('x', NUM_END).attr('y', mid + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 12).attr('font-weight', 700)
      .attr('font-family', 'var(--vq-font-mono)')
      .attr('font-variant-numeric', 'tabular-nums')
      .attr('fill', step.color)
      .text(step.value.toLocaleString());

    // Bar track
    svg.append('rect')
      .attr('x', BAR_X).attr('y', mid - BAR_H / 2)
      .attr('width', BAR_W).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', 'var(--vq-bg)');

    // Bar fill — proportional to % of input
    svg.append('rect')
      .attr('x', BAR_X).attr('y', mid - BAR_H / 2)
      .attr('width', fillW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2)
      .attr('fill', step.color).attr('opacity', 0.8)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(step.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${step.value.toLocaleString()}</span>
          <span class="vq-tooltip__key">of input</span><span>${pctStr}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);

    // Pct label after bar
    svg.append('text')
      .attr('x', PCT_X).attr('y', mid + 4)
      .attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
      .attr('font-variant-numeric', 'tabular-nums')
      .text(pctStr);
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Heuristic score gauge — stacked classification bar + avg score track
// ────────────────────────────────────────────────────────────────────────

function _renderHeuristicGauge(heur) {
  const wrap = document.getElementById('stats-heur-gauge');
  if (!wrap) return;
  wrap.innerHTML = '';

  const total   = heur.scored || 1;
  const known   = heur.viral_known   || 0;
  const unknown = heur.viral_unknown || 0;
  const nonvir  = heur.non_viral     || 0;
  const avg     = heur.avg_score     ?? 0;

  const W    = Math.max(wrap.clientWidth || 0, 180);
  const H    = 62;
  const BAR_Y = 6;
  const BAR_H = 10;

  // Colour-code the avg-score track
  const scoreColor = avg >= 70 ? 'var(--vq-success)'
                   : avg >= 40 ? 'var(--vq-accent)'
                               : 'var(--vq-warning)';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // ── Top: stacked classification bar ──────────────────────────────────
  const segments = [
    { label: 'Viral known',   n: known,   color: 'var(--vq-success)'      },
    { label: 'Viral unknown', n: unknown, color: 'var(--vq-accent)'       },
    { label: 'Non-viral',     n: nonvir,  color: 'var(--vq-warning)'      },
  ].filter(s => s.n > 0);

  // Clip segments to a rounded rect so the stacked fill gets rounded caps
  // without square corners poking out behind an overlay stroke.
  const clipId = `vq-heur-bar-clip-${Math.random().toString(36).slice(2)}`;
  svg.append('defs').append('clipPath')
    .attr('id', clipId)
    .append('rect')
      .attr('x', 0).attr('y', BAR_Y)
      .attr('width', W).attr('height', BAR_H)
      .attr('rx', BAR_H / 2);

  const barG = svg.append('g').attr('clip-path', `url(#${clipId})`);

  let cx = 0;
  segments.forEach(seg => {
    const segW = Math.max((seg.n / total) * W, 2);
    barG.append('rect')
      .attr('x', cx).attr('y', BAR_Y)
      .attr('width', segW).attr('height', BAR_H)
      .attr('fill', seg.color).attr('opacity', 0.85)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(seg.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${seg.n}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(seg.n / total * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
    cx += segW;
  });

  // ── Bottom: avg score track ───────────────────────────────────────────
  const TRACK_Y = BAR_Y + BAR_H + 10;
  const fillW   = Math.max((avg / 100) * W, 3);

  svg.append('text')
    .attr('x', 0).attr('y', TRACK_Y - 2)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)')
    .text('avg score');

  // Track background
  svg.append('rect')
    .attr('x', 0).attr('y', TRACK_Y + 4)
    .attr('width', W).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', 'var(--vq-bg)');

  // Fill
  svg.append('rect')
    .attr('x', 0).attr('y', TRACK_Y + 4)
    .attr('width', fillW).attr('height', BAR_H)
    .attr('rx', BAR_H / 2)
    .attr('fill', scoreColor).attr('opacity', 0.9)
    .on('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">Average heuristic score</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Score</span><span>${avg.toFixed(1)}</span>
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // End labels
  svg.append('text')
    .attr('x', 0).attr('y', TRACK_Y + BAR_H + 14)
    .attr('font-size', 9).attr('fill', 'var(--vq-text-3)').text('0');
  svg.append('text')
    .attr('x', W).attr('y', TRACK_Y + BAR_H + 14)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--vq-text-3)').text('100');

  wrap.appendChild(svg.node());
}


window.vqInitStats = vqInitStats;
})();
