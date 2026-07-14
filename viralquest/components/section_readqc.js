/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_readqc.js — Read QC (fastp) + sequence-quality summary

   Left:  fastp read-QC KPI chips + per-cycle quality / GC curves.
   Right: aggregated sequence-quality signals (dustmask · self-BLAST ·
          jellyfish) across the confirmed viral sequences.
   ============================================================ */

function _chip(label, value, mod, sub) {
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${VQ.esc(label)}</div>
      <div class="vq-stat-chip__value" title="${VQ.esc(String(value))}">${value}</div>
      ${sub ? `<div class="vq-stat-chip__sub">${VQ.esc(sub)}</div>` : ''}
    </div>`;
}

function _pct(x)  { return x == null ? '—' : (x * 100).toFixed(1) + '%'; }
function _num(x)  { return x == null ? '—' : Number(x).toLocaleString(); }

function vqInitReadQc(readQc, rqStats, sqStats) {
  const el = document.getElementById('section-readqc');
  if (!el || !readQc) return;

  const rq = rqStats || {};
  const rt = (readQc.read_type || 'sr').toUpperCase();

  const chips = [
    _chip('Reads (after)',  _num(rq.reads_after),  'accent',
          `${_num(rq.reads_before)} in · ${_pct(rq.pass_rate)} pass`),
    _chip('Q30',            _pct(rq.q30_rate),     (rq.q30_rate ?? 0) >= 0.8 ? 'success' : 'warn'),
    _chip('Q20',            _pct(rq.q20_rate),     'success'),
    _chip('GC',             rq.gc_pct != null ? rq.gc_pct.toFixed(1) + '%' : '—'),
    _chip('Mean length',    rq.mean_length != null ? Math.round(rq.mean_length) + ' bp' : '—'),
    _chip('Duplication',    _pct(rq.dup_rate),     (rq.dup_rate ?? 0) > 0.2 ? 'warn' : ''),
    _chip('Adapter',        _pct(rq.adapter_rate)),
  ].join('');

  const sq = sqStats && sqStats.present ? sqStats : null;
  const sqChips = sq ? [
    _chip('Sequences',       _num(sq.seqs_analyzed)),
    _chip('With repeats',    _num(sq.with_repeats),        sq.with_repeats ? 'warn' : ''),
    _chip('Low complexity',  _num(sq.with_low_complexity), sq.with_low_complexity ? 'warn' : ''),
    _chip('Total repeats',   _num(sq.total_repeats)),
    _chip(`Repeat score (k${sq.kmer_size})`, _pct(sq.mean_repeat_score), '',
          `max ${_pct(sq.max_repeat_score)}`),
  ].join('') : '';

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Read Quality Control</div>
        <div class="vq-section-sub">
          fastp &nbsp;·&nbsp; ${VQ.esc(rt)} reads &nbsp;·&nbsp;
          ${VQ.esc((readQc.reads || []).map(p => p.split('/').pop()).join(', '))}
        </div>
      </div>
    </div>

    <div style="display:flex;flex-wrap:wrap;gap:10px;margin-bottom:18px">${chips}</div>

    <div class="vq-panel" style="padding:14px 16px;margin-bottom:18px">
      <div class="vq-body-label" style="margin-top:0">Per-cycle read metrics</div>
      <div style="display:flex;gap:24px;flex-wrap:wrap">
        <div style="flex:1 1 320px;min-width:280px">
          <div style="font-size:11px;color:var(--vq-text-3);margin-bottom:4px;font-family:var(--vq-font-mono)">Mean base quality</div>
          <div id="rq-quality-chart"></div>
        </div>
        <div style="flex:1 1 320px;min-width:280px">
          <div style="font-size:11px;color:var(--vq-text-3);margin-bottom:4px;font-family:var(--vq-font-mono)">GC content</div>
          <div id="rq-gc-chart"></div>
        </div>
      </div>
    </div>

    ${sq ? `
    <div class="vq-panel" style="padding:14px 16px">
      <div class="vq-body-label" style="margin-top:0">Sequence quality signals</div>
      <div style="display:flex;flex-wrap:wrap;gap:10px">${sqChips}</div>
      <div class="vq-section-sub" style="margin-top:10px">
        Per-sequence dot plots and low-complexity/repeat tracks are shown in the
        <strong>Sequence Viewer</strong>.
      </div>
    </div>` : ''}
  `;

  _lineChart(el.querySelector('#rq-quality-chart'),
    readQc.per_base_quality || [], { color: 'var(--vq-success)', yMax: 42, yTicks: [0, 20, 30, 40] });
  _lineChart(el.querySelector('#rq-gc-chart'),
    (readQc.per_base_gc || []).map(v => v * 100), { color: 'var(--vq-accent)', yMax: 100, yTicks: [0, 25, 50, 75, 100] });
}

// Minimal responsive line chart over per-cycle values.
function _lineChart(container, values, opts) {
  if (!container) return;
  if (!values.length) {
    container.innerHTML = '<div class="vq-empty" style="padding:16px 0;font-size:11px">No per-cycle data.</div>';
    return;
  }
  const o = opts || {};
  const W = 360, H = 150, PAD_L = 32, PAD_B = 22, PAD_T = 8, PAD_R = 8;
  const n = values.length;
  const yMax = o.yMax || Math.max(...values) * 1.1 || 1;

  const x = d3.scaleLinear([0, n - 1], [PAD_L, W - PAD_R]);
  const y = d3.scaleLinear([0, yMax], [H - PAD_B, PAD_T]);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`).attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // Gridlines + y ticks.
  (o.yTicks || y.ticks(4)).forEach(t => {
    svg.append('line')
      .attr('x1', PAD_L).attr('x2', W - PAD_R).attr('y1', y(t)).attr('y2', y(t))
      .attr('stroke', 'var(--vq-border)').attr('stroke-dasharray', '3,3');
    svg.append('text')
      .attr('x', PAD_L - 6).attr('y', y(t) + 3).attr('text-anchor', 'end')
      .attr('font-size', 8).attr('fill', 'var(--vq-text-3)').text(t);
  });

  const area = d3.area().x((d, i) => x(i)).y0(H - PAD_B).y1(d => y(d)).curve(d3.curveMonotoneX);
  svg.append('path').datum(values).attr('d', area).attr('fill', o.color).attr('opacity', 0.14);

  const line = d3.line().x((d, i) => x(i)).y(d => y(d)).curve(d3.curveMonotoneX);
  svg.append('path').datum(values).attr('d', line)
    .attr('fill', 'none').attr('stroke', o.color).attr('stroke-width', 1.6);

  // X axis baseline + end labels (cycle index).
  svg.append('line')
    .attr('x1', PAD_L).attr('x2', W - PAD_R).attr('y1', H - PAD_B).attr('y2', H - PAD_B)
    .attr('stroke', 'var(--vq-border-dark)');
  svg.append('text').attr('x', PAD_L).attr('y', H - 6)
    .attr('font-size', 8).attr('fill', 'var(--vq-text-3)').text('1');
  svg.append('text').attr('x', W - PAD_R).attr('y', H - 6).attr('text-anchor', 'end')
    .attr('font-size', 8).attr('fill', 'var(--vq-text-3)').text('cycle');

  container.appendChild(svg.node());
}

window.vqInitReadQc = vqInitReadQc;

})();
