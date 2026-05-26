/* ============================================================
   section_salmon.js — Section 5: Salmon Quantification
   Left panel: viral sequences — boxplot if clustered, dot if singleton.
   Right panel: housekeeping gene mini-panel per kingdom.
   Toggle: TPM ↔ NumReads.
   ============================================================ */

'use strict';

// Module state
const _SQ = {
  data:    null,   // VQ_REPORT.salmon_quant
  clusters: null,  // VQ_REPORT.clusters
  metric:  'tpm',  // 'tpm' | 'num_reads'
};

function vqInitSalmon(salmonQuant, clusters) {
  const el = document.getElementById('section-salmon');
  if (!el) return;

  _SQ.data     = salmonQuant;
  _SQ.clusters = clusters || [];

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Salmon Quantification</div>
        <div class="vq-section-sub">
          Mapping rate: <strong>${(salmonQuant.mapping_rate ?? 0).toFixed(1)}%</strong>
          &nbsp;·&nbsp; ${(salmonQuant.total_reads ?? 0).toLocaleString()} total reads
        </div>
      </div>
      <div style="display:flex;gap:var(--vq-space-2);align-items:center">
        <button class="vq-btn vq-btn--sm ${_SQ.metric==='tpm'?'vq-btn--primary':'vq-btn--ghost'}"
          id="sq-btn-tpm" onclick="_sqSetMetric('tpm')">TPM</button>
        <button class="vq-btn vq-btn--sm ${_SQ.metric==='num_reads'?'vq-btn--primary':'vq-btn--ghost'}"
          id="sq-btn-reads" onclick="_sqSetMetric('num_reads')">Reads</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportSVG(document.querySelector('#sq-viral-wrap svg'),'salmon_viral.svg')">SVG</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportPNG(document.querySelector('#sq-viral-wrap svg'),'salmon_viral.png')">PNG</button>
      </div>
    </div>

    <div style="display:grid;grid-template-columns:1fr 280px;gap:var(--vq-space-4)">
      <div>
        <div class="vq-card">
          <div class="vq-card__header">
            <div class="vq-card__title">Viral Sequences</div>
          </div>
          <div class="vq-card__body">
            <div id="sq-viral-wrap" class="vq-genome-wrap"></div>
          </div>
        </div>
      </div>
      <div>
        <div class="vq-card" style="height:100%">
          <div class="vq-card__header">
            <div class="vq-card__title">Housekeeping Genes</div>
          </div>
          <div class="vq-card__body" id="sq-cons-wrap"></div>
        </div>
      </div>
    </div>

    ${_renderHostHits(salmonQuant.host_viral_hits || [])}
  `;

  _sqDraw();
}

// ── Metric toggle ────────────────────────────────────────────────────────────

function _sqSetMetric(metric) {
  _SQ.metric = metric;
  document.getElementById('sq-btn-tpm')?.classList.toggle('vq-btn--primary',   metric === 'tpm');
  document.getElementById('sq-btn-tpm')?.classList.toggle('vq-btn--ghost',     metric !== 'tpm');
  document.getElementById('sq-btn-reads')?.classList.toggle('vq-btn--primary', metric === 'num_reads');
  document.getElementById('sq-btn-reads')?.classList.toggle('vq-btn--ghost',   metric !== 'num_reads');
  _sqDraw();
}

// ── Draw both panels ─────────────────────────────────────────────────────────

function _sqDraw() {
  _drawViralPanel();
  _drawConsPanel();
}

// ── Viral panel: boxplot per cluster, dot per singleton ───────────────────────

function _drawViralPanel() {
  const wrap = document.getElementById('sq-viral-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const viral   = (_SQ.data.viral_quant || []);
  const metric  = _SQ.metric;
  const label   = metric === 'tpm' ? 'TPM' : 'Reads';

  if (!viral.length) {
    wrap.innerHTML = '<div class="vq-empty">No viral quantification data.</div>';
    return;
  }

  // Map seq_id → cluster info
  const seqCluster = {};
  (_SQ.clusters || []).forEach(c => {
    c.members.forEach(m => { seqCluster[m.seq_id] = c.cluster_id; });
  });

  // Group by cluster; singletons get their own group keyed by seq_id
  const groups = {};
  viral.forEach(v => {
    const key = seqCluster[v.name] != null ? 'cluster_' + seqCluster[v.name] : v.name;
    if (!groups[key]) groups[key] = { label: key, values: [], seqIds: [] };
    groups[key].values.push(v[metric] ?? 0);
    groups[key].seqIds.push(v.name);
  });

  const keys   = Object.keys(groups);
  const PAD_L  = 180;
  const PAD_R  = 30;
  const PAD_T  = 30;
  const PAD_B  = 60;
  const ROW_H  = 44;
  const W      = Math.max((wrap.clientWidth || 700), 500);
  const H      = PAD_T + keys.length * ROW_H + PAD_B;
  const drawW  = W - PAD_L - PAD_R;

  const allVals = viral.map(v => v[metric] ?? 0);
  const maxVal  = d3.max(allVals) || 1;

  const xScale = d3.scaleLinear([0, maxVal], [0, drawW]);

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('width', W)
    .attr('height', H);

  // X axis
  const axG = svg.append('g')
    .attr('transform', `translate(${PAD_L},${PAD_T})`)
    .call(d3.axisTop(xScale).ticks(5).tickSize(-H + PAD_T + PAD_B));

  axG.select('.domain').remove();
  axG.selectAll('.tick line')
    .attr('stroke', '#dde3ec').attr('stroke-dasharray', '3,3');
  axG.selectAll('.tick text')
    .style('font-size', '10px').style('fill', '#94a3b8');

  // X axis label
  svg.append('text')
    .attr('x', PAD_L + drawW / 2)
    .attr('y', H - 12)
    .attr('text-anchor', 'middle')
    .attr('font-size', 11)
    .attr('fill', 'var(--vq-text-2)')
    .text(label);

  // Draw each group
  keys.forEach((key, i) => {
    const grp  = groups[key];
    const vals = grp.values.slice().sort(d3.ascending);
    const y    = PAD_T + i * ROW_H + ROW_H / 2;
    const isSingleton = vals.length === 1;
    const gEl  = svg.append('g').attr('transform', `translate(${PAD_L},0)`);

    // Label
    svg.append('text')
      .attr('x', PAD_L - 8)
      .attr('y', y + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 10)
      .attr('fill', 'var(--vq-text-2)')
      .attr('font-family', 'var(--vq-font-mono)')
      .text(_sqTruncate(grp.label, 22));

    if (isSingleton) {
      // Single dot
      const cx = xScale(vals[0]);
      gEl.append('circle')
        .attr('cx', cx).attr('cy', y)
        .attr('r', 5)
        .attr('fill', 'var(--vq-accent)')
        .attr('stroke', '#fff').attr('stroke-width', 1)
        .attr('cursor', 'pointer')
        .on('mousemove', evt => vqTooltipShow(`
          <div class="vq-tooltip__title">${grp.seqIds[0]}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">${label}</span><span>${vals[0].toFixed(2)}</span>
          </div>`, evt))
        .on('mouseleave', vqTooltipHide)
        .on('click', () => _jumpToViewer(grp.seqIds[0]));
    } else {
      // Boxplot
      const q1  = d3.quantile(vals, 0.25);
      const med = d3.quantile(vals, 0.5);
      const q3  = d3.quantile(vals, 0.75);
      const iqr = q3 - q1;
      const lo  = Math.max(d3.min(vals), q1 - 1.5 * iqr);
      const hi  = Math.min(d3.max(vals), q3 + 1.5 * iqr);

      const bh = 14;

      // Whiskers
      gEl.append('line')
        .attr('x1', xScale(lo)).attr('x2', xScale(hi))
        .attr('y1', y).attr('y2', y)
        .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5);
      [[lo, lo], [hi, hi]].forEach(([x]) => {
        gEl.append('line')
          .attr('x1', xScale(x)).attr('x2', xScale(x))
          .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
          .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5);
      });

      // Box
      gEl.append('rect')
        .attr('x', xScale(q1)).attr('width', Math.max(xScale(q3) - xScale(q1), 2))
        .attr('y', y - bh / 2).attr('height', bh)
        .attr('rx', 2)
        .attr('fill', 'var(--vq-accent)')
        .attr('opacity', 0.25)
        .attr('stroke', 'var(--vq-accent)').attr('stroke-width', 1.5)
        .on('mousemove', evt => vqTooltipShow(`
          <div class="vq-tooltip__title">${grp.label}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Median ${label}</span><span>${med.toFixed(2)}</span>
            <span class="vq-tooltip__key">Q1–Q3</span><span>${q1.toFixed(2)}–${q3.toFixed(2)}</span>
            <span class="vq-tooltip__key">Members</span><span>${vals.length}</span>
          </div>`, evt))
        .on('mouseleave', vqTooltipHide);

      // Median line
      gEl.append('line')
        .attr('x1', xScale(med)).attr('x2', xScale(med))
        .attr('y1', y - bh / 2).attr('y2', y + bh / 2)
        .attr('stroke', 'var(--vq-primary)').attr('stroke-width', 2);

      // Outlier dots
      const outliers = vals.filter(v => v < lo || v > hi);
      outliers.forEach(v => {
        gEl.append('circle')
          .attr('cx', xScale(v)).attr('cy', y)
          .attr('r', 3)
          .attr('fill', 'var(--vq-danger)')
          .attr('opacity', 0.7);
      });
    }
  });

  wrap.appendChild(svg.node());
}

// ── Conserved (housekeeping) panel ───────────────────────────────────────────

function _drawConsPanel() {
  const wrap = document.getElementById('sq-cons-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const cons   = (_SQ.data.conserved_quant || []);
  const metric = _SQ.metric;
  const label  = metric === 'tpm' ? 'TPM' : 'Reads';

  if (!cons.length) {
    wrap.innerHTML = '<div class="vq-empty" style="font-size:12px">No housekeeping data.</div>';
    return;
  }

  // Group by kingdom
  const byKingdom = {};
  cons.forEach(c => {
    const k = c.kingdom || 'Unknown';
    if (!byKingdom[k]) byKingdom[k] = [];
    byKingdom[k].push(c[metric] ?? 0);
  });

  const kingdoms = Object.keys(byKingdom).sort();
  const maxVal   = d3.max(cons.map(c => c[metric] ?? 0)) || 1;

  const W     = 260;
  const PAD_L = 90;
  const PAD_R = 10;
  const ROW_H = 28;
  const BAR_H = 10;
  const H     = kingdoms.length * ROW_H + 40;

  const xScale = d3.scaleLinear([0, maxVal], [0, W - PAD_L - PAD_R]);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('width', W)
    .attr('height', H);

  svg.append('text')
    .attr('x', PAD_L + (W - PAD_L - PAD_R) / 2)
    .attr('y', 14)
    .attr('text-anchor', 'middle')
    .attr('font-size', 10)
    .attr('fill', 'var(--vq-text-3)')
    .text(`Median ${label} by kingdom`);

  kingdoms.forEach((k, i) => {
    const vals   = byKingdom[k].slice().sort(d3.ascending);
    const median = d3.quantile(vals, 0.5) ?? 0;
    const y      = 24 + i * ROW_H;
    const bw     = Math.max(xScale(median), 2);

    svg.append('text')
      .attr('x', PAD_L - 6).attr('y', y + BAR_H / 2 + 3)
      .attr('text-anchor', 'end')
      .attr('font-size', 9)
      .attr('fill', 'var(--vq-text-2)')
      .text(k);

    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', bw).attr('height', BAR_H)
      .attr('rx', 2)
      .attr('fill', 'var(--vq-success)')
      .attr('opacity', 0.8)
      .on('mousemove', evt => vqTooltipShow(`
        <div class="vq-tooltip__title">${k}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Median ${label}</span><span>${median.toFixed(2)}</span>
          <span class="vq-tooltip__key">Genes</span><span>${vals.length}</span>
        </div>`, evt))
      .on('mouseleave', vqTooltipHide);

    svg.append('text')
      .attr('x', PAD_L + bw + 4).attr('y', y + BAR_H / 2 + 3)
      .attr('font-size', 9).attr('fill', 'var(--vq-text-3)')
      .text(median.toFixed(1));
  });

  wrap.appendChild(svg.node());
}

// ── Host–viral hits table ─────────────────────────────────────────────────────

function _renderHostHits(hits) {
  if (!hits.length) return '';

  const rows = hits.map(h => `
    <tr>
      <td class="vq-td" style="font-family:var(--vq-font-mono);font-size:11px">${h.transcript?.name ?? '—'}</td>
      <td class="vq-td" style="font-family:var(--vq-font-mono);font-size:11px">${h.blast_hits?.[0]?.viral_seq_id ?? '—'}</td>
      <td class="vq-td">${h.blast_hits?.[0]?.pident?.toFixed(1) ?? '—'}%</td>
      <td class="vq-td">${h.blast_hits?.[0]?.qcovhsp ?? '—'}%</td>
      <td class="vq-td">${h.transcript?.tpm?.toFixed(2) ?? '—'}</td>
    </tr>`).join('');

  return `
    <div class="vq-card" style="margin-top:var(--vq-space-4)">
      <div class="vq-card__header">
        <div class="vq-card__title">Host–Viral Transcriptome Hits</div>
        <span style="font-size:var(--vq-text-xs);color:var(--vq-text-3)">${hits.length} transcript${hits.length>1?'s':''} with viral similarity</span>
      </div>
      <div class="vq-card__body" style="overflow-x:auto">
        <table style="width:100%;border-collapse:collapse;font-size:12px">
          <thead>
            <tr>
              ${['Host Transcript','Viral Match','Identity','Coverage','TPM'].map(
                h => `<th style="text-align:left;padding:6px 8px;border-bottom:1px solid var(--vq-border);color:var(--vq-text-3);font-weight:500;font-size:11px">${h}</th>`
              ).join('')}
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </div>`;
}

// ── Helpers ──────────────────────────────────────────────────────────────────

function _sqTruncate(str, n) {
  return str.length > n ? str.slice(0, n - 1) + '…' : str;
}

function _jumpToViewer(seqId) {
  document.querySelector('[data-section="viewer"]')?.click();
  setTimeout(() => {
    const target = document.getElementById('seq-card-' + seqId);
    if (target) {
      target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      target.classList.add('vq-seq-card--highlight');
      setTimeout(() => target.classList.remove('vq-seq-card--highlight'), 2000);
    }
  }, 150);
}
