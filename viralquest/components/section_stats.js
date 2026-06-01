/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_stats.js — Section 1: General Statistics dashboard
   Top row: hero KPIs.
   Charts row: family donut · HMM bars · Salmon panel.
   Detail row: optional cards for input file, CAP3, LLM.
   ============================================================ */

function vqInitStats(report) {
  const el = document.getElementById('section-stats');
  if (!el) return;

  const esc    = VQ.esc;
  const meta   = report.meta         || {};
  const sum    = report.summary      || {};
  const blast  = report.blast_stats  || {};
  const hmm    = report.hmm_stats    || {};
  const llm    = report.llm_stats    || {};
  const salmon = report.salmon_stats || {};
  const seqs   = report.sequences    || [];

  // ── Aggregates ─────────────────────────────────────────────────────────
  const totalSeqs      = sum.total_sequences ?? seqs.length;
  const confirmedViral = sum.confirmed_viral ?? seqs.filter(s => s.is_viral).length;
  const totalOrfs      = sum.total_orfs      ?? seqs.reduce((a, s) => a + (s.orfs || []).length, 0);
  const totalClusters  = sum.total_clusters  ?? (report.clusters || []).length;

  const fmtNum = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();
  const fmtPct = v => (v == null || isNaN(v)) ? '—' : Number(v).toFixed(1) + '%';

  // ── Page scaffold ──────────────────────────────────────────────────────
  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">General Statistics</div>
        <div class="vq-section-sub">
          ${esc(meta.input_file?.name || 'Sample report')}
          &nbsp;·&nbsp; run ${esc(meta.timestamp ? new Date(meta.timestamp).toLocaleString() : '—')}
          &nbsp;·&nbsp; ViralQuest v${esc(meta.viralquest_version || '?')}
        </div>
      </div>
      <div class="vq-section-actions">
        <button class="vq-btn vq-btn--ghost" id="stats-print-pdf" type="button">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
            <polyline points="14 2 14 8 20 8"/>
          </svg>
          Print / Save PDF
        </button>
      </div>
    </div>

    <div class="vq-stats-page">

      <!-- Hero KPI row -->
      <div class="vq-stats-row vq-stats-row--top">
        ${_chip('Total Sequences', fmtNum(totalSeqs), 'accent',
                meta.input_file?.name ? esc(meta.input_file.name) : '')}
        ${_chip('Confirmed Viral', fmtNum(confirmedViral), 'success',
                _pct(confirmedViral, totalSeqs) + ' of input')}
        ${_chip('Total ORFs',      fmtNum(totalOrfs),      '',
                'across confirmed sequences')}
        ${_chip('Clusters',        fmtNum(totalClusters),  'accent',
                'unique viral species')}
      </div>

      <!-- Charts row -->
      <div class="vq-stats-row vq-stats-row--charts">
        <!-- 1. Family / phylum half-donut -->
        <div class="vq-chart-card" id="stats-tax-card">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Viral Family Distribution</div>
              <div class="vq-chart-card__sub" id="stats-tax-sub">across confirmed sequences</div>
            </div>
            <div class="vq-toggle" role="tablist" aria-label="Taxonomy rank">
              <button class="vq-toggle__btn active" type="button" data-rank="family">Family</button>
              <button class="vq-toggle__btn"        type="button" data-rank="phylum">Phylum</button>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;align-items:center;flex:1">
              <div id="stats-tax-svg" style="display:flex;justify-content:center;align-items:center"></div>
              <div id="stats-tax-legend" class="vq-chart-legend"></div>
            </div>
          </div>
        </div>

        <!-- 2. HMM bar chart -->
        <div class="vq-chart-card" id="stats-hmm-card">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">HMM Database Hits</div>
              <div class="vq-chart-card__sub">viral filter databases</div>
            </div>
            <div class="vq-chart-card__big" id="stats-hmm-total">${fmtNum(
              (hmm.rvdb_hits||0) + (hmm.vfam_hits||0) + (hmm.eggnog_hits ?? hmm.eggnong_hits ?? 0) + (hmm.pfam_hits||0)
            )}</div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-hmm-svg" style="width:100%"></div>
          </div>
        </div>

        <!-- 3. Salmon -->
        <div class="vq-chart-card" id="stats-salmon-card">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Salmon Quantification</div>
              <div class="vq-chart-card__sub" id="stats-salmon-sub">
                ${salmon.present
                  ? `${fmtNum(salmon.total_reads)} reads`
                  : 'no quant data'}
              </div>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-salmon-svg" style="width:100%"></div>
          </div>
        </div>
      </div>

      <!-- Pipeline funnel (blue) -->
      <div class="vq-stats-row" style="grid-template-columns:1fr">
        <div class="vq-chart-card" id="stats-funnel-card"
             style="min-height:auto;background:linear-gradient(135deg,#eaf3fb 0%,var(--vq-surface) 55%);border-color:#c2ddf5">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Detection Pipeline</div>
              <div class="vq-chart-card__sub">sequences filtered at each stage</div>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-funnel-svg" style="width:100%"></div>
          </div>
        </div>
      </div>

      <!-- BLASTx identity (green) + sequence length (blue) -->
      <div class="vq-stats-row" style="grid-template-columns:1fr 1fr">
        <div class="vq-chart-card" id="stats-identity-card"
             style="background:linear-gradient(135deg,#e8f5ef 0%,var(--vq-surface) 55%);border-color:#b8dcc8">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">BLASTx Identity Distribution</div>
              <div class="vq-chart-card__sub">best hit % identity · confirmed sequences</div>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-identity-svg" style="width:100%"></div>
          </div>
        </div>
        <div class="vq-chart-card" id="stats-length-card"
             style="background:linear-gradient(135deg,#edf4ff 0%,var(--vq-surface) 55%);border-color:#bdd5f0">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Sequence Length Distribution</div>
              <div class="vq-chart-card__sub">confirmed viral sequences only</div>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-length-svg" style="width:100%"></div>
          </div>
        </div>
      </div>

      <!-- Top detected species -->
      <div class="vq-stats-row" style="grid-template-columns:1fr">
        <div class="vq-chart-card" id="stats-species-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Top Detected Species</div>
              <div class="vq-chart-card__sub" id="stats-species-sub"></div>
            </div>
          </div>
          <div class="vq-chart-card__body" style="padding:0">
            <div id="stats-species-table"></div>
          </div>
        </div>
      </div>

      <!-- Detail row -->
      <div class="vq-stats-row" id="stats-detail-row" style="grid-template-columns:repeat(auto-fit,minmax(280px,1fr))">
        ${_detailCard('BLAST', `
          ${_miniRow('RefSeq candidates', fmtNum(blast.refseq_unique_seqs))}
          ${_miniRow('NR confirmed',      fmtNum(blast.nr_unique_seqs))}
          ${_miniRow('BLASTn hits',       fmtNum(blast.blastn_unique_seqs))}
        `)}
        ${sum.cap3_used ? _detailCard('CAP3 Assembly', `
          ${_miniRow('Contigs',  fmtNum(sum.cap3_contigs))}
          ${_miniRow('Singlets', fmtNum(sum.cap3_singlets))}
        `) : ''}
        ${llm.present ? _detailCard('LLM Scoring', `
          ${_miniRow('Model',         esc(llm.model || '—'))}
          ${_miniRow('Scored',        fmtNum(llm.scored))}
          ${_miniRow('Viral known',   fmtNum(llm.viral_known))}
          ${_miniRow('Viral unknown', fmtNum(llm.viral_unknown))}
          ${_miniRow('Avg VQ score',  llm.avg_score != null ? Number(llm.avg_score).toFixed(1) : '—')}
        `) : ''}
        ${salmon.present ? _detailCard('Salmon Detail', `
          ${_miniRow('Mapping rate',    fmtPct(salmon.mapping_rate))}
          ${_miniRow('Total reads',     fmtNum(salmon.total_reads))}
          ${_miniRow('Host–viral hits', fmtNum(salmon.host_viral_hits))}
        `) : ''}
      </div>
    </div>
  `;

  // ── Wire up ────────────────────────────────────────────────────────────
  document.getElementById('stats-print-pdf')?.addEventListener('click', () => {
    VQ.exportPagePDF('ViralQuest Statistics');
  });

  // Taxonomy donut
  let currentRank = 'family';
  const renderTax = () => _renderTaxDonut(seqs, currentRank);
  document.querySelectorAll('#stats-tax-card [data-rank]').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('#stats-tax-card [data-rank]').forEach(b => {
        b.classList.toggle('active', b === btn);
      });
      currentRank = btn.dataset.rank;
      renderTax();
    });
  });
  renderTax();

  // HMM bars
  _renderHMMBars(hmm);

  // Salmon
  if (salmon.present) _renderSalmonChart(salmon);
  else                _renderSalmonEmpty();

  // New charts
  _renderFunnel(blast);
  _renderIdentityHistogram(seqs);
  _renderLengthHistogram(seqs);
  _renderTopSpecies(seqs);

  // ResizeObserver so charts re-fit when the section is first revealed
  // (children measure 0 while the section is still display:none).
  let lastW = 0;
  const ro = new ResizeObserver(() => {
    const w = el.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    renderTax();
    _renderHMMBars(hmm);
    if (salmon.present) _renderSalmonChart(salmon);
    _renderFunnel(blast);
    _renderIdentityHistogram(seqs);
    _renderLengthHistogram(seqs);
  });
  ro.observe(el);
}

// ────────────────────────────────────────────────────────────────────────
//  Chip / detail helpers
// ────────────────────────────────────────────────────────────────────────

function _chip(label, value, mod = '', sub = '') {
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${VQ.esc(label)}</div>
      <div class="vq-stat-chip__value" title="${VQ.esc(value)}">${value}</div>
      ${sub ? `<div class="vq-stat-chip__sub" title="${VQ.esc(sub)}">${sub}</div>` : ''}
    </div>`;
}

function _detailCard(title, bodyHtml) {
  return `
    <div class="vq-chart-card" style="min-height:auto;padding:14px 16px">
      <div class="vq-chart-card__title" style="margin-bottom:10px">${VQ.esc(title)}</div>
      ${bodyHtml}
    </div>`;
}

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

function _pct(a, b) {
  if (!b || isNaN(a) || isNaN(b)) return '—';
  return ((a / b) * 100).toFixed(1) + '%';
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 1 — Half-donut: family/phylum distribution
// ────────────────────────────────────────────────────────────────────────

const _TAX_DONUT_PALETTE = [
  'var(--vq-tax-flaviviridae)','var(--vq-tax-parvoviridae)',
  'var(--vq-tax-phenuiviridae)','var(--vq-tax-nodaviridae)',
  'var(--vq-tax-rhabdoviridae)','var(--vq-tax-tombusviridae)',
  '#0891b2','#7c3aed','#b45309','#be185d','var(--vq-tax-other)',
];

function _familyColor(family, i) {
  const cssVar = '--vq-tax-' + (family || '').toLowerCase();
  const computed = getComputedStyle(document.documentElement).getPropertyValue(cssVar).trim();
  return computed || _TAX_DONUT_PALETTE[i % _TAX_DONUT_PALETTE.length];
}

function _renderTaxDonut(sequences, rank) {
  const wrap = document.getElementById('stats-tax-svg');
  const lgd  = document.getElementById('stats-tax-legend');
  const sub  = document.getElementById('stats-tax-sub');
  if (!wrap || !lgd) return;
  wrap.innerHTML = '';
  lgd.innerHTML  = '';

  const viral = sequences.filter(s => s.is_viral);
  const counts = {};
  viral.forEach(s => {
    const key = s.taxonomy?.[rank] || 'Unclassified';
    counts[key] = (counts[key] || 0) + 1;
  });

  const entries = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  // Collapse the long tail into "Other" so the donut stays readable.
  const MAX_SLICES = 7;
  let display = entries;
  if (entries.length > MAX_SLICES) {
    const top  = entries.slice(0, MAX_SLICES - 1);
    const rest = entries.slice(MAX_SLICES - 1);
    const otherCount = rest.reduce((a, [, n]) => a + n, 0);
    display = [...top, ['Other (' + rest.length + ')', otherCount]];
  }

  const total = entries.reduce((a, [, n]) => a + n, 0);
  if (sub) sub.textContent = `across ${total} confirmed sequence${total !== 1 ? 's' : ''}`;

  if (!total) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:20px">No viral sequences.</div>';
    return;
  }

  // Geometry
  const W   = 220;
  const H   = 130;
  const cx  = W / 2;
  const cy  = H - 14;
  const rOuter = 110;
  const rInner = 68;

  const arc = d3.arc()
    .innerRadius(rInner).outerRadius(rOuter)
    .cornerRadius(4)
    .padAngle(0.015);

  // Half-donut: π/2 from -π/2 to π/2 won't work for half-pie at bottom;
  // we want a top-opening semicircle: start at -π/2 (top), through 0, to π/2.
  // Actually we want the FLAT side at the bottom (semicircle on top).
  // That means start = -π/2, end = π/2 (sweep clockwise across the top).
  const pie = d3.pie()
    .value(d => d[1])
    .sort(null)
    .startAngle(-Math.PI / 2)
    .endAngle(Math.PI / 2);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('width',  W)
    .attr('height', H)
    .style('overflow', 'visible');

  const g = svg.append('g').attr('transform', `translate(${cx}, ${cy})`);

  const slices = pie(display);
  slices.forEach((s, i) => {
    g.append('path')
      .attr('d', arc(s))
      .attr('fill', _familyColor(s.data[0], i))
      .attr('opacity', 0.92)
      .style('cursor', 'default')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(s.data[0])}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${s.data[1]}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${((s.data[1] / total) * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Center label (big number above the flat side, inside the donut hole)
  const labelG = svg.append('g').attr('transform', `translate(${cx}, ${cy - 20})`);
  labelG.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 22)
    .attr('font-weight', 700)
    .attr('fill', 'var(--vq-primary)')
    .attr('y', 0)
    .text(entries.length);
  labelG.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 10)
    .attr('fill', 'var(--vq-text-3)')
    .attr('y', 14)
    .text(rank === 'family' ? 'families' : 'phyla');

  wrap.appendChild(svg.node());

  // Legend rows
  display.forEach(([name, n], i) => {
    const pct = (n / total) * 100;
    const row = document.createElement('div');
    row.className = 'vq-chart-legend__row';
    row.innerHTML = `
      <span class="vq-chart-legend__dot" style="background:${_familyColor(name, i)}"></span>
      <span style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap"
            title="${VQ.esc(name)}">${VQ.esc(name)}</span>
      <span class="vq-chart-legend__count">${n}</span>
      <span class="vq-chart-legend__pct">${pct.toFixed(0)}%</span>
    `;
    lgd.appendChild(row);
  });
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

  const W   = Math.max(wrap.clientWidth || 0, 240);
  const H   = 200;
  const rate = Math.max(0, Math.min(100, salmon.mapping_rate ?? 0));

  // Geometry: semicircle gauge, big number in the middle.
  const cx     = W / 2;
  const cy     = 130;
  const rOuter = 92;
  const rInner = 70;

  const arcBg = d3.arc()
    .innerRadius(rInner).outerRadius(rOuter)
    .startAngle(-Math.PI / 2).endAngle(Math.PI / 2)
    .cornerRadius(4);
  const arcVal = d3.arc()
    .innerRadius(rInner).outerRadius(rOuter)
    .startAngle(-Math.PI / 2)
    .endAngle(-Math.PI / 2 + (rate / 100) * Math.PI)
    .cornerRadius(4);

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMidYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${cx}, ${cy})`);

  g.append('path').attr('d', arcBg()).attr('fill', 'var(--vq-bg)');

  const rateColor = rate >= 30 ? 'var(--vq-success)'
                  : rate >= 15 ? 'var(--vq-accent)'
                              : 'var(--vq-warning)';
  g.append('path').attr('d', arcVal()).attr('fill', rateColor)
    .on('mousemove', evt => VQ.tooltipShow(`
      <div class="vq-tooltip__title">Mapping rate</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Reads aligned</span><span>${rate.toFixed(1)}%</span>
      </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // Inner number
  g.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 28)
    .attr('font-weight', 700)
    .attr('letter-spacing', '-0.02em')
    .attr('fill', 'var(--vq-primary)')
    .attr('y', -10)
    .text(rate.toFixed(1) + '%');
  g.append('text')
    .attr('text-anchor', 'middle')
    .attr('font-size', 10)
    .attr('fill', 'var(--vq-text-3)')
    .attr('y', 8)
    .text('mapping rate');

  // Gauge end labels (0 / 100)
  svg.append('text').attr('x', cx - rOuter - 4).attr('y', cy + 12)
    .attr('text-anchor', 'end').attr('font-size', 10).attr('fill', 'var(--vq-text-3)').text('0');
  svg.append('text').attr('x', cx + rOuter + 4).attr('y', cy + 12)
    .attr('font-size', 10).attr('fill', 'var(--vq-text-3)').text('100%');

  // Bottom strip — reads / host-viral hits
  const tot  = salmon.total_reads     || 0;
  const host = salmon.host_viral_hits || 0;

  svg.append('text').attr('x', cx).attr('y', cy + 36)
    .attr('text-anchor', 'middle')
    .attr('font-size', 11)
    .attr('fill', 'var(--vq-text-2)')
    .text(`${tot.toLocaleString()} reads · ${host} host–viral hits`);

  wrap.appendChild(svg.node());
}

function _renderSalmonEmpty() {
  const wrap = document.getElementById('stats-salmon-svg');
  if (wrap) wrap.innerHTML =
    '<div class="vq-empty" style="padding:40px 16px;font-size:12px">No Salmon quantification data.</div>';
}

// ────────────────────────────────────────────────────────────────────────
//  Chart 4 — Pipeline funnel
// ────────────────────────────────────────────────────────────────────────

function _renderFunnel(blast) {
  const wrap = document.getElementById('stats-funnel-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const steps = [
    { label: 'Input sequences',    value: blast.total_input         ?? null, color: 'var(--vq-primary)'       },
    { label: 'RefSeq BLASTx hits', value: blast.refseq_unique_seqs  ?? 0,   color: 'var(--vq-accent)'        },
    { label: 'Viral flagged',      value: blast.total_viral_flagged  ?? null, color: 'var(--vq-primary-light)' },
    { label: 'NR confirmed',       value: blast.total_confirmed      ?? 0,   color: 'var(--vq-success)'       },
  ].filter(s => s.value !== null);

  if (!steps.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:20px;font-size:12px">No pipeline data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 400);
  const PAD_L = 178;
  const PAD_R = 130;
  const ROW_H = 46;
  const BAR_H = 22;
  const H     = steps.length * ROW_H + 12;
  const drawW = W - PAD_L - PAD_R;
  const maxV  = steps[0].value || 1;

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  steps.forEach((step, i) => {
    const y    = 6 + i * ROW_H;
    const barW = Math.max(maxV > 0 ? (step.value / maxV) * drawW : 0, step.value > 0 ? 4 : 0);
    const pct  = i === 0 ? '100%' : (steps[0].value > 0
      ? (step.value / steps[0].value * 100).toFixed(1) + '%' : '—');

    svg.append('text')
      .attr('x', PAD_L - 10).attr('y', y + BAR_H / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', 12).attr('font-weight', 500)
      .attr('fill', 'var(--vq-text-2)').text(step.label);

    svg.append('rect')
      .attr('x', PAD_L).attr('y', y)
      .attr('width', drawW).attr('height', BAR_H)
      .attr('rx', BAR_H / 2).attr('fill', 'var(--vq-bg)');

    if (barW > 0) {
      svg.append('rect')
        .attr('x', PAD_L).attr('y', y)
        .attr('width', barW).attr('height', BAR_H)
        .attr('rx', BAR_H / 2).attr('fill', step.color).attr('opacity', 0.86)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${VQ.esc(step.label)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Count</span><span>${step.value.toLocaleString()}</span>
            <span class="vq-tooltip__key">% of input</span><span>${pct}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    }

    svg.append('text')
      .attr('x', PAD_L + barW + 10).attr('y', y + BAR_H / 2 + 4)
      .attr('font-size', 12).attr('font-weight', 600)
      .attr('fill', 'var(--vq-text)').attr('font-variant-numeric', 'tabular-nums')
      .text(`${step.value.toLocaleString()}  ·  ${pct}`);

    if (i < steps.length - 1) {
      const cx = PAD_L + Math.max((steps[i + 1].value / maxV) * drawW / 2, 14);
      svg.append('line')
        .attr('x1', cx).attr('y1', y + BAR_H + 2)
        .attr('x2', cx).attr('y2', y + ROW_H - 2)
        .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '3,2');
    }
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 5 — BLASTx identity histogram (green card)
// ────────────────────────────────────────────────────────────────────────

function _renderIdentityHistogram(sequences) {
  const wrap = document.getElementById('stats-identity-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const values = sequences
    .filter(s => s.is_viral)
    .map(s => {
      const hit = (s.blastx_nr_hits && s.blastx_nr_hits[0]) || (s.blastx_hits && s.blastx_hits[0]);
      return hit ? hit.pct_identity : null;
    })
    .filter(v => v != null);

  if (!values.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No BLASTx data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 240);
  const PAD_L = 34; const PAD_R = 10;
  const PAD_T = 10; const PAD_B = 26;
  const H     = 185;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  const bins = d3.bin().domain([0, 100]).thresholds([0,10,20,30,40,50,60,70,80,90,100])(values);

  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([0, 100], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  function _idColor(x0) {
    if (x0 >= 90) return 'var(--vq-success)';
    if (x0 >= 70) return 'var(--vq-accent)';
    if (x0 >= 50) return 'var(--vq-warning)';
    return 'var(--vq-danger)';
  }

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    g.append('rect')
      .attr('x', x).attr('y', yScale(bin.length))
      .attr('width', bw).attr('height', bh)
      .attr('rx', 3).attr('fill', _idColor(bin.x0)).attr('opacity', 0.83)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${bin.x0}–${bin.x1}% identity</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / values.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-border)').attr('stroke-width', 1);
  [0, 25, 50, 75, 100].forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 4)
      .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 14)
      .attr('text-anchor', 'middle').attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
      .text(tick + '%');
  });

  // Y max label
  g.append('text').attr('x', -4).attr('y', 4)
    .attr('text-anchor', 'end').attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
    .text(yMax);

  // Colour legend
  const legend = [
    { label: '≥90% (known)',    color: 'var(--vq-success)' },
    { label: '70–90%',          color: 'var(--vq-accent)'  },
    { label: '50–70%',          color: 'var(--vq-warning)' },
    { label: '<50% (novel)',     color: 'var(--vq-danger)'  },
  ];
  let lx = PAD_L;
  const ly = H - 2;
  legend.forEach(l => {
    const row = svg.append('g').attr('transform', `translate(${lx},${ly})`);
    row.append('rect').attr('width', 8).attr('height', 8).attr('rx', 2)
      .attr('y', -8).attr('fill', l.color).attr('opacity', 0.85);
    row.append('text').attr('x', 11).attr('y', -1)
      .attr('font-size', 9).attr('fill', 'var(--vq-text-3)').text(l.label);
    lx += 78;
  });

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 6 — Sequence length histogram, viral only (blue card)
// ────────────────────────────────────────────────────────────────────────

function _renderLengthHistogram(sequences) {
  const wrap = document.getElementById('stats-length-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const lengths = sequences
    .filter(s => s.is_viral && s.length > 0)
    .map(s => s.length);

  if (!lengths.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No sequence data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 240);
  const PAD_L = 34; const PAD_R = 10;
  const PAD_T = 10; const PAD_B = 26;
  const H     = 185;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  const maxL = d3.max(lengths);
  const bins = d3.bin().domain([0, maxL * 1.02]).thresholds(9)(lengths);

  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([0, maxL * 1.02], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  function _fmtLen(v) {
    return v >= 1000 ? (v / 1000).toFixed(v % 1000 === 0 ? 0 : 1) + 'k' : String(Math.round(v));
  }

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  bins.forEach((bin, i) => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    const opacity = 0.5 + (i / Math.max(bins.length - 1, 1)) * 0.45;
    g.append('rect')
      .attr('x', x).attr('y', yScale(bin.length))
      .attr('width', bw).attr('height', bh)
      .attr('rx', 3).attr('fill', 'var(--vq-accent)').attr('opacity', opacity)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${_fmtLen(bin.x0)} – ${_fmtLen(bin.x1)} bp</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / lengths.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-border)').attr('stroke-width', 1);
  xScale.ticks(5).forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 4)
      .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 14)
      .attr('text-anchor', 'middle').attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
      .text(_fmtLen(tick));
  });

  // Y max label
  g.append('text').attr('x', -4).attr('y', 4)
    .attr('text-anchor', 'end').attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
    .text(yMax);

  // X axis label
  svg.append('text')
    .attr('x', PAD_L + drawW / 2).attr('y', H)
    .attr('text-anchor', 'middle').attr('font-size', 10).attr('fill', 'var(--vq-text-3)')
    .text('bp');

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 7 — Top detected species table
// ────────────────────────────────────────────────────────────────────────

function _renderTopSpecies(sequences) {
  const tableWrap = document.getElementById('stats-species-table');
  const sub       = document.getElementById('stats-species-sub');
  if (!tableWrap) return;

  const viral = sequences.filter(s => s.is_viral && s.taxonomy && s.taxonomy.species);
  const counts = {};
  viral.forEach(s => {
    const sp = s.taxonomy.species;
    if (!counts[sp]) counts[sp] = { count: 0, family: s.taxonomy.family || '—', order: s.taxonomy.order || '—' };
    counts[sp].count++;
  });

  const entries = Object.entries(counts)
    .sort((a, b) => b[1].count - a[1].count)
    .slice(0, 10);

  if (sub) {
    const uniq = Object.keys(counts).length;
    sub.textContent = `top ${entries.length} of ${uniq} species · ${viral.length} sequence${viral.length !== 1 ? 's' : ''}`;
  }

  if (!entries.length) {
    tableWrap.innerHTML = '<div class="vq-empty" style="padding:20px;font-size:12px">No taxonomy data.</div>';
    return;
  }

  const total = viral.length;
  tableWrap.innerHTML = `
    <table class="vq-table">
      <thead>
        <tr>
          <th style="width:32px">#</th>
          <th>Species</th>
          <th>Family</th>
          <th style="width:80px">Count</th>
          <th style="width:140px">% of viral</th>
        </tr>
      </thead>
      <tbody>
        ${entries.map(([sp, d], i) => {
          const pct = (d.count / total * 100).toFixed(1);
          return `
          <tr>
            <td class="vq-td--num" style="color:var(--vq-text-3)">${i + 1}</td>
            <td style="font-style:italic">${VQ.esc(sp)}</td>
            <td style="color:var(--vq-text-2)">${VQ.esc(d.family)}</td>
            <td class="vq-td--num" style="font-weight:600">${d.count}</td>
            <td>
              <div style="display:flex;align-items:center;gap:6px">
                <div style="flex:1;height:6px;background:var(--vq-bg);border-radius:3px;min-width:60px">
                  <div style="height:100%;width:${pct}%;background:var(--vq-accent);border-radius:3px;opacity:.8"></div>
                </div>
                <span style="font-size:11px;color:var(--vq-text-3);font-variant-numeric:tabular-nums;white-space:nowrap">${pct}%</span>
              </div>
            </td>
          </tr>`;
        }).join('')}
      </tbody>
    </table>
  `;
}


window.vqInitStats = vqInitStats;
})();
