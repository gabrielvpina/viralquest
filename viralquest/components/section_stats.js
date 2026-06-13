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
        ${_chip('Total ORFs', fmtNum(totalOrfs), '',
                'across confirmed sequences')}
        ${_chip('Clusters', fmtNum(totalClusters), 'accent',
                'unique viral species')}
      </div>

      <!-- Uniform 3-column card grid — every card is 1/3 of the row.
           Conditional cards (Salmon, LLM, CAP3) only enter the DOM when
           their data is present; the grid reflows automatically. -->
      <div style="display:grid;grid-template-columns:repeat(3,1fr);gap:var(--vq-space-3)">

        <!-- 1. Detection Pipeline — always -->
        <div class="vq-chart-card" id="stats-funnel-card" style="min-height:auto">
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

        <!-- 2. Viral Family Distribution — always -->
        <div class="vq-chart-card" id="stats-tax-card" style="min-height:auto">
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
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;align-items:center">
              <div id="stats-tax-svg" style="display:flex;justify-content:center;align-items:center"></div>
              <div id="stats-tax-legend" class="vq-chart-legend"></div>
            </div>
          </div>
        </div>

        <!-- 3. HMM Database Hits — always -->
        <div class="vq-chart-card" id="stats-hmm-card" style="min-height:auto">
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

        <!-- 4. BLASTx — always -->
        <div class="vq-chart-card" id="stats-identity-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div class="vq-chart-card__title" id="stats-blastx-title">BLASTx Identity</div>
            <div class="vq-toggle" id="blast-toggle-x" role="tablist">
              <button class="vq-toggle__btn active" type="button" data-blast-metric="identity">Identity</button>
              <button class="vq-toggle__btn"        type="button" data-blast-metric="coverage">Coverage</button>
            </div>
          </div>
          <div class="vq-chart-card__body" style="padding-top:4px">
            <div id="stats-identity-svg" style="width:100%"></div>
          </div>
        </div>

        <!-- 5. BLASTn — always -->
        <div class="vq-chart-card" id="stats-blastn-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div class="vq-chart-card__title" id="stats-blastn-title">BLASTn Identity</div>
            <div class="vq-toggle" id="blast-toggle-n" role="tablist">
              <button class="vq-toggle__btn active" type="button" data-blast-metric="identity">Identity</button>
              <button class="vq-toggle__btn"        type="button" data-blast-metric="coverage">Coverage</button>
            </div>
          </div>
          <div class="vq-chart-card__body" style="padding-top:4px">
            <div id="stats-blastn-svg" style="width:100%"></div>
          </div>
        </div>

        <!-- 6. Viral Length Distribution — always -->
        <div class="vq-chart-card" id="stats-length-card"
             style="min-height:auto;background:var(--vq-c-dark-bg);border-color:#0c2238">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title" style="color:#fff">Viral Length Distribution</div>
              <div class="vq-chart-card__sub" style="color:rgba(255,255,255,.5)">viral · 500 bp bins</div>
            </div>
          </div>
          <div class="vq-chart-card__body">
            <div id="stats-length-svg" style="width:100%"></div>
          </div>
        </div>

        <!-- 7. Top Detected Species — always -->
        <div class="vq-chart-card" id="stats-species-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Top Detected Species</div>
              <div class="vq-chart-card__sub" id="stats-species-sub"></div>
            </div>
          </div>
          <div class="vq-chart-card__body" style="padding:0;justify-content:flex-start">
            <div id="stats-species-table"></div>
          </div>
        </div>

        <!-- 8. NR Classification: Virus vs Phage — always -->
        <div class="vq-chart-card" id="stats-nrclass-card" style="min-height:auto;align-self:start">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title" id="stats-nrclass-title">NR Classification</div>
              <div class="vq-chart-card__sub" id="stats-nrclass-sub">virus · phage breakdown</div>
            </div>
          </div>
          <div class="vq-chart-card__body" style="justify-content:flex-start">
            <div id="stats-nrclass-legend" style="width:100%"></div>
          </div>
        </div>

        <!-- 9. Salmon Quantification — only when salmon was run -->
        ${salmon.present ? `
        <div class="vq-chart-card" id="stats-salmon-card" style="min-height:auto;align-self:start">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Salmon Quantification</div>
              <div class="vq-chart-card__sub">
                ${fmtNum(salmon.total_reads)} reads
                ${salmon.pathway ? ' · ' + esc(salmon.pathway) + ' pathway' : ''}
              </div>
            </div>
          </div>
          <div class="vq-chart-card__body" style="flex-direction:column;justify-content:flex-start;gap:0">
            <div id="stats-salmon-svg" style="width:100%"></div>
            <div id="stats-salmon-rows" style="width:100%;padding:0 4px 4px"></div>
          </div>
        </div>` : ''}

        <!-- 10. LLM Scoring — only when LLM was used -->
        ${llm.present ? `
        <div class="vq-chart-card" id="stats-llm-card" style="min-height:auto;align-self:start">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">LLM Scoring</div>
              <div class="vq-chart-card__sub">${esc(llm.model || '—')}</div>
            </div>
            <div class="vq-chart-card__big">${llm.avg_score != null ? Number(llm.avg_score).toFixed(1) : '—'}</div>
          </div>
          <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
            ${_miniRow('Scored',        fmtNum(llm.scored))}
            ${_miniRow('Viral known',   fmtNum(llm.viral_known))}
            ${_miniRow('Viral unknown', fmtNum(llm.viral_unknown))}
            ${llm.api_error ? _miniRow('API errors', fmtNum(llm.api_error)) : ''}
          </div>
        </div>` : ''}

        <!-- 11. CAP3 Assembly — only when CAP3 was used -->
        ${sum.cap3_used ? `
        <div class="vq-chart-card" id="stats-cap3-card" style="min-height:auto;align-self:start">
          <div class="vq-chart-card__head">
            <div class="vq-chart-card__title">CAP3 Assembly</div>
          </div>
          <div class="vq-chart-card__body" style="padding-top:8px;justify-content:flex-start">
            ${_miniRow('Contigs',  fmtNum(sum.cap3_contigs))}
            ${_miniRow('Singlets', fmtNum(sum.cap3_singlets))}
          </div>
        </div>` : ''}

      </div>
    </div>
  `;

  // ── Wire up ────────────────────────────────────────────────────────────

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

  // Salmon (card only exists in DOM when salmon.present)
  if (salmon.present) _renderSalmonChart(salmon);

  // Pipeline funnel + length + species + NR classification
  _renderFunnel(blast);
  _renderLengthHistogram(seqs);
  _renderTopSpecies(seqs);
  _renderNRClassification(seqs);

  // BLAST toggles — both cards share the same metric state and stay in sync.
  // Either toggle drives both charts simultaneously.
  let blastMetric = 'identity';
  const _allBlastBtns = () => [
    ...document.querySelectorAll('#blast-toggle-x [data-blast-metric]'),
    ...document.querySelectorAll('#blast-toggle-n [data-blast-metric]'),
  ];
  const _syncBlastToggles = () => {
    _allBlastBtns().forEach(b =>
      b.classList.toggle('active', b.dataset.blastMetric === blastMetric));
    const cap = blastMetric === 'coverage' ? 'Coverage' : 'Identity';
    const xt = document.getElementById('stats-blastx-title');
    const nt = document.getElementById('stats-blastn-title');
    if (xt) xt.textContent = `BLASTx ${cap}`;
    if (nt) nt.textContent = `BLASTn ${cap}`;
  };
  const _drawBlast = () => {
    _renderIdentityHistogram(seqs, blastMetric);
    _renderBLAstnHistogram(seqs, blastMetric);
  };
  _allBlastBtns().forEach(btn => {
    btn.addEventListener('click', () => {
      blastMetric = btn.dataset.blastMetric;
      _syncBlastToggles();
      _drawBlast();
    });
  });
  _drawBlast();

  // Print / Save PDF — stats section only
  document.getElementById('stats-print-pdf')?.addEventListener('click', () => {
    document.body.classList.add('vq-print-stats');
    window.addEventListener('afterprint', () => {
      document.body.classList.remove('vq-print-stats');
    }, { once: true });
    VQ.exportPagePDF('ViralQuest – General Statistics');
  });

  // ResizeObserver — re-render charts when section is first revealed
  let lastW = 0;
  const ro = new ResizeObserver(() => {
    const w = el.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    renderTax();
    _renderHMMBars(hmm);
    if (salmon.present) _renderSalmonChart(salmon);
    _renderFunnel(blast);
    _drawBlast();
    _renderLengthHistogram(seqs);
    _renderNRClassification(seqs);
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

function _blobChip(label, value, sub = '') {
  /* Flower-blob KPI chip: 6 overlapping semi-transparent circles form the
     petal ring; a white disc sits on top carrying the number and label. */
  const petals = [
    [104, 85, 'var(--vq-accent)'],
    [95,  68, 'var(--vq-accent-dark)'],
    [75,  68, 'var(--vq-primary-light)'],
    [66,  85, 'var(--vq-accent)'],
    [75, 102, 'var(--vq-accent-dark)'],
    [95, 102, 'var(--vq-primary-light)'],
  ];
  const rings = petals.map(([px, py, col]) =>
    `<circle cx="${px}" cy="${py}" r="54" fill="${col}" opacity="0.28"/>`
  ).join('');
  const subText = sub
    ? `<text x="85" y="108" text-anchor="middle" font-size="8.5"
             font-family="var(--vq-font)" fill="var(--vq-text-3)">${VQ.esc(sub)}</text>`
    : '';
  return `
    <div class="vq-stat-chip vq-stat-chip--blob">
      <svg viewBox="0 0 170 170" width="170" height="170"
           style="overflow:visible;display:block;margin:auto">
        ${rings}
        <circle cx="85" cy="85" r="47" fill="var(--vq-surface)"/>
        <text x="85" y="79" text-anchor="middle"
              font-size="26" font-weight="700" font-family="var(--vq-font)"
              fill="var(--vq-primary)">${VQ.esc(value)}</text>
        <text x="85" y="95" text-anchor="middle"
              font-size="8.5" font-weight="600" letter-spacing=".08em"
              font-family="var(--vq-font)" fill="var(--vq-text-3)"
        >${VQ.esc(label.toUpperCase())}</text>
        ${subText}
      </svg>
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
    ].join('');
  }
}

function _renderSalmonEmpty() {
  const wrap = document.getElementById('stats-salmon-svg');
  if (wrap) wrap.innerHTML =
    '<div class="vq-empty" style="padding:40px 16px;font-size:12px">No Salmon quantification data.</div>';
}

// ── Rounded-top bar path (flat bottom to sit on x-axis) ─────────────────
function _barPath(x, y, w, h, r) {
  if (h <= 0 || w <= 0) return '';
  r = Math.min(r, w / 2, h);
  return `M${x},${y+h} L${x},${y+r} Q${x},${y} ${x+r},${y} L${x+w-r},${y} Q${x+w},${y} ${x+w},${y+r} L${x+w},${y+h} Z`;
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
    { label: 'NR confirmed',       value: blast.total_confirmed       ?? 0,   color: 'var(--vq-success)'       },
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
//  Chart 5 — BLASTx identity histogram (green card)
// ────────────────────────────────────────────────────────────────────────

function _renderIdentityHistogram(sequences, mode = 'identity') {
  const wrap = document.getElementById('stats-identity-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const isCoverage = mode === 'coverage';
  const values = sequences
    .filter(s => s.is_viral)
    .map(s => {
      const hit = (s.blastx_nr_hits && s.blastx_nr_hits[0]) || (s.blastx_hits && s.blastx_hits[0]);
      if (!hit) return null;
      return isCoverage ? (hit.query_coverage ?? null) : hit.pct_identity;
    })
    .filter(v => v != null);

  if (!values.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px 16px;font-size:12px">No BLASTx data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 160);
  const PAD_L = 36; const PAD_R = 6;
  const PAD_T = 6;  const PAD_B = 20;
  const H     = 110;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  const bins   = d3.bin().domain([0, 100]).thresholds([10,20,30,40,50,60,70,80,90])(values);
  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([0, 100], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  const metricLabel = isCoverage ? 'coverage' : 'identity';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + tick labels
  yScale.ticks(4).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.7);
    g.append('text')
      .attr('x', -5).attr('y', yScale(tick) + 3)
      .attr('text-anchor', 'end').attr('font-size', 9)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  // Bars — single accent colour, rounded top / flat bottom
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    const r = Math.min(Math.ceil(bw / 2), 5);
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, r))
      .attr('fill', 'var(--vq-accent)').attr('opacity', 0.82)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${bin.x0}–${bin.x1}% ${metricLabel}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / values.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks (keep x-axis)
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
  [0, 50, 100].forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 3)
      .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 12)
      .attr('text-anchor', 'middle').attr('font-size', 9).attr('fill', 'var(--vq-c-tick)')
      .text(tick + '%');
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
    wrap.innerHTML = '<div style="padding:24px 16px;font-size:12px;color:rgba(255,255,255,.45);text-align:center">No sequence data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 200);
  const PAD_L = 40; const PAD_R = 10;
  const PAD_T = 8;  const PAD_B = 24;
  const H     = 130;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  // 500 bp fixed-width bins
  const maxL       = d3.max(lengths);
  const STEP       = 500;
  const maxBin     = Math.ceil(maxL / STEP) * STEP;
  const thresholds = d3.range(0, maxBin, STEP);
  const bins       = d3.bin().domain([0, maxBin]).thresholds(thresholds)(lengths);

  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([0, maxBin], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  function _fmtLen(v) {
    return v >= 1000 ? (v / 1000).toFixed(v % 1000 === 0 ? 0 : 1) + 'k' : String(Math.round(v));
  }

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // Deep blue background (uses --vq-c-dark-bg = var(--vq-primary))
  svg.append('rect').attr('width', W).attr('height', H).attr('fill', 'var(--vq-c-dark-bg)');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis: horizontal grid lines + tick labels
  yScale.ticks(5).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-dark-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.8);
    g.append('text')
      .attr('x', -8).attr('y', yScale(tick) + 3.5)
      .attr('text-anchor', 'end').attr('font-size', 10)
      .attr('fill', 'var(--vq-c-dark-tick)').text(tick);
  });

  // White bars — rounded top, flat bottom
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 1, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    const r  = Math.min(Math.ceil(bw / 2), 5);
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, r))
      .attr('fill', 'var(--vq-c-dark-bar)')
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
    .attr('stroke', 'var(--vq-c-dark-base)').attr('stroke-width', 1);
  xScale.ticks(Math.min(6, bins.length)).forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 4)
      .attr('stroke', 'var(--vq-c-dark-base)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 14)
      .attr('text-anchor', 'middle').attr('font-size', 10)
      .attr('fill', 'var(--vq-c-dark-tick)').text(_fmtLen(tick));
  });

  // X axis label
  svg.append('text')
    .attr('x', PAD_L + drawW / 2).attr('y', H - 1)
    .attr('text-anchor', 'middle').attr('font-size', 10)
    .attr('fill', 'rgba(255,255,255,.35)').text('bp');

  wrap.appendChild(svg.node());
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 6b — BLASTn identity / coverage histogram (compact, white)
// ────────────────────────────────────────────────────────────────────────

function _renderBLAstnHistogram(sequences, mode) {
  const wrap = document.getElementById('stats-blastn-svg');
  if (!wrap) return;
  wrap.innerHTML = '';

  const isIdentity = mode !== 'coverage';
  const values = sequences
    .filter(s => s.is_viral && s.blastn_hits && s.blastn_hits.length > 0)
    .map(s => {
      const h = s.blastn_hits[0];
      return isIdentity ? (h.pident ?? null) : (h.qcovhsp ?? null);
    })
    .filter(v => v != null);

  if (!values.length) {
    wrap.innerHTML = '<div class="vq-empty" style="padding:24px;font-size:12px">No BLASTn data.</div>';
    return;
  }

  const W     = Math.max(wrap.clientWidth || 0, 200);
  const PAD_L = 40; const PAD_R = 10;
  const PAD_T = 8;  const PAD_B = 28;
  const H     = 130;
  const drawW = W - PAD_L - PAD_R;
  const drawH = H - PAD_T - PAD_B;

  const bins   = d3.bin().domain([0, 100]).thresholds([10,20,30,40,50,60,70,80,90])(values);
  const yMax   = d3.max(bins, b => b.length) || 1;
  const xScale = d3.scaleLinear([0, 100], [0, drawW]);
  const yScale = d3.scaleLinear([0, yMax], [drawH, 0]).nice();

  const barColor = 'var(--vq-primary-light)';
  const metricLabel = isIdentity ? 'identity' : 'coverage';

  const svg = d3.create('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L},${PAD_T})`);

  // Y-axis grid + tick labels
  yScale.ticks(4).forEach(tick => {
    g.append('line')
      .attr('x1', 0).attr('y1', yScale(tick))
      .attr('x2', drawW).attr('y2', yScale(tick))
      .attr('stroke', 'var(--vq-c-grid)').attr('stroke-dasharray', '3,2').attr('stroke-width', 0.8);
    g.append('text')
      .attr('x', -8).attr('y', yScale(tick) + 3.5)
      .attr('text-anchor', 'end').attr('font-size', 9.5)
      .attr('fill', 'var(--vq-c-tick)').text(tick);
  });

  // Bars
  bins.forEach(bin => {
    const x  = xScale(bin.x0);
    const bw = Math.max(xScale(bin.x1) - xScale(bin.x0) - 2, 1);
    const bh = drawH - yScale(bin.length);
    if (bh <= 0) return;
    g.append('path')
      .attr('d', _barPath(x, yScale(bin.length), bw, bh, Math.min(Math.ceil(bw / 2), 5)))
      .attr('fill', barColor).attr('opacity', 0.82)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${bin.x0}–${bin.x1}% ${metricLabel}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Sequences</span><span>${bin.length}</span>
          <span class="vq-tooltip__key">Share</span>
          <span>${(bin.length / values.length * 100).toFixed(1)}%</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
  });

  // Baseline + x ticks
  g.append('line').attr('x1', 0).attr('y1', drawH).attr('x2', drawW).attr('y2', drawH)
    .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
  [0, 25, 50, 75, 100].forEach(tick => {
    const tx = xScale(tick);
    g.append('line').attr('x1', tx).attr('y1', drawH).attr('x2', tx).attr('y2', drawH + 4)
      .attr('stroke', 'var(--vq-c-baseline)').attr('stroke-width', 1);
    g.append('text').attr('x', tx).attr('y', drawH + 14)
      .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('fill', 'var(--vq-c-tick)')
      .text(tick + '%');
  });

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

  const TOP_N  = 8;
  const entries = Object.entries(counts)
    .sort((a, b) => b[1].count - a[1].count)
    .slice(0, TOP_N);

  // Hide card when no species has more than one sequence — data is trivial
  const multiSeqSpecies = Object.values(counts).filter(d => d.count > 1).length;
  const card = document.getElementById('stats-species-card');
  if (multiSeqSpecies <= 1) {
    if (card) card.style.display = 'none';
    return;
  }
  if (card) card.style.display = '';

  if (sub) {
    const uniq = Object.keys(counts).length;
    sub.textContent = `top ${Math.min(TOP_N, uniq)} of ${uniq} species`;
  }

  if (!entries.length) {
    tableWrap.innerHTML = '<div class="vq-empty" style="padding:20px;font-size:12px">No taxonomy data.</div>';
    return;
  }

  const total = viral.length;
  tableWrap.innerHTML = `
    <table class="vq-table" style="font-size:11px">
      <thead>
        <tr>
          <th style="width:24px;padding:4px 8px">#</th>
          <th style="padding:4px 8px">Species</th>
          <th style="padding:4px 8px">Family</th>
          <th style="width:60px;padding:4px 8px">n</th>
          <th style="width:110px;padding:4px 8px">%</th>
        </tr>
      </thead>
      <tbody>
        ${entries.map(([sp, d], i) => {
          const pct = (d.count / total * 100).toFixed(1);
          return `
          <tr>
            <td style="color:var(--vq-text-3);padding:4px 8px">${i + 1}</td>
            <td style="font-style:italic;padding:4px 8px">${VQ.esc(sp)}</td>
            <td style="color:var(--vq-text-2);padding:4px 8px">${VQ.esc(d.family)}</td>
            <td style="font-weight:600;padding:4px 8px;font-variant-numeric:tabular-nums">${d.count}</td>
            <td style="padding:4px 8px">
              <div style="display:flex;align-items:center;gap:5px">
                <div style="flex:1;height:5px;background:var(--vq-bg);border-radius:3px">
                  <div style="height:100%;width:${pct}%;background:var(--vq-accent);border-radius:3px;opacity:.8"></div>
                </div>
                <span style="font-size:10px;color:var(--vq-text-3);white-space:nowrap">${pct}%</span>
              </div>
            </td>
          </tr>`;
        }).join('')}
      </tbody>
    </table>
  `;
}


// ────────────────────────────────────────────────────────────────────────
//  Chart 8 — NR classification: Virus vs Phage (bars only)
// ────────────────────────────────────────────────────────────────────────

const _NR_CLASS_COLORS = {
  'Virus': 'var(--vq-accent)',
  'Phage': 'var(--vq-accent-dark)',
};

function _renderNRClassification(sequences) {
  const lgd   = document.getElementById('stats-nrclass-legend');
  const sub   = document.getElementById('stats-nrclass-sub');
  const title = document.getElementById('stats-nrclass-title');
  if (!lgd) return;
  lgd.innerHTML = '';

  const viral  = sequences.filter(s => s.is_viral);
  const nrUsed = viral.some(s => s.blastx_nr_hits && s.blastx_nr_hits.length > 0);
  if (title) title.textContent = nrUsed ? 'NR Classification' : 'BLASTx Classification';

  let virus = 0, phage = 0;
  viral.forEach(s => {
    const hit = (s.blastx_nr_hits && s.blastx_nr_hits[0]) || (s.blastx_hits && s.blastx_hits[0]);
    if (!hit) return;
    if (/phage|bacteriophage/i.test(hit.subject_title)) phage++;
    else virus++;
  });

  const total   = virus + phage;
  const display = [['Virus', virus], ['Phage', phage]].filter(([, n]) => n > 0);

  if (sub) sub.textContent = `${total} classified sequence${total !== 1 ? 's' : ''}`;

  if (!total) {
    lgd.innerHTML = '<div class="vq-empty" style="padding:20px">No BLASTx hits.</div>';
    return;
  }

  const color = name => _NR_CLASS_COLORS[name] || _TAX_DONUT_PALETTE[0];
  display.forEach(([name, n]) => {
    const pct = (n / total * 100).toFixed(1);
    const col = color(name);
    const row = document.createElement('div');
    row.style.cssText = 'display:flex;align-items:center;gap:10px;padding:6px 0;border-bottom:1px dashed var(--vq-border)';
    row.innerHTML = `
      <span style="width:14px;height:14px;border-radius:50%;background:${col};flex-shrink:0"></span>
      <div style="flex:1;min-width:0">
        <div style="font-size:15px;font-weight:600;color:var(--vq-text)">${VQ.esc(name)}</div>
        <div style="height:4px;background:var(--vq-bg);border-radius:2px;margin-top:5px">
          <div style="height:100%;width:${pct}%;background:${col};border-radius:2px;opacity:.85"></div>
        </div>
      </div>
      <div style="text-align:right;flex-shrink:0">
        <div style="font-size:15px;font-weight:700;font-variant-numeric:tabular-nums;color:var(--vq-text)">${n}</div>
        <div style="font-size:12px;color:var(--vq-text-3)">${pct}%</div>
      </div>
    `;
    lgd.appendChild(row);
  });
}


window.vqInitStats = vqInitStats;
})();
