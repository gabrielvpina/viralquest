/* ============================================================
   section_clusters.js — Section 2: Cluster Viruses
   Linear alignment view: members drawn relative to representative.
   Colour encodes % identity. Click member → jump to viewer.
   ============================================================ */

'use strict';

function vqInitClusters(clusters) {
  const el = document.getElementById('section-clusters');
  if (!el) return;
  if (!clusters.length) {
    el.innerHTML = '<div class="vq-empty">No clusters found.</div>';
    return;
  }

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Cluster Viruses</div>
        <div class="vq-section-sub">
          ${clusters.length} cluster${clusters.length > 1 ? 's' : ''} · members shown as alignments to representative sequence
        </div>
      </div>
      <div style="display:flex;gap:var(--vq-space-2)">
        ${_identityLegend()}
      </div>
    </div>
    <div id="clusters-list"></div>
  `;

  const list = document.getElementById('clusters-list');
  clusters.forEach(c => list.appendChild(_renderCluster(c)));
}

// ── Render one cluster card ──────────────────────────────────────────────────

function _renderCluster(cluster) {
  const card = document.createElement('div');
  card.className = 'vq-card';
  card.style.marginBottom = 'var(--vq-space-4)';

  const rep    = cluster.members.find(m => m.is_representative) || cluster.members[0];
  const repLen = rep?.length || 1;

  card.innerHTML = `
    <div class="vq-card__header">
      <div class="vq-card__title">
        <span class="vq-badge vq-badge--cluster">${cluster.cluster_id}</span>
        ${cluster.species}
      </div>
      <div style="display:flex;gap:var(--vq-space-2);align-items:center">
        <span style="font-size:var(--vq-text-xs);color:var(--vq-text-3)">
          ${cluster.size} member${cluster.size > 1 ? 's' : ''}
        </span>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportPNG(this.closest('.vq-card').querySelector('svg'),'cluster_${cluster.cluster_id}.png')">
          PNG
        </button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportSVG(this.closest('.vq-card').querySelector('svg'),'cluster_${cluster.cluster_id}.svg')">
          SVG
        </button>
      </div>
    </div>
    <div class="vq-card__body">
      <div class="vq-genome-wrap" id="clu-wrap-${cluster.cluster_id}"></div>
    </div>`;

  // Render SVG after insert
  requestAnimationFrame(() => {
    const wrap = document.getElementById('clu-wrap-' + cluster.cluster_id);
    if (wrap) wrap.appendChild(_clusterSVG(cluster, repLen));
  });

  return card;
}

// ── Build the cluster SVG ────────────────────────────────────────────────────

function _clusterSVG(cluster, repLen) {
  const PAD_L  = 160;   // left label area
  const PAD_R  = 20;
  const ROW_H  = 28;
  const GAP    = 6;
  const TRACK  = 14;    // bar height

  const totalH = cluster.members.length * (ROW_H + GAP) + 40;
  const W      = 900;
  const drawW  = W - PAD_L - PAD_R;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('width',   W)
    .attr('height',  totalH);

  const xScale = d3.scaleLinear([0, repLen], [0, drawW]);

  // --- axis ---
  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, 20)`)
    .call(d3.axisTop(xScale).ticks(8).tickSize(-totalH + 36));

  axisG.select('.domain').remove();
  axisG.selectAll('.tick line')
    .attr('stroke', '#dde3ec')
    .attr('stroke-dasharray', '3,3');
  axisG.selectAll('.tick text')
    .style('font-size', '10px')
    .style('fill', '#94a3b8');

  // --- members ---
  cluster.members.forEach((m, i) => {
    const y   = 30 + i * (ROW_H + GAP);
    const col = _identityColor(m.identity);
    const g   = svg.append('g').attr('transform', `translate(${PAD_L}, ${y})`);

    // label
    svg.append('text')
      .attr('x', PAD_L - 8)
      .attr('y', y + TRACK / 2 + 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 11)
      .attr('fill', m.is_representative ? 'var(--vq-primary)' : 'var(--vq-text-2)')
      .attr('font-weight', m.is_representative ? '600' : '400')
      .attr('font-family', 'var(--vq-font-mono)')
      .text(m.seq_id);

    // backbone (full sequence length context)
    g.append('rect')
      .attr('x', 0)
      .attr('width', drawW)
      .attr('y', TRACK / 2 - 1)
      .attr('height', 2)
      .attr('rx', 1)
      .attr('fill', '#dde3ec');

    // alignment bar
    const x1 = xScale(m.aln_start - 1);
    const x2 = xScale(m.aln_end);
    const bw  = Math.max(x2 - x1, 2);

    const bar = g.append('rect')
      .attr('x', x1)
      .attr('width', bw)
      .attr('y', 0)
      .attr('height', TRACK)
      .attr('rx', 3)
      .attr('fill', col)
      .attr('opacity', m.is_representative ? 1 : 0.85)
      .attr('cursor', 'pointer')
      .on('mousemove', evt => {
        vqTooltipShow(`
          <div class="vq-tooltip__title">${m.seq_id}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Identity</span><span>${m.identity.toFixed(1)}%</span>
            <span class="vq-tooltip__key">Coverage</span><span>${m.query_coverage.toFixed(1)}%</span>
            <span class="vq-tooltip__key">Length</span><span>${m.length.toLocaleString()} nt</span>
            <span class="vq-tooltip__key">Aln range</span><span>${m.aln_start}–${m.aln_end}</span>
          </div>`, evt);
      })
      .on('mouseleave', vqTooltipHide)
      .on('click', () => _jumpToViewer(m.seq_id));

    // identity label inside bar if wide enough
    if (bw > 50 && !m.is_representative) {
      g.append('text')
        .attr('x', x1 + bw / 2)
        .attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle')
        .attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.9)')
        .attr('pointer-events', 'none')
        .text(m.identity.toFixed(0) + '%');
    }

    if (m.is_representative) {
      g.append('text')
        .attr('x', x1 + bw / 2)
        .attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle')
        .attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.9)')
        .attr('pointer-events', 'none')
        .text('representative · ' + m.length.toLocaleString() + ' nt');
    }
  });

  return svg.node();
}

// ── Jump to viewer section with sequence highlighted ─────────────────────────

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

// ── Colour scale: red (low) → yellow → green (high identity) ─────────────────

function _identityColor(pct) {
  const t = Math.max(0, Math.min(1, pct / 100));
  return d3.interpolateRdYlGn(t);
}

// ── Legend ───────────────────────────────────────────────────────────────────

function _identityLegend() {
  const stops = [0, 25, 50, 75, 100];
  return `<div class="vq-legend">
    ${stops.map(v =>
      `<div class="vq-legend__item">
         <div class="vq-legend__swatch"
              style="background:${_identityColor(v)};width:18px;height:12px;border-radius:2px"></div>
         <span>${v}% identity</span>
       </div>`
    ).join('')}
  </div>`;
}
