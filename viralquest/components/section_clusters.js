/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_clusters.js — Section 2: Cluster Analysis
   Linear alignment view: members drawn relative to representative.
   Colour encodes % identity. Click member → jump to viewer.
   ============================================================ */

function vqInitClusters(clusters) {
  const el = document.getElementById('section-clusters');
  if (!el) return;
  if (!clusters.length) {
    el.hidden = true;
    document.getElementById('tab-clusters')?.style.setProperty('display', 'none');
    return;
  }

  const esc = VQ.esc;

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Cluster Analysis</div>
        <div class="vq-section-sub">
          ${clusters.length} cluster${clusters.length > 1 ? 's' : ''} ·
          members shown as alignments to representative sequence
        </div>
      </div>
      <div class="vq-section-actions">
        ${_identityLegend()}
      </div>
    </div>
    <div id="clusters-list"></div>
  `;

  const list = document.getElementById('clusters-list');
  clusters.forEach(c => list.appendChild(_renderCluster(c)));
}

// ── Render one cluster card ────────────────────────────────────────────────

function _renderCluster(cluster) {
  const esc   = VQ.esc;
  const safe  = VQ.safeId(cluster.cluster_id);
  const card  = document.createElement('div');
  card.className = 'vq-card';
  card.style.marginBottom = 'var(--vq-space-4)';

  const rep    = cluster.members.find(m => m.is_representative) || cluster.members[0];
  const repLen = rep?.length || 1;

  card.innerHTML = `
    <div class="vq-card__header">
      <div class="vq-card__title">
        <span class="vq-badge vq-badge--cluster">${esc(cluster.cluster_id)}</span>
        ${esc(cluster.species)}
      </div>
      <div class="vq-section-actions">
        <span style="font-size:var(--vq-text-xs);color:var(--vq-text-3)">
          ${cluster.size} member${cluster.size > 1 ? 's' : ''}
        </span>
        ${_exportMenu('clu-menu-' + safe)}
      </div>
    </div>
    <div class="vq-card__body">
      <div class="vq-genome-wrap" id="clu-wrap-${safe}"></div>
    </div>`;

  // Render SVG once the wrap actually has a real width (it's 0 while the
  // section is still display:none on init), and re-render on any width change.
  let lastW = 0;
  let rzTimer = null;
  const wrapEl = card.querySelector('#clu-wrap-' + safe);
  _wireExportMenu(
    card, 'clu-menu-' + safe,
    () => card.querySelector('.vq-genome-wrap svg'),
    `cluster_${cluster.cluster_id}`,
    () => _clusterFasta(cluster),
  );

  const render = () => {
    const w = wrapEl.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    wrapEl.innerHTML = '';
    wrapEl.appendChild(_clusterSVG(cluster, repLen, w));
  };

  const ro = new ResizeObserver(() => {
    clearTimeout(rzTimer);
    rzTimer = setTimeout(render, 60);
  });
  ro.observe(wrapEl);
  // Best-effort first paint (no-op when hidden; ResizeObserver covers reveal)
  requestAnimationFrame(render);

  return card;
}

// ── Reusable export dropdown ────────────────────────────────────────────────

function _exportMenu(id) {
  return `
    <div class="vq-menu" id="${id}">
      <button class="vq-btn vq-btn--sm vq-btn--ghost" data-menu-toggle type="button">
        Export
        <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2.5">
          <polyline points="6 9 12 15 18 9"/>
        </svg>
      </button>
      <div class="vq-menu__panel" role="menu">
        <button class="vq-menu__item" data-fmt="png"   role="menuitem" type="button">PNG image</button>
        <button class="vq-menu__item" data-fmt="svg"   role="menuitem" type="button">SVG vector</button>
        <button class="vq-menu__item" data-fmt="pdf"   role="menuitem" type="button">PDF (print)</button>
        <button class="vq-menu__item" data-fmt="fasta" role="menuitem" type="button">FASTA sequences</button>
      </div>
    </div>`;
}

function _wireExportMenu(root, menuId, getSvg, baseName, getFasta) {
  const menu   = root.querySelector('#' + menuId);
  if (!menu) return;
  const toggle = menu.querySelector('[data-menu-toggle]');
  toggle?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => {
      if (m !== menu) m.classList.remove('open');
    });
    menu.classList.toggle('open');
  });
  menu.querySelectorAll('[data-fmt]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      menu.classList.remove('open');
      const fmt = item.dataset.fmt;

      if (fmt === 'fasta') {
        if (getFasta) VQ.downloadText(getFasta(), baseName + '.fasta');
        return;
      }

      const svg = getSvg();
      if (!svg) return;
      if (fmt === 'png') VQ.exportPNG(svg, baseName + '.png');
      if (fmt === 'svg') VQ.exportSVG(svg, baseName + '.svg');
      if (fmt === 'pdf') VQ.exportPDF(svg, baseName);
    });
  });
}

function _clusterFasta(cluster) {
  const allSeqs  = (typeof VQ_REPORT !== 'undefined' ? VQ_REPORT.sequences : null) || [];
  const memberIds = new Set(cluster.members.map(m => m.seq_id));
  const seqs = allSeqs.filter(s => memberIds.has(s.id));
  return VQ.buildFasta(seqs);
}

// Global click closes any open menu
document.addEventListener('click', () => {
  document.querySelectorAll('.vq-menu.open').forEach(m => m.classList.remove('open'));
});

// ── Build the cluster SVG ───────────────────────────────────────────────────

function _clusterSVG(cluster, repLen, containerWidth) {
  const PAD_L = 170;
  const PAD_R = 20;
  const ROW_H = 28;
  const GAP   = 6;
  const TRACK = 14;

  const W     = Math.max(containerWidth || 0, 720);
  const drawW = W - PAD_L - PAD_R;

  // ── Global coordinate system ─────────────────────────────────────────────
  // offsetNt[i]: where member i's position 0 falls in representative nt-space
  // (0-based). Negative means the member extends to the LEFT of the representative.
  const offsets = cluster.members.map(m => {
    const qStart = m.is_representative ? 1 : (m.query_start ?? 1);
    return (m.aln_start - 1) - (qStart - 1);   // = aln_start - qStart
  });

  // Span covering every sequence's full extent
  const globalMinNt = Math.min(0, ...offsets);
  const globalMaxNt = Math.max(repLen, ...cluster.members.map((m, i) => offsets[i] + m.length));
  const spanNt      = globalMaxNt - globalMinNt;

  // shiftNt: how many nt to add to any representative-space coordinate so that
  // the leftmost sequence starts at pixel 0.
  const shiftNt = -globalMinNt;

  const xScale = d3.scaleLinear([0, spanNt], [0, drawW]);

  const totalH = cluster.members.length * (ROW_H + GAP) + 40;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%')
    .style('height', 'auto')
    .attr('height',  totalH);

  // ── Axis (labels in global nt-space, 0 = leftmost sequence start) ────────
  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, 20)`)
    .call(d3.axisTop(xScale).ticks(8).tickSize(-(totalH - 30)));

  axisG.select('.domain').remove();
  axisG.selectAll('.tick line')
    .attr('stroke', '#dde3ec').attr('stroke-dasharray', '3,3');
  axisG.selectAll('.tick text')
    .style('font-size', '10px').style('fill', '#94a3b8');

  // ── Rows ─────────────────────────────────────────────────────────────────
  cluster.members.forEach((m, i) => {
    const offsetNt = offsets[i];
    const y        = 30 + i * (ROW_H + GAP);
    const col      = _identityColor(m.identity);
    const g        = svg.append('g').attr('transform', `translate(${PAD_L}, ${y})`);

    // Label
    const labelText = m.seq_id;
    const truncated = labelText.length > 18 ? labelText.slice(0, 17) + '…' : labelText;
    const label = g.append('text')
      .attr('x', -8).attr('y', TRACK / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', 11)
      .attr('fill', m.is_representative ? 'var(--vq-primary)' : 'var(--vq-text-2)')
      .attr('font-weight', m.is_representative ? '600' : '400')
      .attr('font-family', 'var(--vq-font-mono)')
      .text(truncated);

    if (truncated !== labelText) {
      label.style('cursor', 'help')
        .on('mousemove', evt => VQ.tooltipShow(
          `<div class="vq-tooltip__title">${VQ.esc(labelText)}</div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    }

    // Full-span backbone (light gray, always the whole drawing width)
    g.append('rect')
      .attr('x', 0).attr('y', TRACK / 2 - 1)
      .attr('width', drawW).attr('height', 2)
      .attr('rx', 1).attr('fill', '#dde3ec');

    // Alignment bar: representative-space coordinates shifted into global space
    const x1 = xScale((m.aln_start - 1) + shiftNt);
    const x2 = xScale(m.aln_end       + shiftNt);
    const bw = Math.max(x2 - x1, 2);

    g.append('rect')
      .attr('x', x1).attr('width', bw)
      .attr('y', 0).attr('height', TRACK)
      .attr('rx', 3)
      .attr('fill', col)
      .attr('opacity', m.is_representative ? 1 : 0.85)
      .attr('cursor', 'pointer')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(m.seq_id)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Identity</span><span>${m.identity.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Coverage</span><span>${m.query_coverage.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Length</span><span>${m.length.toLocaleString()} nt</span>
          <span class="vq-tooltip__key">Aln range</span><span>${m.aln_start}–${m.aln_end}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide)
      .on('click', () => VQ.jumpToViewer(m.seq_id));

    // Member full-length outline drawn on top — no clamping needed because the
    // global coordinate shift guarantees barX >= 0 and barEnd <= drawW.
    if (!m.is_representative) {
      const barX = xScale(offsetNt + shiftNt);
      const barW = xScale(offsetNt + m.length + shiftNt) - barX;
      if (barW > 1) {
        g.append('rect')
          .attr('x', barX).attr('y', 0)
          .attr('width', barW).attr('height', TRACK)
          .attr('rx', 3)
          .attr('fill', 'none')
          .attr('stroke', 'var(--vq-primary)')
          .attr('stroke-width', 1.5)
          .attr('opacity', 0.5)
          .attr('pointer-events', 'none');
      }
    }

    // In-bar label
    if (m.is_representative) {
      g.append('text')
        .attr('x', x1 + bw / 2).attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle').attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.95)').attr('pointer-events', 'none')
        .text('representative · ' + m.length.toLocaleString() + ' nt');
    } else if (bw > 50) {
      g.append('text')
        .attr('x', x1 + bw / 2).attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle').attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.9)').attr('pointer-events', 'none')
        .text(m.identity.toFixed(0) + '%');
    }
  });

  return svg.node();
}

// ── Colour scale (RdYlGn) ──────────────────────────────────────────────────

function _identityColor(pct) {
  const t = Math.max(0, Math.min(1, pct / 100));
  return d3.interpolateRdYlGn(t);
}

// ── Continuous gradient legend ─────────────────────────────────────────────

function _identityLegend() {
  return `
    <div class="vq-gradient-legend" title="Alignment identity (%)">
      <span>% identity</span>
      <div>
        <div class="vq-gradient-legend__bar"></div>
        <div class="vq-gradient-legend__ticks">
          <span>0</span><span>50</span><span>100</span>
        </div>
      </div>
    </div>`;
}


window.vqInitClusters = vqInitClusters;
})();
