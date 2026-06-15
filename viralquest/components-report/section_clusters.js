/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_clusters.js — Section 2: General (cross-sample) Clusters
   Top: bipartite network — each cluster (VQR_CLU_####) is a hub,
        the samples that contain it link to it.
   Bottom: per-cluster alignment view (members vs representative),
        plus a species-homogeneity annotation.
   Dynamic identity / coverage sliders filter everything client-side.
   ============================================================ */

const FLOOR_ID  = 90;   // matches report_clusters.FLOOR_IDENTITY
const FLOOR_COV = 70;   // matches report_clusters.FLOOR_COVERAGE

let _clusters = [];
let _state    = { id: FLOOR_ID, cov: FLOOR_COV };

function vqInitClusters(clusters, samples) {
  const el = document.getElementById('section-clusters');
  if (!el) return;
  _clusters = clusters || [];

  if (!_clusters.length) {
    el.hidden = true;
    document.getElementById('tab-clusters')?.style.setProperty('display', 'none');
    return;
  }

  const esc = VQ.esc;
  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">General Clusters</div>
        <div class="vq-section-sub">
          ${_clusters.length} cross-sample cluster${_clusters.length > 1 ? 's' : ''}
          · sequences shared between samples
        </div>
      </div>
      <div class="vq-section-actions">${_identityLegend()}</div>
    </div>

    <div class="vq-card" style="margin-bottom:var(--vq-space-4)">
      <div class="vq-card__header">
        <div class="vq-card__title">Filters</div>
        <div class="ov-filter-readout" id="clu-readout"></div>
      </div>
      <div class="vq-card__body">
        <div class="ov-filter-row">
          <label>Min identity
            <input type="range" id="clu-id"  min="${FLOOR_ID}"  max="100" step="1" value="${_state.id}">
            <output id="clu-id-out">${_state.id}%</output>
          </label>
          <label>Min coverage
            <input type="range" id="clu-cov" min="${FLOOR_COV}" max="100" step="1" value="${_state.cov}">
            <output id="clu-cov-out">${_state.cov}%</output>
          </label>
        </div>
      </div>
    </div>

    <div class="vq-card" style="margin-bottom:var(--vq-space-4)">
      <div class="vq-card__header">
        <div class="vq-card__title">Cluster Network</div>
        <div class="vq-card__sub" id="clu-net-sub"></div>
      </div>
      <div class="vq-card__body"><div id="clu-graph"></div></div>
    </div>

    <div id="clu-cards"></div>
  `;

  const idIn  = document.getElementById('clu-id');
  const covIn = document.getElementById('clu-cov');
  const apply = () => {
    _state.id  = +idIn.value;
    _state.cov = +covIn.value;
    document.getElementById('clu-id-out').textContent  = _state.id + '%';
    document.getElementById('clu-cov-out').textContent = _state.cov + '%';
    _rerender();
  };
  idIn.addEventListener('input', apply);
  covIn.addEventListener('input', apply);

  _rerender();
}

// ── Filtering ───────────────────────────────────────────────────────────────

function _passing(cluster) {
  // Representative always passes (100/100). Returns members meeting thresholds.
  return cluster.members.filter(m =>
    m.is_representative || (m.identity >= _state.id && m.coverage >= _state.cov));
}

function _activeClusters() {
  // A cluster stays visible only while its passing members span ≥2 samples.
  return _clusters
    .map(c => ({ cluster: c, members: _passing(c) }))
    .filter(({ members }) => new Set(members.map(m => m.sample)).size >= 2);
}

function _rerender() {
  const active = _activeClusters();
  const readout = document.getElementById('clu-readout');
  if (readout)
    readout.textContent = `${active.length} / ${_clusters.length} cluster(s) shown`;
  _renderGraph(document.getElementById('clu-graph'), active);
  _renderCards(document.getElementById('clu-cards'), active);
}

// ── Bipartite network graph ─────────────────────────────────────────────────

function _renderGraph(host, active) {
  if (!host) return;
  host.innerHTML = '';

  if (!active.length) {
    host.innerHTML = `<div class="vq-empty">No clusters pass the current filters.</div>`;
    document.getElementById('clu-net-sub').textContent = '';
    return;
  }

  // Build bipartite nodes/links.
  const nodes = [];
  const links = [];
  const sampleNode = {};
  active.forEach(({ cluster, members }) => {
    const cNode = { id: 'C:' + cluster.gid, type: 'cluster', label: cluster.gid,
                    n: new Set(members.map(m => m.sample)).size };
    nodes.push(cNode);
    [...new Set(members.map(m => m.sample))].forEach(s => {
      if (!sampleNode[s]) {
        sampleNode[s] = { id: 'S:' + s, type: 'sample', label: s };
        nodes.push(sampleNode[s]);
      }
      links.push({ source: cNode.id, target: sampleNode[s].id });
    });
  });

  document.getElementById('clu-net-sub').textContent =
    `${active.length} cluster${active.length > 1 ? 's' : ''} · ` +
    `${Object.keys(sampleNode).length} sample${Object.keys(sampleNode).length > 1 ? 's' : ''}`;

  const W = 720, H = Math.max(320, Math.min(640, 120 + nodes.length * 14));
  const svg = d3.select(host).append('svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMidYMid meet')
    .style('width', '100%').style('height', 'auto').style('display', 'block');

  const link = svg.append('g').attr('stroke', '#cbd5e1').attr('stroke-width', 1.2)
    .selectAll('line').data(links).join('line');

  const node = svg.append('g').selectAll('g').data(nodes).join('g')
    .style('cursor', 'pointer')
    .call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));

  node.each(function (d) {
    const g = d3.select(this);
    if (d.type === 'cluster') {
      g.append('circle').attr('r', 9)
        .attr('fill', 'var(--vq-accent)').attr('stroke', '#fff').attr('stroke-width', 2);
    } else {
      g.append('rect').attr('x', -7).attr('y', -7).attr('width', 14).attr('height', 14)
        .attr('rx', 3).attr('fill', 'var(--vq-primary)').attr('stroke', '#fff').attr('stroke-width', 2);
    }
    g.append('text').attr('x', 0).attr('y', d.type === 'cluster' ? -13 : 19)
      .attr('text-anchor', 'middle').attr('font-size', 10)
      .attr('font-family', d.type === 'cluster' ? 'var(--vq-font-mono)' : 'var(--vq-font)')
      .attr('fill', 'var(--vq-text-2)').text(d.label);
  });

  node.on('mousemove', (e, d) => VQ.tooltipShow(
      d.type === 'cluster'
        ? `<b>${VQ.esc(d.label)}</b><br>${d.n} samples`
        : `<b>${VQ.esc(d.label)}</b><br>sample`, e))
    .on('mouseleave', () => VQ.tooltipHide())
    .on('click', (e, d) => {
      if (d.type === 'cluster') {
        const card = document.getElementById('clu-card-' + VQ.safeId(d.label));
        card?.scrollIntoView({ behavior: 'smooth', block: 'center' });
        card?.classList.add('vq-seq-card--highlight');
        setTimeout(() => card?.classList.remove('vq-seq-card--highlight'), 1600);
      }
    });

  const sim = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(links).id(d => d.id).distance(70).strength(0.6))
    .force('charge', d3.forceManyBody().strength(-220))
    .force('center', d3.forceCenter(W / 2, H / 2))
    .force('collide', d3.forceCollide(24))
    .on('tick', () => {
      link.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
          .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
      node.attr('transform', d => `translate(${d.x},${d.y})`);
    });
}

// ── Cluster cards (alignment view + homogeneity) ────────────────────────────

function _renderCards(host, active) {
  if (!host) return;
  host.innerHTML = '';
  if (!active.length) return;
  active.forEach(({ cluster, members }) => host.appendChild(_renderCluster(cluster, members)));
}

function _renderCluster(cluster, members) {
  const esc  = VQ.esc;
  const safe = VQ.safeId(cluster.gid);
  const card = document.createElement('div');
  card.className = 'vq-card';
  card.id = 'clu-card-' + safe;
  card.style.marginBottom = 'var(--vq-space-4)';

  // Filtered view of the cluster (representative always present).
  const view = Object.assign({}, cluster, { members });
  const rep    = members.find(m => m.is_representative) || members[0];
  const repLen = rep?.length || 1;
  const samples = [...new Set(members.map(m => m.sample))];

  card.innerHTML = `
    <div class="vq-card__header">
      <div class="vq-card__title">
        <span class="vq-badge vq-badge--cluster">${esc(cluster.gid)}</span>
        <span class="clu-species">${esc(cluster.species)}</span>
      </div>
      <div class="vq-section-actions">
        <span style="font-size:var(--vq-text-xs);color:var(--vq-text-3)">
          ${members.length} member${members.length > 1 ? 's' : ''} · ${samples.length} samples
        </span>
        ${_exportMenu('clu-menu-' + safe)}
      </div>
    </div>
    <div class="vq-card__body">
      ${_homogeneityBanner(cluster, members)}
      <div class="vq-genome-wrap" id="clu-wrap-${safe}"></div>
    </div>`;

  const wrapEl = card.querySelector('#clu-wrap-' + safe);
  _wireExportMenu(card, 'clu-menu-' + safe,
    () => card.querySelector('.vq-genome-wrap svg'),
    `cluster_${cluster.gid}`, () => _clusterFasta(view));

  let lastW = 0, rzTimer = null;
  const render = () => {
    const w = wrapEl.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    wrapEl.innerHTML = '';
    wrapEl.appendChild(_clusterSVG(view, repLen, w));
  };
  new ResizeObserver(() => { clearTimeout(rzTimer); rzTimer = setTimeout(render, 60); }).observe(wrapEl);
  requestAnimationFrame(render);

  return card;
}

function _homogeneityBanner(cluster, members) {
  const esc = VQ.esc;
  const named = members.filter(m => m.species_db === 'nr' || m.species_db === 'refseq');
  const dbLabel = cluster.species_db === 'nr' ? 'BLASTx NR'
                : cluster.species_db === 'refseq' ? 'BLASTx RefSeq' : 'best hit';
  const divergent = named.filter(m => m.species && m.species !== cluster.species);

  if (!divergent.length) {
    return `<div class="clu-banner clu-banner--ok">
      ✓ All members align to the same ${esc(dbLabel)} species:
      <b>${esc(cluster.species)}</b></div>`;
  }
  const items = divergent.map(m =>
    `<li><span class="vq-badge">${esc(m.sample)}</span> ${esc(m.seq_id)} → <i>${esc(m.species)}</i></li>`
  ).join('');
  return `<div class="clu-banner clu-banner--warn">
      ⚠ ${divergent.length} member(s) hit a different species than the representative
      (<b>${esc(cluster.species)}</b>, ${esc(dbLabel)}):
      <ul class="clu-divergent">${items}</ul></div>`;
}

// ── Export dropdown (copied from the per-sample report) ─────────────────────

function _exportMenu(id) {
  return `
    <div class="vq-menu" id="${id}">
      <button class="vq-btn vq-btn--sm vq-btn--ghost" data-menu-toggle type="button">
        Export
        <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
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
  const menu = root.querySelector('#' + menuId);
  if (!menu) return;
  const toggle = menu.querySelector('[data-menu-toggle]');
  toggle?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => { if (m !== menu) m.classList.remove('open'); });
    menu.classList.toggle('open');
  });
  menu.querySelectorAll('[data-fmt]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      menu.classList.remove('open');
      const fmt = item.dataset.fmt;
      if (fmt === 'fasta') { if (getFasta) VQ.downloadText(getFasta(), baseName + '.fasta'); return; }
      const svg = getSvg();
      if (!svg) return;
      if (fmt === 'png') VQ.exportPNG(svg, baseName + '.png');
      if (fmt === 'svg') VQ.exportSVG(svg, baseName + '.svg');
      if (fmt === 'pdf') VQ.exportPDF(svg, baseName);
    });
  });
}

function _clusterFasta(cluster) {
  // Cross-sample: match by gid (raw seq_id collides across samples).
  const allSeqs = (typeof VQ_REPORT !== 'undefined' ? VQ_REPORT.sequences : null) || [];
  const ids = new Set(cluster.members.map(m => m.gid));
  const byGid = {};
  allSeqs.forEach(s => { byGid[s.gid] = s; });
  const seqs = cluster.members.map(m => byGid[m.gid]).filter(Boolean);
  return VQ.buildFasta(seqs);
}

document.addEventListener('click', () => {
  document.querySelectorAll('.vq-menu.open').forEach(m => m.classList.remove('open'));
});

// ── Alignment-bar SVG (adapted from the per-sample report) ──────────────────

function _clusterSVG(cluster, repLen, containerWidth) {
  const PAD_L = 170, PAD_R = 20, ROW_H = 28, GAP = 6, TRACK = 14;
  const W = Math.max(containerWidth || 0, 720);
  const drawW = W - PAD_L - PAD_R;

  const offsets = cluster.members.map(m => {
    const qStart = m.is_representative ? 1 : (m.query_start ?? 1);
    return (m.aln_start - 1) - (qStart - 1);
  });
  const globalMinNt = Math.min(0, ...offsets);
  const globalMaxNt = Math.max(repLen, ...cluster.members.map((m, i) => offsets[i] + m.length));
  const spanNt = globalMaxNt - globalMinNt;
  const shiftNt = -globalMinNt;
  const xScale = d3.scaleLinear([0, spanNt], [0, drawW]);
  const totalH = cluster.members.length * (ROW_H + GAP) + 40;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto').attr('height', totalH);

  const axisG = svg.append('g')
    .attr('transform', `translate(${PAD_L}, 20)`)
    .call(d3.axisTop(xScale).ticks(8).tickSize(-(totalH - 30)));
  axisG.select('.domain').remove();
  axisG.selectAll('.tick line').attr('stroke', '#dde3ec').attr('stroke-dasharray', '3,3');
  axisG.selectAll('.tick text').style('font-size', '10px').style('fill', '#94a3b8');

  cluster.members.forEach((m, i) => {
    const offsetNt = offsets[i];
    const y = 30 + i * (ROW_H + GAP);
    const col = _identityColor(m.identity);
    const g = svg.append('g').attr('transform', `translate(${PAD_L}, ${y})`);

    // Label = sample::seq_id (cross-sample needs the sample to disambiguate).
    const labelText = `${m.sample} · ${m.seq_id}`;
    const truncated = labelText.length > 22 ? labelText.slice(0, 21) + '…' : labelText;
    const label = g.append('text')
      .attr('x', -8).attr('y', TRACK / 2 + 4)
      .attr('text-anchor', 'end').attr('font-size', 11)
      .attr('fill', m.is_representative ? 'var(--vq-primary)' : 'var(--vq-text-2)')
      .attr('font-weight', m.is_representative ? '600' : '400')
      .attr('font-family', 'var(--vq-font-mono)').text(truncated);
    if (truncated !== labelText) {
      label.style('cursor', 'help')
        .on('mousemove', evt => VQ.tooltipShow(`<div class="vq-tooltip__title">${VQ.esc(labelText)}</div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    }

    g.append('rect').attr('x', 0).attr('y', TRACK / 2 - 1)
      .attr('width', drawW).attr('height', 2).attr('rx', 1).attr('fill', '#dde3ec');

    const x1 = xScale((m.aln_start - 1) + shiftNt);
    const x2 = xScale(m.aln_end + shiftNt);
    const bw = Math.max(x2 - x1, 2);

    g.append('rect').attr('x', x1).attr('width', bw).attr('y', 0).attr('height', TRACK)
      .attr('rx', 3).attr('fill', col).attr('opacity', m.is_representative ? 1 : 0.85)
      .attr('cursor', 'pointer')
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${VQ.esc(m.sample)} · ${VQ.esc(m.seq_id)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Identity</span><span>${m.identity.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Coverage</span><span>${m.coverage.toFixed(1)}%</span>
          <span class="vq-tooltip__key">Length</span><span>${m.length.toLocaleString()} nt</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide)
      .on('click', () => VQ.jumpToViewer(m.gid));

    if (!m.is_representative) {
      const barX = xScale(offsetNt + shiftNt);
      const barW = xScale(offsetNt + m.length + shiftNt) - barX;
      if (barW > 1) {
        g.append('rect').attr('x', barX).attr('y', 0).attr('width', barW).attr('height', TRACK)
          .attr('rx', 3).attr('fill', 'none').attr('stroke', 'var(--vq-primary)')
          .attr('stroke-width', 1.5).attr('opacity', 0.5).attr('pointer-events', 'none');
      }
    }

    if (m.is_representative) {
      g.append('text').attr('x', x1 + bw / 2).attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle').attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.95)').attr('pointer-events', 'none')
        .text('representative · ' + m.length.toLocaleString() + ' nt');
    } else if (bw > 50) {
      g.append('text').attr('x', x1 + bw / 2).attr('y', TRACK / 2 + 4)
        .attr('text-anchor', 'middle').attr('font-size', 10)
        .attr('fill', 'rgba(255,255,255,.9)').attr('pointer-events', 'none')
        .text(m.identity.toFixed(0) + '%');
    }
  });

  return svg.node();
}

function _identityColor(pct) {
  const t = Math.max(0, Math.min(1, pct / 100));
  return d3.interpolateRdYlGn(t);
}

function _identityLegend() {
  return `
    <div class="vq-gradient-legend" title="Alignment identity (%)">
      <span>% identity</span>
      <div>
        <div class="vq-gradient-legend__bar"></div>
        <div class="vq-gradient-legend__ticks"><span>0</span><span>50</span><span>100</span></div>
      </div>
    </div>`;
}

window.vqInitClusters = vqInitClusters;
})();
