/* ============================================================
   section_taxonomy.js — Section 4: Taxonomy Radial Tree
   D3 cluster layout; nodes coloured by viral family.
   Click leaf → jump to Section 3 viewer.
   Zoom + pan via d3.zoom.
   ============================================================ */

'use strict';

function vqInitTaxonomy(tree, sequences) {
  const el = document.getElementById('section-taxonomy');
  if (!el) return;

  if (!tree && (!sequences || !sequences.length)) {
    el.innerHTML = '<div class="vq-empty">No taxonomy data available.</div>';
    return;
  }

  // Build tree from sequences if no pre-built tree supplied
  const root = tree || _buildTree(sequences);

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Taxonomy</div>
        <div class="vq-section-sub">Radial cluster tree · click leaf to inspect sequence</div>
      </div>
      <div style="display:flex;gap:var(--vq-space-2);align-items:center">
        <button class="vq-btn vq-btn--sm vq-btn--ghost" id="tax-reset-zoom">Reset zoom</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportSVG(document.querySelector('#section-taxonomy svg'),'taxonomy.svg')">SVG</button>
        <button class="vq-btn vq-btn--sm vq-btn--ghost"
          onclick="vqExportPNG(document.querySelector('#section-taxonomy svg'),'taxonomy.png')">PNG</button>
      </div>
    </div>
    <div id="tax-legend" class="vq-legend" style="margin-bottom:var(--vq-space-3)"></div>
    <div id="tax-wrap" class="vq-genome-wrap" style="overflow:hidden;cursor:grab"></div>
  `;

  _renderRadialTree(root);

  document.getElementById('tax-reset-zoom')?.addEventListener('click', () => {
    _resetZoom();
  });
}

// ── Build tree from sequences[].taxonomy ────────────────────────────────────

function _buildTree(sequences) {
  const root = { name: 'Viruses', children: [] };
  const idx  = {};   // path key → node

  sequences.forEach(seq => {
    const tax = seq.taxonomy || {};
    const parts = [
      tax.phylum  || 'Unclassified',
      tax.order   || 'Unclassified',
      tax.family  || 'Unclassified',
      tax.genus   || 'Unclassified',
    ];

    let parent = root;
    let key    = '';
    parts.forEach((name, depth) => {
      key += '/' + name;
      if (!idx[key]) {
        const node = { name, depth, children: [] };
        parent.children.push(node);
        idx[key] = node;
      }
      parent = idx[key];
    });

    // Leaf = sequence
    parent.children.push({
      name:     seq.seq_id,
      seq_id:   seq.seq_id,
      family:   tax.family || 'Unclassified',
      is_leaf:  true,
    });
  });

  return root;
}

// ── Colour by family ─────────────────────────────────────────────────────────

const _TAX_PALETTE = [
  'var(--vq-tax-flaviviridae)',
  'var(--vq-tax-parvoviridae)',
  'var(--vq-tax-phenuiviridae)',
  'var(--vq-tax-nodaviridae)',
  'var(--vq-tax-rhabdoviridae)',
  'var(--vq-tax-other)',
  '#7c3aed','#0891b2','#b45309','#be185d',
];

let _familyColorMap = {};
let _familyIndex    = 0;

function _familyColor(family) {
  if (!family || family === 'Unclassified') return 'var(--vq-text-3)';
  const cssVar = '--vq-tax-' + family.toLowerCase();
  const computed = getComputedStyle(document.documentElement).getPropertyValue(cssVar).trim();
  if (computed) return computed;
  if (!_familyColorMap[family]) {
    _familyColorMap[family] = _TAX_PALETTE[_familyIndex % _TAX_PALETTE.length];
    _familyIndex++;
  }
  return _familyColorMap[family];
}

// ── Radial tree rendering ────────────────────────────────────────────────────

let _taxZoom   = null;
let _taxSvgSel = null;

function _renderRadialTree(treeData) {
  const wrap = document.getElementById('tax-wrap');
  if (!wrap) return;

  const W   = Math.max(wrap.clientWidth || 900, 700);
  const H   = W;
  const R   = W / 2 - 80;

  const hierarchy = d3.hierarchy(treeData)
    .sort((a, b) => d3.ascending(a.data.name, b.data.name));

  const cluster = d3.cluster().size([2 * Math.PI, R]);
  cluster(hierarchy);

  // Collect families for legend
  const families = [...new Set(
    hierarchy.leaves().map(l => l.data.family).filter(Boolean)
  )];

  _renderLegend(families);

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('width', W)
    .attr('height', H);

  const g = svg.append('g')
    .attr('transform', `translate(${W / 2},${H / 2})`);

  // ── Links ─────────────────────────────────────────────────────────────────
  g.append('g')
    .attr('fill', 'none')
    .attr('stroke', 'var(--vq-border)')
    .attr('stroke-width', 1)
    .selectAll('path')
    .data(hierarchy.links())
    .join('path')
    .attr('d', d3.linkRadial()
      .angle(d  => d.x)
      .radius(d => d.y));

  // ── Nodes ─────────────────────────────────────────────────────────────────
  const node = g.append('g')
    .selectAll('g')
    .data(hierarchy.descendants())
    .join('g')
    .attr('transform', d => `rotate(${d.x * 180 / Math.PI - 90}) translate(${d.y},0)`);

  node.append('circle')
    .attr('r', d => d.data.is_leaf ? 4 : 3)
    .attr('fill', d => {
      if (d.data.is_leaf) return _familyColor(d.data.family);
      if (!d.children) return 'var(--vq-text-3)';
      return d.depth === 0 ? 'var(--vq-primary)' : 'var(--vq-accent)';
    })
    .attr('stroke', '#fff')
    .attr('stroke-width', 1)
    .attr('cursor', d => d.data.is_leaf ? 'pointer' : 'default')
    .on('click', (evt, d) => {
      if (d.data.is_leaf && d.data.seq_id) {
        _jumpToViewer(d.data.seq_id);
      }
    })
    .on('mousemove', (evt, d) => {
      const label = d.data.is_leaf
        ? `<div class="vq-tooltip__title">${d.data.seq_id}</div>
           <div class="vq-tooltip__row">
             <span class="vq-tooltip__key">Family</span><span>${d.data.family || '—'}</span>
           </div>`
        : `<div class="vq-tooltip__title">${d.data.name}</div>`;
      vqTooltipShow(label, evt);
    })
    .on('mouseleave', vqTooltipHide);

  // ── Labels ────────────────────────────────────────────────────────────────
  node.append('text')
    .attr('dy', '0.31em')
    .attr('x', d => d.x < Math.PI === !d.children ? 6 : -6)
    .attr('text-anchor', d => d.x < Math.PI === !d.children ? 'start' : 'end')
    .attr('transform', d => d.x >= Math.PI ? 'rotate(180)' : null)
    .attr('font-size', d => d.data.is_leaf ? 9 : 10)
    .attr('font-weight', d => d.data.is_leaf ? '400' : '600')
    .attr('fill', d => d.data.is_leaf
      ? _familyColor(d.data.family)
      : d.depth === 0 ? 'var(--vq-primary)' : 'var(--vq-text-2)')
    .attr('font-family', d => d.data.is_leaf ? 'var(--vq-font-mono)' : 'var(--vq-font)')
    .attr('cursor', d => d.data.is_leaf ? 'pointer' : 'default')
    .text(d => d.data.name)
    .on('click', (evt, d) => {
      if (d.data.is_leaf && d.data.seq_id) _jumpToViewer(d.data.seq_id);
    });

  // ── Zoom + pan ────────────────────────────────────────────────────────────
  _taxZoom = d3.zoom()
    .scaleExtent([0.3, 5])
    .on('zoom', evt => g.attr('transform', evt.transform.translate(W / 2, H / 2) + ''));

  // Simpler: let zoom manage full transform
  _taxZoom = d3.zoom()
    .scaleExtent([0.3, 5])
    .on('zoom', evt => {
      g.attr('transform',
        `translate(${evt.transform.x + W / 2},${evt.transform.y + H / 2}) scale(${evt.transform.k})`
      );
    });

  _taxSvgSel = svg;
  svg.call(_taxZoom);

  wrap.innerHTML = '';
  wrap.appendChild(svg.node());
}

function _resetZoom() {
  if (_taxSvgSel && _taxZoom) {
    _taxSvgSel.transition().duration(400).call(_taxZoom.transform, d3.zoomIdentity);
  }
}

// ── Legend ───────────────────────────────────────────────────────────────────

function _renderLegend(families) {
  const el = document.getElementById('tax-legend');
  if (!el || !families.length) return;
  el.innerHTML = families.map(f => `
    <div class="vq-legend__item">
      <div class="vq-legend__swatch"
           style="background:${_familyColor(f)};width:12px;height:12px;border-radius:50%"></div>
      <span>${f}</span>
    </div>`).join('');
}

// ── Cross-section navigation ──────────────────────────────────────────────────

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
