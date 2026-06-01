/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_taxonomy.js — Section 4: Taxonomy (Linear Tree)
   Horizontal phylogeny: root at left, leaves at right.
   In-card layout, full width, family filter, PNG/SVG export.
   ============================================================ */

const _TAX = {
  tree:     null,
  families: new Set(),   // currently selected family filters (empty = show all)
  resizeT:  null,
};

const _TAX_PALETTE = [
  'var(--vq-tax-flaviviridae)','var(--vq-tax-parvoviridae)',
  'var(--vq-tax-phenuiviridae)','var(--vq-tax-nodaviridae)',
  'var(--vq-tax-rhabdoviridae)','var(--vq-tax-tombusviridae)',
  'var(--vq-tax-other)',
  '#7c3aed','#0891b2','#b45309','#be185d',
];

let _familyColorMap = {};
let _familyIndex    = 0;

function vqInitTaxonomy(tree, sequences) {
  const el = document.getElementById('section-taxonomy');
  if (!el) return;

  if (!tree && (!sequences || !sequences.length)) {
    el.innerHTML = '<div class="vq-empty">No taxonomy data available.</div>';
    return;
  }

  _TAX.tree = tree || _buildTree(sequences);
  const families = _allFamilies(_TAX.tree);

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">Taxonomy</div>
        <div class="vq-section-sub">
          Linear phylogeny · click a leaf to inspect the sequence
        </div>
      </div>
      <div class="vq-section-actions">
        <div class="vq-menu" id="tax-export-menu">
          <button class="vq-btn vq-btn--sm vq-btn--ghost" data-menu-toggle type="button">
            Export
            <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="2.5">
              <polyline points="6 9 12 15 18 9"/>
            </svg>
          </button>
          <div class="vq-menu__panel" role="menu">
            <button class="vq-menu__item" data-fmt="png" role="menuitem" type="button">PNG image</button>
            <button class="vq-menu__item" data-fmt="svg" role="menuitem" type="button">SVG vector</button>
          </div>
        </div>
      </div>
    </div>

    <div class="vq-toolbar__filter-row" style="margin-bottom:var(--vq-space-3)">
      <span class="vq-filter-label">Families</span>
      ${families.map(f => `
        <button class="vq-tax-pill" data-fam="${VQ.esc(f)}" type="button"
                style="--pill-color:${_familyColor(f)}">
          <span class="vq-tax-pill__dot" style="background:${_familyColor(f)}"></span>
          ${VQ.esc(f)}
        </button>`).join('')}
      <button class="vq-btn vq-btn--sm vq-btn--ghost" id="tax-clear-filter"
              type="button" style="display:none;margin-left:auto">
        Clear filter
      </button>
    </div>

    <div class="vq-card">
      <div class="vq-card__header">
        <div class="vq-card__title">Viral Taxonomy</div>
      </div>
      <div class="vq-card__body" style="padding:0;min-height:520px;display:flex;">
        <div id="tax-wrap" style="width:100%;overflow:auto;display:flex;align-items:stretch"></div>
      </div>
    </div>
  `;

  // Family pill toggle
  el.querySelectorAll('.vq-tax-pill[data-fam]').forEach(btn => {
    btn.addEventListener('click', () => _toggleFamily(btn.dataset.fam));
  });

  // Clear all
  document.getElementById('tax-clear-filter')?.addEventListener('click', () => {
    _TAX.families.clear();
    _updatePills();
    _draw();
  });

  // Export menu
  const exportMenu = document.getElementById('tax-export-menu');
  const exportTog  = exportMenu?.querySelector('[data-menu-toggle]');
  exportTog?.addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(m => {
      if (m !== exportMenu) m.classList.remove('open');
    });
    exportMenu.classList.toggle('open');
  });
  exportMenu?.querySelectorAll('[data-fmt]').forEach(item => {
    item.addEventListener('click', e => {
      e.stopPropagation();
      exportMenu.classList.remove('open');
      const svg = document.querySelector('#tax-wrap svg');
      if (!svg) return;
      if (item.dataset.fmt === 'png') VQ.exportPNG(svg, 'taxonomy.png');
      else                            VQ.exportSVG(svg, 'taxonomy.svg');
    });
  });

  // Render + resize handling. ResizeObserver covers the 0→real width
  // transition that happens when the section is first revealed.
  let lastW = 0;
  const wrap = document.getElementById('tax-wrap');
  const ro = new ResizeObserver(() => {
    const w = wrap.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    clearTimeout(_TAX.resizeT);
    _TAX.resizeT = setTimeout(_draw, 60);
  });
  if (wrap) ro.observe(wrap);
  _draw();
}

// ── Build the tree from sequences ──────────────────────────────────────────

function _buildTree(sequences) {
  const root = { name: 'Viruses', children: [] };
  const idx  = {};

  sequences.forEach(seq => {
    const tax = seq.taxonomy || {};
    const parts = [
      tax.phylum || 'Unclassified',
      tax.order  || 'Unclassified',
      tax.family || 'Unclassified',
      tax.genus  || 'Unclassified',
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

    parent.children.push({
      name:    seq.id,
      seq_id:  seq.id,
      family:  tax.family || 'Unclassified',
      is_leaf: true,
    });
  });

  return root;
}

function _allFamilies(tree) {
  const set = new Set();
  function walk(node) {
    if (node.is_leaf && node.family) set.add(node.family);
    (node.children || []).forEach(walk);
  }
  walk(tree);
  return [...set].sort();
}

// ── Filter ─────────────────────────────────────────────────────────────────

function _filteredTree() {
  if (!_TAX.families.size) return _TAX.tree;
  function prune(node) {
    if (node.is_leaf) {
      return _TAX.families.has(node.family) ? { ...node } : null;
    }
    const kids = (node.children || []).map(prune).filter(Boolean);
    if (!kids.length && node !== _TAX.tree) return null;
    return { ...node, children: kids };
  }
  const pruned = prune(_TAX.tree);
  return pruned || { name: 'Viruses', children: [] };
}

function _toggleFamily(fam) {
  if (_TAX.families.has(fam)) _TAX.families.delete(fam);
  else                         _TAX.families.add(fam);
  _updatePills();
  _draw();
}

function _updatePills() {
  document.querySelectorAll('.vq-tax-pill[data-fam]').forEach(btn => {
    btn.classList.toggle('active', _TAX.families.has(btn.dataset.fam));
  });
  const clearBtn = document.getElementById('tax-clear-filter');
  if (clearBtn) clearBtn.style.display = _TAX.families.size ? '' : 'none';
}

// ── Family colour ──────────────────────────────────────────────────────────

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

// ── Draw linear tree ───────────────────────────────────────────────────────

function _draw() {
  try {
    _drawInner();
  } catch (e) {
    console.error('taxonomy draw failed:', e.message, e.stack);
  }
}

function _drawInner() {
  const wrap = document.getElementById('tax-wrap');
  if (!wrap) return;
  wrap.innerHTML = '';

  const data = _filteredTree();
  const leaves = _countLeaves(data);
  if (!leaves) {
    wrap.innerHTML = '<div class="vq-empty">No leaves match the filter.</div>';
    return;
  }

  // Geometry — fill BOTH the available width and the card body height
  const W = Math.max(wrap.clientWidth || 0, 720);

  const PAD_T = 28;
  const PAD_B = 28;
  const PAD_L = 90;
  const RIGHT_LABEL_W = 260;
  const MIN_ROW = 26;
  const MAX_ROW = 36;

  // Pick a row height that uses the available vertical space, capped.
  const availH = Math.max(wrap.clientHeight || 0, 520) - PAD_T - PAD_B;
  const rowH   = Math.max(MIN_ROW, Math.min(MAX_ROW, Math.floor(availH / leaves)));
  const treeH  = leaves * rowH;
  const H      = PAD_T + PAD_B + treeH;
  const drawW  = W - PAD_L - RIGHT_LABEL_W;

  const hierarchy = d3.hierarchy(data)
    .sort((a, b) => d3.ascending(a.data.name, b.data.name));
  d3.cluster().size([treeH, drawW])(hierarchy);

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${H}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  const g = svg.append('g').attr('transform', `translate(${PAD_L}, ${PAD_T})`);

  // Step links (Cartesian elbow connectors)
  g.append('g')
    .attr('fill', 'none')
    .attr('stroke', 'var(--vq-border-dark)')
    .attr('stroke-width', 1.2)
    .selectAll('path')
    .data(hierarchy.links())
    .join('path')
    .attr('d', d => {
      const sx = d.source.y, sy = d.source.x;
      const tx = d.target.y, ty = d.target.x;
      return `M${sx},${sy}V${ty}H${tx}`;
    });

  // Render nodes — internal nodes get a pill chip with the taxon name inside,
  // leaves get a coloured circle + monospace sequence id label.
  const node = g.append('g')
    .selectAll('g')
    .data(hierarchy.descendants())
    .join('g')
    .attr('transform', d => `translate(${d.y},${d.x})`);

  // ── Leaves ──────────────────────────────────────────────────────────────
  const leafSel = node.filter(d => d.data.is_leaf);

  leafSel.append('circle')
    .attr('r', 4.5)
    .attr('fill', d => _familyColor(d.data.family))
    .attr('stroke', '#fff').attr('stroke-width', 1.5)
    .attr('cursor', 'pointer')
    .on('click', (evt, d) => { if (d.data.seq_id) VQ.jumpToViewer(d.data.seq_id); })
    .on('mousemove', (evt, d) => VQ.tooltipShow(
      `<div class="vq-tooltip__title">${VQ.esc(d.data.seq_id)}</div>
       <div class="vq-tooltip__row">
         <span class="vq-tooltip__key">Family</span><span>${VQ.esc(d.data.family || '—')}</span>
       </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  leafSel.append('text')
    .attr('x', 9).attr('y', 4)
    .attr('font-size', 12)
    .attr('font-family', 'var(--vq-font-mono)')
    .attr('fill', d => _familyColor(d.data.family))
    .attr('cursor', 'pointer')
    .text(d => d.data.name)
    .on('click', (evt, d) => { if (d.data.seq_id) VQ.jumpToViewer(d.data.seq_id); });

  // ── Internal nodes — pill chip centered on the horizontal incoming branch.
  const internalSel = node.filter(d => !d.data.is_leaf && d.depth > 0);

  internalSel.each(function (d) {
    const sel       = d3.select(this);
    const text      = d.data.name;
    // Measure roughly: 6.4px per char @ 12px font, plus padding.
    const w         = Math.max(24, text.length * 6.6 + 14);
    const h         = 18;
    // Center horizontally on the half-branch leading INTO this node.
    const parentY   = d.parent ? d.parent.y : 0;
    const segLen    = d.y - parentY;
    const cx        = -segLen / 2;   // relative to node (which is at d.y)
    const cy        = -h / 2 - 1;    // sit ABOVE the horizontal link

    sel.append('rect')
      .attr('x', cx - w / 2).attr('y', cy)
      .attr('width', w).attr('height', h)
      .attr('rx', h / 2).attr('ry', h / 2)
      .attr('fill', 'var(--vq-surface)')
      .attr('stroke', 'var(--vq-border-dark)')
      .attr('stroke-width', 1);
    sel.append('text')
      .attr('x', cx).attr('y', cy + h / 2 + 4)
      .attr('text-anchor', 'middle')
      .attr('font-size', 11)
      .attr('font-weight', 500)
      .attr('fill', 'var(--vq-text-2)')
      .attr('font-family', 'var(--vq-font)')
      .text(text);

    // Small node dot at the actual branch junction
    sel.append('circle')
      .attr('r', 3.5)
      .attr('fill', d.depth === 1 ? 'var(--vq-accent)' : 'var(--vq-accent-dark)')
      .attr('stroke', '#fff').attr('stroke-width', 1.5);
  });

  internalSel
    .on('mousemove', (evt, d) => VQ.tooltipShow(
      `<div class="vq-tooltip__title">${VQ.esc(d.data.name)}</div>
       <div class="vq-tooltip__row">
         <span class="vq-tooltip__key">Depth</span><span>${d.depth}</span>
         <span class="vq-tooltip__key">Leaves</span><span>${d.leaves().length}</span>
       </div>`, evt))
    .on('mouseleave', VQ.tooltipHide);

  // Root label (depth 0)
  node.filter(d => d.depth === 0)
    .append('circle')
      .attr('r', 5)
      .attr('fill', 'var(--vq-primary)')
      .attr('stroke', '#fff').attr('stroke-width', 2);
  node.filter(d => d.depth === 0)
    .append('text')
      .attr('x', -10).attr('y', 4)
      .attr('text-anchor', 'end')
      .attr('font-size', 12)
      .attr('font-weight', 600)
      .attr('fill', 'var(--vq-primary)')
      .text(d => d.data.name);

  wrap.appendChild(svg.node());
}

function _countLeaves(node) {
  if (node.is_leaf) return 1;
  return (node.children || []).reduce((a, c) => a + _countLeaves(c), 0);
}


window.vqInitTaxonomy = vqInitTaxonomy;
})();
