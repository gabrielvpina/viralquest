/* This file's private helpers live in an IIFE so they don't collide across sections. */
(function () {
'use strict';

/* ============================================================
   section_synteny.js — Synteny with the BLASTn reference genome
   One synteny cluster at a time (picked from a list):
     • reference genome on top — its GenBank CDS (+ mat_peptides)
     • every member sequence below, placed where it aligns on the
       reference, with ONLY the ORFs that DIAMOND blastp matched to a
       reference CDS (the canonical ones) and their HMM domains
     • ribbons joining each aligned ORF stretch to its CDS stretch,
       shaded by amino-acid identity
     • a table of the ORF ↔ CDS pairs
   Data: report.synteny = { references: {acc: RefGenome}, clusters:
   [{reference, members}] } and sequence.synteny (see viralquest/synteny.py).
   Drawing blocks (ORF arrows, domain lanes, domain colours) come from the
   viewer via VQ.genome, so both tabs look alike.
   ============================================================ */

const _SY = {
  refs:     {},
  clusters: [],
  seqById:  {},
  idx:      0,
  showAll:  false,   // also draw the ORFs with no reference match (greyed)
  el:       null,
};

const PAD_L = 16, PAD_R = 16;
const ORF_H = 14, DOM_H = 6, DOM_G = 1.5, LANE_G = 6;
const PEP_H = 6, RIB_H = 64, LABEL_H = 16;

// ── Init ───────────────────────────────────────────────────────────────────

function vqInitSynteny(report, mountId) {
  const el = document.getElementById(mountId || 'section-synteny');
  if (!el) return;
  const syn = (report && report.synteny) || {};
  _SY.el       = el;
  _SY.refs     = syn.references || {};
  _SY.clusters = (syn.clusters || []).filter(c => _SY.refs[c.reference]);
  _SY.seqById  = {};
  (report.sequences || []).forEach(s => { _SY.seqById[s.id] = s; });
  _SY.idx = 0;

  if (!_SY.clusters.length) {
    el.innerHTML = `
      <div class="vq-section-header"><div class="vq-section-title">Synteny</div></div>
      <div class="vq-empty">No synteny data — the step needs BLASTn hits with a viral title
        (and NCBI access to download their GenBank records).</div>`;
    return;
  }

  const n = _SY.clusters.length;
  el.innerHTML = `
    <div class="vq-section-header">
      <div class="vq-section-title">
        Synteny
        <span class="vq-count-label">${n} reference${n !== 1 ? 's' : ''}</span>
      </div>
      <div class="vq-section-actions">
        <label class="vq-filter-check">
          <input type="checkbox" id="syn-show-all"> Show non-canonical ORFs
        </label>
        <div class="vq-menu" id="syn-export-menu">
          <button class="vq-btn vq-btn--ghost vq-btn--sm" data-menu-toggle type="button">
            Export
            <svg width="10" height="10" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="2.5"><polyline points="6 9 12 15 18 9"/></svg>
          </button>
          <div class="vq-menu__panel" role="menu">
            <button class="vq-menu__item" data-fmt="png" role="menuitem" type="button">PNG image</button>
            <button class="vq-menu__item" data-fmt="svg" role="menuitem" type="button">SVG vector</button>
            <button class="vq-menu__item" data-fmt="pdf" role="menuitem" type="button">PDF (print)</button>
            <button class="vq-menu__item" data-fmt="tsv" role="menuitem" type="button">ORF ↔ CDS pairs (TSV)</button>
          </div>
        </div>
      </div>
    </div>

    <div class="vq-toolbar vq-syn-picker">
      <button class="vq-btn vq-btn--ghost vq-btn--sm" id="syn-prev" type="button" aria-label="Previous reference">◀</button>
      <div class="vq-search-wrap">
        <svg width="14" height="14" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="2" aria-hidden="true">
          <circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/>
        </svg>
        <input class="vq-input" id="syn-search" type="search"
               placeholder="Filter by virus, accession or sequence ID…" aria-label="Filter references">
      </div>
      <select class="vq-select vq-syn-select" id="syn-select" aria-label="Synteny cluster"></select>
      <button class="vq-btn vq-btn--ghost vq-btn--sm" id="syn-next" type="button" aria-label="Next reference">▶</button>
      <span class="vq-count-label" id="syn-pos"></span>
    </div>

    <div id="syn-body"></div>`;

  _fillSelect('');
  document.getElementById('syn-search').addEventListener('input', e => _fillSelect(e.target.value));
  document.getElementById('syn-select').addEventListener('change', e => _show(+e.target.value));
  document.getElementById('syn-prev').addEventListener('click', () => _step(-1));
  document.getElementById('syn-next').addEventListener('click', () => _step(1));
  document.getElementById('syn-show-all').addEventListener('change', e => {
    _SY.showAll = e.target.checked;
    _render();
  });
  _wireExport();
  _show(0);
  VQ.redrawOnResize(el, () => _drawMap());
}

function _clusterLabel(c) {
  const ref = _SY.refs[c.reference];
  const n   = c.members.length;
  return `${ref.organism || ref.title} — ${ref.accession} (${n} seq${n !== 1 ? 's' : ''})`;
}

/* The list keeps every cluster index as its value; the filter only hides. */
function _fillSelect(query) {
  const sel = document.getElementById('syn-select');
  const q   = query.trim().toLowerCase();
  const opts = _SY.clusters.map((c, i) => {
    const ref = _SY.refs[c.reference];
    const hay = [ref.organism, ref.title, ref.accession, ...c.members].join(' ').toLowerCase();
    return (!q || hay.includes(q)) ? `<option value="${i}">${VQ.esc(_clusterLabel(c))}</option>` : '';
  }).join('');
  sel.innerHTML = opts || '<option disabled>No reference matches the filter</option>';
  if (opts && !sel.querySelector(`option[value="${_SY.idx}"]`)) {
    _show(+sel.options[0].value);
  } else {
    sel.value = String(_SY.idx);
  }
}

function _step(dir) {
  const sel  = document.getElementById('syn-select');
  const vals = [...sel.options].filter(o => !o.disabled).map(o => +o.value);
  if (!vals.length) return;
  const at = vals.indexOf(_SY.idx);
  _show(vals[(at + dir + vals.length) % vals.length]);
}

function _show(i) {
  _SY.idx = i;
  const sel = document.getElementById('syn-select');
  if (sel) sel.value = String(i);
  const pos = document.getElementById('syn-pos');
  if (pos) pos.textContent = `${i + 1} / ${_SY.clusters.length}`;
  _render();
}

// ── Cluster model ──────────────────────────────────────────────────────────

/* Members of the current cluster, with their ORFs split into matched / not. */
function _current() {
  const c   = _SY.clusters[_SY.idx];
  const ref = _SY.refs[c.reference];
  const members = c.members.map(id => _SY.seqById[id]).filter(s => s && s.synteny).map(seq => {
    const syn  = seq.synteny;
    const hits = syn.hits || [];
    const matched = new Set(hits.map(h => h.orf_name));
    return {
      seq, syn, hits,
      orfs:  (seq.orfs || []).filter(o => matched.has(o.name)),
      other: (seq.orfs || []).filter(o => !matched.has(o.name)),
    };
  });
  return { c, ref, members };
}

/* Contig position → reference position, from the placement in synteny.py. */
function _toRef(m, x) {
  return m.syn.flipped ? m.syn.offset - x : x + m.syn.offset;
}
function _ivToRef(m, a, b) {
  const p = _toRef(m, a), q = _toRef(m, b);
  return [Math.min(p, q), Math.max(p, q)];
}
function _shownStrand(m, strand) {
  return m.syn.flipped ? (strand === '-' ? '+' : '-') : strand;
}

/* aa range on an ORF → [start, end) on the contig (strand-aware). */
function _orfAaToNt(orf, a, b) {
  if (orf.strand === '-') return [orf.stop_position - 3 * b, orf.stop_position - 3 * (a - 1)];
  return [orf.start_position + 3 * (a - 1), orf.start_position + 3 * b];
}

/* Greedy lane packing of [start, end) items. */
function _lanes(items) {
  const ends = [];
  return [...items].sort((a, b) => a.start - b.start).map(it => {
    let lane = ends.findIndex(e => e <= it.start);
    if (lane < 0) { lane = ends.length; ends.push(0); }
    ends[lane] = it.end;
    return { ...it, lane };
  });
}

// ── Render ─────────────────────────────────────────────────────────────────

function _render() {
  const body = document.getElementById('syn-body');
  if (!body) return;
  const { ref, members } = _current();
  const esc  = VQ.esc;
  const matchedCds = new Set(members.flatMap(m => m.hits.map(h => h.cds_index)));
  const noMatch    = members.filter(m => !m.hits.length).map(m => m.seq.id);
  const refUrl     = VQ.ncbiUrl(ref.accession, 'nuccore');
  const accHtml    = refUrl
    ? `<a class="vq-acc-link" href="${esc(refUrl)}" target="_blank" rel="noopener noreferrer">${esc(ref.accession)}</a>`
    : esc(ref.accession);

  body.innerHTML = `
    <div class="vq-card vq-syn-summary">
      <div class="vq-syn-summary__title">${esc(ref.organism || ref.title)}</div>
      <div class="vq-syn-summary__sub">${esc(ref.title)}</div>
      <div class="vq-syn-kv">
        <span>Reference</span><span>${accHtml}</span>
        <span>Length</span><span>${ref.length.toLocaleString()} nt</span>
        <span>CDS matched</span><span>${matchedCds.size} / ${ref.cds.length}</span>
        <span>Sequences</span><span>${members.length}</span>
      </div>
      ${noMatch.length ? `<div class="vq-syn-note">No ORF matched the reference proteins in:
        ${noMatch.map(id => `<code>${esc(id)}</code>`).join(', ')}</div>` : ''}
    </div>
    <div class="vq-card vq-syn-map"><div class="vq-genome-wrap" id="syn-map"></div></div>
    <div class="vq-card vq-syn-pairs">${_pairsTable(ref, members)}</div>`;

  body.querySelectorAll('[data-jump]').forEach(a => a.addEventListener('click', e => {
    e.preventDefault();
    VQ.jumpToViewer(a.dataset.jump);
  }));
  _drawMap();
}

function _drawMap() {
  const wrap = document.getElementById('syn-map');
  if (!wrap || !_SY.clusters.length) return;
  wrap.innerHTML = '';
  wrap.appendChild(_syntenySVG(_current(), wrap.clientWidth));
}

// ── SVG ────────────────────────────────────────────────────────────────────

function _syntenySVG({ ref, members }, containerWidth) {
  const G   = VQ.genome;
  const esc = VQ.esc;
  const W   = Math.max(containerWidth || 0, 720);

  // Common scale: the reference plus every member's placed extent.
  let lo = 0, hi = ref.length;
  members.forEach(m => {
    const [a, b] = _ivToRef(m, 0, m.seq.length || 0);
    lo = Math.min(lo, a); hi = Math.max(hi, b);
  });
  const X = d3.scaleLinear([lo, hi], [PAD_L, W - PAD_R]);

  // Reference CDS lanes and the peptide row under each lane.
  const cdsItems = _lanes(ref.cds.map((c, i) => ({ start: c.start, end: c.end, c, i })));
  const nCdsLanes = Math.max(1, ...cdsItems.map(it => it.lane + 1));
  const pepByCds  = {};
  (ref.peptides || []).forEach(p => {
    if (p.parent != null) (pepByCds[p.parent] = pepByCds[p.parent] || []).push(p);
  });
  const REF_LANE_H = ORF_H + PEP_H + 3 + LANE_G;
  const matchedCds = new Set(members.flatMap(m => m.hits.map(h => h.cds_index)));

  // Member layouts.
  const AXIS_Y = 24;
  let y = AXIS_Y + 10 + LABEL_H;
  const refY = y;
  y += nCdsLanes * REF_LANE_H;

  const layouts = members.map(m => {
    const item = (o, matched) => {
      const [a, b] = _ivToRef(m, o.start_position, o.stop_position);
      const doms = G.assignDomainLanes(G.bestDomainPerDatabase(o.domains || []));
      return { start: a, end: b, orf: o, doms, matched };
    };
    // Matched ORFs take the top lanes, where the ribbons land; non-canonical
    // ones, when shown, are packed in the lanes below.
    const top   = _lanes(m.orfs.map(o => item(o, true)));
    const nTop  = Math.max(0, ...top.map(it => it.lane + 1));
    const below = _SY.showAll
      ? _lanes(m.other.map(o => item(o, false))).map(it => ({ ...it, lane: it.lane + nTop }))
      : [];
    const items = top.concat(below);
    const laneDom = [];
    items.forEach(it => {
      const nd = Math.max(0, ...it.doms.map(d => d.lane + 1));
      laneDom[it.lane] = Math.max(laneDom[it.lane] || 0, nd);
    });
    const nLanes = Math.max(1, laneDom.length);
    const laneY = [];
    const ribTop = y;
    let ly = y + RIB_H + LABEL_H;
    for (let l = 0; l < nLanes; l++) {
      laneY.push(ly);
      ly += ORF_H + (laneDom[l] || 0) * (DOM_H + DOM_G) + LANE_G;
    }
    const layout = { m, items, laneY, ribTop, labelY: ribTop + RIB_H + 11 };
    y = ly + 4;
    return layout;
  });
  const totalH = y + 30;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');

  // Axis (reference coordinates).
  const axisG = svg.append('g').attr('transform', `translate(0, ${AXIS_Y})`)
    .call(d3.axisTop(X).ticks(Math.min(14, Math.max(4, Math.round(W / 120)))).tickSize(0)
      .tickFormat(hi - lo > 50000 ? d3.format(',') : null));
  axisG.selectAll('.tick line').attr('y1', 0).attr('y2', totalH - AXIS_Y - 26)
    .attr('stroke', 'var(--vq-border)').attr('stroke-dasharray', '3,3');
  axisG.select('.domain').attr('stroke', 'var(--vq-border-dark)');
  axisG.selectAll('.tick text').style('font-size', '9.5px').style('fill', 'var(--vq-text-3)');

  const label = (x, yy, text, bold) => svg.append('text')
    .attr('x', x).attr('y', yy).attr('font-size', 10.5)
    .attr('font-weight', bold ? 600 : 400).attr('fill', 'var(--vq-text-2)').text(text);

  // ── Ribbons (behind the tracks) ──
  const ribG = svg.append('g').attr('class', 'vq-syn-ribbons');
  const refLaneOf = {};
  cdsItems.forEach(it => { refLaneOf[it.i] = it.lane; });
  const refBottom = refY + nCdsLanes * REF_LANE_H - LANE_G;

  layouts.forEach(L => {
    const laneOfOrf = {};
    L.items.forEach(it => { laneOfOrf[it.orf.name] = it.lane; });
    L.m.hits.forEach(h => {
      if (laneOfOrf[h.orf_name] == null) return;
      const yTop = refBottom;
      const yBot = L.laneY[0] - 2;
      const cds  = ref.cds[h.cds_index];
      const op   = 0.12 + 0.55 * Math.max(0, Math.min(1, (h.pident - 30) / 70));
      // One ribbon per CDS part (exon) the alignment spans; reports without
      // segments fall back to the envelope.
      const segs = (h.segments && h.segments.length)
        ? h.segments : [[h.ref_nt[0], h.ref_nt[1], h.orf_nt[0], h.orf_nt[1]]];
      segs.forEach(([r0, r1, o0, o1]) => {
      const [c0, c1] = _ivToRef(L.m, o0, o1);
      ribG.append('polygon')
        .attr('points', `${X(r0)},${yTop} ${X(r1)},${yTop} ${X(c1)},${yBot} ${X(c0)},${yBot}`)
        .attr('fill', 'var(--vq-accent)').attr('opacity', op)
        .attr('stroke', 'var(--vq-accent-dark)').attr('stroke-width', 0.4).attr('stroke-opacity', 0.5)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(h.orf_name)} ↔ ${esc(cds.product)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Identity</span><span>${h.pident.toFixed(1)}%</span>
            <span class="vq-tooltip__key">E-value</span><span>${h.evalue.toExponential(1)}</span>
            <span class="vq-tooltip__key">ORF cov</span><span>${h.q_cov}% (aa ${h.q_start}–${h.q_end})</span>
            <span class="vq-tooltip__key">CDS cov</span><span>${h.s_cov}% (aa ${h.s_start}–${h.s_end})</span>
            <span class="vq-tooltip__key">Protein</span><span>${esc(cds.protein_id || '—')}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
      });
    });
  });

  // ── Reference track ──
  label(PAD_L, refY - 6, `Reference · ${ref.accession}`, true);
  svg.append('line').attr('x1', X(0)).attr('x2', X(ref.length))
    .attr('y1', refY + ORF_H / 2).attr('y2', refY + ORF_H / 2)
    .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1.5);
  cdsItems.forEach(it => {
    const c  = it.c;
    const yy = refY + it.lane * REF_LANE_H;
    const x1 = X(c.start), x2 = Math.max(X(c.end), x1 + 4);
    const hit = matchedCds.has(it.i);
    const g  = svg.append('g');
    g.append('polygon')
      .attr('points', G.orfArrowPoints(x1, x2, yy, ORF_H, c.strand))
      .attr('fill', hit ? 'var(--vq-primary-light)' : 'var(--vq-surface-2)')
      .attr('stroke', hit ? 'none' : 'var(--vq-border-dark)')
      .attr('stroke-dasharray', hit ? null : '3,2')
      .attr('opacity', hit ? 0.9 : 1)
      .on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(c.product)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Protein</span><span>${esc(c.protein_id || '—')}</span>
          <span class="vq-tooltip__key">Position</span><span>${c.start + 1}–${c.end} (${c.strand})</span>
          <span class="vq-tooltip__key">Length</span><span>${c.length_aa} aa</span>
          <span class="vq-tooltip__key">Matched</span><span>${hit ? 'yes' : 'no ORF of these sequences'}</span>
        </div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
    if (x2 - x1 > c.product.length * 5.6 + 14) {
      g.append('text').attr('x', (x1 + x2) / 2).attr('y', yy + ORF_H / 2 + 3.5)
        .attr('text-anchor', 'middle').attr('font-size', 9.5).attr('pointer-events', 'none')
        .attr('fill', hit ? 'rgba(255,255,255,.95)' : 'var(--vq-text-2)').text(c.product);
    }
    (pepByCds[it.i] || []).forEach((p, k) => {
      const px1 = X(p.start), px2 = Math.max(X(p.end), px1 + 2);
      svg.append('rect').attr('x', px1).attr('width', px2 - px1)
        .attr('y', yy + ORF_H + 2).attr('height', PEP_H).attr('rx', 1)
        .attr('fill', k % 2 ? 'var(--vq-text-3)' : 'var(--vq-border-dark)').attr('opacity', 0.75)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(p.product)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Mature peptide</span><span>${p.start + 1}–${p.end}</span>
            <span class="vq-tooltip__key">Length</span><span>~${p.length_aa} aa</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
    });
  });

  // ── Member tracks ──
  layouts.forEach(L => {
    const { m } = L;
    const [e0, e1] = _ivToRef(m, 0, m.seq.length || 0);
    label(PAD_L, L.labelY,
      `${m.seq.id} · ${(m.seq.length || 0).toLocaleString()} nt${m.syn.flipped ? ' · reverse complement' : ''}`
      + (m.hits.length ? '' : ' · no matching ORF'), true)
      .style('cursor', 'pointer').on('click', () => VQ.jumpToViewer(m.seq.id));
    svg.append('line').attr('x1', X(e0)).attr('x2', X(e1))
      .attr('y1', L.laneY[0] + ORF_H / 2).attr('y2', L.laneY[0] + ORF_H / 2)
      .attr('stroke', 'var(--vq-border-dark)').attr('stroke-width', 1.5).attr('stroke-dasharray', '2,3');

    L.items.forEach(it => {
      const o  = it.orf;
      const yy = L.laneY[it.lane];
      const x1 = X(it.start), x2 = Math.max(X(it.end), x1 + 4);
      const st = _shownStrand(m, o.strand);
      const col = !it.matched ? 'var(--vq-muted)' : (st === '-' ? 'var(--vq-warning)' : 'var(--vq-accent)');
      const g = svg.append('g').attr('opacity', it.matched ? 1 : 0.45);
      g.append('polygon')
        .attr('points', G.orfArrowPoints(x1, x2, yy, ORF_H, st))
        .attr('fill', col).attr('opacity', .88).style('cursor', 'pointer')
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(o.name)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Frame</span><span>${esc(o.frame)}</span>
            <span class="vq-tooltip__key">Length</span><span>${o.length_aa} aa</span>
            <span class="vq-tooltip__key">Position</span><span>${o.start_position}–${o.stop_position} (${esc(o.strand)})</span>
            <span class="vq-tooltip__key">Reference</span><span>${it.matched ? 'matched' : 'no match (non-canonical)'}</span>
          </div>
          <div style="margin-top:5px;font-size:10px;color:rgba(255,255,255,.55)">Click to open in the Sequence Viewer</div>`, evt))
        .on('mouseleave', VQ.tooltipHide)
        .on('click', () => VQ.jumpToViewer(m.seq.id));

      it.doms.forEach(d => {
        const [a, b] = _orfAaToNt(o, d.dom.start, d.dom.stop);
        const [ra, rb] = _ivToRef(m, a, b);
        const dx1 = X(ra), dw = Math.max(X(rb) - dx1, 3);
        g.append('rect').attr('x', dx1).attr('width', dw)
          .attr('y', yy + ORF_H + d.lane * (DOM_H + DOM_G) + 1).attr('height', DOM_H).attr('rx', 2)
          .attr('fill', G.domainColor(d.dom.target)).attr('opacity', .9)
          .on('mousemove', evt => VQ.tooltipShow(`
            <div class="vq-tooltip__title">${esc(d.dom.target)}</div>
            <div class="vq-tooltip__row">
              <span class="vq-tooltip__key">Database</span><span>${esc(d.dom.database)}</span>
              <span class="vq-tooltip__key">Score</span><span>${d.dom.score}</span>
              <span class="vq-tooltip__key">E-value</span><span>${d.dom.e_value?.toExponential(2) ?? '—'}</span>
              <span class="vq-tooltip__key">aa range</span><span>${d.dom.start}–${d.dom.stop}</span>
            </div>`, evt))
          .on('mouseleave', VQ.tooltipHide);
      });
    });
  });

  // ── Legend ──
  const lg = svg.append('g').attr('transform', `translate(${PAD_L}, ${totalH - 18})`);
  const txt = (x, t) => lg.append('text').attr('x', x).attr('y', 8.5)
    .style('font-size', '9.5px').style('fill', 'var(--vq-text-3)').text(t);
  lg.append('polygon').attr('points', G.orfArrowPoints(0, 18, 0, 9, '+')).attr('fill', 'var(--vq-primary-light)');
  txt(24, 'Reference CDS');
  lg.append('polygon').attr('points', G.orfArrowPoints(104, 122, 0, 9, '+'))
    .attr('fill', 'var(--vq-surface-2)').attr('stroke', 'var(--vq-border-dark)').attr('stroke-dasharray', '3,2');
  txt(128, 'CDS without match');
  lg.append('polygon').attr('points', G.orfArrowPoints(226, 244, 0, 9, '+')).attr('fill', 'var(--vq-accent)');
  txt(250, 'Matched ORF');
  lg.append('rect').attr('x', 322).attr('width', 12).attr('height', 5).attr('y', 2).attr('rx', 1)
    .attr('fill', 'var(--vq-dom-1)');
  txt(340, 'HMM domain');
  [0.15, 0.4, 0.67].forEach((op, k) => lg.append('rect').attr('x', 416 + k * 12).attr('width', 12)
    .attr('height', 9).attr('fill', 'var(--vq-accent)').attr('opacity', op));
  txt(456, 'aa identity 30% → 100%');

  return svg.node();
}

// ── Pairs table ────────────────────────────────────────────────────────────

function _pairRows(ref, members) {
  return members.flatMap(m => m.hits.map(h => {
    const orf = (m.seq.orfs || []).find(o => o.name === h.orf_name) || {};
    const cds = ref.cds[h.cds_index];
    const doms = _orfDomainNames(orf);
    return { m, h, orf, cds, doms };
  }));
}

function _orfDomainNames(orf) {
  const seen = [];
  VQ.genome.bestDomainPerDatabase(orf.domains || []).forEach(d => {
    if (!seen.includes(d.target)) seen.push(d.target);
  });
  return seen;
}

function _pairsTable(ref, members) {
  const esc  = VQ.esc;
  const rows = _pairRows(ref, members);
  if (!rows.length) return '<div class="vq-empty">No ORF ↔ CDS pairs for this reference.</div>';
  const body = rows.map(({ m, h, orf, cds, doms }) => {
    const url = VQ.ncbiUrl(cds.protein_id, 'protein');
    return `
    <tr>
      <td class="vq-td--mono"><a href="#" data-jump="${esc(m.seq.id)}">${esc(m.seq.id)}</a></td>
      <td class="vq-td--mono">${esc(h.orf_name)}</td>
      <td class="vq-td--num">${esc(orf.frame ?? '—')}</td>
      <td class="vq-td--num">${orf.length_aa ?? '—'}</td>
      <td>${esc(cds.product)}</td>
      <td class="vq-td--mono">${url
        ? `<a class="vq-acc-link" href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(cds.protein_id)}</a>`
        : esc(cds.protein_id || '—')}</td>
      <td class="vq-td--num">${h.pident.toFixed(1)}%</td>
      <td class="vq-td--num">${h.q_cov}%</td>
      <td class="vq-td--num">${h.s_cov}%</td>
      <td class="vq-td--num">${h.evalue.toExponential(1)}</td>
      <td>${doms.map(esc).join(', ') || '—'}</td>
    </tr>`;
  }).join('');
  return `
    <table class="vq-table">
      <thead><tr>
        <th>Sequence</th><th>ORF</th><th>Frame</th><th>aa</th><th>Reference CDS</th>
        <th>Protein</th><th>Identity</th><th>ORF cov</th><th>CDS cov</th><th>E-value</th><th>Domains</th>
      </tr></thead>
      <tbody>${body}</tbody>
    </table>`;
}

function _pairsTSV(ref, members) {
  const cols = ['sequence', 'orf', 'frame', 'orf_aa', 'reference', 'cds_product', 'protein_id',
                'pident', 'orf_cov', 'cds_cov', 'evalue', 'bitscore', 'domains'];
  const lines = _pairRows(ref, members).map(({ m, h, orf, cds, doms }) => [
    m.seq.id, h.orf_name, orf.frame ?? '', orf.length_aa ?? '', ref.accession, cds.product,
    cds.protein_id || '', h.pident, h.q_cov, h.s_cov, h.evalue, h.bit_score, doms.join('; '),
  ].map(v => String(v).replace(/[\t\r\n]+/g, ' ')).join('\t'));
  return [cols.join('\t')].concat(lines).join('\n') + '\n';
}

// ── Export ─────────────────────────────────────────────────────────────────

function _wireExport() {
  const menu = document.getElementById('syn-export-menu');
  if (!menu) return;
  menu.querySelector('[data-menu-toggle]').addEventListener('click', e => {
    e.stopPropagation();
    document.querySelectorAll('.vq-menu.open').forEach(x => { if (x !== menu) x.classList.remove('open'); });
    menu.classList.toggle('open');
  });
  menu.querySelectorAll('[data-fmt]').forEach(item => item.addEventListener('click', e => {
    e.stopPropagation();
    menu.classList.remove('open');
    const cur  = _current();
    const name = `synteny_${cur.ref.accession}`;
    const svg  = document.querySelector('#syn-map svg');
    const fmt  = item.dataset.fmt;
    if (fmt === 'tsv') { VQ.downloadText(_pairsTSV(cur.ref, cur.members), name + '.tsv'); return; }
    if (!svg) return;
    if (fmt === 'png') VQ.exportPNG(svg, name + '.png');
    if (fmt === 'svg') VQ.exportSVG(svg, name + '.svg');
    if (fmt === 'pdf') VQ.exportPDF(svg, `Synteny — ${cur.ref.organism || cur.ref.accession}`);
  }));
}

window.vqInitSynteny = vqInitSynteny;
})();
