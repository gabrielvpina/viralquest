/* This file's private helpers live in an IIFE so they don't collide across sections. */
(function () {
'use strict';

/* ============================================================
   section_synteny.js — Synteny with the BLASTn reference genome
   One synteny cluster (one reference accession) at a time, picked from a
   list grouped by organism, so the segments of a segmented virus sit
   together:
     • summary card — reference, taxonomy tags, KPIs (genome covered, CDS
       matched, canonical ORFs, identity) and one row per member sequence
       with its span on the reference; links to the other segments
     • map (gggenomes-style) — the reference genome and the member
       sequences as stacked tracks on the reference axis, genes coloured by
       reference gene family (CDS, or mature peptide), straight link bands
       between tracks: reference ↔ member (ORF ↔ CDS hits) and member ↔
       member (contig ↔ contig ORF blastp). Only canonical ORFs — those that
       matched a reference CDS — are drawn unless asked otherwise
     • a table of the ORF ↔ CDS pairs
   Data: report.synteny = { references: {acc: RefGenome}, clusters:
   [{reference, members}] } and sequence.synteny (see viralquest/synteny.py).
   ORF arrows, domain lanes and domain colours come from the viewer via
   VQ.genome, so both tabs look alike.
   ============================================================ */

const _SY = {
  refs:       {},
  clusters:   [],
  seqById:    {},
  idx:        0,
  showAll:    false,   // also draw the ORFs with no reference match (greyed)
  identShade: false,   // links shaded by aa identity instead of a flat tint
  refAxis:    false,   // tracks on the reference axis; default: each on its own, left-aligned
};


// HMM banks used only to call a sequence viral; the pairs table lists Pfam.
const _FILTER_DBS = new Set(['RVDB', 'Vfam', 'EggNOG']);

const _fmt = n => (n == null || isNaN(n)) ? '—' : Math.round(n).toLocaleString('en-US');

// ── Init ───────────────────────────────────────────────────────────────────

function vqInitSynteny(report, mountId) {
  const el = document.getElementById(mountId || 'section-synteny');
  if (!el) return;
  const syn = (report && report.synteny) || {};
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
  const nOrg = new Set(_SY.clusters.map(c => _orgKey(_SY.refs[c.reference]))).size;
  el.innerHTML = `
    <div class="vq-section-header">
      <div class="vq-section-title">
        Synteny
        <span class="vq-count-label">${n} reference${n !== 1 ? 's' : ''} · ${nOrg} virus${nOrg !== 1 ? 'es' : ''}</span>
      </div>
      <div class="vq-section-actions">
        <label class="vq-filter-check"><input type="checkbox" id="syn-show-all"> Non-canonical ORFs</label>
        <label class="vq-filter-check"><input type="checkbox" id="syn-refaxis"> Align tracks to reference</label>
        <label class="vq-filter-check"><input type="checkbox" id="syn-ident"> Shade links by identity</label>
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
    _SY.showAll = e.target.checked; _render();
  });
  document.getElementById('syn-ident').addEventListener('change', e => {
    _SY.identShade = e.target.checked; _drawMap();
  });
  document.getElementById('syn-refaxis').addEventListener('change', e => {
    _SY.refAxis = e.target.checked; _drawMap();
  });
  _wireExport();
  _show(0);
  VQ.redrawOnResize(el, () => _drawMap());
}

const _orgKey = ref => (ref.organism || ref.title || ref.accession).toLowerCase();

/* Options grouped by organism (segments of one virus together); the value is
   always the cluster index, the filter only hides. */
function _fillSelect(query) {
  const sel = document.getElementById('syn-select');
  const q   = query.trim().toLowerCase();
  const groups = new Map();
  _SY.clusters.forEach((c, i) => {
    const ref = _SY.refs[c.reference];
    const hay = [ref.organism, ref.title, ref.accession, ...c.members].join(' ').toLowerCase();
    if (q && !hay.includes(q)) return;
    const key = _orgKey(ref);
    if (!groups.has(key)) groups.set(key, { label: ref.organism || ref.title, opts: [] });
    const n = c.members.length;
    groups.get(key).opts.push(`<option value="${i}">${VQ.esc(
      `${ref.accession} · ${_fmt(ref.length)} nt · ${n} seq${n !== 1 ? 's' : ''}`)}</option>`);
  });
  sel.innerHTML = groups.size
    ? [...groups.values()].map(g =>
        `<optgroup label="${VQ.esc(g.label)}">${g.opts.join('')}</optgroup>`).join('')
    : '<option disabled>No reference matches the filter</option>';
  if (groups.size && !sel.querySelector(`option[value="${_SY.idx}"]`)) {
    _show(+sel.querySelector('option').value);
  } else {
    sel.value = String(_SY.idx);
  }
}

function _step(dir) {
  const vals = [...document.querySelectorAll('#syn-select option')]
    .filter(o => !o.disabled).map(o => +o.value);
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

/* Contig position → reference position, from the placement in synteny.py. */
function _toRef(syn, x) { return syn.flipped ? syn.offset - x : x + syn.offset; }
function _ivToRef(syn, a, b) {
  const p = _toRef(syn, a), q = _toRef(syn, b);
  return [Math.min(p, q), Math.max(p, q)];
}

/* Members of the current cluster, left to right along the reference. */
function _current() {
  const c   = _SY.clusters[_SY.idx];
  const ref = _SY.refs[c.reference];
  const members = c.members.map(id => _SY.seqById[id]).filter(s => s && s.synteny).map(seq => {
    const syn  = seq.synteny;
    const hits = syn.hits || [];
    const names = new Set(hits.map(h => h.orf_name));
    const bits  = hits.reduce((a, h) => a + h.bit_score, 0);
    return {
      seq, syn, hits,
      orfs:   (seq.orfs || []).filter(o => names.has(o.name)),
      other:  (seq.orfs || []).filter(o => !names.has(o.name)),
      span:   _ivToRef(syn, 0, seq.length || 0),
      ident:  bits ? hits.reduce((a, h) => a + h.pident * h.bit_score, 0) / bits : null,
    };
  }).sort((a, b) => a.span[0] - b.span[0] || (b.seq.length || 0) - (a.seq.length || 0));
  return { c, ref, members };
}

/* Total length of the union of [start, end) intervals. */
function _unionLen(ivs) {
  let total = 0, cur = null;
  [...ivs].sort((a, b) => a[0] - b[0]).forEach(([a, b]) => {
    if (!cur || a > cur[1]) { if (cur) total += cur[1] - cur[0]; cur = [a, b]; }
    else cur[1] = Math.max(cur[1], b);
  });
  return cur ? total + cur[1] - cur[0] : 0;
}

/* Domains listed for an ORF in the pairs table: Pfam (architecture) only. */
function _shownDomains(orf) {
  return VQ.genome.bestDomainPerDatabase(orf.domains || []).filter(d => !_FILTER_DBS.has(d.database));
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

/* Ribbon / identity shading: 50 % → faint, 100 % → strong. */
const _identOpacity = p => 0.14 + 0.56 * Math.max(0, Math.min(1, (p - 50) / 50));

// ── Render ─────────────────────────────────────────────────────────────────

function _render() {
  const body = document.getElementById('syn-body');
  if (!body) return;
  const cur = _current();
  body.innerHTML = `
    ${_summaryCard(cur)}
    <div class="vq-card vq-syn-map"><div class="vq-genome-wrap" id="syn-map"></div></div>
    <div class="vq-card vq-syn-pairs">${_pairsTable(cur.ref, cur.members)}</div>`;

  body.querySelectorAll('[data-jump]').forEach(a => a.addEventListener('click', e => {
    e.preventDefault();
    VQ.jumpToViewer(a.dataset.jump);
  }));
  body.querySelectorAll('[data-cluster]').forEach(a => a.addEventListener('click', e => {
    e.preventDefault();
    _show(+a.dataset.cluster);
  }));
  _drawMap();
}

function _drawMap() {
  const wrap = document.getElementById('syn-map');
  if (!wrap || !_SY.clusters.length) return;
  wrap.innerHTML = '';
  wrap.appendChild(_syntenySVG(_current(), wrap.clientWidth));
}

// ── Summary card ───────────────────────────────────────────────────────────

function _summaryCard({ ref, members }) {
  const esc = VQ.esc;
  const allHits    = members.flatMap(m => m.hits);
  const matchedCds = new Set(allHits.map(h => h.cds_index));
  const covered    = _unionLen(allHits.map(h => h.ref_nt));
  const bits       = allHits.reduce((a, h) => a + h.bit_score, 0);
  const ident      = bits ? allHits.reduce((a, h) => a + h.pident * h.bit_score, 0) / bits : null;
  const nCanon     = members.reduce((a, m) => a + m.orfs.length, 0);
  const nOrfs      = members.reduce((a, m) => a + m.orfs.length + m.other.length, 0);

  // Taxonomy of the reference: the members' BLASTn-hit taxonomy names it.
  const tax = members.map(m => (m.seq.blastn_taxonomy || {}).taxonomy || m.seq.taxonomy)
    .find(t => t && (t.family || t.genus)) || {};
  const refUrl = VQ.ncbiUrl(ref.accession, 'nuccore');

  // Other references of the same organism (segments of a segmented virus).
  const siblings = _SY.clusters
    .map((c, i) => ({ c, i, r: _SY.refs[c.reference] }))
    .filter(x => x.i !== _SY.idx && _orgKey(x.r) === _orgKey(ref));

  const tags = [
    refUrl ? `<a class="vq-tag" href="${esc(refUrl)}" target="_blank" rel="noopener noreferrer"
                 title="Open the GenBank record at NCBI">${esc(ref.accession)} ↗</a>`
           : `<span class="vq-tag">${esc(ref.accession)}</span>`,
    tax.family ? `<span class="vq-tag vq-tag--plain">${esc(tax.family)}</span>` : '',
    tax.genus  ? `<span class="vq-tag vq-tag--plain"><i>${esc(tax.genus)}</i></span>` : '',
    tax.genome ? `<span class="vq-tag vq-tag--plain">${esc(tax.genome)}</span>` : '',
  ].join('');

  const chip = VQ.statChip;
  const kpis = [
    chip('Reference', `${_fmt(ref.length)} nt`, 'accent',
         `${ref.cds.length} CDS${ref.peptides.length ? ` · ${ref.peptides.length} mature peptides` : ''}`),
    chip('Genome covered', ref.length ? `${(100 * covered / ref.length).toFixed(1)}%` : '—', 'success',
         `${_fmt(covered)} nt aligned by ORFs`),
    chip('CDS matched', `${matchedCds.size} / ${ref.cds.length}`, '',
         matchedCds.size < ref.cds.length ? `${ref.cds.length - matchedCds.size} without a matching ORF` : 'every CDS matched'),
    chip('Canonical ORFs', `${nCanon} / ${nOrfs}`, '',
         `${nOrfs - nCanon} non-canonical ${_SY.showAll ? 'shown in grey' : 'hidden'}`),
    chip('aa identity', ident != null ? `${ident.toFixed(1)}%` : '—', '', 'bit-score weighted mean'),
  ].join('');

  const rows = members.map(m => {
    const L  = ref.length || 1;
    const pc = v => `${(100 * Math.max(0, Math.min(L, v)) / L).toFixed(2)}%`;
    const span  = `left:${pc(m.span[0])};width:calc(${pc(m.span[1])} - ${pc(m.span[0])})`;
    const marks = m.hits.map(h =>
      `<span class="vq-syn-span__hit" style="left:${pc(h.ref_nt[0])};width:calc(${pc(h.ref_nt[1])} - ${pc(h.ref_nt[0])});opacity:${_identOpacity(h.pident) + 0.25}"></span>`
    ).join('');
    return `
      <div class="vq-syn-member">
        <a class="vq-syn-member__id" href="#" data-jump="${esc(m.seq.id)}" title="Open in the Sequence Viewer">${esc(m.seq.id)}</a>
        <span class="vq-syn-member__meta">${_fmt(m.seq.length)} nt</span>
        <span class="vq-badge ${m.syn.flipped ? 'vq-syn-badge--rc' : 'vq-badge--cluster'}"
              title="${m.syn.flipped ? 'Reverse-complemented relative to the reference (drawn flipped)' : 'Same orientation as the reference'}">
          ${m.syn.flipped ? '⇄ rev. comp.' : '→ forward'}</span>
        <span class="vq-syn-span" title="Placed at ${_fmt(Math.max(0, m.span[0]) + 1)}–${_fmt(Math.min(L, m.span[1]))} on the reference">
          <span class="vq-syn-span__extent" style="${span}"></span>${marks}
        </span>
        <span class="vq-syn-member__meta">${m.hits.length ? `${m.orfs.length} ORF${m.orfs.length !== 1 ? 's' : ''} · ${m.ident.toFixed(1)}%` : 'no matching ORF'}</span>
      </div>`;
  }).join('');

  return `
    <div class="vq-card vq-syn-head">
      <div class="vq-syn-head__top">
        <div class="vq-syn-head__id">
          <div class="vq-syn-head__eyebrow">Reference genome · best viral BLASTn hit</div>
          <div class="vq-syn-head__title">${esc(ref.organism || ref.title)}</div>
          <div class="vq-syn-head__sub">${esc(ref.title)}</div>
          <div class="vq-tag-row">${tags}</div>
        </div>
        ${siblings.length ? `
        <div class="vq-syn-head__sib">
          <div class="vq-syn-head__eyebrow">Other references of this virus</div>
          <div class="vq-tag-row">${siblings.map(x => `
            <a class="vq-tag" href="#" data-cluster="${x.i}"
               title="${esc(x.r.title)}">${esc(x.r.accession)} · ${_fmt(x.r.length)} nt</a>`).join('')}</div>
        </div>` : ''}
      </div>
      <div class="vq-stats-row vq-stats-row--top vq-syn-kpis">${kpis}</div>
      <div class="vq-syn-members">
        <div class="vq-syn-head__eyebrow">Sequences on this reference (${members.length})</div>
        ${rows}
      </div>
    </div>`;
}

// ── SVG (gggenomes-style: stacked tracks + curved links) ──────────────────
//
// Every track (the reference, then the member sequences) is drawn on the
// reference coordinate axis — members flipped / shifted by their placement.
// A link is drawn for EVERY blastp HSP between two drawn genes: reference
// CDS ↔ member ORF, and member ORF ↔ member ORF (the contig blastp of
// synteny.py; older reports derive these from the shared reference). Each
// link leaves its upper gene exactly at the HSP's start/end, runs straight
// behind any track in between, bends in the gaps and lands on the lower
// gene at that HSP's start/end there. Links are near-invisible until hovered.
// Genes are coloured by "gene family": the reference CDS, or its mature
// peptides when it has them, so a polyprotein reads as capsid · E · NS1…

const GENE_PALETTE = [
  '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02',
  '#1f78b4', '#a6761d', '#b2abd2', '#fb8072', '#80b1d3', '#b3de69',
];
const NA_FILL     = '#d9d2c3';   // part of a gene with no family (gap between peptides, …)
const GENE_STROKE = '#3b3b3b';
const LINK_OPACITY = 0.1;        // resting opacity — links show on hover

/* Gene families of a reference: one per mature peptide of a CDS that has them,
   else one per CDS. Coloured in genome order. */
function _families(ref) {
  const fams = [];
  ref.cds.forEach((c, i) => {
    const peps = (ref.peptides || []).filter(p => p.parent === i).sort((a, b) => a.start - b.start);
    if (peps.length) {
      peps.forEach(p => fams.push({ cds: i, start: p.start, end: p.end,
        name: _shortPeptide(p.product), full: p.product }));
    } else {
      fams.push({ cds: i, start: c.start, end: c.end, name: c.product, full: c.product });
    }
  });
  fams.sort((a, b) => a.start - b.start);
  // Same name twice (e.g. three "hypothetical protein" CDS) → number them.
  const seen = {};
  fams.forEach(f => { seen[f.name] = (seen[f.name] || 0) + 1; });
  const k = {};
  fams.forEach((f, i) => {
    if (seen[f.name] > 1) { k[f.name] = (k[f.name] || 0) + 1; f.name = `${f.name} ${k[f.name]}`; }
    f.color = GENE_PALETTE[i % GENE_PALETTE.length];
  });
  return fams;
}

/* Family sections of the reference range [r0, r1] of CDS `cds`, carried onto
   the displayed range [d0, d1] (linear). Returns [{f, r:[a,b], d:[a,b]}]. */
function _sections(fams, cds, r0, r1, d0, d1) {
  const span = r1 - r0 || 1;
  const at = x => d0 + (x - r0) / span * (d1 - d0);
  return fams.filter(f => f.cds === cds && f.end > r0 && f.start < r1).map(f => {
    const a = Math.max(r0, f.start), b = Math.min(r1, f.end);
    return { f, r: [a, b], d: [at(a), at(b)] };
  });
}

/* Order the members so consecutive tracks share as much reference as possible. */
function _trackOrder(ref, members) {
  const ov = (a, b) => Math.max(0, Math.min(a[1], b[1]) - Math.max(a[0], b[0]));
  const left = [...members];
  const out  = [];
  let prev = [0, ref.length];
  while (left.length) {
    let best = 0, score = -1;
    left.forEach((m, i) => {
      const s = ov(prev, m.span);
      if (s > score || (s === score && m.span[0] < left[best].span[0])) { best = i; score = s; }
    });
    const m = left.splice(best, 1)[0];
    out.push(m);
    prev = m.span;
  }
  return out;
}

const GUT = 172, PAD_R2 = 24;
const GENE_H = 16, LANE_GAP = 5;
const LINK_H = 58, REF_LABEL_H = 48;

function _syntenySVG({ ref, members }, containerWidth) {
  const esc = VQ.esc;
  const W   = Math.max(containerWidth || 0, 760);
  const fams = _families(ref);
  const order = _trackOrder(ref, members);

  // Two layouts. Reference axis: every track at its placement on the
  // reference. Own coordinates (default): every track starts at the left
  // edge, on one shared nt scale, so links stretch and bend between the
  // aligned stretches. Either way a member keeps the reference's orientation.
  // Gene / HSP positions stay in reference coordinates; each track's px()
  // turns them into pixels.
  let lo = 0, hi = ref.length;
  order.forEach(m => { lo = Math.min(lo, m.span[0]); hi = Math.max(hi, m.span[1]); });
  const maxLen = Math.max(ref.length, ...order.map(m => m.span[1] - m.span[0]));
  const X = _SY.refAxis
    ? d3.scaleLinear([lo, hi], [GUT, W - PAD_R2])
    : d3.scaleLinear([0, maxLen], [GUT, W - PAD_R2]);
  const pxOf = span => _SY.refAxis ? (v => X(v)) : (v => X(v - span[0]));

  // ── Tracks: genes (display coords on the reference axis) packed in lanes ──
  const refGenes = ref.cds.map((c, i) => ({
    start: c.start, end: c.end, strand: c.strand, key: `ref__${i}`, c, i,
    secs: _sections(fams, i, c.start, c.end, c.start, c.end),
  }));
  const tracks = [{ kind: 'ref', ref, genes: _lanes(refGenes), span: [0, ref.length],
                   px: pxOf([0, ref.length]) }];

  order.forEach(m => {
    const hitsBy = {};
    m.hits.forEach(h => { (hitsBy[h.orf_name] = hitsBy[h.orf_name] || []).push(h); });
    const gene = (o, canon) => {
      const [a, b] = _ivToRef(m.syn, o.start_position, o.stop_position);
      const secs = canon ? (hitsBy[o.name] || []).flatMap(h => {
        const d = _ivToRef(m.syn, h.orf_nt[0], h.orf_nt[1]);
        return _sections(fams, h.cds_index, h.ref_nt[0], h.ref_nt[1], d[0], d[1]);
      }) : [];
      return {
        start: a, end: b, orf: o, canon, key: _linkKey(m, o.name),
        strand: m.syn.flipped ? (o.strand === '-' ? '+' : '-') : o.strand, secs,
      };
    };
    const top   = _lanes(m.orfs.map(o => gene(o, true)));
    const nTop  = Math.max(0, ...top.map(g => g.lane + 1));
    const below = _SY.showAll
      ? _lanes(m.other.map(o => gene(o, false))).map(g => ({ ...g, lane: g.lane + nTop }))
      : [];
    tracks.push({ kind: 'seq', m, genes: top.concat(below), span: m.span, px: pxOf(m.span) });
  });

  // ── Vertical layout ──
  let y = 22 + REF_LABEL_H;
  tracks.forEach((t, i) => {
    if (i) y += LINK_H;
    const nL = Math.max(1, ...t.genes.map(g => g.lane + 1));
    t.top = y;
    t.laneY = [];
    for (let l = 0; l < nL; l++) { t.laneY.push(y); y += GENE_H + (l < nL - 1 ? LANE_GAP : 0); }
    t.bottom = y;
  });
  const axisY   = y + 26;
  const legendY = axisY + 34;
  const legend  = _legendItems(fams);
  const colW    = Math.min(300, Math.max(140, 40 + 6 * Math.max(...legend.map(it => it.label.length))));
  const perRow  = Math.max(1, Math.floor((W - GUT - PAD_R2) / colW));
  const totalH  = legendY + Math.ceil(legend.length / perRow) * 18 + 16;

  const svg = d3.create('svg')
    .attr('class', 'vq-genome-svg vq-genome-svg--fluid vq-syn-svg')
    .attr('viewBox', `0 0 ${W} ${totalH}`)
    .attr('preserveAspectRatio', 'xMinYMin meet')
    .style('width', '100%').style('height', 'auto');
  svg.append('rect').attr('width', W).attr('height', totalH).attr('fill', 'var(--vq-surface)');
  const defs = svg.append('defs');
  const uid  = `syn${Math.random().toString(36).slice(2, 8)}`;

  // Links first, so the tracks paint over the stretches that run behind them.
  const linkG = svg.append('g').attr('class', 'vq-syn-links');
  _collectLinks(tracks, fams, ref, esc).forEach(l => _drawLink(linkG, svg, tracks, l));

  tracks.forEach((t, i) => _drawTrack(svg, defs, `${uid}_${i}`, t, esc));
  _drawRefLabels(svg, tracks[0], fams);

  // ── Axis ──
  const kb = (_SY.refAxis ? hi - lo : maxLen) > 5000;
  const axisG = svg.append('g').attr('transform', `translate(0, ${axisY})`)
    .call(d3.axisBottom(X).ticks(Math.min(12, Math.max(4, Math.round((W - GUT) / 110))))
      .tickFormat(d => kb ? `${d3.format(',')(d / 1000)}k` : d3.format(',')(d)));
  axisG.select('.domain').attr('stroke', 'var(--vq-text-2)');
  axisG.selectAll('.tick line').attr('stroke', 'var(--vq-text-2)');
  axisG.selectAll('.tick text').style('font-size', '10px').style('fill', 'var(--vq-text-2)');
  svg.append('text').attr('x', GUT - 12).attr('y', axisY + 14).attr('text-anchor', 'end')
    .attr('font-size', 9.5).attr('fill', 'var(--vq-text-3)')
    .text(_SY.refAxis ? 'reference position' : 'position from track start');

  _drawLegend(svg, legend, legendY, perRow, colW);
  return svg.node();
}

function _drawTrack(svg, defs, id, t, esc) {
  const X = t.px;
  const G = VQ.genome;
  const isRef = t.kind === 'ref';
  const g = svg.append('g').attr('class', 'vq-syn-track');
  const yMid = t.laneY[0] + GENE_H / 2;

  // Name in the left gutter.
  const name = isRef ? t.ref.accession : t.m.seq.id;
  const sub  = isRef ? `reference · ${_fmt(t.ref.length)} nt`
    : `${_fmt(t.m.seq.length)} nt${t.m.syn.flipped ? ' · rev. comp.' : ''}`;
  const lbl = g.append('text').attr('x', GUT - 14).attr('y', yMid - 1).attr('text-anchor', 'end');
  lbl.append('tspan').attr('font-size', 11.5).attr('font-weight', 600)
    .attr('fill', isRef ? 'var(--vq-primary)' : 'var(--vq-text)')
    .text(name.length > 24 ? name.slice(0, 23) + '…' : name);
  lbl.append('tspan').attr('x', GUT - 14).attr('dy', 13).attr('font-size', 9.5)
    .attr('fill', 'var(--vq-text-3)').text(sub);
  if (!isRef) lbl.style('cursor', 'pointer').on('click', () => VQ.jumpToViewer(t.m.seq.id))
    .append('title').text('Open in the Sequence Viewer');

  // Backbone with end caps.
  const x0 = X(t.span[0]), x1 = X(t.span[1]);
  g.append('line').attr('x1', x0).attr('x2', x1).attr('y1', yMid).attr('y2', yMid)
    .attr('stroke', GENE_STROKE).attr('stroke-width', 1.2);
  [x0, x1].forEach((x, k) => g.append('rect')
    .attr('x', k ? x : x - 5).attr('width', 5).attr('y', yMid - GENE_H / 2 - 3).attr('height', GENE_H + 6)
    .attr('fill', isRef ? 'var(--vq-primary-light)' : 'var(--vq-text-3)').attr('opacity', 0.8));

  t.genes.forEach((gn, n) => {
    const gy  = t.laneY[gn.lane];
    const gx1 = X(gn.start), gx2 = Math.max(X(gn.end), gx1 + 4);
    const pts = G.orfArrowPoints(gx1, gx2, gy, GENE_H, gn.strand);
    const cid = `${id}_${n}`;
    const gg  = g.append('g').attr('opacity', isRef || gn.canon ? 1 : 0.55);
    defs.append('clipPath').attr('id', cid).append('polygon').attr('points', pts);
    gg.append('polygon').attr('points', pts).attr('fill', NA_FILL);
    const body = gg.append('g').attr('clip-path', `url(#${cid})`);
    gn.secs.forEach(s => body.append('rect')
      .attr('x', X(Math.min(s.d[0], s.d[1]))).attr('width', Math.abs(X(s.d[1]) - X(s.d[0])))
      .attr('y', gy).attr('height', GENE_H).attr('fill', s.f.color));
    gg.append('polygon').attr('points', pts)
      .attr('class', 'vq-syn-gene').attr('data-link', gn.key)
      .attr('fill', 'transparent').attr('stroke', GENE_STROKE).attr('stroke-width', 0.9)
      .style('cursor', isRef ? 'default' : 'pointer')
      .on('mouseenter', () => _highlight(svg, [gn.key], true))
      .on('mousemove', evt => VQ.tooltipShow(_geneTip(t, gn, esc), evt))
      .on('mouseleave', () => { _highlight(svg, [gn.key], false); VQ.tooltipHide(); })
      .on('click', () => { if (!isRef) VQ.jumpToViewer(t.m.seq.id); });
  });
}

function _geneTip(t, gn, esc) {
  const fams = [...new Set(gn.secs.map(s => s.f.full))];
  if (t.kind === 'ref') {
    const c = gn.c;
    return `
      <div class="vq-tooltip__title">${esc(c.product)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Protein</span><span>${esc(c.protein_id || '—')}</span>
        <span class="vq-tooltip__key">Position</span><span>${_fmt(c.start + 1)}–${_fmt(c.end)} (${c.strand})</span>
        <span class="vq-tooltip__key">Length</span><span>${_fmt(c.length_aa)} aa</span>
      </div>
      ${fams.length > 1 ? `<div style="margin-top:5px;font-size:10px;color:rgba(255,255,255,.72)">${fams.map(esc).join(' · ')}</div>` : ''}`;
  }
  const o = gn.orf;
  return `
    <div class="vq-tooltip__title">${esc(o.name)}</div>
    <div class="vq-tooltip__row">
      <span class="vq-tooltip__key">Frame</span><span>${esc(o.frame)}</span>
      <span class="vq-tooltip__key">Length</span><span>${_fmt(o.length_aa)} aa</span>
      <span class="vq-tooltip__key">Position</span><span>${_fmt(o.start_position)}–${_fmt(o.stop_position)} (${esc(o.strand)})</span>
      <span class="vq-tooltip__key">Reference</span><span>${gn.canon ? esc(fams.join(' · ') || 'matched') : 'no match (non-canonical)'}</span>
    </div>
    <div style="margin-top:5px;font-size:10px;color:rgba(255,255,255,.55)">Click to open in the Sequence Viewer</div>`;
}

/* Gene-family names above the reference: level when the name fits over its
   gene, else rotated (skipping ones that would collide). Overlapping genes
   sit in lower lanes, so their names are left to the legend. */
function _drawRefLabels(svg, refTrack, fams) {
  const X = refTrack.px;
  const y = refTrack.top - 6;
  const lane0 = new Set(refTrack.genes.filter(g => g.lane === 0).map(g => g.i));
  let lastX = -Infinity;
  fams.filter(f => lane0.has(f.cds)).forEach(f => {
    const x1 = X(f.start), x2 = X(f.end), xm = (x1 + x2) / 2;
    if (x2 - x1 < 4) return;
    const t = svg.append('text').attr('font-size', 10).attr('fill', 'var(--vq-text)');
    if (x2 - x1 > f.name.length * 6 + 12) {
      t.attr('x', xm).attr('y', y).attr('text-anchor', 'middle').text(f.name);
      lastX = x2;
    } else {
      if (xm - lastX < 12) { t.remove(); return; }
      t.attr('transform', `translate(${xm},${y}) rotate(-35)`)
        .text(f.name.length > 18 ? f.name.slice(0, 17) + '…' : f.name);
      lastX = xm;
    }
    t.append('title').text(f.full);
  });
}

// ── Links ──────────────────────────────────────────────────────────────────

/* Every link to draw: {up:{t, key, iv}, down:{t, key, iv}, color, pident, tip}
   where iv is the HSP's [start, end) on that gene in display coordinates. */
function _collectLinks(tracks, fams, ref, esc) {
  const links = [];
  const where = new Map();             // gene key → track index, for drawn genes
  tracks.forEach((t, ti) => t.genes.forEach(g => {
    if (t.kind === 'ref' || g.canon) where.set(g.key, ti);
  }));
  const add = (ka, ia, kb, ib, color, pident, tip) => {
    const ta = where.get(ka), tb = where.get(kb);
    if (ta == null || tb == null || ta === tb) return;
    const [u, d] = ta < tb ? [[ta, ka, ia], [tb, kb, ib]] : [[tb, kb, ib], [ta, ka, ia]];
    // `r`: the link's stretch in reference coordinates, used to route it
    // through the homologous stretch of any track it crosses.
    links.push({ up: { t: u[0], key: u[1], iv: u[2] }, down: { t: d[0], key: d[1], iv: d[2] },
                 r: [Math.min(...ia), Math.max(...ia)], color, pident, tip });
  };

  const members = tracks.slice(1).map(t => t.m);

  // Reference CDS ↔ member ORF: every HSP, per CDS part, split by family.
  members.forEach(m => m.hits.forEach(h => {
    const cds  = ref.cds[h.cds_index];
    const tip  = `
      <div class="vq-tooltip__title">${esc(h.orf_name)} ↔ ${esc(cds.product)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Identity</span><span>${h.pident.toFixed(1)}%</span>
        <span class="vq-tooltip__key">E-value</span><span>${h.evalue.toExponential(1)}</span>
        <span class="vq-tooltip__key">ORF aa</span><span>${h.q_start}–${h.q_end} (${h.q_cov}%)</span>
        <span class="vq-tooltip__key">CDS aa</span><span>${h.s_start}–${h.s_end} (${h.s_cov}%)</span>
      </div>`;
    const segs = (h.segments && h.segments.length)
      ? h.segments : [[h.ref_nt[0], h.ref_nt[1], h.orf_nt[0], h.orf_nt[1]]];
    segs.forEach(([r0, r1, o0, o1]) => {
      const d = _ivToRef(m.syn, o0, o1);
      const secs = _sections(fams, h.cds_index, r0, r1, d[0], d[1]);
      (secs.length ? secs : [{ f: { color: NA_FILL }, r: [r0, r1], d }]).forEach(s =>
        add(`ref__${h.cds_index}`, s.r, _linkKey(m, h.orf_name), s.d, s.f.color, h.pident, tip));
    });
  }));

  // Member ORF ↔ member ORF: every contig blastp HSP, for every pair of members.
  const hasLinks = members.some(m => Array.isArray(m.syn.links));
  const bestHit = (m, orf) => m.hits.filter(h => h.orf_name === orf)
    .sort((p, q) => q.bit_score - p.bit_score)[0];
  for (let i = 0; i < members.length; i++) {
    for (let j = i + 1; j < members.length; j++) {
      const a = members[i], b = members[j];
      const pairs = [];
      if (hasLinks) {
        (a.syn.links || []).filter(l => l.other_seq === b.seq.id).forEach(l => pairs.push(
          { a, b, oa: l.orf_name, ob: l.other_orf, na: l.orf_nt, nb: l.other_nt, l }));
        (b.syn.links || []).filter(l => l.other_seq === a.seq.id).forEach(l => pairs.push(
          { a, b, oa: l.other_orf, ob: l.orf_name, na: l.other_nt, nb: l.orf_nt, l,
            swap: true }));
      } else {
        // Older reports: the reference stretch two hits on one CDS share.
        a.hits.forEach(ha => b.hits.forEach(hb => {
          if (ha.cds_index !== hb.cds_index) return;
          const s = Math.max(ha.ref_nt[0], hb.ref_nt[0]), e = Math.min(ha.ref_nt[1], hb.ref_nt[1]);
          if (e <= s) return;
          const sub = (h, x) => h.orf_nt[0] + (x - h.ref_nt[0]) / ((h.ref_nt[1] - h.ref_nt[0]) || 1) * (h.orf_nt[1] - h.orf_nt[0]);
          pairs.push({ a, b, oa: ha.orf_name, ob: hb.orf_name, derived: true,
            na: [sub(ha, s), sub(ha, e)].sort((p, q) => p - q),
            nb: [sub(hb, s), sub(hb, e)].sort((p, q) => p - q),
            l: { pident: Math.min(ha.pident, hb.pident) } });
        }));
      }
      pairs.forEach(p => {
        const da = _ivToRef(a.syn, p.na[0], p.na[1]);
        const db = _ivToRef(b.syn, p.nb[0], p.nb[1]);
        const ha = bestHit(a, p.oa);
        let secs = ha ? _sections(fams, ha.cds_index, da[0], da[1], db[0], db[1]) : [];
        if (!secs.length) secs = [{ f: { color: NA_FILL }, r: da, d: db }];
        const l = p.l;
        const aa = p.derived ? '' : (p.swap
          ? `<span class="vq-tooltip__key">aa</span><span>${l.s_start ?? '?'}–${l.s_end ?? '?'} ↔ ${l.q_start ?? '?'}–${l.q_end ?? '?'}</span>`
          : `<span class="vq-tooltip__key">aa</span><span>${l.q_start ?? '?'}–${l.q_end ?? '?'} ↔ ${l.s_start ?? '?'}–${l.s_end ?? '?'}</span>`);
        const tip = `
          <div class="vq-tooltip__title">${esc(p.oa)} ↔ ${esc(p.ob)}</div>
          <div class="vq-tooltip__row">
            ${p.derived
              ? `<span class="vq-tooltip__key">Via</span><span>shared reference stretch</span>`
              : `<span class="vq-tooltip__key">Identity</span><span>${l.pident.toFixed(1)}%</span>
                 <span class="vq-tooltip__key">E-value</span><span>${l.evalue.toExponential(1)}</span>${aa}`}
          </div>`;
        secs.forEach(s => add(_linkKey(a, p.oa), s.r, _linkKey(b, p.ob), s.d,
          s.f.color, l.pident ?? 100, tip));
      });
    }
  }
  return links;
}

/* One link as a band: it leaves the upper gene at the HSP's start/end, runs
   straight behind every track in between (at a position eased from one end
   to the other), bends in each gap with vertical tangents, and lands on the
   lower gene at that gene's HSP start/end. */
function _drawLink(g, svg, tracks, l) {
  const U = tracks[l.up.t], D = tracks[l.down.t];
  const gU = U.genes.find(x => x.key === l.up.key);
  const gD = D.genes.find(x => x.key === l.down.key);
  if (!gU || !gD) return;
  const ends = (t, iv) => { const p = [t.px(iv[0]), t.px(iv[1])]; return [Math.min(...p), Math.max(...p)]; };
  const a = ends(U, l.up.iv);
  const b = ends(D, l.down.iv);

  // Vertical levels the band passes straight through.
  const yStart = U.laneY[gU.lane] + GENE_H;      // bottom of the upper gene
  const yEnd   = D.laneY[gD.lane];               // top of the lower gene
  const levels = [{ y0: yStart, y1: U.bottom + 3 }];
  for (let k = l.up.t + 1; k < l.down.t; k++) {
    levels.push({ y0: tracks[k].top - 3, y1: tracks[k].bottom + 3 });
  }
  levels.push({ y0: D.top - 3, y1: yEnd });
  const span = levels[levels.length - 1].y0 - levels[0].y1 || 1;
  const ease = t => t * t * (3 - 2 * t);
  levels.forEach((lv, i) => {
    // A track in between that also covers this reference stretch: pass
    // through its homologous stretch. Otherwise ease from one end to the other.
    const T = tracks[l.up.t + i];
    const ov = T ? Math.min(l.r[1], T.span[1]) - Math.max(l.r[0], T.span[0]) : 0;
    if (i > 0 && i < levels.length - 1 && ov >= 0.5 * (l.r[1] - l.r[0])) {
      [lv.l, lv.r] = ends(T, l.r);
      return;
    }
    const t = ease(Math.max(0, Math.min(1, ((lv.y0 + lv.y1) / 2 - levels[0].y1) / span)));
    lv.l = a[0] + (b[0] - a[0]) * t;
    lv.r = a[1] + (b[1] - a[1]) * t;
  });
  levels[0].l = a[0]; levels[0].r = a[1];
  const last = levels[levels.length - 1];
  last.l = b[0]; last.r = b[1];

  let d = `M${levels[0].l},${levels[0].y0} L${levels[0].l},${levels[0].y1}`;
  for (let i = 1; i < levels.length; i++) {
    const p = levels[i - 1], q = levels[i], ym = (p.y1 + q.y0) / 2;
    d += ` C${p.l},${ym} ${q.l},${ym} ${q.l},${q.y0} L${q.l},${q.y1}`;
  }
  d += ` L${last.r},${last.y1} L${last.r},${last.y0}`;
  for (let i = levels.length - 2; i >= 0; i--) {
    const p = levels[i], q = levels[i + 1], ym = (p.y1 + q.y0) / 2;
    d += ` C${q.r},${ym} ${p.r},${ym} ${p.r},${p.y1} L${p.r},${p.y0}`;
  }
  d += ' Z';

  const keys = [l.up.key, l.down.key];
  const op = _SY.identShade ? 0.03 + 0.17 * Math.max(0, Math.min(1, (l.pident - 50) / 50)) : LINK_OPACITY;
  g.append('path').attr('d', d)
    .attr('class', 'vq-syn-link').attr('data-a', keys[0]).attr('data-b', keys[1])
    .attr('fill', l.color).attr('fill-opacity', op)
    .on('mouseenter', function () { d3.select(this).classed('vq-syn-hl', true); _highlight(svg, keys, true); })
    .on('mousemove', evt => VQ.tooltipShow(l.tip, evt))
    .on('mouseleave', function () {
      d3.select(this).classed('vq-syn-hl', false);
      _highlight(svg, keys, false);
      VQ.tooltipHide();
    });
}

function _legendItems(fams) {
  const items = fams.map(f => ({ color: f.color, label: f.name, title: f.full }));
  items.push({ color: NA_FILL, label: 'no reference family' });
  return items;
}

function _drawLegend(svg, items, y0, perRow, colW) {
  items.forEach((it, i) => {
    const x = GUT + (i % perRow) * colW;
    const y = y0 + Math.floor(i / perRow) * 18;
    svg.append('polygon').attr('points', VQ.genome.orfArrowPoints(x, x + 20, y, 11, '+'))
      .attr('fill', it.color).attr('stroke', GENE_STROKE).attr('stroke-width', 0.8);
    svg.append('text').attr('x', x + 26).attr('y', y + 9.5).attr('font-size', 10.5)
      .attr('fill', 'var(--vq-text-2)')
      .text(it.label.length > 42 ? it.label.slice(0, 41) + '…' : it.label)
      .append('title').text(it.title || it.label);
  });
  svg.append('text').attr('x', GUT - 14).attr('y', y0 + 9.5).attr('text-anchor', 'end')
    .attr('font-size', 10.5).attr('font-weight', 600).attr('fill', 'var(--vq-text-2)').text('Genes');
}

/* Short gene-family name: "nonstructural protein NS5" → "NS5",
   "envelope protein E" → "E", "2K protein" → "2K", "capsid protein" → "capsid". */
const _GENERIC_WORDS = /^(protein|peptide|polypeptide|polyprotein|precursor|glycoprotein|mature|putative)$/i;
function _shortPeptide(product) {
  const words = String(product || '').trim().split(/\s+/).filter(Boolean);
  while (words.length > 1 && _GENERIC_WORDS.test(words[words.length - 1])) words.pop();
  const last = words[words.length - 1] || '';
  if (/^[A-Za-z0-9'_-]{1,8}$/.test(last)) return last;
  return words.join(' ') || product || '';
}

/* One key per (sequence, ORF): ties an ORF to its links. */
const _linkKey = (m, orfName) => VQ.safeId(`${m.seq.id}__${orfName}`);

/* Hovering a gene lights up it, every link touching it and the genes at the
   other ends; hovering a link lights up that link and its two genes. */
function _highlight(svg, keys, on) {
  const genes = new Set(keys);
  if (keys.length === 1) {
    svg.selectAll(`[data-a="${keys[0]}"],[data-b="${keys[0]}"]`).each(function () {
      d3.select(this).classed('vq-syn-hl', on);
      genes.add(this.dataset.a); genes.add(this.dataset.b);
    });
  }
  genes.forEach(k => svg.selectAll(`[data-link="${k}"]`).classed('vq-syn-hl', on));
  svg.classed('vq-syn-dim', on);
}

// ── Pairs table ────────────────────────────────────────────────────────────

function _orfDomainNames(orf) {
  return [...new Set(_shownDomains(orf).map(d => d.target))];
}

function _pairRows(ref, members) {
  return members.flatMap(m => m.hits.map(h => {
    const orf = (m.seq.orfs || []).find(o => o.name === h.orf_name) || {};
    return { m, h, orf, cds: ref.cds[h.cds_index], doms: _orfDomainNames(orf) };
  }));
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
      <td class="vq-td--num">${_fmt(orf.length_aa)}</td>
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
        <th>Protein</th><th>Identity</th><th>ORF cov</th><th>CDS cov</th><th>E-value</th>
        <th>Pfam domains</th>
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
