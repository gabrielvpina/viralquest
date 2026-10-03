/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_quant.js — Section 4: RNA Quantification
   Target-centred: pick one species, family or cross-sample cluster at
   the top; every card below describes that target only — its load in
   each sample, its share of the sample's virome, its level against the
   host housekeeping baseline, what it is made of, its contigs and the
   viruses that co-occur with it.

   Metrics.  Each sample's Salmon index holds that sample's own contigs
   (plus the housekeeping panels / transcriptome), so TPM shares a
   denominator only within one sample.  Across samples the tab defaults
   to reads per million input reads (RPM, from NumReads and the library
   size) and offers TPM relative to the sample's housekeeping median
   (× HK) — both comparable between samples.  Raw TPM stays available,
   flagged as within-sample only.

   Fragments of one virus are summed per sample: Salmon splits reads
   between near-identical contigs, so their counts are parts of one whole.
   ============================================================ */

const _METRICS = {
  rpm: { label: 'RPM',  long: 'viral reads per million input reads — comparable across samples' },
  hk:  { label: '× HK', long: 'TPM relative to the sample’s housekeeping-gene median — comparable across samples' },
  tpm: { label: 'TPM',  long: 'TPM within each sample’s own index — compare within a sample only' },
};
const _GROUP_LABEL = { species: 'species', family: 'family', cluster: 'cluster' };
const _PART_COLORS = [
  '#2563eb', '#16a34a', '#db2777', '#d97706', '#7c3aed', '#0891b2',
  '#dc2626', '#65a30d', '#c026d3', '#0d9488', '#ea580c', '#4f46e5',
];
const _OTHER_COLOR = '#cbd5e1';
const _TOP_PARTS = 10;
const _NO_HIT = 'No BLASTx hit';

let _S = { group: 'species', metric: 'rpm', target: null, q: '', pendingCluster: null };
let _samples  = [];          // [{name, input, mapped, rate, hkMed, viralReads, salmon}]
let _contigs  = [];          // one per quantified viral contig
let _clusterSpecies = new Map();
let _ents     = [];          // entities of the current grouping
let _order    = [];          // sample display order (indices into _samples)

function vqInitQuant(clusters, samples) {
  const el = document.getElementById('section-quant');
  if (!el) return;

  const qs = (samples || []).filter(s => s.salmon);
  if (!qs.length) {
    el.hidden = true;
    document.getElementById('tab-quant')?.style.setProperty('display', 'none');
    return;
  }

  _samples = qs.map(s => ({
    name:   s.sample,
    salmon: s.salmon,
    rate:   s.salmon.mapping_rate ?? null,
    mapped: s.salmon.total_reads ?? null,
    input:  s.salmon.input_reads
            ?? (s.salmon.total_reads && s.salmon.mapping_rate
                ? s.salmon.total_reads / (s.salmon.mapping_rate / 100) : null),
    hkMed:  s.salmon.ref_hk?.median || null,
    viralReads: s.salmon.viral_reads ?? null,
  }));
  const inQuant = new Set(_samples.map(s => s.name));

  const clusterOf = new Map();
  (clusters || []).forEach(c => {
    _clusterSpecies.set(c.gid, c.species || c.gid);
    c.members.forEach(m => clusterOf.set(m.gid, c.gid));
  });

  const NOV = window.vqNovelty;
  const seqs = (typeof VQ_REPORT !== 'undefined' ? VQ_REPORT.sequences : null) || [];
  _contigs = seqs
    .filter(q => inQuant.has(q.sample) && (q.reads != null || q.tpm != null))
    .map(q => {
      const hits = (q.blastx_nr_hits || []).length ? q.blastx_nr_hits : (q.blastx_hits || []);
      return {
        sample:  q.sample, gid: q.gid, id: q._orig_id || q.id, length: q.length || 0,
        reads:   q.reads ?? null, tpm: q.tpm ?? null,
        species: (hits[0] && (hits[0].species || _titleSpecies(hits[0].subject_title))) || _NO_HIT,
        family:  q.taxonomy?.family || 'Unclassified',
        nov:     NOV ? NOV.tier(q) : null,
        cluster: clusterOf.get(q.gid) || null,
      };
    });

  const hasHk  = _samples.some(s => s.hkMed);
  const hasClu = _contigs.some(c => c.cluster);

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">RNA Quantification</div>
        <div class="vq-section-sub">
          Salmon read support in ${_samples.length} sample${_samples.length > 1 ? 's' : ''}
          · choose a species, family or cluster — every card below follows it
        </div>
      </div>
    </div>

    <div class="vq-stats-page">
      <div class="vq-chart-card qt-target" style="min-height:auto">
        <div class="qt-target__row">
          <span class="qt-controls__label">Target</span>
          <div class="vq-toggle" id="qt-group" role="tablist" aria-label="Target type">
            <button class="vq-toggle__btn active" type="button" data-group="species">Species</button>
            <button class="vq-toggle__btn" type="button" data-group="family">Family</button>
            ${hasClu ? `<button class="vq-toggle__btn" type="button" data-group="cluster">Cluster</button>` : ''}
          </div>
          <input class="vq-input vq-input--sm qt-target__search" type="search" id="qt-search"
                 placeholder="Search…" aria-label="Search targets">
          <select class="vq-select qt-target__select" id="qt-select" aria-label="Target"></select>
          <button class="vq-btn vq-btn--sm vq-btn--ghost" type="button" id="qt-prev" title="Previous target">‹</button>
          <button class="vq-btn vq-btn--sm vq-btn--ghost" type="button" id="qt-next" title="Next target">›</button>
        </div>
        <div class="qt-target__row">
          <span class="qt-controls__label">Metric</span>
          <div class="vq-toggle" id="qt-metric" role="tablist" aria-label="Abundance metric">
            <button class="vq-toggle__btn active" type="button" data-metric="rpm">RPM</button>
            ${hasHk ? `<button class="vq-toggle__btn" type="button" data-metric="hk">× HK</button>` : ''}
            <button class="vq-toggle__btn" type="button" data-metric="tpm">TPM</button>
          </div>
          <span class="qt-controls__hint" id="qt-metric-hint"></span>
        </div>
        <div class="qt-target__head" id="qt-head"></div>
      </div>

      <div class="vq-stats-row vq-stats-row--top" id="qt-kpis"></div>

      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Abundance across samples</h3>
          <span class="vq-group__hint">how much of the target each library carries</span>
        </div>
        <div class="vq-chart-card" id="qt-load-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Load per Sample</div>
              <div class="vq-chart-card__sub" id="qt-load-sub"></div>
            </div>
            <button class="vq-btn vq-btn--sm" id="qt-load-csv" type="button"
                    title="Download the per-sample values of the target as CSV">Export CSV</button>
          </div>
          <div class="qt-focus__grid">
            <div class="qt-scroll" id="qt-load"></div>
            <div class="vq-vt-wrap qt-focus__table" id="qt-load-table"></div>
          </div>
        </div>
      </section>

      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Within each sample</h3>
          <span class="vq-group__hint">the target against the rest of the virome and against the host</span>
        </div>
        <div class="qt-pair">
          <div class="vq-chart-card" id="qt-share-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title">Share of the Viral Reads</div>
                <div class="vq-chart-card__sub" id="qt-share-sub"></div>
              </div>
            </div>
            <div class="vq-chart-card__body qt-scroll" id="qt-share"></div>
          </div>
          ${hasHk ? `
          <div class="vq-chart-card" id="qt-hk-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title">Target vs Host Housekeeping</div>
                <div class="vq-chart-card__sub">log₂(target TPM / housekeeping median) per sample · above 0 = out-transcribes the host baseline · grey = rest of the virome</div>
              </div>
            </div>
            <div class="vq-chart-card__body qt-scroll" id="qt-hk"></div>
          </div>` : ''}
        </div>
      </section>

      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Inside the target</h3>
          <span class="vq-group__hint">what the target’s reads are made of, and its contigs</span>
        </div>
        <div class="qt-stack">
          <div class="vq-chart-card" id="qt-parts-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title" id="qt-parts-title">Composition</div>
                <div class="vq-chart-card__sub" id="qt-parts-sub"></div>
              </div>
            </div>
            <div class="vq-chart-card__body qt-scroll" id="qt-parts"></div>
            <div class="ov-legend qt-legend" id="qt-parts-legend"></div>
          </div>
          <div class="vq-chart-card" id="qt-contigs-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title">Contigs</div>
                <div class="vq-chart-card__sub" id="qt-contigs-sub"></div>
              </div>
              <button class="vq-btn vq-btn--sm" id="qt-contigs-csv" type="button"
                      title="Download the target's contigs as CSV">Export CSV</button>
            </div>
            <div class="vq-chart-card__body" style="justify-content:flex-start">
              <div class="vq-vt-wrap qt-contigs-wrap" id="qt-contigs"></div>
            </div>
          </div>
        </div>
      </section>

      ${_samples.length >= 3 ? `
      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Co-occurrence</h3>
          <span class="vq-group__hint">viruses whose abundance rises and falls with the target across samples</span>
        </div>
        <div class="vq-chart-card" id="qt-co-card" style="min-height:auto">
          <div class="vq-chart-card__head">
            <div>
              <div class="vq-chart-card__title">Co-occurring Viruses</div>
              <div class="vq-chart-card__sub" id="qt-co-sub"></div>
            </div>
          </div>
          <div class="vq-chart-card__body qt-scroll" id="qt-co"></div>
        </div>
      </section>` : ''}
    </div>
  `;

  // ── Controls ──
  el.querySelectorAll('#qt-group [data-group]').forEach(b => b.addEventListener('click', () => {
    if (_S.group === b.dataset.group) return;
    _S.group = b.dataset.group;
    el.querySelectorAll('#qt-group [data-group]').forEach(x => x.classList.toggle('active', x === b));
    _S.target = (_S.group === 'cluster' && _S.pendingCluster) || null;
    _S.q = ''; document.getElementById('qt-search').value = '';
    _rebuild();
  }));
  el.querySelectorAll('#qt-metric [data-metric]').forEach(b => b.addEventListener('click', () => {
    _S.metric = b.dataset.metric;
    el.querySelectorAll('#qt-metric [data-metric]').forEach(x => x.classList.toggle('active', x === b));
    VQ.tooltipHide();
    _fillSelect();
    _renderTarget();
  }));
  document.getElementById('qt-search').addEventListener('input', e => {
    _S.q = e.target.value.trim().toLowerCase(); _fillSelect();
  });
  document.getElementById('qt-select').addEventListener('change', e => _setTarget(e.target.value));
  document.getElementById('qt-prev').addEventListener('click', () => _step(-1));
  document.getElementById('qt-next').addEventListener('click', () => _step(1));
  document.getElementById('qt-load-csv').addEventListener('click', () =>
    VQ.downloadText(_loadCsv(), `viralquest_${_S.group}_${_safe(_target()?.label)}_per_sample.csv`));
  document.getElementById('qt-contigs-csv').addEventListener('click', () =>
    VQ.downloadText(_contigsCsv(), `viralquest_${_S.group}_${_safe(_target()?.label)}_contigs.csv`));
  document.getElementById('qt-contigs').addEventListener('click', e => {
    const tr = e.target.closest('tr[data-gid]');
    if (tr) VQ.jumpToViewer(tr.dataset.gid);
  });

  // A cluster picked in General Clusters becomes the target (Cluster grouping).
  document.addEventListener('vq:cluster-select', e => {
    _S.pendingCluster = e.detail;
    if (_S.group === 'cluster' && _ents.some(x => x.key === e.detail)) _setTarget(e.detail);
  });

  _rebuild();
  VQ.redrawOnResize(document.getElementById('qt-load'), _renderTarget);
}

// ── Data ────────────────────────────────────────────────────────────────────

function _titleSpecies(title) {
  const m = /\[([^\]]+)\]\s*$/.exec(title || '');
  return m ? m[1] : null;
}

function _mostCommon(counts) {
  let best = null, n = -1;
  counts.forEach((k, v) => { if (k > n) { best = v; n = k; } });
  return best;
}

function _keyOf(c, group = _S.group) {
  if (group === 'family')  return c.family;
  if (group === 'cluster') return c.cluster;
  return c.species;
}

/* Group the contigs into entities (species / family / cluster). */
function _buildEntities() {
  const map = new Map();
  _contigs.forEach(c => {
    const key = _keyOf(c);
    if (!key) return;                                   // unclustered contig in Cluster mode
    if (!map.has(key)) map.set(key, {
      key,
      label: _S.group === 'cluster' ? (_clusterSpecies.get(key) || key) : key,
      sub:   _S.group === 'cluster' ? key : '',
      fam: new Map(), nov: new Map(), bySample: new Map(), contigs: [],
    });
    const e = map.get(key);
    e.contigs.push(c);
    e.fam.set(c.family, (e.fam.get(c.family) || 0) + 1);
    if (c.nov) e.nov.set(c.nov, (e.nov.get(c.nov) || 0) + 1);
    if (!e.bySample.has(c.sample)) e.bySample.set(c.sample, { reads: 0, tpm: 0, n: 0 });
    const b = e.bySample.get(c.sample);
    b.reads += c.reads || 0; b.tpm += c.tpm || 0; b.n += 1;
  });
  return [...map.values()].map(e => {
    e.family  = _mostCommon(e.fam) || '';
    e.novelty = _mostCommon(e.nov) || null;
    e.reads   = d3.sum([...e.bySample.values()], b => b.reads);
    e.prev    = e.bySample.size;
    return e;
  });
}

/* Value of an entity in a sample for a metric; null = not detected there. */
function _val(e, s, metric = _S.metric) {
  const b = e.bySample.get(s.name);
  if (!b) return null;
  if (metric === 'rpm') return s.input ? b.reads / s.input * 1e6 : null;
  if (metric === 'hk')  return s.hkMed ? b.tpm / s.hkMed : null;
  return b.tpm;
}

function _meanVal(e, metric = _S.metric) {
  return d3.mean(_samples, s => _val(e, s, metric) || 0) || 0;
}

const _target = () => _ents.find(e => e.key === _S.target) || null;
const _ranked = () => _ents.slice().sort((a, b) => _meanVal(b) - _meanVal(a) || b.reads - a.reads);

function _rebuild() {
  _ents  = _buildEntities();
  _order = _samples.map((_, i) => i);
  if (!_target()) _S.target = _ranked()[0]?.key || null;
  _fillSelect();
  _renderTarget();
}

function _fillSelect() {
  const sel = document.getElementById('qt-select');
  if (!sel) return;
  const esc = VQ.esc, unit = _METRICS[_S.metric].label;
  const q = _S.q;
  let list = _ranked().filter(e => !q || `${e.label} ${e.sub} ${e.family}`.toLowerCase().includes(q));
  const cur = _target();
  if (cur && !list.includes(cur)) list = [cur].concat(list);
  sel.innerHTML = list.map(e => `<option value="${esc(e.key)}"${e.key === _S.target ? ' selected' : ''}>${
    esc(`${_trunc(e.label, 48)}${e.sub ? ' · ' + e.sub : ''} — ${e.prev}/${_samples.length} samples · mean ${_fmtVal(_meanVal(e))} ${unit}`)}</option>`).join('')
    || '<option value="">no match</option>';
}

function _setTarget(key) {
  if (!key || key === _S.target || !_ents.some(e => e.key === key)) return;
  _S.target = key;
  VQ.tooltipHide();
  _fillSelect();
  _renderTarget();
}

function _step(dir) {
  const list = [...document.getElementById('qt-select').options].map(o => o.value).filter(Boolean);
  const i = list.indexOf(_S.target);
  if (i >= 0 && list[i + dir]) _setTarget(list[i + dir]);
}

// ── Target rendering ────────────────────────────────────────────────────────

function _renderTarget() {
  const hint = document.getElementById('qt-metric-hint');
  if (hint) {
    hint.textContent = _METRICS[_S.metric].long;
    hint.classList.toggle('is-warn', _S.metric === 'tpm');
  }
  const e = _target();
  _renderHead(e);
  _renderKpis(e);
  _drawLoad(e);
  _drawShare(e);
  _drawHk(e);
  _drawParts(e);
  _renderContigs(e);
  _drawCooccurrence(e);
}

function _renderHead(e) {
  const host = document.getElementById('qt-head');
  if (!host) return;
  const esc = VQ.esc, NOV = window.vqNovelty;
  if (!e) { host.innerHTML = '<div class="vq-empty">No target available for this grouping.</div>'; return; }
  host.innerHTML = `
    <div class="qt-target__name">${esc(e.label)}</div>
    <div class="clu-detail__pills">
      ${e.sub ? `<span class="vq-badge vq-badge--cluster">${esc(e.sub)}</span>` : ''}
      <span class="vq-wf-opt"><span class="vq-wf-opt__k">Type</span>${_GROUP_LABEL[_S.group]}</span>
      ${e.family && _S.group !== 'family' ? `<span class="vq-wf-opt"><span class="vq-wf-opt__k">Family</span>${esc(e.family)}</span>` : ''}
      ${e.novelty && NOV ? `<span class="vq-nov vq-nov--${e.novelty}" title="dominant novelty tier of its contigs">${esc(NOV.label[e.novelty])}</span>` : ''}
    </div>`;
}

function _shareOf(e, s) {
  const tot = d3.sum(_ents, x => x.bySample.get(s.name)?.reads || 0);
  const mine = e.bySample.get(s.name)?.reads || 0;
  const rank = mine > 0
    ? 1 + _ents.filter(x => (x.bySample.get(s.name)?.reads || 0) > mine).length : null;
  const of = _ents.filter(x => (x.bySample.get(s.name)?.reads || 0) > 0).length;
  return { tot, mine, share: tot > 0 ? mine / tot : 0, rank, of };
}

function _renderKpis(e) {
  const host = document.getElementById('qt-kpis');
  if (!host) return;
  if (!e) { host.innerHTML = ''; return; }
  const fmt = VQ.fmtNum, unit = _METRICS[_S.metric].label;
  const vals = _samples.map(s => ({ s, v: _val(e, s) })).filter(d => d.v != null);
  const peak = vals.reduce((a, d) => (!a || d.v > a.v ? d : a), null);
  const hk = _samples.map(s => _val(e, s, 'hk')).filter(v => v != null && v > 0).sort((a, b) => a - b);
  const hkMed = hk.length ? d3.quantileSorted(hk, 0.5) : null;
  const shares = _samples.filter(s => e.bySample.has(s.name)).map(s => _shareOf(e, s).share).sort((a, b) => a - b);
  const shMed = shares.length ? d3.quantileSorted(shares, 0.5) : null;
  host.innerHTML = [
    VQ.statChip('Detected In', `${e.prev}/${_samples.length}`, 'accent', 'samples with Salmon reads'),
    VQ.statChip('Total Reads', fmt(Math.round(e.reads)), '', `${fmt(e.contigs.length)} contig${e.contigs.length === 1 ? '' : 's'}`),
    VQ.statChip(`Peak ${unit}`, peak ? _fmtVal(peak.v) : '—', 'success', peak ? `in ${peak.s.name}` : ''),
    VQ.statChip('Share of Viral Reads', shMed == null ? '—' : _fmtPct(shMed), '', 'median where detected'),
    ...(hkMed != null ? [VQ.statChip('vs Host Baseline', _fmtVal(hkMed, 'hk'), hkMed >= 1 ? 'warn' : '',
      hkMed >= 1 ? 'median · above the housekeeping median' : 'median · below the housekeeping median')] : []),
  ].join('');
}

function _fmtVal(v, metric = _S.metric) {
  if (v == null || isNaN(v)) return '—';
  if (metric === 'hk') return v >= 10 ? v.toFixed(0) + '×' : v >= 0.1 ? v.toFixed(2) + '×' : v.toExponential(1) + '×';
  if (v === 0) return '0';
  if (v >= 1000) return Math.round(v).toLocaleString();
  if (v >= 10) return v.toFixed(1);
  if (v >= 0.1) return v.toFixed(2);
  return v.toExponential(1);
}

function _fmtPct(f) {
  if (f == null) return '—';
  const p = 100 * f;
  return p >= 10 ? p.toFixed(1) + '%' : p >= 0.1 ? p.toFixed(2) + '%' : p > 0 ? p.toExponential(1) + '%' : '0%';
}

function _fmtShort(v) {
  if (v == null || isNaN(v)) return '';
  if (Math.abs(v) < 1e-9) return '0';
  if (v >= 1e6) return (v / 1e6).toFixed(1) + 'M';
  if (v >= 1e3) return (v / 1e3).toFixed(v >= 1e4 ? 0 : 1) + 'k';
  if (v >= 10) return v.toFixed(0);
  if (v >= 1) return v.toFixed(1);
  return v.toFixed(2);
}

function _fmtTpm(v) {
  if (v == null || isNaN(v)) return '—';
  if (v === 0)   return '0';
  if (v >= 1000) return v.toFixed(0);
  if (v >= 10)   return v.toFixed(1);
  return v.toFixed(2);
}

function _trunc(s, n) {
  s = String(s ?? '');
  return s.length > n ? s.slice(0, n - 1) + '…' : s;
}

const _safe = s => String(s || 'target').replace(/[^A-Za-z0-9._-]+/g, '_').slice(0, 60);

function _pxSvg(host, minW, H) {
  const W = Math.max(host.clientWidth || 0, minW);
  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`)
    .style('display', 'block');
  return { svg, W };
}

// ── Load per sample (lollipop + table) ──────────────────────────────────────

function _drawLoad(e) {
  const plot = document.getElementById('qt-load'), table = document.getElementById('qt-load-table');
  if (!plot) return;
  const esc = VQ.esc;
  const sub = document.getElementById('qt-load-sub');
  if (!e) { plot.innerHTML = ''; if (table) table.innerHTML = ''; return; }
  if (sub) sub.textContent = `${_METRICS[_S.metric].label} per sample, log scale · — = not detected in that sample`;
  const cols = _order.map(i => _samples[i]);
  const hasHk = _samples.some(s => s.hkMed);

  table.innerHTML = `
    <table class="vq-table vq-vt">
      <thead><tr><th>Sample</th><th class="vq-vt-num">Contigs</th><th class="vq-vt-num">Reads</th>
        <th class="vq-vt-num">RPM</th><th class="vq-vt-num">TPM</th>${hasHk ? '<th class="vq-vt-num">× HK</th>' : ''}</tr></thead>
      <tbody>${cols.map(s => {
        const b = e.bySample.get(s.name);
        if (!b) return `<tr class="qt-absent"><td>${esc(s.name)}</td><td class="vq-vt-num" colspan="${hasHk ? 5 : 4}"><span class="vq-vt-na">not detected</span></td></tr>`;
        return `<tr><td>${esc(s.name)}</td><td class="vq-vt-num">${b.n}</td>
          <td class="vq-vt-num">${Math.round(b.reads).toLocaleString()}</td>
          <td class="vq-vt-num">${_fmtVal(_val(e, s, 'rpm'), 'rpm')}</td>
          <td class="vq-vt-num">${_fmtTpm(b.tpm)}</td>
          ${hasHk ? `<td class="vq-vt-num">${_fmtVal(_val(e, s, 'hk'), 'hk')}</td>` : ''}</tr>`;
      }).join('')}</tbody>
    </table>`;

  const vals = cols.map(s => ({ s, v: _val(e, s) }));
  const present = vals.filter(d => d.v != null && d.v > 0);
  const padL = 58, padR = 14, padT = 16, padB = 74, H = 290;
  const { svg, W } = _pxSvg(plot, 320, H);
  const x = d3.scaleBand().domain(cols.map(s => s.name)).range([padL, W - padR]).padding(0.35);
  const isHk = _S.metric === 'hk';
  const lo = present.length ? d3.min(present, d => d.v) : 1, hi = present.length ? d3.max(present, d => d.v) : 10;
  const y = d3.scaleLog().domain(isHk ? [Math.min(lo, 0.5) / 1.5, Math.max(hi, 2) * 1.5]
                                      : [Math.max(1e-3, lo / 3), hi * 2]).range([H - padB, padT]).nice();
  const ticks = y.ticks(6).filter(t => Math.abs(Math.log10(t) % 1) < 1e-9);
  (ticks.length >= 2 ? ticks : y.ticks(4)).forEach(t => {
    svg.append('line').attr('class', isHk && Math.abs(t - 1) < 1e-9 ? 'qt-zero' : 'ov-grid')
      .attr('x1', padL).attr('x2', W - padR).attr('y1', y(t)).attr('y2', y(t));
    svg.append('text').attr('class', 'qt-axis').attr('x', padL - 6).attr('y', y(t))
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central').text(_fmtShort(t) + (isHk ? '×' : ''));
  });
  svg.append('text').attr('class', 'qt-axis').attr('transform', `translate(12,${(padT + H - padB) / 2}) rotate(-90)`)
    .attr('text-anchor', 'middle').text(`${_METRICS[_S.metric].label} (log scale)`);
  const base = isHk ? y(1) : H - padB;
  vals.forEach(({ s, v }) => {
    const cx = x(s.name) + x.bandwidth() / 2;
    svg.append('text').attr('class', 'qt-axis')
      .attr('transform', `translate(${cx + 3},${H - padB + 12}) rotate(-40)`).attr('text-anchor', 'end').text(_trunc(s.name, 16));
    if (v == null || v <= 0) {
      svg.append('text').attr('class', 'qt-axis').attr('x', cx).attr('y', H - padB - 6).attr('text-anchor', 'middle').text('—');
      return;
    }
    svg.append('line').attr('class', 'qt-lollipop').attr('x1', cx).attr('x2', cx).attr('y1', base).attr('y2', y(v));
    svg.append('circle').attr('class', 'qt-lollipop__dot').attr('cx', cx).attr('cy', y(v)).attr('r', 6.5)
      .on('mousemove', evt => VQ.tooltipShow(`<div class="vq-tooltip__title">${esc(s.name)}</div>
        <div class="vq-tooltip__row"><span class="vq-tooltip__key">${esc(_METRICS[_S.metric].label)}</span><span>${_fmtVal(v)}</span>
        <span class="vq-tooltip__key">Reads</span><span>${Math.round(e.bySample.get(s.name).reads).toLocaleString()}</span></div>`, evt))
      .on('mouseleave', VQ.tooltipHide);
    svg.append('text').attr('class', 'qt-num').attr('x', cx).attr('y', y(v) - 12).attr('text-anchor', 'middle').text(_fmtShort(v));
  });
}

// ── Share of the viral reads ────────────────────────────────────────────────

function _drawShare(e) {
  const host = document.getElementById('qt-share');
  if (!host) return;
  const esc = VQ.esc;
  const sub = document.getElementById('qt-share-sub');
  const plural = { species: 'species', family: 'families', cluster: 'clusters' }[_S.group];
  if (sub) sub.textContent = `fraction of each sample’s viral reads that belong to the target · rank among the sample’s ${plural}`;
  if (!e) { host.innerHTML = ''; return; }
  const cols = _order.map(i => _samples[i]);
  const rows = cols.map(s => ({ s, ..._shareOf(e, s) }));

  const padL = 130, padR = 120, rowH = 28;
  const { svg, W } = _pxSvg(host, 380, rows.length * rowH + 8);
  const x = d3.scaleLinear().domain([0, 1]).range([padL, W - padR]);
  rows.forEach((r, j) => {
    const y = 4 + j * rowH;
    const g = svg.append('g').attr('class', 'qt-row');
    g.append('rect').attr('class', 'qt-row__bg').attr('x', 0).attr('y', y).attr('width', W).attr('height', rowH);
    g.append('text').attr('class', 'qt-axis-strong').attr('x', 6).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(_trunc(r.s.name, 18));
    g.append('rect').attr('class', 'qt-track').attr('x', padL).attr('y', y + 8).attr('width', W - padR - padL)
      .attr('height', rowH - 16).attr('rx', 5);
    if (r.mine > 0) g.append('rect').attr('class', 'qt-share__bar').attr('x', padL).attr('y', y + 8)
      .attr('width', Math.max(3, x(r.share) - padL)).attr('height', rowH - 16).attr('rx', 5);
    g.append('text').attr('class', 'qt-num').attr('x', W - padR + 10).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(r.mine > 0 ? _fmtPct(r.share) : '—');
    g.append('text').attr('class', 'qt-num qt-num--muted').attr('x', W - 6).attr('y', y + rowH / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(r.rank ? `#${r.rank} of ${r.of}` : '');
    g.on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(r.s.name)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Target reads</span><span>${Math.round(r.mine).toLocaleString()}</span>
          <span class="vq-tooltip__key">Viral reads</span><span>${Math.round(r.tot).toLocaleString()}</span>
          <span class="vq-tooltip__key">Share</span><span>${_fmtPct(r.share)}</span>
          ${r.rank ? `<span class="vq-tooltip__key">Rank</span><span>${r.rank} of ${r.of} ${plural}</span>` : ''}
        </div>`, evt))
     .on('mouseleave', VQ.tooltipHide);
  });
}

// ── Target vs host housekeeping ─────────────────────────────────────────────

function _drawHk(e) {
  const host = document.getElementById('qt-hk');
  if (!host) return;
  const esc = VQ.esc;
  if (!e) { host.innerHTML = ''; return; }
  const cols = _order.map(i => _samples[i]);
  const pts = [];
  _ents.forEach(x => cols.forEach(s => {
    const v = _val(x, s, 'hk');
    if (v != null && v > 0) pts.push({ e: x, s, l: Math.log2(v) });
  }));
  const mine = pts.filter(p => p.e === e);
  if (!mine.length) { host.innerHTML = '<div class="vq-empty">The target has no TPM against a housekeeping baseline.</div>'; return; }

  const padL = 44, padR = 14, padT = 12, padB = 70, H = 300;
  const { svg, W } = _pxSvg(host, 340, H);
  const lo = Math.min(-2, d3.min(pts, p => p.l)), hi = Math.max(2, d3.max(pts, p => p.l));
  const y = d3.scaleLinear().domain([lo, hi]).nice().range([H - padB, padT]);
  const x = d3.scalePoint().domain(cols.map(s => s.name)).range([padL, W - padR]).padding(0.5);
  const step = x.step();

  y.ticks(6).forEach(t => {
    svg.append('line').attr('class', t === 0 ? 'qt-zero' : 'ov-grid').attr('x1', padL).attr('x2', W - padR).attr('y1', y(t)).attr('y2', y(t));
    svg.append('text').attr('class', 'qt-axis').attr('x', padL - 6).attr('y', y(t))
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central').text(t);
  });
  svg.append('text').attr('class', 'qt-axis qt-zero__lbl').attr('x', W - padR - 2).attr('y', y(0) - 5)
    .attr('text-anchor', 'end').text('= host HK median');
  cols.forEach(s => svg.append('text').attr('class', 'qt-axis')
    .attr('transform', `translate(${x(s.name) + 3},${H - padB + 12}) rotate(-40)`)
    .attr('text-anchor', 'end').text(_trunc(s.name, 16)));

  const jit = (k, s) => ((Math.sin((k.length * 31 + s.length) * 12.9898 + k.charCodeAt(0)) * 43758.5453) % 1 + 1) % 1 - 0.5;
  pts.filter(p => p.e !== e).forEach(p => svg.append('circle').attr('class', 'qt-hk-pt qt-hk-pt--rest')
    .attr('cx', x(p.s.name) + jit(p.e.key, p.s.name) * Math.min(step * 0.55, 30)).attr('cy', y(p.l)).attr('r', 2.6)
    .on('mousemove', evt => tip(evt, p)).on('mouseleave', VQ.tooltipHide)
    .on('click', () => _setTarget(p.e.key)));
  if (mine.length > 1) svg.append('path').attr('class', 'qt-hk-line')
    .attr('d', d3.line().x(p => x(p.s.name)).y(p => y(p.l))(mine));
  mine.forEach(p => svg.append('circle').attr('class', 'qt-hk-focus')
    .attr('cx', x(p.s.name)).attr('cy', y(p.l)).attr('r', 6)
    .on('mousemove', evt => tip(evt, p)).on('mouseleave', VQ.tooltipHide));

  function tip(evt, p) {
    VQ.tooltipShow(`
      <div class="vq-tooltip__title">${esc(p.e.label)}</div>
      <div class="vq-tooltip__row">
        <span class="vq-tooltip__key">Sample</span><span>${esc(p.s.name)}</span>
        <span class="vq-tooltip__key">TPM</span><span>${_fmtTpm(p.e.bySample.get(p.s.name).tpm)}</span>
        <span class="vq-tooltip__key">HK median</span><span>${_fmtTpm(p.s.hkMed)} TPM</span>
        <span class="vq-tooltip__key">Ratio</span><span>${_fmtVal(Math.pow(2, p.l), 'hk')} (log₂ ${p.l.toFixed(2)})</span>
        ${p.e !== e ? '<span class="vq-tooltip__key"></span><span>click to make it the target</span>' : ''}
      </div>`, evt);
  }
}

// ── Inside the target: what its reads are made of ───────────────────────────
// family → its species · species → the clusters its contigs fall in ·
// cluster → the species each member resolved to (species agreement).

function _partOf(c) {
  if (_S.group === 'family')  return c.species;
  if (_S.group === 'cluster') return c.species;
  return c.cluster ? `${c.cluster} · ${_clusterSpecies.get(c.cluster) || ''}` : 'not in a cross-sample cluster';
}

function _drawParts(e) {
  const host = document.getElementById('qt-parts');
  if (!host) return;
  const esc = VQ.esc;
  const title = document.getElementById('qt-parts-title');
  const sub = document.getElementById('qt-parts-sub');
  const what = { family: 'Species in the Family', species: 'Clusters of the Species', cluster: 'Species of the Cluster Members' }[_S.group];
  if (title) title.textContent = what;
  if (sub) sub.textContent = {
    family:  'share of the family’s reads per species in each sample',
    species: 'share of the species’ reads per cross-sample cluster in each sample — the same colour in two samples is the same genotype',
    cluster: 'share of the cluster’s reads by each member’s own best-hit species — one colour everywhere means the members agree',
  }[_S.group];
  const legend = document.getElementById('qt-parts-legend');
  if (!e) { host.innerHTML = ''; if (legend) legend.innerHTML = ''; return; }

  const totals = new Map();
  e.contigs.forEach(c => totals.set(_partOf(c), (totals.get(_partOf(c)) || 0) + (c.reads || 0)));
  const top = [...totals.entries()].sort((a, b) => b[1] - a[1]).slice(0, _TOP_PARTS).map(([k]) => k);
  const color = new Map(top.map((k, i) => [k, _PART_COLORS[i % _PART_COLORS.length]]));
  const cols = _order.map(i => _samples[i]).filter(s => e.bySample.has(s.name));

  const padL = 140, padR = 100, rowH = 26;
  const { svg, W } = _pxSvg(host, 520, cols.length * rowH + 8);
  const x = d3.scaleLinear().domain([0, 1]).range([padL, W - padR]);
  cols.forEach((s, j) => {
    const y = 4 + j * rowH;
    const parts = new Map();
    e.contigs.filter(c => c.sample === s.name).forEach(c => {
      const k = top.includes(_partOf(c)) ? _partOf(c) : null;
      parts.set(k, (parts.get(k) || 0) + (c.reads || 0));
    });
    const tot = d3.sum([...parts.values()]);
    const g = svg.append('g').attr('class', 'qt-row');
    g.append('rect').attr('class', 'qt-row__bg').attr('x', 0).attr('y', y).attr('width', W).attr('height', rowH);
    g.append('text').attr('class', 'qt-axis-strong').attr('x', 6).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(_trunc(s.name, 18));
    g.append('rect').attr('class', 'qt-track').attr('x', padL).attr('y', y + 5).attr('width', W - padR - padL)
      .attr('height', rowH - 10).attr('rx', 4);
    let acc = 0;
    if (tot > 0) top.concat([null]).forEach(k => {
      const r = parts.get(k) || 0;
      if (r <= 0) return;
      const f = r / tot;
      g.append('rect').attr('class', 'qt-seg').attr('x', x(acc)).attr('y', y + 5)
        .attr('width', Math.max(0.5, x(acc + f) - x(acc))).attr('height', rowH - 10)
        .attr('fill', k ? color.get(k) : _OTHER_COLOR)
        .on('mousemove', evt => VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(k || 'Other')}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Sample</span><span>${esc(s.name)}</span>
            <span class="vq-tooltip__key">Share</span><span>${_fmtPct(f)} of the target’s reads</span>
            <span class="vq-tooltip__key">Reads</span><span>${Math.round(r).toLocaleString()}</span>
          </div>`, evt))
        .on('mouseleave', VQ.tooltipHide);
      acc += f;
    });
    g.append('text').attr('class', 'qt-num').attr('x', W - 6).attr('y', y + rowH / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(`${Math.round(tot).toLocaleString()} reads`);
  });
  if (legend) legend.innerHTML = top.map(k =>
    `<span class="ov-legend__item"><span class="ov-legend__dot" style="background:${color.get(k)}"></span>${esc(_trunc(k, 48))}</span>`).join('')
    + (totals.size > top.length ? `<span class="ov-legend__item"><span class="ov-legend__dot" style="background:${_OTHER_COLOR}"></span>Other (${totals.size - top.length})</span>` : '');
}

// ── Contigs of the target ───────────────────────────────────────────────────

function _renderContigs(e) {
  const host = document.getElementById('qt-contigs');
  if (!host) return;
  const esc = VQ.esc, NOV = window.vqNovelty;
  const sub = document.getElementById('qt-contigs-sub');
  if (!e) { host.innerHTML = ''; return; }
  const byName = new Map(_samples.map(s => [s.name, s]));
  const rows = e.contigs.slice().sort((a, b) => (b.reads || 0) - (a.reads || 0));
  const hasHk = _samples.some(s => s.hkMed);
  const showSp = _S.group !== 'species', showClu = _S.group !== 'cluster' && rows.some(c => c.cluster);
  if (sub) sub.textContent = `${rows.length} contig${rows.length === 1 ? '' : 's'} · sorted by reads · click a row to open it in the Sequence Viewer`;
  host.innerHTML = `
    <table class="vq-table vq-vt">
      <thead><tr>
        <th>Sample</th><th>Contig</th>${showSp ? '<th>Species</th>' : ''}${showClu ? '<th>Cluster</th>' : ''}
        <th>Novelty</th><th class="vq-vt-num">Length (bp)</th><th class="vq-vt-num">Reads</th>
        <th class="vq-vt-num">RPM</th><th class="vq-vt-num">TPM</th>${hasHk ? '<th class="vq-vt-num">× HK</th>' : ''}
      </tr></thead>
      <tbody>${rows.map(c => {
        const s = byName.get(c.sample);
        const rpm = c.reads != null && s.input ? c.reads / s.input * 1e6 : null;
        const hk = c.tpm != null && s.hkMed ? c.tpm / s.hkMed : null;
        return `<tr data-gid="${esc(c.gid)}" title="Open ${esc(c.id)} in the Sequence Viewer">
          <td>${esc(c.sample)}</td>
          <td class="vq-td--mono">${esc(c.id)}</td>
          ${showSp ? `<td class="vq-vt-species">${esc(c.species)}</td>` : ''}
          ${showClu ? `<td class="vq-td--mono">${c.cluster ? esc(c.cluster) : '<span class="vq-vt-na">—</span>'}</td>` : ''}
          <td>${c.nov && NOV ? `<span class="vq-nov vq-nov--${c.nov}">${esc(NOV.label[c.nov])}</span>` : '<span class="vq-vt-na">—</span>'}</td>
          <td class="vq-vt-num">${c.length.toLocaleString()}</td>
          <td class="vq-vt-num">${c.reads == null ? '—' : Math.round(c.reads).toLocaleString()}</td>
          <td class="vq-vt-num">${_fmtVal(rpm, 'rpm')}</td>
          <td class="vq-vt-num">${_fmtTpm(c.tpm)}</td>
          ${hasHk ? `<td class="vq-vt-num">${_fmtVal(hk, 'hk')}</td>` : ''}
        </tr>`;
      }).join('')}</tbody>
    </table>`;
}

// ── Co-occurrence ───────────────────────────────────────────────────────────
// Pearson correlation of log abundance across samples (absent = 0) between
// the target and every other entity seen in ≥ 2 samples.

function _drawCooccurrence(e) {
  const host = document.getElementById('qt-co');
  if (!host) return;
  const esc = VQ.esc;
  const sub = document.getElementById('qt-co-sub');
  if (!e) { host.innerHTML = ''; return; }
  if (sub) sub.textContent = `Pearson r of log(1 + ${_METRICS[_S.metric].label}) across ${_samples.length} samples (absent = 0) · strongest positive and negative · click a bar to make it the target`;
  if (e.prev < 2) {
    host.innerHTML = '<div class="vq-empty">The target occurs in a single sample — co-occurrence needs at least two.</div>';
    return;
  }
  const vec = x => _samples.map(s => Math.log10(1 + (_val(x, s) || 0)));
  const t = vec(e);
  const corr = (a, b) => {
    const ma = d3.mean(a), mb = d3.mean(b);
    let n = 0, da = 0, db = 0;
    a.forEach((v, i) => { n += (v - ma) * (b[i] - mb); da += (v - ma) ** 2; db += (b[i] - mb) ** 2; });
    return da && db ? n / Math.sqrt(da * db) : null;
  };
  const all = _ents.filter(x => x !== e && x.prev >= 2)
    .map(x => ({ x, r: corr(t, vec(x)), both: _samples.filter(s => x.bySample.has(s.name) && e.bySample.has(s.name)).length }))
    .filter(d => d.r != null);
  const pos = all.filter(d => d.r > 0).sort((a, b) => b.r - a.r).slice(0, 8);
  const neg = all.filter(d => d.r < 0).sort((a, b) => a.r - b.r).slice(0, 4);
  const rows = pos.concat(neg);
  if (!rows.length) { host.innerHTML = '<div class="vq-empty">No other virus varies with the target.</div>'; return; }

  const padL = 300, padR = 180, rowH = 24;
  const { svg, W } = _pxSvg(host, 640, rows.length * rowH + 34 + (neg.length && pos.length ? 10 : 0));
  const x = d3.scaleLinear().domain([-1, 1]).range([padL, W - padR]);
  [-1, -0.5, 0, 0.5, 1].forEach(v => {
    svg.append('line').attr('class', v === 0 ? 'qt-zero' : 'ov-grid').attr('x1', x(v)).attr('x2', x(v))
      .attr('y1', 4).attr('y2', rows.length * rowH + 14);
    svg.append('text').attr('class', 'qt-axis').attr('x', x(v)).attr('y', rows.length * rowH + 28 + (neg.length && pos.length ? 10 : 0))
      .attr('text-anchor', 'middle').text(v);
  });
  rows.forEach((d, j) => {
    const y = 6 + j * rowH + (j >= pos.length && pos.length ? 10 : 0);
    const g = svg.append('g').attr('class', 'qt-row').style('cursor', 'pointer').on('click', () => _setTarget(d.x.key));
    g.append('rect').attr('class', 'qt-row__bg').attr('x', 0).attr('y', y).attr('width', W).attr('height', rowH);
    g.append('text').attr('class', 'qt-co__lbl').attr('x', padL - 10).attr('y', y + rowH / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central')
      .text(_trunc(d.x.label, 40) + (d.x.sub ? ` · ${d.x.sub}` : ''));
    g.append('rect').attr('class', d.r > 0 ? 'qt-co__pos' : 'qt-co__neg')
      .attr('x', Math.min(x(0), x(d.r))).attr('y', y + 5).attr('width', Math.abs(x(d.r) - x(0))).attr('height', rowH - 10).attr('rx', 3);
    g.append('text').attr('class', 'qt-num').attr('x', W - padR + 10).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(`r = ${d.r.toFixed(2)}`);
    g.append('text').attr('class', 'qt-num qt-num--muted').attr('x', W - 6).attr('y', y + rowH / 2)
      .attr('text-anchor', 'end').attr('dominant-baseline', 'central').text(`${d.both} shared`);
    g.on('mousemove', evt => VQ.tooltipShow(`
        <div class="vq-tooltip__title">${esc(d.x.label)}</div>
        <div class="vq-tooltip__row">
          <span class="vq-tooltip__key">Pearson r</span><span>${d.r.toFixed(3)}</span>
          <span class="vq-tooltip__key">Samples with both</span><span>${d.both} of ${_samples.length}</span>
          <span class="vq-tooltip__key">Its samples</span><span>${d.x.prev}</span>
        </div>`, evt))
     .on('mouseleave', VQ.tooltipHide);
  });
}

// ── CSV exports (plain numbers, RFC 4180 quoting) ───────────────────────────

function _csv(head, rows) {
  const cell = v => {
    if (v == null || (typeof v === 'number' && isNaN(v))) return '';
    const s = String(v);
    return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [head.map(cell).join(','), ...rows.map(r => r.map(cell).join(','))].join('\r\n') + '\r\n';
}

const _num = (v, p = 6) => (v == null || isNaN(v)) ? null : +(+v).toPrecision(p);

function _loadCsv() {
  const e = _target();
  if (!e) return _csv(['sample'], []);
  return _csv(
    ['target_type', 'target', 'sample', 'detected', 'contigs', 'reads', 'rpm', 'tpm', 'hk_ratio', 'share_of_viral_reads', 'rank_in_sample'],
    _order.map(i => _samples[i]).map(s => {
      const b = e.bySample.get(s.name), sh = _shareOf(e, s);
      return [_S.group, e.sub || e.label, s.name, b ? 1 : 0, b ? b.n : 0, b ? _num(b.reads) : 0,
              _num(_val(e, s, 'rpm')), b ? _num(b.tpm) : null, _num(_val(e, s, 'hk')), _num(sh.share), sh.rank];
    }));
}

function _contigsCsv() {
  const e = _target();
  if (!e) return _csv(['contig_id'], []);
  const byName = new Map(_samples.map(s => [s.name, s]));
  return _csv(
    ['sample', 'contig_id', 'species', 'family', 'cluster', 'novelty', 'length_bp', 'reads', 'rpm', 'tpm', 'hk_ratio'],
    e.contigs.slice().sort((a, b) => (b.reads || 0) - (a.reads || 0)).map(c => {
      const s = byName.get(c.sample);
      return [c.sample, c.id, c.species, c.family, c.cluster, c.nov, c.length, _num(c.reads),
              c.reads != null && s.input ? _num(c.reads / s.input * 1e6) : null, _num(c.tpm),
              c.tpm != null && s.hkMed ? _num(c.tpm / s.hkMed) : null];
    }));
}

window.vqInitQuant = vqInitQuant;
})();
