/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_clusters.js — Section 2: General (cross-sample) Clusters
   Cluster-centred page: pick a cluster in the explorer table and the
   detail card shows only that cluster — sample presence, members,
   species agreement and the alignment against the representative.
   Below, two across-sample views of every cluster passing the filters:
   a cluster × sample presence matrix and a sample-similarity heatmap
   (UPGMA-ordered), both highlighting the selected cluster.
   Identity / coverage sliders and the other filters act client-side.
   ============================================================ */

const FLOOR_ID  = 90;   // matches report_clusters.FLOOR_IDENTITY
const FLOOR_COV = 70;   // matches report_clusters.FLOOR_COVERAGE

let _clusters = [];
let _samples  = [];          // every loaded sample name, report order
let _seqByGid = new Map();
let _hasTpm   = false;
let _hasRpm   = false;
let _inputReads = new Map();   // sample → library size (reads per million)
let _state = {
  id: FLOOR_ID, cov: FLOOR_COV, minSamples: 2, agree: 'all', family: 'all', novelty: 'all',
  q: '', sort: 'nSamples', dir: -1, sel: null, matrixMode: 'presence',
};
let _view = [];              // clusters passing the filters, table order

function vqInitClusters(clusters, samples) {
  const el = document.getElementById('section-clusters');
  if (!el) return;
  _clusters = clusters || [];

  if (!_clusters.length) {
    el.hidden = true;
    document.getElementById('tab-clusters')?.style.setProperty('display', 'none');
    return;
  }

  _samples = (samples || []).map(s => s.sample).filter(Boolean);
  _clusters.forEach(c => c.members.forEach(m => {
    if (!_samples.includes(m.sample)) _samples.push(m.sample);
  }));
  const seqs = (typeof VQ_REPORT !== 'undefined' ? VQ_REPORT.sequences : null) || [];
  _seqByGid = new Map(seqs.map(s => [s.gid, s]));
  _hasTpm = _clusters.some(c => c.members.some(m => m.tpm != null));
  (samples || []).forEach(s => { if (s.salmon?.input_reads) _inputReads.set(s.sample, s.salmon.input_reads); });
  _hasRpm = _inputReads.size > 0 && _clusters.some(c => c.members.some(m => m.reads != null));
  // RPM compares across samples; TPM only within one, so it is the fallback.
  _state.matrixMode = _hasRpm ? 'rpm' : _hasTpm ? 'tpm' : 'presence';

  // Family / novelty of each cluster come from its representative contig.
  const NOV = window.vqNovelty;
  _clusters.forEach(c => {
    const rep = _seqByGid.get(c.representative);
    c._family  = rep?.taxonomy?.family || '';
    c._novelty = rep && NOV ? NOV.tier(rep) : '';
  });

  const families  = [...new Set(_clusters.map(c => c._family).filter(Boolean))].sort();
  const novelties = NOV ? NOV.order.filter(t => _clusters.some(c => c._novelty === t)) : [];
  const esc = VQ.esc;

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">General Clusters</div>
        <div class="vq-section-sub">
          ${_clusters.length} cross-sample cluster${_clusters.length > 1 ? 's' : ''}
          · sequences shared between ${_samples.length} samples · pick a cluster to explore it
        </div>
      </div>
    </div>

    <div class="vq-stats-page">
      <div class="vq-stats-row vq-stats-row--top" id="clu-kpis"></div>

      <div class="vq-chart-card clu-filters" style="min-height:auto">
        <div class="clu-filters__row">
          <label>Min identity
            <input type="range" id="clu-id"  min="${FLOOR_ID}"  max="100" step="1" value="${_state.id}">
            <output id="clu-id-out">${_state.id}%</output>
          </label>
          <label>Min coverage
            <input type="range" id="clu-cov" min="${FLOOR_COV}" max="100" step="1" value="${_state.cov}">
            <output id="clu-cov-out">${_state.cov}%</output>
          </label>
          <label>Min samples
            <select class="vq-select" id="clu-minsamples">
              ${_samples.slice(1).map((_, i) => `<option value="${i + 2}">${i + 2}</option>`).join('')}
            </select>
          </label>
          <label>Species
            <select class="vq-select" id="clu-agree">
              <option value="all">all clusters</option>
              <option value="ok">members agree</option>
              <option value="warn">members disagree</option>
            </select>
          </label>
          ${families.length ? `
          <label>Family
            <select class="vq-select" id="clu-family">
              <option value="all">all</option>
              ${families.map(f => `<option value="${esc(f)}">${esc(f)}</option>`).join('')}
            </select>
          </label>` : ''}
          ${novelties.length ? `
          <label>Novelty
            <select class="vq-select" id="clu-novelty">
              <option value="all">all</option>
              ${novelties.map(t => `<option value="${t}">${esc(NOV.label[t])}</option>`).join('')}
            </select>
          </label>` : ''}
          <span class="clu-filters__readout" id="clu-readout"></span>
        </div>
      </div>

      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Cluster explorer</h3>
          <span class="vq-group__hint">pick a cluster — the card below shows only that cluster</span>
        </div>
        <div class="clu-stack">
          <div class="vq-chart-card" id="clu-table-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title">Clusters</div>
                <div class="vq-chart-card__sub" id="clu-table-sub"></div>
              </div>
              <div class="vq-vt-tools">
                <input class="vq-input vq-input--sm" type="search" id="clu-search"
                       placeholder="Filter cluster, species, sample…" aria-label="Filter clusters">
                <button class="vq-btn vq-btn--sm" type="button" id="clu-csv"
                        title="Download the clusters shown (filters and sort applied) as CSV">Export CSV</button>
              </div>
            </div>
            <div class="vq-chart-card__body" style="justify-content:flex-start">
              <div class="vq-vt-wrap clu-table-wrap" id="clu-table"></div>
            </div>
          </div>
          <div class="vq-chart-card clu-detail" id="clu-detail" style="min-height:auto"></div>
        </div>
      </section>

      <section class="vq-group">
        <div class="vq-group__head">
          <h3 class="vq-group__title">Across samples</h3>
          <span class="vq-group__hint">every cluster passing the filters · the selected cluster is highlighted</span>
        </div>
        <div class="vq-chart-card" id="clu-matrix-card" style="min-height:auto">
            <div class="vq-chart-card__head">
              <div>
                <div class="vq-chart-card__title">Cluster Presence</div>
                <div class="vq-chart-card__sub" id="clu-matrix-sub"></div>
              </div>
              <div class="vq-toggle" id="clu-matrix-mode" role="tablist" aria-label="Cell value">
                <button class="vq-toggle__btn${_state.matrixMode === 'presence' ? ' active' : ''}" type="button" data-mode="presence">Presence</button>
                <button class="vq-toggle__btn${_state.matrixMode === 'members' ? ' active' : ''}" type="button" data-mode="members">Members</button>
                ${_hasRpm ? `<button class="vq-toggle__btn${_state.matrixMode === 'rpm' ? ' active' : ''}" type="button" data-mode="rpm" title="reads per million input reads — comparable across samples">RPM</button>` : ''}
                ${_hasTpm ? `<button class="vq-toggle__btn${_state.matrixMode === 'tpm' ? ' active' : ''}" type="button" data-mode="tpm" title="TPM within each sample's own index — compare within a sample only">TPM</button>` : ''}
              </div>
            </div>
            <div class="vq-chart-card__body clu-scroll" id="clu-matrix"></div>
        </div>
      </section>
    </div>
  `;

  // ── Filters ──
  const idIn = document.getElementById('clu-id'), covIn = document.getElementById('clu-cov');
  const apply = () => {
    _state.id  = +idIn.value;
    _state.cov = +covIn.value;
    document.getElementById('clu-id-out').textContent  = _state.id + '%';
    document.getElementById('clu-cov-out').textContent = _state.cov + '%';
    _rerender();
  };
  idIn.addEventListener('input', apply);
  covIn.addEventListener('input', apply);
  [['clu-minsamples', 'minSamples', Number], ['clu-agree', 'agree', String],
   ['clu-family', 'family', String], ['clu-novelty', 'novelty', String]].forEach(([id, key, cast]) =>
    document.getElementById(id)?.addEventListener('change', e => { _state[key] = cast(e.target.value); _rerender(); }));

  // ── Table ──
  document.getElementById('clu-search').addEventListener('input', e => {
    _state.q = e.target.value.trim().toLowerCase(); _rerender();
  });
  document.getElementById('clu-table').addEventListener('click', e => {
    const th = e.target.closest('th[data-sort]');
    if (th) {
      const key = th.dataset.sort;
      const asc = ['gid', 'species', 'family', 'novRank', 'agreeRank'];
      _state.dir = _state.sort === key ? -_state.dir : (asc.includes(key) ? 1 : -1);
      _state.sort = key;
      _rerender();
      return;
    }
    const tr = e.target.closest('tr[data-gid]');
    if (tr) _select(tr.dataset.gid);
  });
  document.getElementById('clu-csv').addEventListener('click', () =>
    VQ.downloadText(_clustersCsv(), 'viralquest_clusters.csv'));

  // ── Across-sample toggles ──
  document.querySelectorAll('#clu-matrix-mode [data-mode]').forEach(b => b.addEventListener('click', () => {
    _state.matrixMode = b.dataset.mode;
    document.querySelectorAll('#clu-matrix-mode [data-mode]').forEach(x => x.classList.toggle('active', x === b));
    VQ.tooltipHide();
    _drawMatrix();
  }));

  _rerender();
  VQ.redrawOnResize(document.getElementById('clu-matrix'), _drawAcross);
}

// ── Filtering & per-cluster summary ─────────────────────────────────────────

function _passing(cluster) {
  // Representative always passes (100/100). Returns members meeting thresholds.
  return cluster.members.filter(m =>
    m.is_representative || (m.identity >= _state.id && m.coverage >= _state.cov));
}

/* Summary of a cluster restricted to its passing members. */
function _summary(c) {
  const members = _passing(c);
  const rep     = members.find(m => m.is_representative) || members[0];
  const others  = members.filter(m => !m.is_representative);
  const bySample = new Map();
  members.forEach(m => {
    if (!bySample.has(m.sample)) bySample.set(m.sample, { n: 0, tpm: 0, reads: 0, hasTpm: false, hasReads: false });
    const b = bySample.get(m.sample);
    b.n += 1;
    if (m.reads != null) { b.reads += m.reads; b.hasReads = true; }
    if (m.tpm != null) { b.tpm += m.tpm; b.hasTpm = true; }
  });
  const named = members.filter(m => m.species_db === 'nr' || m.species_db === 'refseq');
  const divergent = named.filter(m => m.species && m.species !== c.species).length;
  const NOV = window.vqNovelty;
  return {
    c, gid: c.gid, members, rep, bySample,
    species:   c.species || '',
    family:    c._family,
    novelty:   c._novelty,
    novRank:   NOV && c._novelty ? NOV.order.indexOf(c._novelty) : 99,
    size:      members.length,
    nSamples:  bySample.size,
    minId:     others.length ? d3.min(others, m => m.identity) : null,
    meanId:    others.length ? d3.mean(others, m => m.identity) : null,
    minCov:    others.length ? d3.min(others, m => m.coverage) : null,
    repLen:    rep?.length || 0,
    tpm:       _hasTpm ? d3.sum([...bySample.values()], b => b.tpm) : null,
    divergent,
    agreeRank: divergent ? 1 : 0,
  };
}

function _filtered() {
  const q = _state.q;
  let rows = _clusters.map(_summary).filter(r =>
    r.nSamples >= Math.max(2, _state.minSamples) &&
    (_state.agree === 'all' || (_state.agree === 'ok' ? !r.divergent : r.divergent > 0)) &&
    (_state.family === 'all' || r.family === _state.family) &&
    (_state.novelty === 'all' || r.novelty === _state.novelty) &&
    (!q || `${r.gid} ${r.species} ${r.family} ${[...r.bySample.keys()].join(' ')}`.toLowerCase().includes(q)));
  const key = _state.sort, dir = _state.dir;
  return rows.sort((a, b) => {
    const x = a[key], y = b[key];
    let d;
    if (typeof x === 'string' || typeof y === 'string') {
      if (!x !== !y) return !x ? 1 : -1;
      d = dir * String(x).localeCompare(String(y));
    } else {
      if ((x == null) !== (y == null)) return x == null ? 1 : -1;
      d = dir * ((x ?? 0) - (y ?? 0));
    }
    return d || b.nSamples - a.nSamples || b.size - a.size || a.gid.localeCompare(b.gid);
  });
}

function _rerender() {
  _view = _filtered();
  if (!_view.some(r => r.gid === _state.sel)) _state.sel = _view.length ? _view[0].gid : null;
  const readout = document.getElementById('clu-readout');
  if (readout) readout.textContent = `${_view.length} / ${_clusters.length} clusters shown`;
  _renderKpis();
  _renderTable();
  _renderDetail();
  _drawAcross();
}

function _select(gid) {
  if (!gid || gid === _state.sel) return;
  _state.sel = gid;
  // The RNA Quantification tab follows the cluster picked here.
  document.dispatchEvent(new CustomEvent('vq:cluster-select', { detail: gid }));
  _renderTable();
  _renderDetail();
  _drawAcross();
}

// ── KPIs ────────────────────────────────────────────────────────────────────

function _renderKpis() {
  const host = document.getElementById('clu-kpis');
  if (!host) return;
  const fmt = VQ.fmtNum;
  const connected = new Set(_view.flatMap(r => [...r.bySample.keys()]));
  const core = _view.filter(r => r.nSamples === _samples.length).length;
  const hetero = _view.filter(r => r.divergent > 0).length;
  const largest = _view.reduce((a, r) => (!a || r.size > a.size ? r : a), null);
  host.innerHTML = [
    VQ.statChip('Clusters', fmt(_view.length), 'accent',
      _view.length === _clusters.length ? 'shared by ≥2 samples' : `of ${fmt(_clusters.length)} · filters applied`),
    VQ.statChip('Samples Connected', `${connected.size}/${_samples.length}`, '',
      'share at least one cluster'),
    VQ.statChip('Shared by All', fmt(core), 'success', `present in all ${_samples.length} samples`),
    VQ.statChip('Species Disagree', fmt(hetero), hetero ? 'warn' : '',
      'members hit a different species'),
    VQ.statChip('Largest Cluster', largest ? fmt(largest.size) : '—', '',
      largest ? `${largest.gid} · ${largest.nSamples} samples` : ''),
  ].join('');
}

// ── Cluster table (the selector) ────────────────────────────────────────────

const _fmt1 = v => (v == null || isNaN(v)) ? '—' : v.toFixed(1);

function _novPill(t) {
  const NOV = window.vqNovelty;
  if (!t || !NOV) return '<span class="vq-vt-na">—</span>';
  return `<span class="vq-nov vq-nov--${t}" title="${VQ.esc(NOV.long[t])}">${VQ.esc(NOV.label[t])}</span>`;
}

function _agreePill(r) {
  return r.divergent
    ? `<span class="vq-wf-pill vq-wf-pill--partial" title="${r.divergent} member(s) hit a different species">⚠ ${r.divergent} disagree</span>`
    : `<span class="vq-wf-pill vq-wf-pill--done" title="every named member hits the representative's species">✓ agree</span>`;
}

/* Identity vs the representative, coloured like the virus table (nt cut-offs). */
function _idCell(v) {
  if (v == null || isNaN(v)) return '<td class="vq-vt-num"><span class="vq-vt-na">—</span></td>';
  const q = v >= 95 ? ['hi', '≥95% nt · same species'] : v >= 85 ? ['ok', '85–95% nt · variant']
          : v >= 70 ? ['mid', '70–85% nt · distant'] : ['lo', '<70% nt'];
  return `<td class="vq-vt-num"><span class="vq-q vq-q--${q[0]}" title="${q[1]}">${v.toFixed(1)}</span></td>`;
}

function _covCell(v) {
  if (v == null || isNaN(v)) return '<td class="vq-vt-num"><span class="vq-vt-na">—</span></td>';
  const q = v >= 70 ? 'hi' : v >= 40 ? 'mid' : 'lo';
  return `<td class="vq-vt-num"><span class="vq-q vq-q--${q}">${v.toFixed(1)}</span></td>`;
}

function _renderTable() {
  const wrap = document.getElementById('clu-table');
  if (!wrap) return;
  const esc = VQ.esc;
  const sub = document.getElementById('clu-table-sub');
  if (sub) sub.textContent =
    `${_view.length} cluster${_view.length === 1 ? '' : 's'} · identity vs the representative · click a row to select`;

  const th = (k, label, cls = '') => {
    const on = _state.sort === k;
    return `<th data-sort="${k}" class="${cls}${on ? ' vq-vt-sorted' : ''}" scope="col">
      ${label}<span class="vq-vt-arrow">${on ? (_state.dir > 0 ? '▲' : '▼') : ''}</span></th>`;
  };

  wrap.innerHTML = _view.length ? `
    <table class="vq-table vq-vt">
      <thead><tr>
        ${th('gid', 'Cluster')}
        ${th('species', 'Species')}
        ${th('family', 'Family')}
        ${th('novRank', 'Novelty')}
        ${th('nSamples', 'Samples', 'vq-vt-num')}
        ${th('size', 'Members', 'vq-vt-num')}
        ${th('minId', 'Min id %', 'vq-vt-num')}
        ${th('meanId', 'Mean id %', 'vq-vt-num')}
        ${th('minCov', 'Min cov %', 'vq-vt-num')}
        ${th('repLen', 'Rep. length', 'vq-vt-num')}
        ${_hasTpm ? th('tpm', 'TPM', 'vq-vt-num') : ''}
        ${th('agreeRank', 'Species agreement')}
      </tr></thead>
      <tbody>
        ${_view.map(r => `
          <tr data-gid="${esc(r.gid)}" class="${r.gid === _state.sel ? 'is-selected' : ''}">
            <td class="vq-td--mono">${esc(r.gid)}</td>
            <td class="vq-vt-species">${r.species ? esc(r.species) : '<span class="vq-vt-na">—</span>'}</td>
            <td>${r.family ? esc(r.family) : '<span class="vq-vt-na">—</span>'}</td>
            <td>${_novPill(r.novelty)}</td>
            <td class="vq-vt-num" title="${esc([...r.bySample.keys()].join(', '))}">${r.nSamples}/${_samples.length}</td>
            <td class="vq-vt-num">${r.size}</td>
            ${_idCell(r.minId)}
            ${_idCell(r.meanId)}
            ${_covCell(r.minCov)}
            <td class="vq-vt-num">${r.repLen.toLocaleString()}</td>
            ${_hasTpm ? `<td class="vq-vt-num">${_fmtTpm(r.tpm)}</td>` : ''}
            <td>${_agreePill(r)}</td>
          </tr>`).join('')}
      </tbody>
    </table>` : '<div class="vq-empty" style="padding:24px;font-size:12px">No cluster passes the filters.</div>';

  // Keep the selected row visible inside the table only (never scroll the page).
  const tr = wrap.querySelector('tr.is-selected');
  if (tr) {
    const head = wrap.querySelector('thead')?.offsetHeight || 0;
    const top = tr.getBoundingClientRect().top - wrap.getBoundingClientRect().top + wrap.scrollTop;
    if (top - head < wrap.scrollTop || top + tr.offsetHeight > wrap.scrollTop + wrap.clientHeight)
      wrap.scrollTop = Math.max(0, top - head - (wrap.clientHeight - head) / 2);
  }
}

function _fmtTpm(v) {
  if (v == null || isNaN(v)) return '—';
  if (v === 0) return '0';
  if (v >= 1000) return Math.round(v).toLocaleString();
  if (v >= 10) return v.toFixed(1);
  return v.toFixed(2);
}

/* CSV of the clusters shown: plain numbers, RFC 4180 quoting, lists joined by ";". */
function _clustersCsv() {
  const cols = [
    ['cluster_id',        r => r.gid],
    ['species',           r => r.species],
    ['family',            r => r.family],
    ['novelty',           r => r.novelty],
    ['n_samples',         r => r.nSamples],
    ['samples',           r => [...r.bySample.keys()].join(';')],
    ['n_members',         r => r.size],
    ['min_identity',      r => r.minId == null ? null : +r.minId.toFixed(2)],
    ['mean_identity',     r => r.meanId == null ? null : +r.meanId.toFixed(2)],
    ['min_coverage',      r => r.minCov == null ? null : +r.minCov.toFixed(2)],
    ['representative',    r => r.c.representative],
    ['representative_length', r => r.repLen],
    ...(_hasTpm ? [['total_tpm', r => +r.tpm.toFixed(3)]] : []),
    ['species_agreement', r => r.divergent ? 'disagree' : 'agree'],
    ['divergent_members', r => r.divergent],
  ];
  return _csv(cols, _view);
}

function _membersCsv(r) {
  const cols = [
    ['cluster_id',     () => r.gid],
    ['sample',         m => m.sample],
    ['contig_id',      m => m.seq_id],
    ['representative', m => m.is_representative ? 1 : 0],
    ['length_bp',      m => m.length],
    ['identity',       m => m.identity],
    ['coverage',       m => m.coverage],
    ['species',        m => m.species],
    ['species_source', m => m.species_db],
    ...(_hasTpm ? [['tpm', m => m.tpm]] : []),
  ];
  return _csv(cols, r.members);
}

function _csv(cols, rows) {
  const cell = v => {
    if (v == null || (typeof v === 'number' && isNaN(v))) return '';
    const s = String(v);
    return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [cols.map(([h]) => h).join(','),
          ...rows.map(r => cols.map(([, f]) => cell(f(r))).join(','))].join('\r\n') + '\r\n';
}

// ── Selected cluster detail ─────────────────────────────────────────────────

function _renderDetail() {
  const host = document.getElementById('clu-detail');
  if (!host) return;
  const esc = VQ.esc;
  const r = _view.find(x => x.gid === _state.sel);
  if (!r) {
    host.innerHTML = '<div class="vq-empty" style="padding:24px">No cluster selected — relax the filters.</div>';
    return;
  }
  const idx = _view.indexOf(r);
  const safe = VQ.safeId(r.gid);
  // Reads per million when Salmon reads + library size are known (comparable
  // across samples); TPM otherwise.
  const rpmOf = (s, b) => b.hasReads && _inputReads.get(s) ? b.reads / _inputReads.get(s) * 1e6 : null;
  const maxRpm = d3.max(_samples, s => { const b = r.bySample.get(s); return b ? rpmOf(s, b) : null; }) || 0;
  const maxTpm = d3.max([...r.bySample.values()], b => b.tpm) || 0;

  const chips = _samples.map(s => {
    const b = r.bySample.get(s);
    if (!b) return `<span class="clu-presence clu-presence--off" title="${esc(s)}: not in this cluster">
                      <span class="clu-presence__name">${esc(s)}</span><span class="clu-presence__val">—</span></span>`;
    const rpm = rpmOf(s, b);
    const qTxt = rpm != null ? ` · ${_fmtTpm(rpm)} RPM` : b.hasTpm ? ` · ${_fmtTpm(b.tpm)} TPM` : '';
    const w = rpm != null && maxRpm > 0 ? Math.max(4, 100 * rpm / maxRpm)
            : b.hasTpm && maxTpm > 0 ? Math.max(4, 100 * b.tpm / maxTpm) : 100;
    return `<span class="clu-presence" title="${esc(s)}: ${b.n} member${b.n > 1 ? 's' : ''}${qTxt}">
              <span class="clu-presence__name">${esc(s)}</span>
              <span class="clu-presence__val">${b.n}×${qTxt}</span>
              ${rpm != null || b.hasTpm ? `<span class="clu-presence__bar"><span style="width:${w.toFixed(1)}%"></span></span>` : ''}
            </span>`;
  }).join('');

  const stat = (k, v) => `<div class="clu-stat"><span>${esc(k)}</span><strong>${v}</strong></div>`;

  host.innerHTML = `
    <div class="vq-chart-card__head">
      <div>
        <div class="clu-detail__title">
          <span class="vq-badge vq-badge--cluster">${esc(r.gid)}</span>
          <span class="clu-species">${esc(r.species || 'no species label')}</span>
        </div>
        <div class="clu-detail__pills">
          ${r.family ? `<span class="vq-wf-opt"><span class="vq-wf-opt__k">Family</span>${esc(r.family)}</span>` : ''}
          ${_novPill(r.novelty)}
          ${_agreePill(r)}
        </div>
      </div>
      <div class="vq-section-actions">
        <button class="vq-btn vq-btn--sm vq-btn--ghost" type="button" id="clu-prev"
                ${idx <= 0 ? 'disabled' : ''} title="Previous cluster in the table">‹</button>
        <span class="clu-detail__pos">${idx + 1} / ${_view.length}</span>
        <button class="vq-btn vq-btn--sm vq-btn--ghost" type="button" id="clu-next"
                ${idx >= _view.length - 1 ? 'disabled' : ''} title="Next cluster in the table">›</button>
        ${_exportMenu('clu-menu-' + safe)}
      </div>
    </div>

    <div class="clu-stats">
      ${stat('Members', r.size)}
      ${stat('Samples', `${r.nSamples}/${_samples.length}`)}
      ${stat('Identity', r.minId == null ? '—' : `${_fmt1(r.minId)}–${_fmt1(d3.max(r.members.filter(m => !m.is_representative), m => m.identity))}%`)}
      ${stat('Min coverage', r.minCov == null ? '—' : _fmt1(r.minCov) + '%')}
      ${stat('Representative', `${r.repLen.toLocaleString()} nt`)}
      ${_hasTpm ? stat('Total TPM', _fmtTpm(r.tpm)) : ''}
    </div>

    <div class="clu-block">
      <div class="clu-block__title">Presence across samples <span>${_inputReads.size ? 'members · reads per million (Salmon)' : 'members per sample'}</span></div>
      <div class="clu-presence-row">${chips}</div>
    </div>

    <div class="clu-block">${_homogeneityBanner(r.c, r.members)}</div>

    <div class="clu-block">
      <div class="clu-block__title">Members <span>click a contig to open it in the Sequence Viewer</span></div>
      <div class="vq-vt-wrap clu-members-wrap">
        <table class="vq-table vq-vt">
          <thead><tr>
            <th>Sample</th><th>Contig</th><th class="vq-vt-num">Length (bp)</th>
            <th class="vq-vt-num">Identity %</th><th class="vq-vt-num">Coverage %</th>
            <th>Species (own best hit)</th>${_hasTpm ? '<th class="vq-vt-num">TPM</th>' : ''}
          </tr></thead>
          <tbody>
            ${r.members.map(m => `
              <tr data-gid="${esc(m.gid)}"${m.is_representative ? ' class="clu-rep"' : ''}>
                <td>${esc(m.sample)}</td>
                <td class="vq-td--mono">${esc(m.seq_id)}${m.is_representative ? ' <span class="clu-rep-tag">rep.</span>' : ''}</td>
                <td class="vq-vt-num">${(m.length || 0).toLocaleString()}</td>
                ${m.is_representative ? '<td class="vq-vt-num"><span class="vq-vt-na">ref</span></td><td class="vq-vt-num"><span class="vq-vt-na">ref</span></td>'
                                      : _idCell(m.identity) + _covCell(m.coverage)}
                <td class="vq-vt-species${m.species && r.species && m.species !== r.species && (m.species_db === 'nr' || m.species_db === 'refseq') ? ' clu-species--off' : ''}">
                  ${m.species ? esc(m.species) : '<span class="vq-vt-na">—</span>'}
                  ${m.species_db ? `<span class="clu-db">${esc(m.species_db)}</span>` : ''}</td>
                ${_hasTpm ? `<td class="vq-vt-num">${_fmtTpm(m.tpm)}</td>` : ''}
              </tr>`).join('')}
          </tbody>
        </table>
      </div>
    </div>

    <div class="clu-block">
      <div class="clu-block__title">Alignment against the representative
        <span>bar = aligned region · outline = full contig</span></div>
      <div class="vq-genome-wrap" id="clu-wrap-${safe}"></div>
    </div>`;

  host.querySelector('#clu-prev')?.addEventListener('click', () => _select(_view[idx - 1]?.gid));
  host.querySelector('#clu-next')?.addEventListener('click', () => _select(_view[idx + 1]?.gid));
  host.querySelectorAll('.clu-members-wrap tr[data-gid]').forEach(tr =>
    tr.addEventListener('click', () => VQ.jumpToViewer(tr.dataset.gid)));

  const view = Object.assign({}, r.c, { members: r.members });
  _wireExportMenu(host, 'clu-menu-' + safe,
    () => host.querySelector('.vq-genome-wrap svg'),
    `cluster_${r.gid}`, () => _clusterFasta(view), () => _membersCsv(r));

  const wrapEl = host.querySelector('#clu-wrap-' + safe);
  const draw = () => {
    if (!wrapEl.clientWidth) return;
    wrapEl.innerHTML = '';
    wrapEl.appendChild(_clusterSVG(view, r.repLen || 1, wrapEl.clientWidth));
  };
  requestAnimationFrame(draw);
  VQ.redrawOnResize(wrapEl, draw);
}

// ── Across samples: presence matrix + sample similarity ─────────────────────

/* Per-sample value of every filtered cluster (members or summed TPM). */
function _sampleVectors() {
  const vec = new Map(_samples.map(s => [s, new Map()]));
  _view.forEach(r => r.bySample.forEach((b, s) => vec.get(s)?.set(r.gid, b)));
  return vec;
}

/* Jaccard similarity of the samples' cluster sets — orders the matrix columns
   so samples sharing the same clusters sit next to each other. */
function _similarity(vec) {
  const n = _samples.length;
  const M = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    const a = vec.get(_samples[i]), b = vec.get(_samples[j]);
    const inter = [...a.keys()].filter(k => b.has(k)).length;
    const uni = new Set([...a.keys(), ...b.keys()]).size;
    M[i][j] = M[j][i] = uni ? inter / uni : (i === j ? 1 : 0);
  }
  return { M };
}

function _drawAcross() {
  _drawMatrix();
}

function _drawMatrix() {
  const host = document.getElementById('clu-matrix');
  if (!host) return;
  const esc = VQ.esc;
  if (!_view.length) { host.innerHTML = '<div class="vq-empty">No cluster passes the filters.</div>'; return; }

  const order = _simOrder();
  const cols = order.map(i => _samples[i]);
  // Rows: shared core first (prevalence), then size.
  const rows = _view.slice().sort((a, b) => b.nSamples - a.nSamples || b.size - a.size || a.gid.localeCompare(b.gid));
  const mode = _state.matrixMode;
  const selS = _selSamples();

  const sub = document.getElementById('clu-matrix-sub');
  if (sub) sub.textContent = 'clusters × samples · shared core on top · columns grouped by shared clusters · '
    + { presence: 'filled = cluster present', members: 'colour = member contigs in the sample',
        rpm: 'colour = reads per million input reads (log)', tpm: 'colour = TPM (log) · compare within a sample only' }[mode]
    + ' · click a row to select';

  // Cell value per mode (null = absent).
  const valOf = (b, s) => {
    if (!b) return null;
    if (mode === 'members') return b.n;
    if (mode === 'rpm') return b.hasReads && _inputReads.get(s) ? b.reads / _inputReads.get(s) * 1e6 : null;
    if (mode === 'tpm') return b.hasTpm ? b.tpm : null;
    return 1;
  };
  const vals = rows.flatMap(r => cols.map(s => valOf(r.bySample.get(s), s))).filter(v => v != null);
  const maxV = d3.max(vals) || 1;
  const t01 = v => mode === 'members' ? Math.sqrt(v / maxV)
                 : mode === 'presence' ? 0.7
                 : Math.log10(1 + v) / Math.log10(1 + maxV);
  const color = v => d3.interpolateBlues(0.18 + 0.77 * t01(v));
  const fmtCell = v => mode === 'members' ? String(v)
    : v >= 1e4 ? (v / 1e3).toFixed(0) + 'k' : v >= 1e3 ? (v / 1e3).toFixed(1) + 'k'
    : v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v >= 0.01 ? v.toFixed(2) : '<0.01';

  // Layout: cells share the card width (the matrix is the card's only content).
  const avail = host.clientWidth || 900;
  const idW = 96, padR = 96, rowH = 24;
  const spW = Math.min(300, Math.max(150, avail * 0.2));
  const padL = idW + spW;
  const longest = d3.max(cols, s => s.length) || 4;
  const cell = Math.max(22, Math.min(140, (avail - padL - padR) / cols.length));
  const rotate = cell < longest * 7 + 10;
  const headH = rotate ? Math.min(120, longest * 5.2 + 20) : 26;
  const legendH = mode === 'presence' ? 0 : 40;
  const W = Math.max(avail, padL + padR + cell * cols.length);
  const H = headH + rows.length * rowH + 8 + legendH;
  const spChars = Math.floor((spW - 14) / 6.4);
  const showVals = mode !== 'presence' && cell >= 40;

  host.innerHTML = '';
  const svg = d3.select(host).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`).style('display', 'block');

  cols.forEach((s, i) => {
    const x = padL + i * cell + cell / 2;
    const t = svg.append('text').attr('class', 'clu-axis' + (selS.has(s) ? ' is-on' : '')).text(s);
    if (rotate) t.attr('transform', `translate(${x + 3},${headH - 6}) rotate(-45)`);
    else t.attr('x', x).attr('y', headH - 9).attr('text-anchor', 'middle');
  });
  svg.append('text').attr('class', 'clu-axis').attr('x', padL + cols.length * cell + 12).attr('y', headH - 9)
    .text('samples');

  rows.forEach((r, j) => {
    const y = headH + j * rowH;
    const sel = r.gid === _state.sel;
    const g = svg.append('g').attr('class', 'clu-mrow' + (sel ? ' is-selected' : '')).style('cursor', 'pointer')
      .on('click', () => { _select(r.gid); document.getElementById('clu-detail')?.scrollIntoView({ behavior: 'smooth', block: 'start' }); });
    g.append('rect').attr('class', 'clu-mrow__bg').attr('x', 0).attr('y', y).attr('width', W).attr('height', rowH);
    g.append('text').attr('class', 'clu-mrow__id').attr('x', 8).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(r.gid);
    g.append('text').attr('class', 'clu-mrow__sp').attr('x', idW + 4).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(_trunc(r.species, spChars)).append('title').text(r.species);
    cols.forEach((s, i) => {
      const b = r.bySample.get(s);
      const x = padL + i * cell;
      if (!b) {
        g.append('rect').attr('class', 'clu-cell--off')
          .attr('x', x + cell / 2 - 2).attr('y', y + rowH / 2 - 2).attr('width', 4).attr('height', 4).attr('rx', 2);
        return;
      }
      const v = valOf(b, s);
      const fill = v == null ? 'var(--vq-border-dark)' : color(v);
      g.append('rect').attr('class', 'clu-cell')
        .attr('x', x + 2).attr('y', y + 2.5).attr('width', cell - 4).attr('height', rowH - 5).attr('rx', 4)
        .attr('fill', fill)
        .on('mousemove', evt => {
          const rpm = b.hasReads && _inputReads.get(s) ? b.reads / _inputReads.get(s) * 1e6 : null;
          VQ.tooltipShow(`
          <div class="vq-tooltip__title">${esc(r.gid)} · ${esc(s)}</div>
          <div class="vq-tooltip__row">
            <span class="vq-tooltip__key">Species</span><span>${esc(r.species || '—')}</span>
            <span class="vq-tooltip__key">Members</span><span>${b.n}</span>
            ${rpm != null ? `<span class="vq-tooltip__key">RPM</span><span>${_fmtTpm(rpm)}</span>` : ''}
            ${b.hasTpm ? `<span class="vq-tooltip__key">TPM</span><span>${_fmtTpm(b.tpm)}</span>` : ''}
          </div>`, evt); })
        .on('mouseleave', VQ.tooltipHide);
      if (showVals && v != null) g.append('text').attr('class', 'clu-cell__v')
        .attr('x', x + cell / 2).attr('y', y + rowH / 2).attr('text-anchor', 'middle').attr('dominant-baseline', 'central')
        .attr('fill', t01(v) > 0.55 ? '#fff' : 'var(--vq-text)').text(fmtCell(v));
    });
    // Prevalence: bar + k/N
    const px = padL + cols.length * cell + 12, pw = 40;
    g.append('rect').attr('class', 'clu-prev__track').attr('x', px).attr('y', y + rowH / 2 - 3).attr('width', pw).attr('height', 6).attr('rx', 3);
    g.append('rect').attr('x', px).attr('y', y + rowH / 2 - 3).attr('width', pw * r.nSamples / _samples.length).attr('height', 6).attr('rx', 3)
      .attr('fill', r.nSamples === _samples.length ? 'var(--vq-success)' : 'var(--vq-accent)');
    g.append('text').attr('class', 'clu-mrow__prev').attr('x', px + pw + 6).attr('y', y + rowH / 2)
      .attr('dominant-baseline', 'central').text(`${r.nSamples}/${_samples.length}`);
  });

  // Colour legend
  if (legendH) {
    const ly = headH + rows.length * rowH + 20, lw = Math.min(220, cols.length * cell);
    const grad = svg.append('defs').append('linearGradient').attr('id', 'clu-mx-grad');
    d3.range(0, 1.01, 0.1).forEach(f => grad.append('stop').attr('offset', f)
      .attr('stop-color', d3.interpolateBlues(0.18 + 0.77 * f)));
    svg.append('rect').attr('x', padL).attr('y', ly).attr('width', lw).attr('height', 8).attr('rx', 4).attr('fill', 'url(#clu-mx-grad)');
    const lo = mode === 'members' ? '1' : '0';
    svg.append('text').attr('class', 'clu-axis').attr('x', padL).attr('y', ly + 20).text(lo);
    svg.append('text').attr('class', 'clu-axis').attr('x', padL + lw).attr('y', ly + 20).attr('text-anchor', 'end')
      .text(`${fmtCell(maxV)} ${mode === 'members' ? 'members' : mode.toUpperCase()}`);
    svg.append('text').attr('class', 'clu-axis').attr('x', padL + lw + 14).attr('y', ly + 8)
      .text(mode === 'members' ? '√ scale' : 'log scale');
  }
}

function _selSamples() {
  const r = _view.find(x => x.gid === _state.sel);
  return new Set(r ? r.bySample.keys() : []);
}

let _simCache = null;
function _simOrder() {
  return _simData().tree.leaves;
}

function _simData() {
  const key = _view.map(r => r.gid).join(',') + '|' + _state.id + '|' + _state.cov;
  if (_simCache && _simCache.key === key) return _simCache;
  const sim = _similarity(_sampleVectors());
  _simCache = { key, ...sim, tree: VQ.upgma(sim.M) };
  return _simCache;
}

function _trunc(s, n) {
  s = String(s ?? '');
  return s.length > n ? s.slice(0, n - 1) + '…' : s;
}

// ── Species agreement banner ───────────────────────────────────────────────

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
        <button class="vq-menu__item" data-fmt="csv"   role="menuitem" type="button">Members CSV</button>
      </div>
    </div>`;
}

function _wireExportMenu(root, menuId, getSvg, baseName, getFasta, getCsv) {
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
      if (fmt === 'csv')   { if (getCsv) VQ.downloadText(getCsv(), baseName + '_members.csv'); return; }
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
    const col = _barColor(m);
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

/* Alignment bar colour: one colour, not an identity gradient — members are
   already filtered by minimum identity and coverage, so a 0–100 scale would
   never reach its red end. Identity is in the bar label and the tooltip. */
function _barColor(m) {
  return m.is_representative ? 'var(--vq-primary-light)' : 'var(--vq-accent)';
}

window.vqInitClusters = vqInitClusters;
})();
