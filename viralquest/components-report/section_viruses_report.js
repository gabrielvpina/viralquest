/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_viruses_report.js — Virome tab (multi-sample)
   Wraps the per-sample vqInitViruses (components/section_viruses.js — the
   same file the per-sample report uses), adding a sample checkbox filter on
   top. Sequences are keyed by their global id (sample::id), as in the
   Sequence Viewer, so table rows and scatter points open the right card.
   ============================================================ */

let _allSeqs  = [];
let _samples  = [];
let _selected = new Set();

function vqInitVirusesReport(sequences) {
  const el = document.getElementById('section-virome');
  if (!el) return;
  const esc = VQ.esc;

  _allSeqs = (sequences || []).map(s => Object.assign({}, s, {
    id: s.gid || s.id, _orig_id: s.id,
  }));
  _samples  = [...new Set(_allSeqs.map(s => s.sample).filter(Boolean))].sort();
  _selected = new Set(_samples);   // all samples on by default

  const checks = _samples.map(sm => `
    <label class="vr-sample">
      <input type="checkbox" value="${esc(sm)}" checked> ${esc(sm)}
      <span class="vr-sample__n">${_allSeqs.filter(s => s.sample === sm).length}</span>
    </label>`).join('');

  el.innerHTML = `
    <div class="vr-sample-bar">
      <span class="vr-sample-bar__label">Samples</span>
      <div class="vr-sample-bar__items" id="vv-sample-checks">${checks}</div>
      <div class="vr-sample-bar__actions">
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="vv-sample-all" type="button">All</button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="vv-sample-none" type="button">None</button>
      </div>
    </div>
    <div id="virome-host"></div>
  `;

  document.getElementById('vv-sample-checks').addEventListener('change', e => {
    if (e.target.matches('input[type="checkbox"]')) {
      if (e.target.checked) _selected.add(e.target.value);
      else _selected.delete(e.target.value);
      _render();
    }
  });
  document.getElementById('vv-sample-all')?.addEventListener('click', () => _setAll(true));
  document.getElementById('vv-sample-none')?.addEventListener('click', () => _setAll(false));

  _render();
}

function _setAll(on) {
  _selected = on ? new Set(_samples) : new Set();
  document.querySelectorAll('#vv-sample-checks input[type="checkbox"]')
    .forEach(c => { c.checked = on; });
  _render();
}

function _render() {
  const subset = _allSeqs.filter(s => _selected.has(s.sample));
  VQ.tooltipHide();
  vqInitViruses({ sequences: subset }, 'virome-host');
}

window.vqInitVirusesReport = vqInitVirusesReport;
})();
