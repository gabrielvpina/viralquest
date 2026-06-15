/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_viewer_report.js — Section 3: Sequence Viewer (multi-sample)
   Wraps the per-sample vqInitViewer, adding a sample checkbox filter
   on top.  Sequences are keyed by their global id (sample::id) so cards
   are unique across samples and cluster "jump to viewer" links resolve.
   ============================================================ */

let _allSeqs  = [];   // viewer-input copies (id = gid)
let _samples  = [];
let _selected = new Set();

function vqInitViewerReport(sequences) {
  const el = document.getElementById('section-viewer');
  if (!el) return;
  const esc = VQ.esc;

  // Use the global id as the display/lookup id so nothing collides across
  // samples; keep the original contig id available for reference.
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
      <div class="vr-sample-bar__items" id="vr-sample-checks">${checks}</div>
      <div class="vr-sample-bar__actions">
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="vr-sample-all" type="button">All</button>
        <button class="vq-btn vq-btn--ghost vq-btn--sm" id="vr-sample-none" type="button">None</button>
      </div>
    </div>
    <div id="viewer-host"></div>
  `;

  const checksWrap = document.getElementById('vr-sample-checks');
  checksWrap.addEventListener('change', e => {
    if (e.target.matches('input[type="checkbox"]')) {
      if (e.target.checked) _selected.add(e.target.value);
      else _selected.delete(e.target.value);
      _render();
    }
  });
  document.getElementById('vr-sample-all')?.addEventListener('click', () => _setAll(true));
  document.getElementById('vr-sample-none')?.addEventListener('click', () => _setAll(false));

  _render();
}

function _setAll(on) {
  _selected = on ? new Set(_samples) : new Set();
  document.querySelectorAll('#vr-sample-checks input[type="checkbox"]')
    .forEach(c => { c.checked = on; });
  _render();
}

function _render() {
  const subset = _allSeqs.filter(s => _selected.has(s.sample));
  // Delegate the heavy lifting to the per-sample viewer, mounted in our host.
  vqInitViewer(subset, 'viewer-host');
}

window.vqInitViewerReport = vqInitViewerReport;
})();
