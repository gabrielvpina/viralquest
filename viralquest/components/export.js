/* ============================================================
   export.js — shared PNG / SVG / PDF (print) export utilities
   All functions are pure; they receive the SVG element directly.
   ============================================================ */

'use strict';

// ── PNG export ──────────────────────────────────────────────────────────────
// Serialises the target SVG, renders it to a canvas, then triggers download.

function vqExportPNG(svgEl, filename = 'viralquest_export.png', scale = 2) {
  const serialiser = new XMLSerializer();
  const bbox       = svgEl.getBoundingClientRect();
  const w = Math.round(bbox.width  * scale) || 800;
  const h = Math.round(bbox.height * scale) || 400;

  // Inline computed styles so the PNG looks right outside the document
  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width',  w);
  clone.setAttribute('height', h);

  // Copy CSS custom property values onto the clone
  const styleEl = document.createElement('style');
  styleEl.textContent = _inlineCSSVars();
  clone.prepend(styleEl);

  const svgStr  = serialiser.serializeToString(clone);
  const blob    = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  const url     = URL.createObjectURL(blob);

  const img = new Image();
  img.onload = () => {
    const canvas = document.createElement('canvas');
    canvas.width  = w;
    canvas.height = h;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);
    URL.revokeObjectURL(url);
    _download(canvas.toDataURL('image/png'), filename);
  };
  img.onerror = () => { URL.revokeObjectURL(url); console.error('PNG export failed'); };
  img.src = url;
}

// ── SVG export ──────────────────────────────────────────────────────────────

function vqExportSVG(svgEl, filename = 'viralquest_export.svg') {
  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

  const styleEl = document.createElement('style');
  styleEl.textContent = _inlineCSSVars();
  clone.prepend(styleEl);

  const svgStr = new XMLSerializer().serializeToString(clone);
  const blob   = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  _download(URL.createObjectURL(blob), filename);
}

// ── PDF export (browser print dialog) ───────────────────────────────────────
// Opens the browser print dialog with the target SVG isolated.

function vqExportPDF(svgEl, title = 'ViralQuest Export') {
  const w = window.open('', '_blank', 'width=900,height=700');
  if (!w) { alert('Allow pop-ups for PDF export.'); return; }

  const clone    = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.style.cssText = 'max-width:100%;height:auto;';

  const styleText = _inlineCSSVars() + `
    body { margin: 1cm; font-family: sans-serif; }
    h3   { margin-bottom: 12px; font-size: 14px; color: #1a3a5c; }
    svg  { max-width: 100%; }
  `;

  w.document.write(`
    <!DOCTYPE html>
    <html><head>
      <title>${title}</title>
      <style>${styleText}</style>
    </head><body>
      <h3>${title}</h3>
      ${new XMLSerializer().serializeToString(clone)}
      <script>window.onload = () => { window.print(); window.close(); }<\/script>
    </body></html>
  `);
  w.document.close();
}

// ── Batch PNG export (multiple SVGs → individual PNGs) ──────────────────────
// Returns a Promise that resolves when all PNGs have been queued for download.

function vqExportBatchPNG(items, delay = 300) {
  return items.reduce((chain, { svgEl, filename }, i) =>
    chain.then(() => new Promise(resolve => {
      setTimeout(() => { vqExportPNG(svgEl, filename); resolve(); }, i * delay);
    })),
    Promise.resolve()
  );
}

// ── Tooltip helpers ──────────────────────────────────────────────────────────

const _tooltip = document.getElementById('vq-tooltip');

function vqTooltipShow(html, event) {
  if (!_tooltip) return;
  _tooltip.innerHTML = html;
  _tooltip.classList.add('visible');
  vqTooltipMove(event);
}

function vqTooltipMove(event) {
  if (!_tooltip) return;
  const x = event.clientX + 14;
  const y = event.clientY + 14;
  const tw = _tooltip.offsetWidth;
  const th = _tooltip.offsetHeight;
  _tooltip.style.left = (x + tw > window.innerWidth  ? x - tw - 28 : x) + 'px';
  _tooltip.style.top  = (y + th > window.innerHeight ? y - th - 28 : y) + 'px';
}

function vqTooltipHide() {
  if (_tooltip) _tooltip.classList.remove('visible');
}

// ── Internal helpers ─────────────────────────────────────────────────────────

function _inlineCSSVars() {
  const style = getComputedStyle(document.documentElement);
  const vars  = [
    '--vq-primary','--vq-accent','--vq-success','--vq-warning','--vq-danger',
    '--vq-bg','--vq-surface','--vq-surface-2','--vq-border','--vq-text','--vq-text-2','--vq-text-3',
    '--vq-font','--vq-font-mono',
    '--vq-dom-1','--vq-dom-2','--vq-dom-3','--vq-dom-4',
    '--vq-dom-5','--vq-dom-6','--vq-dom-7','--vq-dom-8',
    '--vq-tax-flaviviridae','--vq-tax-parvoviridae','--vq-tax-phenuiviridae',
    '--vq-tax-nodaviridae','--vq-tax-rhabdoviridae','--vq-tax-other',
  ];
  const decls = vars.map(v => `${v}:${style.getPropertyValue(v)}`).join(';');
  return `:root{${decls}}`;
}

function _download(url, filename) {
  const a  = document.createElement('a');
  a.href   = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}
