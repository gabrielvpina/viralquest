/* ============================================================
   export.js — shared export utilities + helpers
   Functions are pure: they receive the SVG element directly.
   Exposes everything on window.VQ so other section scripts can
   call them without relying on bare global names.
   ============================================================ */

'use strict';

// ── Small helpers ──────────────────────────────────────────────────────────

/** Escape a value for safe HTML interpolation. */
function vqEsc(v) {
  if (v == null) return '';
  return String(v)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

/** Replace special chars in a string before using it as a DOM id. */
function vqSafeId(s) {
  return String(s ?? '').replace(/[^A-Za-z0-9_\-]/g, '_');
}

// ── PNG export ──────────────────────────────────────────────────────────────

function vqExportPNG(svgEl, filename = 'viralquest_export.png', scale = 2) {
  if (!svgEl) return;
  const serialiser = new XMLSerializer();
  const bbox       = svgEl.getBoundingClientRect();
  const w = Math.round(bbox.width  * scale) || 800;
  const h = Math.round(bbox.height * scale) || 400;

  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width',  w);
  clone.setAttribute('height', h);

  _applyInlineVars(clone);

  const svgStr = serialiser.serializeToString(clone);
  const blob   = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  const url    = URL.createObjectURL(blob);

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
  if (!svgEl) return;
  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  _applyInlineVars(clone);

  const svgStr = new XMLSerializer().serializeToString(clone);
  const blob   = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  _download(URL.createObjectURL(blob), filename);
}

// ── PDF export (browser print dialog) ───────────────────────────────────────

function vqExportPDF(svgEl, title = 'ViralQuest Export') {
  if (!svgEl) return;
  const w = window.open('', '_blank', 'width=900,height=700');
  if (!w) { alert('Allow pop-ups for PDF export.'); return; }

  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.style.cssText = 'max-width:100%;height:auto;';
  _applyInlineVars(clone);

  const styleText = `
    body { margin: 1cm; font-family: -apple-system, 'Segoe UI', system-ui, sans-serif; }
    h3   { margin-bottom: 12px; font-size: 14px; color: #15314b; }
    svg  { max-width: 100%; }
  `;

  w.document.write(`
    <!DOCTYPE html>
    <html><head>
      <meta charset="UTF-8">
      <title>${vqEsc(title)}</title>
      <style>${styleText}</style>
    </head><body>
      <h3>${vqEsc(title)}</h3>
      ${new XMLSerializer().serializeToString(clone)}
      <script>window.onload = () => { window.print(); window.close(); }<\/script>
    </body></html>
  `);
  w.document.close();
}

/** Print the whole report (or a target node) via the browser's print dialog. */
function vqExportPagePDF(title = 'ViralQuest Report') {
  const prevTitle = document.title;
  document.title  = title;
  window.print();
  setTimeout(() => { document.title = prevTitle; }, 1000);
}

// ── Batch PNG export ────────────────────────────────────────────────────────

function vqExportBatchPNG(items, delay = 300) {
  return items.reduce((chain, { svgEl, filename }, i) =>
    chain.then(() => new Promise(resolve => {
      setTimeout(() => { vqExportPNG(svgEl, filename); resolve(); }, i * delay);
    })),
    Promise.resolve()
  );
}

// ── Composite export — stitch multiple SVGs into ONE file ───────────────────
// Each entry: { svgEl, label }.  We measure each SVG's intrinsic size from
// its viewBox or bounding rect and stack vertically, keeping every chart's
// horizontal proportions intact.

function _measureSvg(svg) {
  // Prefer viewBox so we preserve the original aspect ratio.
  const vb = svg.viewBox?.baseVal;
  if (vb && vb.width) return { w: vb.width, h: vb.height };
  const r = svg.getBoundingClientRect();
  return { w: r.width || 800, h: r.height || 400 };
}

/**
 * Build a single SVG that contains each input SVG as a <g> stacked vertically.
 * All inputs are scaled to a common width so proportions are preserved.
 */
function vqBuildCompositeSVG(items, opts = {}) {
  const W       = opts.width      || 1200;
  const labelH  = opts.labelHeight || 24;
  const gap     = opts.gap         || 16;
  const padding = opts.padding     || 24;

  // Plan: measure, compute scaled height per item, sum.
  const placed = items.map(it => {
    const { w, h } = _measureSvg(it.svgEl);
    const scaledH  = (W - padding * 2) * (h / w);
    return { ...it, w, h, scaledH };
  });

  const totalH = padding * 2
    + placed.reduce((a, p) => a + labelH + p.scaledH + gap, 0)
    - gap;

  const NS  = 'http://www.w3.org/2000/svg';
  const out = document.createElementNS(NS, 'svg');
  out.setAttribute('xmlns', NS);
  out.setAttribute('width',  W);
  out.setAttribute('height', totalH);
  out.setAttribute('viewBox', `0 0 ${W} ${totalH}`);

  // Inline CSS vars so colours render outside the document.
  _applyInlineVars(out);

  // Background
  const bg = document.createElementNS(NS, 'rect');
  bg.setAttribute('width',  W);
  bg.setAttribute('height', totalH);
  bg.setAttribute('fill',   '#ffffff');
  out.appendChild(bg);

  let y = padding;
  placed.forEach(p => {
    // Label
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', padding);
    t.setAttribute('y', y + labelH - 8);
    t.setAttribute('font-family', 'var(--vq-font), system-ui, sans-serif');
    t.setAttribute('font-size', '13');
    t.setAttribute('font-weight', '600');
    t.setAttribute('fill', '#15314b');
    t.textContent = p.label || '';
    out.appendChild(t);
    y += labelH;

    // Inner SVG cloned + scaled into a <g>
    const inner = p.svgEl.cloneNode(true);
    // Make sure the clone has explicit width/height equal to its intrinsic.
    inner.removeAttribute('width');
    inner.removeAttribute('height');
    inner.setAttribute('viewBox', `0 0 ${p.w} ${p.h}`);
    inner.setAttribute('preserveAspectRatio', 'xMinYMin meet');
    const scale = (W - padding * 2) / p.w;
    const g = document.createElementNS(NS, 'g');
    g.setAttribute('transform', `translate(${padding}, ${y}) scale(${scale})`);
    // Move all child nodes into g
    while (inner.firstChild) g.appendChild(inner.firstChild);
    out.appendChild(g);

    y += p.scaledH + gap;
  });

  return out;
}

function vqExportCompositeSVG(items, filename = 'viralquest_composite.svg', opts = {}) {
  if (!items.length) return;
  const composite = vqBuildCompositeSVG(items, opts);
  const svgStr    = new XMLSerializer().serializeToString(composite);
  const blob      = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  _download(URL.createObjectURL(blob), filename);
}

function vqExportCompositePNG(items, filename = 'viralquest_composite.png', opts = {}) {
  if (!items.length) return;
  const scale     = opts.scale || 2;
  const composite = vqBuildCompositeSVG(items, opts);
  const vb        = composite.viewBox.baseVal;
  const W         = vb.width  * scale;
  const H         = vb.height * scale;
  composite.setAttribute('width',  W);
  composite.setAttribute('height', H);

  const svgStr = new XMLSerializer().serializeToString(composite);
  const blob   = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  const url    = URL.createObjectURL(blob);

  const img = new Image();
  img.onload = () => {
    const canvas = document.createElement('canvas');
    canvas.width  = W;
    canvas.height = H;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);
    ctx.drawImage(img, 0, 0, W, H);
    URL.revokeObjectURL(url);
    _download(canvas.toDataURL('image/png'), filename);
  };
  img.onerror = () => { URL.revokeObjectURL(url); console.error('Composite PNG export failed'); };
  img.src = url;
}

// ── Text / FASTA download ───────────────────────────────────────────────────

function vqDownloadText(text, filename) {
  const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
  _download(URL.createObjectURL(blob), filename);
}

function vqBuildFasta(seqs) {
  return seqs.map(s => {
    const body = (s.sequence || '').replace(/^>.*\n/, '');
    const header = `>${s.id} ${s.taxonomy?.species || ''}`.trim();
    return body.startsWith('>') ? body : header + '\n' + body;
  }).join('\n') + '\n';
}

// ── Tooltip ─────────────────────────────────────────────────────────────────

let _tooltip = null;
function _getTooltip() {
  if (!_tooltip) _tooltip = document.getElementById('vq-tooltip');
  return _tooltip;
}

function vqTooltipShow(html, event) {
  const t = _getTooltip();
  if (!t) return;
  t.innerHTML = html;
  t.classList.add('visible');
  vqTooltipMove(event);
}

function vqTooltipMove(event) {
  const t = _getTooltip();
  if (!t || !event) return;
  const x  = event.clientX + 14;
  const y  = event.clientY + 14;
  const tw = t.offsetWidth;
  const th = t.offsetHeight;
  t.style.left = (x + tw > window.innerWidth  ? x - tw - 28 : x) + 'px';
  t.style.top  = (y + th > window.innerHeight ? y - th - 28 : y) + 'px';
}

function vqTooltipHide() {
  const t = _getTooltip();
  if (t) t.classList.remove('visible');
}

// ── Cross-section navigation (used by clusters + taxonomy) ──────────────────

function vqJumpToViewer(seqId) {
  const tab = document.querySelector('[data-section="viewer"]');
  if (tab) tab.click();
  setTimeout(() => {
    const target = document.getElementById('seq-card-' + vqSafeId(seqId));
    if (target) {
      target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      target.classList.add('vq-seq-card--highlight');
      setTimeout(() => target.classList.remove('vq-seq-card--highlight'), 2000);
    }
  }, 150);
}

// ── Internal helpers ────────────────────────────────────────────────────────

/**
 * Inline computed values for CSS custom properties used by SVG export.
 * Two-pronged: (a) drop a <style> block that sets the vars on the root,
 * and (b) write each var as an inline style on the SVG root, so any
 * descendant `fill="var(--vq-…)"` style resolves when the SVG is rendered
 * standalone (inside an <img> for PNG, or saved as .svg).
 */
function _applyInlineVars(svg) {
  const style = getComputedStyle(document.documentElement);
  const vars  = [
    '--vq-primary','--vq-primary-light','--vq-accent','--vq-accent-dark',
    '--vq-accent-subtle','--vq-success','--vq-warning','--vq-danger',
    '--vq-muted','--vq-bg','--vq-surface','--vq-surface-2','--vq-border',
    '--vq-border-dark','--vq-text','--vq-text-2','--vq-text-3',
    '--vq-font','--vq-font-mono',
    '--vq-dom-1','--vq-dom-2','--vq-dom-3','--vq-dom-4',
    '--vq-dom-5','--vq-dom-6','--vq-dom-7','--vq-dom-8',
    '--vq-tax-flaviviridae','--vq-tax-parvoviridae','--vq-tax-phenuiviridae',
    '--vq-tax-nodaviridae','--vq-tax-rhabdoviridae','--vq-tax-tombusviridae',
    '--vq-tax-other',
  ];
  const inline = vars
    .map(v => `${v}:${style.getPropertyValue(v).trim()}`)
    .filter(s => !s.endsWith(':'))
    .join(';');

  const styleEl = document.createElement('style');
  styleEl.textContent = `:root{${inline}}svg{${inline}}`;
  svg.prepend(styleEl);
  // Also apply directly to <svg> as a fallback for renderers that don't
  // resolve `:root` inside an SVG document.
  const existing = svg.getAttribute('style') || '';
  svg.setAttribute('style', (existing ? existing + ';' : '') + inline);
}

function _download(url, filename) {
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}

// ── Expose on window.VQ ─────────────────────────────────────────────────────

window.VQ = Object.assign(window.VQ || {}, {
  esc:          vqEsc,
  safeId:       vqSafeId,
  exportPNG:    vqExportPNG,
  exportSVG:    vqExportSVG,
  exportPDF:    vqExportPDF,
  exportPagePDF: vqExportPagePDF,
  exportBatchPNG: vqExportBatchPNG,
  exportCompositeSVG: vqExportCompositeSVG,
  exportCompositePNG: vqExportCompositePNG,
  buildCompositeSVG:  vqBuildCompositeSVG,
  downloadText: vqDownloadText,
  buildFasta:   vqBuildFasta,
  tooltipShow:  vqTooltipShow,
  tooltipMove:  vqTooltipMove,
  tooltipHide:  vqTooltipHide,
  jumpToViewer: vqJumpToViewer,
});
