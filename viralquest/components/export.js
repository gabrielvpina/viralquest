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
    const body = (s.sequence || s.sequence_nt || '').replace(/^>.*\n/, '');
    const header = `>${s.id}`;   // original contig name only (no species suffix)
    return body.startsWith('>') ? body : header + '\n' + body;
  }).join('\n') + '\n';
}

// ── Sequence table (TSV / XLSX) ─────────────────────────────────────────────
// One row per sequence. Mirrors viralquest/seq_table.py (the table written to
// the run's output folder) — keep the columns and rules in sync.

const VQ_TABLE_COLUMNS = [
  'Sample', 'Order', 'Family', 'Genus', 'Species', 'Genome type',
  'Sequence (nt)', 'Length', 'Domains',
  'BLASTn hit', 'BLASTn identity (%)', 'BLASTn coverage (%)',
  'BLASTn hit ID', 'BLASTn hit sequence', 'BLASTn e-value',
  'BLASTx hit', 'BLASTx identity (%)', 'BLASTx coverage (%)',
  'BLASTx hit ID', 'BLASTx hit sequence', 'BLASTx e-value',
];
// Accession columns → the NCBI database their record lives in.
const _TABLE_LINKS = { 'BLASTn hit ID': 'nuccore', 'BLASTx hit ID': 'protein' };
// Excel refuses longer cells; such values stay complete in the TSV.
const _XLSX_CELL_LIMIT = 32767;

/** NCBI record URL for a bare accession; null for pipe-packed / local ids. */
function vqNcbiUrl(accession, db) {
  if (!accession || !/^[A-Za-z0-9_]+(\.[0-9]+)?$/.test(accession)) return null;
  return `https://www.ncbi.nlm.nih.gov/${db}/${accession}`;
}

function _tableBest(hits) {
  return hits && hits.length
    ? hits.reduce((a, b) => ((b.bit_score || 0) > (a.bit_score || 0) ? b : a))
    : null;
}

function _tableSample(s) {
  if (s.sample || s.sample_name) return s.sample || s.sample_name;
  const meta = (typeof VQ_REPORT !== 'undefined' && VQ_REPORT.meta) || {};
  const name = (meta.input_file && meta.input_file.name) || '';
  return name.replace(/\.[^.]*$/, '');
}

function _tableDomains(s) {
  const seen = [];
  (s.orfs || []).forEach(o => (o.domains || []).forEach(d => {
    const label = `${d.target || ''} (${d.database || ''})`;
    if (!seen.includes(label)) seen.push(label);
  }));
  return seen.join('; ');
}

/** Table rows (column name → value) for exported sequence records. */
function vqSeqTableRows(seqs) {
  return seqs.map(s => {
    const bt  = s.blastn_taxonomy || {};
    const t   = s.taxonomy || bt.taxonomy || {};
    const bn  = _tableBest(s.blastn_hits);
    const bx  = _tableBest(s.blastx_nr_hits && s.blastx_nr_hits.length ? s.blastx_nr_hits : s.blastx_hits);
    const g   = (h, k) => (h && h[k] != null ? h[k] : null);
    return {
      'Sample':              _tableSample(s),
      'Order':               t.order || '',
      'Family':              t.family || '',
      'Genus':               t.genus || '',
      'Species':             t.species || t.scientific_name || '',
      'Genome type':         t.genome || '',
      'Sequence (nt)':       s.sequence || '',
      'Length':              s.length ?? null,
      'Domains':             _tableDomains(s),
      'BLASTn hit':          g(bn, 'stitle') || '',
      'BLASTn identity (%)': g(bn, 'pident'),
      'BLASTn coverage (%)': g(bn, 'qcovhsp'),
      'BLASTn hit ID':       g(bn, 'accession') || '',
      'BLASTn hit sequence': g(bn, 'subject_seq') || '',
      'BLASTn e-value':      g(bn, 'evalue'),
      'BLASTx hit':          g(bx, 'subject_title') || '',
      'BLASTx identity (%)': g(bx, 'pct_identity'),
      'BLASTx coverage (%)': g(bx, 'query_coverage'),
      'BLASTx hit ID':       g(bx, 'subject_id') || '',
      'BLASTx hit sequence': g(bx, 'subject_seq') || '',
      'BLASTx e-value':      g(bx, 'e_value'),
    };
  });
}

function vqBuildTSV(rows) {
  const cell = v => (v == null ? '' : String(v).replace(/[\t\r\n]+/g, ' '));
  return [VQ_TABLE_COLUMNS.join('\t')]
    .concat(rows.map(r => VQ_TABLE_COLUMNS.map(c => cell(r[c])).join('\t')))
    .join('\n') + '\n';
}

// ── XLSX: minimal SpreadsheetML in an uncompressed ZIP (no library) ─────────

const _xmlEsc = s => String(s)
  .replace(/[\x00-\x08\x0b\x0c\x0e-\x1f]/g, '')
  .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
  .replace(/"/g, '&quot;');

function _xlsxCol(i) {
  let s = '';
  for (i += 1; i; i = Math.floor((i - 1) / 26)) s = String.fromCharCode(65 + (i - 1) % 26) + s;
  return s;
}

function _xlsxCell(ref, v, style = 0) {
  const st = style ? ` s="${style}"` : '';
  if (v == null || v === '') return style ? `<c r="${ref}"${st}/>` : '';
  if (typeof v === 'number' && isFinite(v)) return `<c r="${ref}"${st}><v>${v}</v></c>`;
  let text = String(v);
  if (text.length > _XLSX_CELL_LIMIT) {
    text = `[${text.length.toLocaleString('en-US')} characters: too long for an Excel cell, see the TSV]`;
  }
  return `<c r="${ref}"${st} t="inlineStr"><is><t xml:space="preserve">${_xmlEsc(text)}</t></is></c>`;
}

const _XLSX_WIDTHS = {
  'Sample': 18, 'Species': 30, 'BLASTn hit': 40, 'BLASTx hit': 40, 'Domains': 30,
  'Sequence (nt)': 30, 'BLASTn hit sequence': 30, 'BLASTx hit sequence': 30,
  'BLASTn hit ID': 16, 'BLASTx hit ID': 16,
};

function _xlsxParts(rows, sheet = 'Sequences') {
  const X   = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
  const ns  = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"';
  const rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"';
  const rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
  const pkg = 'http://schemas.openxmlformats.org/package/2006/relationships';
  const ct  = 'application/vnd.openxmlformats-officedocument.spreadsheetml';
  const last = _xlsxCol(VQ_TABLE_COLUMNS.length - 1);
  const lastRow = Math.max(1, rows.length + 1);

  const xmlRows = ['<row r="1">' + VQ_TABLE_COLUMNS.map((c, i) => _xlsxCell(`${_xlsxCol(i)}1`, c, 1)).join('') + '</row>'];
  const links = [];
  rows.forEach((r, k) => {
    const n = k + 2;
    const cells = VQ_TABLE_COLUMNS.map((c, i) => {
      const ref = `${_xlsxCol(i)}${n}`;
      const url = _TABLE_LINKS[c] ? vqNcbiUrl(r[c], _TABLE_LINKS[c]) : null;
      if (url) links.push([ref, url]);
      return _xlsxCell(ref, r[c], url ? 2 : 0);
    });
    xmlRows.push(`<row r="${n}">${cells.join('')}</row>`);
  });

  const cols = VQ_TABLE_COLUMNS.map((c, i) =>
    `<col min="${i + 1}" max="${i + 1}" width="${_XLSX_WIDTHS[c] || 12}" customWidth="1"/>`).join('');
  const hyper = links.length
    ? '<hyperlinks>' + links.map(([ref], k) => `<hyperlink ref="${ref}" r:id="rId${k + 1}"/>`).join('') + '</hyperlinks>'
    : '';

  const parts = {
    '[Content_Types].xml': X +
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' +
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' +
      '<Default Extension="xml" ContentType="application/xml"/>' +
      `<Override PartName="/xl/workbook.xml" ContentType="${ct}.sheet.main+xml"/>` +
      `<Override PartName="/xl/worksheets/sheet1.xml" ContentType="${ct}.worksheet+xml"/>` +
      `<Override PartName="/xl/styles.xml" ContentType="${ct}.styles+xml"/>` +
      '</Types>',
    '_rels/.rels': X + `<Relationships xmlns="${pkg}">` +
      `<Relationship Id="rId1" Type="${rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>`,
    'xl/workbook.xml': X + `<workbook ${ns} ${rns}>` +
      `<sheets><sheet name="${_xmlEsc(sheet)}" sheetId="1" r:id="rId1"/></sheets>` +
      '<definedNames><definedName name="_xlnm._FilterDatabase" localSheetId="0" hidden="1">' +
      `'${_xmlEsc(sheet)}'!$A$1:$${last}$${lastRow}</definedName></definedNames></workbook>`,
    'xl/_rels/workbook.xml.rels': X + `<Relationships xmlns="${pkg}">` +
      `<Relationship Id="rId1" Type="${rel}/worksheet" Target="worksheets/sheet1.xml"/>` +
      `<Relationship Id="rId2" Type="${rel}/styles" Target="styles.xml"/></Relationships>`,
    'xl/styles.xml': X +
      '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' +
      '<fonts count="3"><font><sz val="11"/><name val="Calibri"/></font>' +
      '<font><b/><sz val="11"/><name val="Calibri"/></font>' +
      '<font><u/><sz val="11"/><color rgb="FF0563C1"/><name val="Calibri"/></font></fonts>' +
      '<fills count="2"><fill><patternFill patternType="none"/></fill>' +
      '<fill><patternFill patternType="gray125"/></fill></fills>' +
      '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>' +
      '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' +
      '<cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>' +
      '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>' +
      '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>' +
      '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' +
      '</styleSheet>',
    'xl/worksheets/sheet1.xml': X + `<worksheet ${ns} ${rns}>` +
      '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" ' +
      'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>' +
      `<cols>${cols}</cols><sheetData>${xmlRows.join('')}</sheetData>` +
      `<autoFilter ref="A1:${last}${lastRow}"/>${hyper}</worksheet>`,
  };
  if (links.length) {
    parts['xl/worksheets/_rels/sheet1.xml.rels'] = X + `<Relationships xmlns="${pkg}">` +
      links.map(([, url], k) =>
        `<Relationship Id="rId${k + 1}" Type="${rel}/hyperlink" Target="${_xmlEsc(url)}" TargetMode="External"/>`
      ).join('') + '</Relationships>';
  }
  return parts;
}

let _crcTable = null;
function _crc32(bytes) {
  if (!_crcTable) {
    _crcTable = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1;
      _crcTable[n] = c >>> 0;
    }
  }
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) crc = _crcTable[(crc ^ bytes[i]) & 0xFF] ^ (crc >>> 8);
  return (crc ^ 0xFFFFFFFF) >>> 0;
}

/** Uncompressed (stored) ZIP of { name: text } → Uint8Array. */
function _zipStore(files) {
  const enc = new TextEncoder();
  const chunks = [], central = [];
  let offset = 0;
  const u16 = (v, o, d) => d.setUint16(o, v, true);
  const u32 = (v, o, d) => d.setUint32(o, v >>> 0, true);
  Object.entries(files).forEach(([name, text]) => {
    const nameB = enc.encode(name), data = enc.encode(text), crc = _crc32(data);
    const lh = new DataView(new ArrayBuffer(30));
    u32(0x04034b50, 0, lh); u16(20, 4, lh); u16(0x0800, 6, lh); u16(0, 8, lh);
    u16(0, 10, lh); u16(0x21, 12, lh); u32(crc, 14, lh);
    u32(data.length, 18, lh); u32(data.length, 22, lh); u16(nameB.length, 26, lh); u16(0, 28, lh);
    const ch = new DataView(new ArrayBuffer(46));
    u32(0x02014b50, 0, ch); u16(20, 4, ch); u16(20, 6, ch); u16(0x0800, 8, ch); u16(0, 10, ch);
    u16(0, 12, ch); u16(0x21, 14, ch); u32(crc, 16, ch); u32(data.length, 20, ch);
    u32(data.length, 24, ch); u16(nameB.length, 28, ch); u32(offset, 42, ch);
    chunks.push(new Uint8Array(lh.buffer), nameB, data);
    central.push(new Uint8Array(ch.buffer), nameB);
    offset += 30 + nameB.length + data.length;
  });
  const cdSize = central.reduce((a, c) => a + c.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  const n = Object.keys(files).length;
  u32(0x06054b50, 0, end); u16(n, 8, end); u16(n, 10, end); u32(cdSize, 12, end); u32(offset, 16, end);
  const all = chunks.concat(central, [new Uint8Array(end.buffer)]);
  const out = new Uint8Array(all.reduce((a, c) => a + c.length, 0));
  let p = 0;
  all.forEach(c => { out.set(c, p); p += c.length; });
  return out;
}

function vqBuildXLSX(rows) {
  return _zipStore(_xlsxParts(rows));
}

/** Download the sequence table for `seqs` as `<basename>.tsv` / `.xlsx`. */
function vqExportSeqTable(seqs, basename, fmt = 'xlsx') {
  if (!seqs || !seqs.length) return;
  const rows = vqSeqTableRows(seqs);
  if (fmt === 'tsv') {
    const blob = new Blob([vqBuildTSV(rows)], { type: 'text/tab-separated-values;charset=utf-8' });
    _download(URL.createObjectURL(blob), basename + '.tsv');
    return;
  }
  const blob = new Blob([vqBuildXLSX(rows)], {
    type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
  });
  _download(URL.createObjectURL(blob), basename + '.xlsx');
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
  // Two frames: the viewer panel must be shown (laid out) before measuring.
  requestAnimationFrame(() => requestAnimationFrame(() => {
    // The viewer opens the card (and lifts filters hiding it) when it can.
    const card = window.VQ?.viewerFocus
      ? window.VQ.viewerFocus(seqId)
      : document.getElementById('seq-card-' + vqSafeId(seqId));
    if (!card) return;

    // Land the card just below the fixed nav bar, which scrollIntoView ignores.
    const gap  = (document.querySelector('.vq-nav')?.offsetHeight || 0) + 16;
    const goTo = () => window.scrollTo({
      top: Math.max(0, card.getBoundingClientRect().top + window.scrollY - gap),
      behavior: 'smooth',
    });
    goTo();
    // Content above may still settle while scrolling (lazy charts): correct once.
    setTimeout(() => {
      if (Math.abs(card.getBoundingClientRect().top - gap) > 8) goTo();
    }, 750);

    // Highlight the whole card for 6 s (restart if it is already highlighted).
    card.classList.remove('vq-seq-card--highlight');
    void card.offsetWidth;
    card.classList.add('vq-seq-card--highlight');
    clearTimeout(card._vqHighlight);
    card._vqHighlight = setTimeout(() => card.classList.remove('vq-seq-card--highlight'), 6000);
  }));
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

// ── Overview tabs (Run · Virome) — card / layout helpers ─────────────────────

/** KPI chip: label, big value and an optional sub-line. */
function vqStatChip(label, value, mod = '', sub = '') {
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${vqEsc(label)}</div>
      <div class="vq-stat-chip__value" title="${vqEsc(value)}">${value}</div>
      ${sub ? `<div class="vq-stat-chip__sub" title="${vqEsc(sub)}">${sub}</div>` : ''}
    </div>`;
}

/** Locale-formatted number, or an em dash when missing. */
function vqFmtNum(n) {
  return (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();
}

/** a / b as a one-decimal percentage string ('—' when b is 0 or missing). */
function vqPct(a, b) {
  if (!b || isNaN(a) || isNaN(b)) return '—';
  return ((a / b) * 100).toFixed(1) + '%';
}

/** Titled group of cards on a grid. `slots` = grid columns the cards fill
    (a span-2 card counts as 2): 1–2 slots use a 2-column grid so a short group
    does not leave an empty third column, and 4 slots a 2 × 2 grid rather than
    3 + 1 (a lone card on its own row). Empty groups vanish. */
function vqCardGroup(title, hint, cardsHtml, slots = 3) {
  if (!cardsHtml.trim() || !slots) return '';
  return `
    <section class="vq-group">
      <div class="vq-group__head">
        <h3 class="vq-group__title">${vqEsc(title)}</h3>
        ${hint ? `<span class="vq-group__hint">${vqEsc(hint)}</span>` : ''}
      </div>
      <div class="vq-grid${slots <= 2 || slots === 4 ? ' vq-grid--2' : ''}">${cardsHtml}</div>
    </section>`;
}

/** Re-render a panel's charts when it is first shown or resized: a hidden tab
    has zero width, so its charts are drawn at a fallback size until revealed. */
function vqRedrawOnResize(el, redraw) {
  // One observer per element: a tab re-mounted with new data (e.g. the
  // multi-sample report's sample filter) replaces the previous observer, whose
  // closure would otherwise redraw the old data on the next resize.
  if (el.__vqResizeObserver) el.__vqResizeObserver.disconnect();
  // Starts at 0 so the observer's first callback (right after mounting) redraws
  // at the settled width — e.g. once a scrollbar has appeared.
  let lastW = 0;
  const ro = new ResizeObserver(() => {
    const w = el.clientWidth;
    if (!w || w === lastW) return;
    lastW = w;
    redraw();
  });
  ro.observe(el);
  el.__vqResizeObserver = ro;
}

/** Average-linkage (UPGMA) clustering of n items from a similarity matrix
    (1 = identical).  Returns the root {leaves, h, kids}; root.leaves is the
    display order.  O(n³) — fine for the tens of samples / rows a report has. */
function vqUpgma(sim) {
  let nodes = sim.map((_, i) => ({ leaves: [i], h: 0, kids: null }));
  const dist = (A, B) => {
    let s = 0;
    A.leaves.forEach(i => B.leaves.forEach(j => { s += 1 - sim[i][j]; }));
    return s / (A.leaves.length * B.leaves.length);
  };
  while (nodes.length > 1) {
    let best = null;
    for (let a = 0; a < nodes.length; a++) for (let b = a + 1; b < nodes.length; b++) {
      const d = dist(nodes[a], nodes[b]);
      if (!best || d < best.d) best = { a, b, d };
    }
    const A = nodes[best.a], B = nodes[best.b];
    nodes = nodes.filter((_, k) => k !== best.a && k !== best.b)
      .concat([{ leaves: A.leaves.concat(B.leaves), h: best.d, kids: [A, B] }]);
  }
  return nodes[0] || { leaves: [], h: 0, kids: null };
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
  ncbiUrl:        vqNcbiUrl,
  seqTableRows:   vqSeqTableRows,
  buildTSV:       vqBuildTSV,
  buildXLSX:      vqBuildXLSX,
  exportSeqTable: vqExportSeqTable,
  tooltipShow:  vqTooltipShow,
  tooltipMove:  vqTooltipMove,
  tooltipHide:  vqTooltipHide,
  jumpToViewer: vqJumpToViewer,
  statChip:       vqStatChip,
  fmtNum:         vqFmtNum,
  pct:            vqPct,
  cardGroup:      vqCardGroup,
  redrawOnResize: vqRedrawOnResize,
  upgma:          vqUpgma,
});
