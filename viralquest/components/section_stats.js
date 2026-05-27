/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_stats.js — Section 1: General Statistics
   Dashboard of cards, one per pipeline step.
   ============================================================ */

function vqInitStats(report) {
  const el = document.getElementById('section-stats');
  if (!el) return;

  const esc    = VQ.esc;
  const meta   = report.meta         || {};
  const sum    = report.summary      || {};
  const blast  = report.blast_stats  || {};
  const hmm    = report.hmm_stats    || {};
  const llm    = report.llm_stats    || {};
  const salmon = report.salmon_stats || {};

  // ── Formatters ─────────────────────────────────────────────────────────
  const fmtBytes = b => {
    if (b == null || isNaN(b)) return '—';
    return b > 1e6 ? (b / 1e6).toFixed(1) + ' MB'
                   : (b / 1e3).toFixed(1) + ' KB';
  };
  const fmtNum = n => (n == null || isNaN(n)) ? '—' : Number(n).toLocaleString();
  const fmtPct = v => (v == null || isNaN(v)) ? '—' : Number(v).toFixed(1) + '%';
  const pct    = (a, b) => (b && !isNaN(a) && !isNaN(b))
                            ? ((a / b) * 100).toFixed(1) + '%'
                            : '—';

  // ── Aggregate from sequences/clusters when summary is missing ──────────
  const seqs     = report.sequences || [];
  const clusters = report.clusters  || [];

  const totalSeqs      = sum.total_sequences ?? seqs.length;
  const confirmedViral = sum.confirmed_viral ?? seqs.filter(s => s.is_viral).length;
  const totalOrfs      = sum.total_orfs ?? seqs.reduce((a, s) => a + (s.orfs || []).length, 0);
  const totalClusters  = sum.total_clusters ?? clusters.length;

  // Tolerate both `eggnog_hits` and the legacy typo `eggnong_hits`.
  const eggnogHits = hmm.eggnog_hits ?? hmm.eggnong_hits;

  el.innerHTML = `
    <div class="vq-section-header">
      <div>
        <div class="vq-section-title">General Statistics</div>
        <div class="vq-section-sub">
          Run on ${esc(meta.timestamp ? new Date(meta.timestamp).toLocaleString() : '—')}
          &nbsp;·&nbsp; ViralQuest v${esc(meta.viralquest_version || '?')}
        </div>
      </div>
      <div class="vq-section-actions">
        <button class="vq-btn vq-btn--ghost" id="stats-print-pdf" type="button">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
            <polyline points="14 2 14 8 20 8"/>
          </svg>
          Print / Save PDF
        </button>
      </div>
    </div>

    ${_card('Input File', `
      <div class="vq-stats-grid">
        ${_chip('Filename',        esc(meta.input_file?.name || '—'), '', '', true)}
        ${_chip('Total Sequences', fmtNum(totalSeqs), 'accent')}
        ${_chip('File Size',       fmtBytes(meta.input_file?.size_bytes))}
        ${_chip('CAP3 Assembly',   sum.cap3_used ? 'Yes' : 'No', sum.cap3_used ? 'success' : '')}
        ${sum.cap3_used ? _chip('CAP3 Contigs',  fmtNum(sum.cap3_contigs))  : ''}
        ${sum.cap3_used ? _chip('CAP3 Singlets', fmtNum(sum.cap3_singlets)) : ''}
      </div>
    `)}

    ${_card('BLAST Results', `
      <div class="vq-stats-grid">
        ${_chip('RefSeq Candidates',  fmtNum(blast.refseq_unique_seqs), 'accent',
                pct(blast.refseq_unique_seqs, totalSeqs) + ' of input')}
        ${_chip('NR Confirmed Viral', fmtNum(blast.nr_unique_seqs), 'success',
                pct(blast.nr_unique_seqs, blast.refseq_unique_seqs) + ' of candidates')}
        ${_chip('BLASTn Hits',        fmtNum(blast.blastn_unique_seqs), '',
                'on confirmed sequences')}
        ${_chip('Total Confirmed',    fmtNum(confirmedViral), 'success',
                pct(confirmedViral, totalSeqs) + ' of input')}
      </div>
    `)}

    ${_card('HMM Domain Annotation', `
      <div class="vq-stats-grid">
        ${_chip('RVDB Hits',   fmtNum(hmm.rvdb_hits),  'accent',  'viral filter DB')}
        ${_chip('Vfam Hits',   fmtNum(hmm.vfam_hits),  'accent',  'viral filter DB')}
        ${_chip('EggNOG Hits', fmtNum(eggnogHits),     '',        'viral filter DB')}
        ${_chip('Pfam Hits',   fmtNum(hmm.pfam_hits),  'success', 'functional characterisation')}
        ${_chip('Total ORFs',  fmtNum(totalOrfs),      '',        'across confirmed sequences')}
        ${_chip('Clusters',    fmtNum(totalClusters),  'accent',  'unique viral species')}
      </div>
    `)}

    ${llm.present ? _card('LLM Scoring', `
      <div class="vq-stats-grid">
        ${_chip('Model',            esc(llm.model || '—'), '', '', true)}
        ${_chip('Mode',             esc(llm.mode  || '—'), '', '', true)}
        ${_chip('Sequences Scored', fmtNum(llm.scored),       'accent')}
        ${_chip('Viral Known',      fmtNum(llm.viral_known),  'success',
                pct(llm.viral_known, llm.scored))}
        ${_chip('Viral Unknown',    fmtNum(llm.viral_unknown), 'warn',
                pct(llm.viral_unknown, llm.scored))}
        ${_chip('Avg VQ Score',
                llm.avg_score != null ? Number(llm.avg_score).toFixed(1) : '—',
                'accent')}
      </div>
    `) : ''}

    ${salmon.present ? _card('Salmon Quantification', `
      <div class="vq-stats-grid">
        ${_chip('Mapping Rate', fmtPct(salmon.mapping_rate),
                (salmon.mapping_rate != null && salmon.mapping_rate > 30) ? 'success' : 'warn')}
        ${_chip('Total Reads',     fmtNum(salmon.total_reads))}
        ${_chip('Host–Viral Hits', fmtNum(salmon.host_viral_hits), 'accent',
                'transcriptome seqs with viral similarity')}
      </div>
    `) : ''}
  `;

  document.getElementById('stats-print-pdf')?.addEventListener('click', () => {
    VQ.exportPagePDF('ViralQuest Statistics');
  });
}

// ── Private helpers ────────────────────────────────────────────────────────

function _card(title, bodyHtml) {
  return `
    <div class="vq-card" style="margin-bottom:var(--vq-space-4)">
      <div class="vq-card__header">
        <div class="vq-card__title">${VQ.esc(title)}</div>
      </div>
      <div class="vq-card__body">${bodyHtml}</div>
    </div>`;
}

function _chip(label, value, mod = '', sub = '', isText = false) {
  return `
    <div class="vq-stat-chip${mod ? ' vq-stat-chip--' + mod : ''}">
      <div class="vq-stat-chip__label">${VQ.esc(label)}</div>
      <div class="vq-stat-chip__value${isText ? ' vq-stat-chip__value--text' : ''}"
           title="${VQ.esc(value)}">${value}</div>
      ${sub ? `<div class="vq-stat-chip__sub" title="${VQ.esc(sub)}">${VQ.esc(sub)}</div>` : ''}
    </div>`;
}


window.vqInitStats = vqInitStats;
})();
