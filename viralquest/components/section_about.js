/* This file's private helpers live in an IIFE so they don't collide across sections. */
// __VQ_WRAPPED__
(function () {
'use strict';

/* ============================================================
   section_about.js — About ViralQuest
   Logo · description · paper citation · GitHub · PyPI
   ============================================================ */

const _PAPER_URL  = 'https://doi.org/10.1186/s12859-026-06391-6';
const _GITHUB_URL = 'https://github.com/gabrielvpina/viralquest';
const _PYPI_URL   = 'https://pypi.org/project/viralquest/';

const _CITATION =
  'Rodrigues, G.V.P., Ferreira, L.Y.M. & Aguiar, E.R.G.R. ' +
  'ViralQuest: a user-friendly interactive pipeline for viral-sequences analysis and curation. ' +
  'BMC Bioinformatics 27, 64 (2026). ' +
  'https://doi.org/10.1186/s12859-026-06391-6';

function vqInitAbout() {
  const el = document.getElementById('section-about');
  if (!el) return;

  const meta    = (typeof VQ_REPORT !== 'undefined' ? VQ_REPORT.meta : null) || {};
  const version = meta.viralquest_version || '';
  const logoSrc = document.querySelector('.vq-nav__logo img')?.src || '';

  el.innerHTML = `
    <div class="vq-about">

      <!-- ── Hero ──────────────────────────────────────────── -->
      <div class="vq-about__hero">
        ${logoSrc ? `<img class="vq-about__logo" src="${logoSrc}" alt="ViralQuest">` : ''}
        <h1 class="vq-about__title">ViralQuest</h1>
        <p class="vq-about__tagline">
          A user-friendly pipeline for viral diversity discovery and characterization
        </p>
        ${version ? `<span class="vq-about__version">v${version}</span>` : ''}
      </div>

      <!-- ── Link cards ────────────────────────────────────── -->
      <div class="vq-about__links">

        <a class="vq-about__card" href="${_PAPER_URL}" target="_blank" rel="noopener noreferrer">
          <div class="vq-about__card-icon vq-about__card-icon--paper">
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">
              <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
              <polyline points="14 2 14 8 20 8"/>
              <line x1="16" y1="13" x2="8" y2="13"/>
              <line x1="16" y1="17" x2="8" y2="17"/>
              <polyline points="10 9 9 9 8 9"/>
            </svg>
          </div>
          <div class="vq-about__card-body">
            <div class="vq-about__card-title">Research Paper</div>
            <div class="vq-about__card-sub">BMC Bioinformatics · 2026</div>
          </div>
          <svg class="vq-about__card-arrow" width="16" height="16" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="7" y1="17" x2="17" y2="7"/>
            <polyline points="7 7 17 7 17 17"/>
          </svg>
        </a>

        <a class="vq-about__card" href="${_GITHUB_URL}" target="_blank" rel="noopener noreferrer">
          <div class="vq-about__card-icon vq-about__card-icon--github">
            <svg width="28" height="28" viewBox="0 0 24 24" fill="currentColor">
              <path d="M12 .297c-6.63 0-12 5.373-12 12 0 5.303 3.438 9.8 8.205 11.385.6.113.82-.258.82-.577
                       0-.285-.01-1.04-.015-2.04-3.338.724-4.042-1.61-4.042-1.61C4.422 18.07 3.633 17.7
                       3.633 17.7c-1.087-.744.084-.729.084-.729 1.205.084 1.838 1.236 1.838 1.236
                       1.07 1.835 2.809 1.305 3.495.998.108-.776.417-1.305.76-1.605-2.665-.3-5.466-1.332-5.466-5.93
                       0-1.31.465-2.38 1.235-3.22-.135-.303-.54-1.523.105-3.176 0 0 1.005-.322 3.3 1.23
                       .96-.267 1.98-.399 3-.405 1.02.006 2.04.138 3 .405 2.28-1.552 3.285-1.23 3.285-1.23
                       .645 1.653.24 2.873.12 3.176.765.84 1.23 1.91 1.23 3.22 0 4.61-2.805 5.625-5.475
                       5.92.42.36.81 1.096.81 2.22 0 1.606-.015 2.896-.015 3.286 0 .315.21.69.825.57
                       C20.565 22.092 24 17.592 24 12.297c0-6.627-5.373-12-12-12"/>
            </svg>
          </div>
          <div class="vq-about__card-body">
            <div class="vq-about__card-title">GitHub</div>
            <div class="vq-about__card-sub">Source code &amp; documentation</div>
          </div>
          <svg class="vq-about__card-arrow" width="16" height="16" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="7" y1="17" x2="17" y2="7"/>
            <polyline points="7 7 17 7 17 17"/>
          </svg>
        </a>

        <a class="vq-about__card" href="${_PYPI_URL}" target="_blank" rel="noopener noreferrer">
          <div class="vq-about__card-icon vq-about__card-icon--pypi">
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none"
                 stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">
              <path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8
                       a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/>
              <polyline points="3.27 6.96 12 12.01 20.73 6.96"/>
              <line x1="12" y1="22.08" x2="12" y2="12"/>
            </svg>
          </div>
          <div class="vq-about__card-body">
            <div class="vq-about__card-title">PyPI</div>
            <div class="vq-about__card-sub"><code class="vq-about__code">pip install viralquest</code></div>
          </div>
          <svg class="vq-about__card-arrow" width="16" height="16" viewBox="0 0 24 24" fill="none"
               stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="7" y1="17" x2="17" y2="7"/>
            <polyline points="7 7 17 7 17 17"/>
          </svg>
        </a>

      </div>

      <!-- ── Citation ───────────────────────────────────────── -->
      <div class="vq-about__cite-wrap">
        <div class="vq-card">
          <div class="vq-card__header">
            <div class="vq-card__title">Cite this work</div>
            <button class="vq-btn vq-btn--ghost vq-btn--sm" id="about-copy-btn" type="button">
              Copy citation
            </button>
          </div>
          <div class="vq-card__body">
            <blockquote class="vq-about__citation">
              Rodrigues, G.V.P., Ferreira, L.Y.M. &amp; Aguiar, E.R.G.R.
              ViralQuest: a user-friendly interactive pipeline for viral-sequences analysis and curation.
              <em>BMC Bioinformatics</em> 27, 64 (2026).
              <a class="vq-about__doi" href="${_PAPER_URL}" target="_blank" rel="noopener noreferrer">
                https://doi.org/10.1186/s12859-026-06391-6
              </a>
            </blockquote>
          </div>
        </div>
      </div>

    </div>
  `;

  // Copy citation to clipboard
  document.getElementById('about-copy-btn')?.addEventListener('click', function () {
    navigator.clipboard?.writeText(_CITATION).then(() => {
      this.textContent = 'Copied!';
      setTimeout(() => { this.textContent = 'Copy citation'; }, 2000);
    }).catch(() => {
      this.textContent = 'Copy failed';
      setTimeout(() => { this.textContent = 'Copy citation'; }, 2000);
    });
  });
}

window.vqInitAbout = vqInitAbout;
})();
