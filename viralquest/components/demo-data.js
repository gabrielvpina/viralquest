/* ============================================================
   demo-data.js — mock VQ_REPORT for the in-browser preview.
   The Python report builder will substitute this with real JSON.
   ============================================================ */

const FAMILIES = [
  'Flaviviridae','Parvoviridae','Phenuiviridae','Nodaviridae',
  'Rhabdoviridae','Tombusviridae','Other',
];
const GENERA = {
  Flaviviridae:   ['Flavivirus','Pegivirus','Pestivirus'],
  Parvoviridae:   ['Densovirus','Iteradensovirus'],
  Phenuiviridae:  ['Phlebovirus','Phasivirus'],
  Nodaviridae:    ['Alphanodavirus','Betanodavirus'],
  Rhabdoviridae:  ['Vesiculovirus','Sigmavirus'],
  Tombusviridae:  ['Tombusvirus','Carmovirus'],
  Other:          ['Unclassified'],
};

const REF_NAMES = {
  Flaviviridae:  ['Dengue virus','West Nile virus','Zika virus','Hepatitis C virus'],
  Parvoviridae:  ['Aedes aegypti densovirus','Culex pipiens densovirus'],
  Phenuiviridae: ['Phasi Charoen-like virus','Wuhan Mosquito Virus 1'],
  Nodaviridae:   ['Nodamura virus','Wuhan insect virus 12'],
  Rhabdoviridae: ['Vesicular stomatitis virus','Sigma virus'],
  Tombusviridae: ['Tomato bushy stunt virus','Carnation mottle virus'],
  Other:         ['Unclassified RNA virus','Unclassified DNA virus'],
};

function _pick(a) { return a[Math.floor(Math.random()*a.length)]; }
function _rand(lo,hi){ return lo+Math.random()*(hi-lo); }
function _randi(lo,hi){ return Math.floor(_rand(lo,hi+1)); }

function _orfs(seqLen, allowOverlap = true) {
  const AA = 'ACDEFGHIKLMNPQRSTVWY';
  const n = _randi(2, 6);
  const orfs = [];
  let cursor = 80;
  for (let i = 0; i < n; i++) {
    const orfLen = _randi(450, Math.min(4200, seqLen - cursor - 120));
    if (orfLen < 250 || cursor + orfLen > seqLen) break;
    // Some sequences have an overlapping next-ORF
    const overlap = allowOverlap && i > 0 && Math.random() < 0.45;
    const start  = overlap ? Math.max(80, cursor - _randi(50, 300)) : cursor + _randi(0, 180);
    const stop   = Math.min(seqLen - 10, start + orfLen);
    if (stop <= start + 100) break;
    const strand = Math.random() > 0.35 ? '+' : '-';
    const aaLen  = Math.floor((stop - start) / 3);

    // Domains — sometimes overlapping inside the same ORF
    const nDom = _randi(0, 4);
    const domains = [];
    let dpos = 5;
    for (let j = 0; j < nDom; j++) {
      const overlapDom = j > 0 && Math.random() < 0.35;
      const dLen = _randi(40, 160);
      const dStart = overlapDom
        ? Math.max(5, dpos - _randi(20, 80))
        : dpos + _randi(8, 50);
      const dStop = Math.min(aaLen - 5, dStart + dLen);
      if (dStop <= dStart + 20) break;
      domains.push({
        target: _pick(['Helicase','RdRp','Methyltransferase','Capsid','Protease',
                       'Glycoprotein','Polymerase','Envelope','RNase','Zinc finger']),
        database: _pick(['Pfam','RVDB','Vfam']),
        description: 'mock domain hit',
        start: dStart, stop: dStop,
        score: +(_rand(40, 220)).toFixed(1),
        e_value: Math.pow(10, -_randi(5, 60)),
      });
      dpos = dStop + 8;
    }

    // Synthetic AA sequence (length_aa residues, capped at 600 for the preview)
    const previewLen = Math.min(aaLen, 600);
    let aa = 'M';
    for (let k = 1; k < previewLen - 1; k++) aa += AA[Math.floor(Math.random() * AA.length)];
    aa += '*';

    orfs.push({
      name: `ORF${i + 1}`,
      strand,
      frame: (strand === '+' ? '+' : '-') + ((i % 3) + 1),
      length_aa: aaLen,
      start_position: start,
      stop_position:  stop,
      orf_type: domains.length ? 'putative protein' : 'hypothetical',
      domains,
      aa_sequence: aa,
    });
    cursor = stop + _randi(10, 100);
  }
  return orfs;
}

function _fasta(id, len) {
  const bases = 'ACGT';
  let s = '';
  for (let i = 0; i < len; i++) s += bases[Math.floor(Math.random() * 4)];
  // wrap at 60 chars
  const lines = [`>${id}`];
  for (let i = 0; i < s.length; i += 60) lines.push(s.slice(i, i + 60));
  return lines.join('\n');
}

function _blastHits(seqLen, kind, family) {
  const n = _randi(3, 7);
  const hits = [];
  const refs = REF_NAMES[family] || REF_NAMES.Other;
  for (let i = 0; i < n; i++) {
    const pid = _rand(60, 99) - i * _rand(1, 4);
    const cov = _rand(35, 95) - i * _rand(0, 3);
    const accession = (kind === 'blastn'
      ? ['NC_', 'KY', 'MN', 'MK', 'OL', 'MT']
      : ['YP_', 'NP_', 'QHR', 'AAB', 'AGW'])[_randi(0,5)]
      + _randi(100000, 999999);
    hits.push({
      accession,
      title: `${_pick(refs)} ${kind === 'blastn' ? 'genome' : 'polyprotein'}`,
      pident:    +pid.toFixed(1),
      qcovhsp:   +cov.toFixed(0),
      evalue:    Math.pow(10, -_randi(10, 180)),
      bitscore:  _randi(120, 2500) - i * 50,
      length:    _randi(50, 2500),
      qstart:    _randi(1, Math.max(2, Math.floor(seqLen * 0.1))),
      qend:      _randi(Math.floor(seqLen * 0.4), seqLen - 1),
      sstart:    _randi(1, 200),
      send:      _randi(200, 5000),
    });
  }
  return hits.sort((a, b) => b.bitscore - a.bitscore);
}

function _makeSequences() {
  const seqs = [];
  const N = 18;
  for (let i = 0; i < N; i++) {
    const family = _pick(FAMILIES);
    const genus  = _pick(GENERA[family]);
    const seqLen = _randi(2500, 12500);
    const id     = `contig_${String(i + 1).padStart(3, '0')}`;
    const cls    = ['viral-known','viral-unknown','non-viral'][_randi(0, 2)];
    const score  = cls === 'viral-known' ? _randi(70, 99)
                 : cls === 'viral-unknown' ? _randi(40, 75)
                 : _randi(0, 35);

    seqs.push({
      id,
      length: seqLen,
      gc_content: _rand(28, 62),
      orfs: _orfs(seqLen),
      taxonomy: {
        phylum:  'Riboviria',
        order:   _pick(['Mononegavirales','Picornavirales','Nidovirales','Tymovirales']),
        family,
        genus,
        species: `${genus} sp. ${i + 1}`,
      },
      llm_output: {
        classification: cls,
        vq_score: score,
        analysis: `Sequence aligns to multiple ${family} reference genomes with ` +
                  `${cls === 'non-viral' ? 'borderline' : 'strong'} support. Predicted ORFs ` +
                  `include canonical replication-associated proteins; domain coverage ` +
                  `is consistent with the expected polyprotein layout.`,
      },
      is_viral:     cls !== 'non-viral',
      sequence_nt:  _fasta(id, Math.min(seqLen, 600)),  // truncated preview
      blastn_hits:  _blastHits(seqLen, 'blastn', family),
      blastx_refseq_hits: _blastHits(seqLen, 'blastx', family),
      blastx_nr_hits:     _blastHits(seqLen, 'blastx', family),
    });
  }

  // Assign cluster_ids
  seqs[0].cluster_id = 'CL1'; seqs[1].cluster_id = 'CL1'; seqs[2].cluster_id = 'CL1';
  seqs[3].cluster_id = 'CL2'; seqs[4].cluster_id = 'CL2';
  seqs[5].cluster_id = 'CL2'; seqs[6].cluster_id = 'CL2';
  seqs[7].cluster_id = 'CL3'; seqs[8].cluster_id = 'CL3';

  // Force one big genome to demonstrate the > 20 kb horizontal scroller.
  const big = seqs[9];
  big.length = 26_400;
  big.orfs   = _orfs(big.length);
  big.sequence_nt = _fasta(big.id, 600);

  return seqs;
}

function _makeClusters(seqs) {
  const byCid = {};
  seqs.forEach(s => {
    if (!s.cluster_id) return;
    (byCid[s.cluster_id] ??= []).push(s);
  });
  return Object.keys(byCid).map(cid => {
    const members = byCid[cid];
    return {
      cluster_id: cid,
      species: members[0].taxonomy.species.replace(/sp\. \d+/, 'sp.'),
      size: members.length,
      members: members.map((m, i) => ({
        seq_id: m.id,
        length: m.length,
        identity: i === 0 ? 100 : _rand(72, 99),
        query_coverage: _rand(60, 99),
        aln_start: _randi(1, 200),
        aln_end:   m.length - _randi(0, 200),
        is_representative: i === 0,
      })),
    };
  });
}

function _makeSalmon(seqs) {
  const viral_quant = seqs.map(s => ({
    name: s.id,
    tpm: _rand(0.1, 800),
    num_reads: _randi(10, 50000),
  }));
  const kingdoms = ['Bacteria','Archaea','Eukaryota','Viruses'];
  const conserved_quant = [];
  kingdoms.forEach(k => {
    const n = _randi(4, 9);
    for (let i = 0; i < n; i++) {
      conserved_quant.push({
        name: `${k.slice(0,3).toLowerCase()}_${_pick(['actb','gapdh','rpl19','tubA','hsp70','elf1a','rps18'])}_${i+1}`,
        kingdom: k,
        tpm: _rand(5, 1200),
        num_reads: _randi(100, 200000),
      });
    }
  });
  const host_viral_hits = seqs.slice(0, 4).map((s, i) => ({
    transcript: { name: `host_tx_${i + 1}`, tpm: _rand(2, 80) },
    blast_hits: [{
      viral_seq_id: s.id,
      pident: _rand(72, 99),
      qcovhsp: _randi(40, 95),
    }],
  }));
  return {
    mapping_rate: _rand(35, 78),
    total_reads:  _randi(2_000_000, 14_000_000),
    viral_quant,
    conserved_quant,
    host_viral_hits,
  };
}

const _SEQS     = _makeSequences();
const _CLUSTERS = _makeClusters(_SEQS);
const _SALMON   = _makeSalmon(_SEQS);

const confirmedViral = _SEQS.filter(s => s.is_viral).length;
const totalOrfs      = _SEQS.reduce((a, s) => a + s.orfs.length, 0);

window.VQ_REPORT = {
  meta: {
    timestamp: '2026-05-20T14:32:00Z',
    viralquest_version: '0.8.2',
    input_file: { name: 'mosquito_pool_A.fasta', size_bytes: 18_400_000 },
    sample_name: 'mosquito_pool_A',
  },
  summary: {
    total_sequences: 124,                 // pretend the raw input had 124 reads/contigs
    confirmed_viral: confirmedViral,
    total_orfs:      totalOrfs,
    total_clusters:  _CLUSTERS.length,
    cap3_used: true,
    cap3_contigs:  142,
    cap3_singlets: 387,
  },
  blast_stats: {
    refseq_unique_seqs: 24,
    nr_unique_seqs:     18,
    blastn_unique_seqs: 16,
  },
  hmm_stats: {
    rvdb_hits: 41,
    vfam_hits: 28,
    eggnog_hits: 17,
    pfam_hits: 53,
  },
  llm_stats: {
    present: true,
    model: 'claude-haiku-4-5',
    mode:  'batch',
    scored: _SEQS.length,
    viral_known:   _SEQS.filter(s => s.llm_output.classification === 'viral-known').length,
    viral_unknown: _SEQS.filter(s => s.llm_output.classification === 'viral-unknown').length,
    avg_score: +(_SEQS.reduce((a,s)=>a+s.llm_output.vq_score,0)/_SEQS.length).toFixed(1),
  },
  salmon_stats: {
    present: true,
    mapping_rate: _SALMON.mapping_rate,
    total_reads:  _SALMON.total_reads,
    host_viral_hits: _SALMON.host_viral_hits.length,
  },
  sequences:    _SEQS,
  clusters:     _CLUSTERS,
  salmon_quant: _SALMON,
};
