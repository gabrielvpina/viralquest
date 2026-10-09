import re
import uuid
import pathlib as Path
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Input FASTA file
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class InputFasta:
    """FASTA file used in input"""
    name: str
    num_seqs: int
    path: str
    size: int # bytes



# ---------------------------------------------------------------------------
# CAP3 DATA
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class Cap3:
    input_fasta: Path
    contigs: Path
    singlets: Path
    info: Path
    ace: Path
    log: Path
    
    @property
    def is_successful(self) -> bool:
        """verify all files"""
        return self.contigs.exists() and self.singlets.exists()

    def get_combined_fasta(self, output_path: Path) -> Path:
        """combines the contigs and singlets files"""
        with open(output_path, 'w') as outfile:
            if self.contigs.exists():
                outfile.write(self.contigs.read_text())
            if self.singlets.exists():
                outfile.write(self.singlets.read_text())
        return output_path



# ---------------------------------------------------------------------------
# HMM domain hit
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class HmmDomain:
    """One pyhmmer domain hit attached to an ORF."""
    database: str
    target: str
    score: float
    e_value: float
    start: int
    stop: int
    length: int
    description: str
    type: str
    details: str


# ---------------------------------------------------------------------------
# Open Reading Frame
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class Orf:
    """Stores a single Open Reading Frame and its annotation results."""
    start_codon: str
    stop_codon: str
    start_position: int
    stop_position: int
    strand: str
    frame: int
    length_aa: int
    length_nt: int
    bigger_than_50: bool
    aa_sequence: str
    nuc_sequence: str
    orf_type: str
    # ids
    name: str
    uid: uuid.UUID = field(default_factory=uuid.uuid4, init=False)
    # hmm data - composition
    domains: list[HmmDomain] = field(default_factory=list, init=False)
    # raw HMM hit count per database (score >= threshold, before positional
    # de-dup) — quantitative signal for the heuristic scorer. Keyed by database
    # name ("RVDB", "Vfam", "EggNOG", "Pfam").
    raw_hmm_counts: dict[str, int] = field(default_factory=dict, init=False)


# ---------------------------------------------------------------------------
# BLASTx result  (Diamond)
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class BlastxResult:
    """One Diamond BLASTx hit."""
    query_id: str
    subject_id: str
    subject_title: str          # fetched separately or from --outfmt stitle
    pct_identity: float
    aln_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    e_value: float
    bit_score: float
    subject_seq: str = ""       # full subject protein (DIAMOND full_sseq); "" when absent
    query_coverage: float = field(init=False, default=0.0)
    species: str = field(init=False, default="")

    def __post_init__(self):
        m = re.search(r'\[([^\[\]]+)\]\s*$', self.subject_title.strip())
        self.species = m.group(1) if m else ""

    def compute_coverage(self, query_length: int) -> None:
        """
        Single-HSP coverage from this hit's query coordinates.
        For multi-HSP merged coverage use merge_query_intervals().
        """
        if query_length > 0:
            aligned = abs(self.query_end - self.query_start)
            self.query_coverage = round((aligned / query_length) * 100, 2)
 




# ---------------------------------------------------------------------------
# Multi-HSP coverage helper
# ---------------------------------------------------------------------------
 
def merge_query_intervals(hits: list['BlastxResult']) -> int:
    """
    Merges all query-coordinate intervals across multiple HSPs and returns
    the total number of non-overlapping aligned bases.
 
    Example:
        hits covering query [10-50], [30-80], [200-300]
        merged → [10-80] + [200-300] = 70 + 100 = 170 bases
 
    Intervals are normalised (start < end) to handle hits on the minus strand
    where Diamond may report query_end < query_start.
    """
    if not hits:
        return 0
 
    intervals = sorted(
        (min(h.query_start, h.query_end), max(h.query_start, h.query_end))
        for h in hits
    )
 
    merged_bases = 0
    cur_start, cur_end = intervals[0]
 
    for start, end in intervals[1:]:
        if start <= cur_end:          # overlapping or adjacent — extend
            cur_end = max(cur_end, end)
        else:                         # gap — commit and start new interval
            merged_bases += cur_end - cur_start
            cur_start, cur_end = start, end
 
    merged_bases += cur_end - cur_start   # commit last interval
    return merged_bases






# ---------------------------------------------------------------------------
# BLASTn result  (local binary or NCBI qblast)
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class BlastnResult:
    """One BLASTn hit — compatible with both local and online search."""
    qseqid: str
    qlen: int
    slen: int       # may be 0 on DBs not indexed with -parse_seqids
    qcovhsp: int    # query coverage per HSP, integer percentage 0-100
    pident: float
    evalue: float
    bit_score: float
    stitle: str
    accession: str | None = None
    query_start: int | None = None
    query_end: int | None = None
    subject_seq: str | None = None   # full subject nt sequence (viral-titled hits only)


# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class Taxonomy:
    """NCBI taxonomy record resolved from a BLASTx species name."""
    tax_id: int
    scientific_name: str
    no_rank: str | None
    clade: str | None
    kingdom: str | None
    phylum: str | None
    class_: str | None          # 'class' is a Python keyword
    order: str | None
    family: str | None
    subfamily: str | None
    genus: str | None
    species: str | None
    genome: str | None


@dataclass(slots=True)
class BlastnTaxonomy:
    """
    Taxonomy of the best BLASTn hit's subject.

    status  — "resolved" or "unresolved".
    method  — "title" (virus name scanned out of the hit title against the
              bundled viral taxonomy), or None when unresolved.
    """
    status: str
    method: str | None
    accession: str | None
    subject: str
    taxid: int | None = None         # taxid of the matched viral taxon
    matched_name: str | None = None  # the name matched in the title
    synonym_of: str | None = None    # current name, when matched_name is an NCBI synonym
    taxonomy: Taxonomy | None = None


# ---------------------------------------------------------------------------
# Viral family / order / genus descriptive info
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ViralFamilyInfo:
    """Merged descriptive record from viral-family-info JSON files."""
    source: str          # "ViralZone" or "ICTV"
    type: str            # "Family", "Order", or "Genus"
    name: str            # taxon name, e.g. "Orthomyxoviridae"
    info_high: str       # full-text / high-token description
    info_low: str        # standardized / low-token description


# ---------------------------------------------------------------------------
# LLM scoring output
# ---------------------------------------------------------------------------

# Novelty tiers — how far a contig sits from known viruses (finer than the
# classification). Ordered from closest to farthest; "non-viral" mirrors the
# classification. Criteria live in score_heuristic (CALIBRATION).
NOVELTY_TIERS = (
    "known",             # same species, close strain   (>=95% nt)
    "variant",           # known species, divergent variant (>=90% aa or >=85% nt)
    "novel-species",     # putative new species          (70–90% aa)
    "divergent",         # putative new genus or above   (40–70% aa)
    "highly-divergent",  # remote homology / HMM only    (<40% aa)
    "non-viral",
)
# Qualifiers attached to a novelty tier.
NOVELTY_FLAGS = ("low-coverage", "nt-only", "hmm-only")


# LlmOutput.classification values that mean "no assessment was made": the
# request failed (api-error) or the model's reply could not be interpreted
# (parse-error). Their vq_score is a placeholder 0 and must not be read as a score.
LLM_ERROR_CLASSES = frozenset({"api-error", "parse-error"})


@dataclass(slots=True)
class LlmOutput:
    """Result produced by score_ai.SequenceScorer for one NucSequence."""
    seq_id:         str
    model:          str
    mode:           str        # "high" or "low"
    vq_score:       int        # 0–100
    classification: str        # "viral-known" | "viral-unknown" | "non-viral"
    analysis:       str        # plain-text summary (max ~200 words)
    blastn_species: str        = field(default="")
    error:          str | None = field(default=None)
    novelty:        str        = field(default="")   # one of NOVELTY_TIERS, "" if not given


# ---------------------------------------------------------------------------
# Heuristic (non-LLM) scoring output
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ScoreComponents:
    """
    Per-component breakdown of a HeuristicScore.

    Each component is on the same 0-100 scale as the final vq_score. A value of
    None means the component was absent (e.g. BLASTn not run) and was therefore
    excluded from the weighted average — distinct from a present-but-zero value.
    """
    blastn:  float | None       # None when no viral BLASTn evidence is present
    blastx:  float | None       # None when no BLASTx hit is present
    hmm:     float | None       # None when no FILTER-database domain is present
    weights: dict[str, float]   # weights actually applied after renormalisation
    penalty: float = 1.0        # multiplicative factor applied after the weighted
                                # sum (1.0 = none; < 1.0 = false-positive damping)


@dataclass(slots=True)
class HeuristicScore:
    """
    Deterministic, rule-based score produced by score_heuristic.HeuristicScorer.

    Mirrors the shape of LlmOutput so the report, stats, and viewer can treat
    the two scoring engines symmetrically — the numeric metric is named
    `vq_score` in both.
    """
    seq_id:         str
    vq_score:       int          # 0-100
    classification: str          # "viral-known" | "viral-unknown" | "non-viral"
    blastn_species: str = ""
    analysis:       str = ""     # plain-text, template-generated breakdown
    components:     ScoreComponents | None = None
    novelty:        str = ""     # one of NOVELTY_TIERS
    novelty_flags:  list[str] = field(default_factory=list)   # subset of NOVELTY_FLAGS


# ---------------------------------------------------------------------------
# Sequence clustering
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ClusterMember:
    """Alignment record for one sequence within a ViralCluster."""
    seq_id:            str
    length:            int
    is_representative: bool
    identity:          float   # % identity vs representative (100.0 for rep itself)
    query_coverage:    float   # % of this sequence covered by the alignment
    aln_start:         int     # 1-based start on representative (0 = no alignment)
    aln_end:           int     # 1-based end on representative   (0 = no alignment)
    query_start:       int     # 1-based start within this sequence (0 = no alignment)
    query_end:         int     # 1-based end within this sequence   (0 = no alignment)


@dataclass(slots=True)
class ViralCluster:
    """Group of NucSequences sharing the same best BLASTx species hit."""
    cluster_id:        str
    species:           str
    representative_id: str
    members:           list[ClusterMember]

    @property
    def size(self) -> int:
        return len(self.members)


# ---------------------------------------------------------------------------
# Salmon quantification
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class PfamHkEntry:
    """One assembled contig identified as a housekeeping gene via Pfam domain hit."""
    seq_id:       str     # original contig ID
    pfam_target:  str     # Pfam_TargetID  (plot label)
    pfam_acc:     str     # e.g. PF00022
    pfam_desc:    str     # Pfam_Description
    pfam_details: str     # Pfam_Details (HTML, stripped for tooltips)
    tpm:          float
    num_reads:    float
    score:        float   # HMM bit score
    e_value:      float


@dataclass(slots=True)
class SalmonEntry:
    """One row from quant.sf, tagged by sequence origin."""
    name:       str
    length:     int
    eff_length: float
    tpm:        float
    num_reads:  float
    seq_type:   str    # "viral" | "host" | "conserved"
    kingdom:    str    # e.g. "MAMMALS" — empty string for viral / host entries


@dataclass(slots=True)
class HostViralHit:
    """BLASTn alignment from one host transcript to one viralquest viral sequence."""
    host_transcript_id: str
    viral_seq_id:       str
    pident:             float
    qcovhsp:            int
    evalue:             float
    bit_score:          float


@dataclass(slots=True)
class HostViralRecord:
    """Host transcript with its Salmon quant data and all viral BLASTn hits."""
    transcript: SalmonEntry
    blast_hits: list[HostViralHit]


@dataclass(slots=True)
class CoverageProfile:
    """
    Per-base read-coverage profile of one viral sequence.

    Produced by aligning the input reads to the confirmed viral sequences
    (minimap2) and computing per-base depth (samtools depth).  Used by the
    genome viewer to draw a coverage track aligned with the ORF/HMM lanes;
    a coverage discontinuity or internal drop flags a possible mis-assembly
    (chimeric contig).
    """
    seq_id:          str
    length:          int
    mean_depth:      float
    max_depth:       float
    cv:              float              # stdev/mean — coverage uniformity (lower = more uniform)
    breadth_1x:      float              # fraction of bases with depth >= 1
    bins:            list[float]        # downsampled mean depth per bin (<= N_BINS points)
    low_cov_regions: list[list[int]]    # internal [start, end] runs below the drop threshold
    quality_bins:    list[float] = field(default_factory=list)  # mean base-quality (Phred) per bin, aligned to `bins`
    mean_quality:    float = 0.0        # mean per-base read quality (Phred) across the sequence
    # Strand-split depth (sense = forward-mapping reads, antisense = reverse),
    # derived from the case of the mpileup read-bases column; aligned to `bins`.
    bins_fwd:        list[float] = field(default_factory=list)  # forward (sense) depth per bin
    bins_rev:        list[float] = field(default_factory=list)  # reverse (antisense) depth per bin
    mean_depth_fwd:  float = 0.0        # mean forward-strand depth across the sequence
    mean_depth_rev:  float = 0.0        # mean reverse-strand depth across the sequence


@dataclass(slots=True)
class SelfRepeat:
    """
    One internal repeat of a sequence against itself (self-BLASTn HSP).

    The full-length self-diagonal hit is filtered out upstream, so every record
    here is a *duplicated* or *inverted* region — a direct-repeat pair when
    ``strand == "plus"`` and an inverted repeat when ``strand == "minus"``.
    """
    q_start: int      # 0-based, inclusive — first copy
    q_end:   int
    s_start: int      # 0-based, inclusive — second copy
    s_end:   int
    strand:  str      # "plus" (direct repeat) | "minus" (inverted repeat)
    pct_id:  float    # alignment percent identity
    length:  int      # alignment length (bp)


@dataclass(slots=True)
class DotPlot:
    """
    Self-similarity dot plot, precomputed for rendering as a native <svg>.

    Each segment is an extended k-mer match drawn as a line: matches on the main
    diagonal are the sequence itself; off-diagonal segments mark direct repeats,
    anti-diagonal segments mark inverted repeats.  Coordinates are in base pairs.
    """
    word:     int                    # k-mer word size used to seed matches
    length:   int                    # sequence length (both axes)
    segments: list[list[int]]        # each [x1, y1, x2, y2] (bp); anti-diagonal ⇒ inverted


@dataclass(slots=True)
class SequenceQuality:
    """
    Structural quality signals for one viral sequence, complementary to the
    read-coverage / base-quality track.  Flags low-complexity stretches,
    internal repeats and overall k-mer repetitiveness that may indicate
    mis-assembly or sequencing artefacts.
    """
    seq_id:                 str
    length:                 int
    low_complexity_regions: list[list[int]]  # dustmasker [start, end] runs (0-based, inclusive)
    self_repeats:           list[SelfRepeat] # self-BLASTn internal repeats
    kmer_size:              int              # jellyfish k
    kmer_distinct:          int              # number of distinct k-mers
    kmer_total:             int              # total k-mers counted
    kmer_repeat_score:      float            # fraction of k-mers with multiplicity > 1 (0–1)
    low_complexity_frac:    float            # fraction of bases inside low-complexity runs (0–1)
    dotplot:                "DotPlot | None" = None


@dataclass(slots=True)
class SalmonQuantReport:
    """Full Salmon quantification result, attached as an optional pipeline section."""
    reads:           list[str]             # read file(s) passed to salmon quant
    mapping_rate:    float                 # overall mapping rate (%)
    total_reads:     int                   # total reads counted across all entries
    pathway:         str                   # "reference" | "de_novo"
    viral_quant:     list[SalmonEntry]     # TPM for each viralquest viral sequence
    conserved_quant: list[SalmonEntry]     # bundled HK (reference) or HK-matched contigs (de_novo)
    ref_hk_quant:    list[SalmonEntry]     # user-provided reference HK genes (reference pathway only)
    host_viral_hits: list[HostViralRecord] # host transcripts with similarity to viral seqs
    pfam_hk_quant:   list[PfamHkEntry]    # contigs identified as HK genes via Pfam domains (de_novo)


# ---------------------------------------------------------------------------
# Nucleotide sequence  (central object)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Synteny  (input ORFs ↔ CDS of the best viral BLASTn reference)
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class RefFeature:
    """A CDS or mat_peptide of a reference genome (GenBank), genome coordinates."""
    kind: str                 # "CDS" | "mat_peptide"
    product: str
    protein_id: str | None
    start: int                # 0-based, inclusive (envelope of all parts)
    end: int                  # exclusive
    strand: str               # "+" | "-"
    length_aa: int
    parent: int | None = None # mat_peptide → index of the CDS it lies in


@dataclass(slots=True)
class RefGenome:
    """Reference virus genome: the subject of a sequence's best viral BLASTn hit."""
    accession: str            # with version, as in the GenBank record
    title: str
    organism: str
    length: int
    cds: list[RefFeature] = field(default_factory=list)
    peptides: list[RefFeature] = field(default_factory=list)


@dataclass(slots=True)
class SyntenyHit:
    """One DIAMOND blastp HSP: an input ORF against a reference CDS (every HSP
    is kept, so one ORF/CDS pair may have several, each with its own range)."""
    orf_name: str
    cds_index: int            # index into RefGenome.cds
    pident: float
    evalue: float
    bit_score: float
    q_start: int              # 1-based aa on the ORF
    q_end: int
    s_start: int              # 1-based aa on the CDS
    s_end: int
    q_cov: float              # % of the ORF aligned
    s_cov: float              # % of the CDS aligned
    orf_nt: list[int]         # [start, end) of the aligned stretch on the contig
    ref_nt: list[int]         # [start, end) envelope of the aligned stretch on the reference
    # One [ref_start, ref_end, orf_start, orf_end] per CDS part the alignment
    # spans — several for a spliced / joined CDS, so each exon gets its ribbon.
    segments: list[list[int]] = field(default_factory=list)


@dataclass(slots=True)
class SyntenyLink:
    """One DIAMOND blastp HSP between canonical ORFs of two sequences that share
    a reference (contig ↔ contig correspondence). Every HSP is kept; each pair
    is stored once, on the sequence whose id sorts first."""
    orf_name: str
    other_seq: str
    other_orf: str
    pident: float
    evalue: float
    bit_score: float
    q_start: int              # 1-based aa on this ORF
    q_end: int
    s_start: int              # 1-based aa on the other ORF
    s_end: int
    orf_nt: list[int]         # [start, end) of the aligned stretch on this contig
    other_nt: list[int]       # [start, end) of the aligned stretch on the other contig


@dataclass(slots=True)
class SyntenyResult:
    """A sequence's synteny with its reference (see synteny.py)."""
    reference: str            # RefGenome.accession
    flipped: bool             # contig is reverse-complemented relative to the reference
    offset: float             # reference position ≈ offset ± contig position (see synteny.py)
    hits: list[SyntenyHit] = field(default_factory=list)
    links: list[SyntenyLink] = field(default_factory=list)


@dataclass(slots=True)
class NucSequence:
    """A parsed nucleotide sequence with all downstream annotation results."""
    id: str
    sequence: str
    sample_name: str = ""           # set by the pipeline from the input filename

    uid: uuid.UUID = field(default_factory=uuid.uuid4, init=False)

    # annotation results
    orfs:          list[Orf]         = field(default_factory=list, init=False)
    blastx_hits:   list[BlastxResult] = field(default_factory=list, init=False)  # phase-1 refseq hits
    blastx_nr_hits: list[BlastxResult] = field(default_factory=list, init=False) # phase-2 NR hits
    blastn_hits:   list[BlastnResult] = field(default_factory=list, init=False)

    # viral selection flag — set to True after phase-1 filter
    is_viral: bool = field(default=False, init=False)

    # taxonomy resolved from best BLASTx species
    taxonomy: Taxonomy | None = field(default=None, init=False)

    # taxonomy of the best BLASTn hit's subject (accession → taxid, else name)
    blastn_taxonomy: BlastnTaxonomy | None = field(default=None, init=False)

    # viral family / order / genus annotation
    viral_family_info: ViralFamilyInfo | None = field(default=None, init=False)

    # LLM scoring result
    llm_output: LlmOutput | None = field(default=None, init=False)

    # heuristic (non-LLM) scoring result
    heuristic_output: HeuristicScore | None = field(default=None, init=False)

    # cluster membership
    cluster_id: str | None = field(default=None, init=False)

    # Salmon quantification (populated before LLM scoring when --reads is given)
    salmon_tpm:   float | None = field(default=None, init=False)
    salmon_reads: float | None = field(default=None, init=False)

    # Read-coverage profile (populated when --reads is given)
    coverage: "CoverageProfile | None" = field(default=None, init=False)

    # Structural quality signals — dustmask / self-BLAST / jellyfish / dot plot
    # (populated when --reads is given)
    seq_quality: "SequenceQuality | None" = field(default=None, init=False)

    # Synteny with the best viral BLASTn reference (populated unless --no-synteny)
    synteny: SyntenyResult | None = field(default=None, init=False)

    # computed
    length:    int   = field(init=False)
    n_count:   int   = field(init=False)
    gc_content: float = field(init=False)

    def __post_init__(self):
        self.length = len(self.sequence)
        self.n_count = self.sequence.count('N')
        if self.length > 0:
            g = self.sequence.count('G')
            c = self.sequence.count('C')
            self.gc_content = round(((g + c) / self.length) * 100, 2)
        else:
            self.gc_content = 0.0

    @property
    def best_blastx(self) -> BlastxResult | None:
        """Top NR BLASTx hit (highest bit_score), or refseq hit as fallback."""
        pool = self.blastx_nr_hits or self.blastx_hits
        return max(pool, key=lambda h: h.bit_score) if pool else None

    @property
    def best_blastn(self) -> BlastnResult | None:
        return self.blastn_hits[0] if self.blastn_hits else None

    @property
    def hmm_domains(self) -> list[HmmDomain]:
        """Flat list of all domains across all ORFs."""
        return [d for orf in self.orfs for d in orf.domains]

    @property
    def hmm_raw_counts(self) -> dict[str, int]:
        """
        Raw HMM hit counts (score >= threshold, before positional de-dup)
        aggregated across all ORFs, keyed by database name. Quantitative
        complement to hmm_domains for the heuristic scorer.
        """
        out: dict[str, int] = {}
        for orf in self.orfs:
            for db, n in orf.raw_hmm_counts.items():
                out[db] = out.get(db, 0) + n
        return out