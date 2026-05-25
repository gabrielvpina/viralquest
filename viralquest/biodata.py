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
    slen: int     # slen — may be 0 on dbs not indexed with -parse_seqids
    qcovhsp: int     # qcovs — blastn reports integer percentage (0-100)
    pident: float
    evalue: float
    stitle: str


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


# ---------------------------------------------------------------------------
# Nucleotide sequence  (central object)
# ---------------------------------------------------------------------------

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