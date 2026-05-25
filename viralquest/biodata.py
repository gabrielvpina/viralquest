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
    qseqid: str
    qlen: int
    slen: int
    qcovhsp: float
    pident: float
    evalue: float
    stitle: str

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