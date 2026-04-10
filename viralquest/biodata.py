import uuid
from dataclasses import dataclass, field

@dataclass(slots=True)
class Orf:
    """Stores a single Open Reading Frame metadata"""
    start_codon: str
    stop_codon: str
    start_position: int
    stop_position: int
    length_aa: int
    length_nt: int
    bigger_than_50: bool

@dataclass(slots=True)
class BlastnResult:
    """Stores metadata from a single BLASTn hit"""
    hit_id: str
    hit_order: int
    e_value: float
    query_coverage: float
    identity_perc: float
    hsp: float


@dataclass(slots=True)
class NucSequence:
    """Valid biological sequences. Mutable fields."""
    id: str
    sequence: str

    # unique code id 
    uid: uuid.UUID = field(default_factory=uuid.uuid4, init=False)
    
    # blastn metadata
    blast_hits: list[BlastnResult] = field(default_factory=list, init=False)

    # orfs metadata
    orfs: list[Orf] = field(default_factory=list, init=False)

    # other parameters
    length: int = field(init=False)
    n_count: int = field(init=False)
    gc_content: float = field(init=False)

    def __post_init__(self):

        # length
        self.length = len(self.sequence)
        
        # count N bases 
        self.n_count = self.sequence.count('N')
       
        # gc_content
        if self.length > 0:
            g_count = self.sequence.count('G')
            c_count = self.sequence.count('C')
            self.gc_content = ((g_count + c_count) / self.length) * 100
        else:
            self.gc_content = 0.0

    def to_fasta_format(self) -> str:
        """returns sequence in fasta"""
        return f">{self.id}\n{self.sequence}"
