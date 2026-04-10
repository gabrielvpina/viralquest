from dataclasses import dataclass, field

@dataclass
class NucSequence:
    """Valid biological sequences."""
    id: str
    sequence: str
    
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
