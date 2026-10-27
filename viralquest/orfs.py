import orfipy_core
from Bio.Seq import Seq
from loguru import logger
from .biodata import NucSequence, Orf

class OrfAnalyzer:
    """Find, filter and manipulate ORFs in valid sequences."""
    

    def __init__(self, min_len_nt: int = 150):
        self.min_len_nt = min_len_nt
        logger.info(f"Starting ORF search (Minimal length: {self.min_len_nt} nt)")



    def _parse_orfipy_description(self, description: str) -> dict:

        desc_dict = {}
        for item in description.split(';'):
            if '=' in item:
                k, v = item.split('=')
                desc_dict[k] = v
            elif ':' in item:
                k, v = item.split(':')
                desc_dict[k] = v
        return desc_dict



    def find_n_save_orfs(self, nuc_seq: NucSequence) -> None:
        """uses orfipy to find and index the ORFs in NucSequence"""
        
        logger.debug(f"Searching ORFs in sequence {nuc_seq.id}...")
        
        seq_str = nuc_seq.sequence
        
        # partial5 / partial3: ORFs that run off either end of the contig (no
        # start or no stop codon inside it).  Assembled viral contigs are often
        # fragments of a gene or of a polyprotein, so complete-only ORF calling
        # would miss them and hide their domains from the HMM searches.
        # The ORF_type field ("complete", "5-prime-partial", "3-prime-partial",
        # "no_start_no_stop") records which kind each ORF is.
        results = orfipy_core.orfs(seq_str, minlen=self.min_len_nt,
                                   partial5=True, partial3=True)
        
        orfs_number = 0
        
        for start, stop, strand, description in results:

            metadata = self._parse_orfipy_description(description)
            
            bio_seq = Seq(seq_str[start:stop])
            
            # reverse complement to negative strands 
            if strand == '-':
                bio_seq = bio_seq.reverse_complement()
                
            nuc_sequence = str(bio_seq)
            
            # translate to aminoacids  
            try:
                aa_sequence = str(bio_seq.translate())
            except Exception as e:
                logger.warning(f"Translate ORF from {start}-{stop}: {e}")
                aa_sequence = ""
            
            # Orf dataclass
            nova_orf = Orf(
                start_codon=metadata.get('Start', ''),
                stop_codon=metadata.get('Stop', ''),
                orf_type=metadata.get('ORF_type', ''),
                start_position=start,
                stop_position=stop,
                strand=strand,
                frame=int(metadata.get('ORF_frame', 0)),
                length_nt=len(nuc_sequence),
                length_aa=len(aa_sequence),
                bigger_than_50=len(aa_sequence) > 50,
                aa_sequence=aa_sequence,
                nuc_sequence=nuc_sequence,
                name=f"{nuc_seq.id}_ORF_{start}_{stop}{strand}"
            )
            

            nuc_seq.orfs.append(nova_orf)
            orfs_number += 1
            
        logger.info(f"Success: {orfs_number} ORFs found in {nuc_seq.id}.")

    # some future methods
    #
    #
    #
    
    def get_longest_orf(self, nuc_seq: NucSequence) -> Orf | None:
        """returns biggest ORF"""
        if not nuc_seq.orfs:
            return None
        return max(nuc_seq.orfs, key=lambda orf: orf.length_aa)
        
    def filter_complete_orfs(self, nuc_seq: NucSequence) -> list[Orf]:
        """returns ORFs with valid start and stop codons"""
        return [orf for orf in nuc_seq.orfs if orf.start_codon and orf.stop_codon]


# ── Older reports (hidden --reload) ─────────────────────────────────────────

# Partial ORFs (no start / no stop codon inside the contig) are called since
# 3.0.9; earlier reports carry complete ORFs only, which misses most gene
# fragments of assembled viral contigs.
PARTIAL_ORFS_SINCE = (3, 0, 9)


def _version(v: str | None) -> tuple[int, ...]:
    try:
        return tuple(int(x) for x in str(v).split(".")[:3])
    except ValueError:
        return (0,)


def needs_orf_recall(report: dict) -> bool:
    """True for a report whose ORFs predate partial-ORF calling (and were not re-called yet)."""
    meta = report.get("meta") or {}
    return (_version(meta.get("viralquest_version")) < PARTIAL_ORFS_SINCE
            and not meta.get("orfs_recalled_with"))


def recall_orfs(seq: dict, analyzer: "OrfAnalyzer") -> int:
    """
    Re-call the ORFs of an exported sequence (dict, in place) with the current
    caller. A new ORF sharing an old one's strand and 3' end is the same gene
    extended upstream: it inherits the old ORF's HMM domains, shifted by the
    extension (domain coordinates are aa from the ORF's N-terminus). Old ORFs
    with no such match are kept. Returns how many ORFs are new or longer.
    """
    from .exporter import ReportExporter

    tmp = NucSequence(id=seq["id"], sequence=seq.get("sequence") or "")
    analyzer.find_n_save_orfs(tmp)
    three = lambda o: (o["strand"], o["stop_position"] if o["strand"] == "+" else o["start_position"])
    old = {three(o): o for o in seq.get("orfs") or []}
    out, changed = [], 0
    for orf in tmp.orfs:
        d = ReportExporter._orf_to_dict(orf)
        prev = old.pop(three(d), None)
        if prev is None:
            changed += 1
        else:
            shift = ((prev["start_position"] - d["start_position"]) if d["strand"] == "+"
                     else (d["stop_position"] - prev["stop_position"])) // 3
            changed += shift > 0
            d["domains"] = [{**dom, "start": dom["start"] + shift, "stop": dom["stop"] + shift}
                            for dom in prev.get("domains") or []]
        out.append(d)
    out += old.values()
    seq["orfs"] = sorted(out, key=lambda o: (o["start_position"], o["stop_position"]))
    return changed
