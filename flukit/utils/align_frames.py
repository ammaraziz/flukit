from numpy import inf
from Bio import SeqRecord
from typing import Tuple
from .utils import safe_translate, load_features, get_reference
from .codon_align import get_cds, codon_align, safe_translate


# alignment code from nextstrain influenza build

scoring_params = {"score_match": 3, "score_mismatch": -
                  1, "score_gapext": -1, "score_gapopen": -10}

def align_pairwise(seq1, seq2):
    from Bio import pairwise2
    aln = pairwise2.align.globalms(seq1, seq2,
                                   scoring_params['score_match'], scoring_params['score_mismatch'],
                                   scoring_params['score_gapopen'], scoring_params['score_gapext'],
                                   penalize_end_gaps=False, one_alignment_only=True)[0]
    return aln[2], aln[0], aln[1]

def premature_stop(seq, refstr, refAA):
    '''
    Check to see if aa sequence has premature stop compared to aa ref sequence
    
    returns tuple(boolan, position of prem-stop)
    '''
    seqstr = str(seq.seq).upper()
    score, refaln, seqaln = align_pairwise(refstr, seqstr)
    
    if score < 0:  # did not align
        return True, 0

    seqAA = safe_translate(seqaln)
    scoreAA, refalnAA, seqalnAA = align_pairwise(refAA, seqAA)
    return(refalnAA, seqalnAA)

def align(lineage: str, input_record: SeqRecord.SeqRecord) -> Tuple[str, str]:
    '''
    align gene to reference
    returns seqAA and refAA
    '''

    refname, ref = get_reference(lineage, input_record.gene)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(
        ref=ref, refname=refname, input_gene=input_record.gene)

    try:
        seq_aln = codon_align(input_record, refCDS, refAA, 0, cds_end)
        if seq_aln is None:
            raise ValueError(f"didn't translate properly - {input_record}")
    except Exception:
        raise ValueError("Sequence didn't align.")
        

    seqAA = safe_translate(seq_aln)

    return(seqAA, refAA)