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


def get_cds(ref, refname=None, input_gene=None):
    '''
    assuming there is one contiguous coding region which might be
    split into multiple sub-proteins like HA1 and HA2.
    loop over all features, pull out min and max of their union

    Parameters
        ref     :   Seq.Record of reference 
        gene    :   str
            'pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns'
    Return
        refstr, refCDS, refAA, cds_start, cds_end 
    '''
    # setup and vars
    feature_dict = {'mp': 'M2', 'ns': 'NS1', 'pb2': 'PB2', 'pb1': 'PB1',
                    'pa': 'PA', 'np': 'NP', 'na': 'NA', 'ha': 'HA'}

    cds_start, cds_end = inf, 0
    input_gene = input_gene.lower()
    feature_name = feature_dict[input_gene]

    if input_gene is None or input_gene == 'ha':
        for feature in ref.features:
            if feature.type == 'CDS':
                # skip over Sigpep so we can get the correct positions
                if 'SigPep' not in feature.qualifiers['gene']:
                    if feature.location.start < cds_start:
                        cds_start = feature.location.start
                    if feature.location.end > cds_end:
                        cds_end = feature.location.end

        refstr = str(ref.seq).upper()
        refCDS = refstr[cds_start:cds_end]
        refAA = safe_translate(refstr[cds_start:cds_end])

    elif input_gene in feature_dict.keys():
        feat = load_features(refname, feature_name)
        refstr = str(feat[feature_name].extract(ref).seq).upper()
        refCDS = refstr  # why?
        refAA = safe_translate(refstr)
        cds_start, cds_end = 0, len(refCDS)

    else:
        raise ValueError("Error gene input or reference is wrong.")

    return refstr, refCDS, refAA, cds_start, cds_end

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

def align(lineage: str, input_record: SeqRecord) -> Tuple[str, str]:
    '''
    align gene to reference
    returns seqAA and refAA
    '''

    refname, ref = get_reference(lineage, input_record.gene)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(
        ref=ref, refname=refname, input_gene=input_record.gene)

    try:
        seq_aln = codon_align(input_record, refCDS, refAA, 0, cds_end)
    except Exception:
        raise ValueError("Sequence didn't align.")
    except seq_aln is None:
        raise ValueError(f"didn't translate properly - {input_record}")

    seqAA = safe_translate(seq_aln)

    return(seqAA, refAA)