import sys
from numpy import array, inf
from flukit.utils.utils import safe_translate, load_features

# from nextstrain influenza build

scoring_params = {
    "score_match": 3, 
    "score_mismatch": -1, 
    "score_gapext": -1, 
    "score_gapopen": -10
    }


def align_pairwise(seq1: str, seq2: str) -> tuple:
    '''
    Align sequences
    returns: score, refaln, seqaln
    '''
    
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = scoring_params['score_match']
    aligner.mismatch_score = scoring_params['score_mismatch']
    aligner.open_gap_score = scoring_params['score_gapopen']
    aligner.extend_gap_score = scoring_params['score_gapext']
    aln = aligner.align(seq1, seq2)[0]

    return(aln.score, aln[0], aln[1])

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
    feature_dict = {
        'mp': 'M2', 
        'ns': 'NS1', 
        'pb2': 'PB2', 
        'pb1': 'PB1',
        'pa': 'PA', 
        'np': 'NP', 
        'na': 'NA', 
        'ha': 'HA'
        }

    cds_start, cds_end = inf, 0
    input_gene = input_gene.lower()
    gene_feat = feature_dict[input_gene]

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
        feat = load_features(refname, gene_feat)
        refstr = str(feat[gene_feat].extract(ref).seq).upper()
        refCDS = refstr
        refAA = safe_translate(refstr)
        cds_start, cds_end = 0, len(refCDS)

    else:
        raise ValueError("Error gene input or reference is wrong.")

    return refstr, refCDS, refAA, cds_start, cds_end


def codon_align(seq, refstr, refAA, cds_start, cds_end):
    seqstr = str(seq.seq).upper()
    score, refaln, seqaln = align_pairwise(refstr, seqstr)
    if score < 0:  # did not align
        return None
    ref_aln_array = array(list(refaln))
    seq_aln_array = array(list(seqaln))

    # stip gaps
    ungapped = ref_aln_array != '-'
    ref_aln_array_ungapped = ref_aln_array[ungapped]
    seq_aln_array_ungapped = seq_aln_array[ungapped]

    seq5pUTR = "".join(seq_aln_array_ungapped[:cds_start])
    seq3pUTR = "".join(seq_aln_array_ungapped[cds_end:])
    seqCDS = "".join(seq_aln_array_ungapped[cds_start:cds_end])
    seqCDS_ungapped = seqCDS.replace('-', '')
    seqAA = safe_translate(seqCDS_ungapped)

    scoreAA, refalnAA, seqalnAA = align_pairwise(refAA, seqAA)
    if scoreAA < 0 or sum(seqAA.count(x) for x in ['*', 'X']) > 10 or refalnAA.count('-') > 10:
        print(seq.id, "didn't translate properly", file=sys.stderr)
        return('')

    seqCDS_aln = seq5pUTR
    pos = 0
    for aa_ref, aa_seq in zip(refalnAA, seqalnAA):
        if aa_seq == '-':
            seqCDS_aln += '---'
            # if the nucleotide sequence is gapped
            # (i.e. because of missing data at the 5p and 3p end, advance pos)
            if seqCDS_ungapped[pos:pos+3] == '---':
                pos += 3
        else:
            if len(seqCDS_ungapped) >= pos+3:
                seqCDS_aln += seqCDS_ungapped[pos:pos+3]
            else:
                seqCDS_aln += '---'
            pos += 3

    return ''.join(seqCDS_aln)+seq3pUTR