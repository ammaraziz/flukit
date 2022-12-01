import argparse
from numpy import inf
from utils.utils import safe_translate, load_features
from Bio import SeqIO, Seq, AlignIO
from Bio.Align import MultipleSeqAlignment

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




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--reference", required=True, help="annotated genbank file")
    parser.add_argument("--output", required=True, help="FASTA file of extracted sample sequences")
    parser.add_argument("--gene", required=True, help="Gene")
    args = parser.parse_args()

    aln = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(args.reference, 'genbank')

    # get sequence as string, CDS seq, amino acid sequence, and start/end pos
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(ref, refname = args.reference, input_gene=args.gene)

    alignment = []
    for seq in aln:
        seq_aln = codon_align(seq,  refstr, refAA, cds_start, cds_end)
        if seq_aln:
            # if len(seq_aln)!=len(refstr):
            #     print(seq.name, seq_aln, refstr)
            # else:
                seq.seq=Seq.Seq(seq_aln)
                alignment.append(seq)

    # output
    AlignIO.write(MultipleSeqAlignment(alignment), args.output, 'fasta')
