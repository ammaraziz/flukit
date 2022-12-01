#!/usr/bin/env python3
'''
Variant calling for influenza gene segement. 
'''
from .codon_align import get_cds, codon_align, safe_translate
from flukit.utils.utils import get_reference

def set_gene(seqdict):
    '''
    Modify dict with SeqRecords adding gene info
    '''
    seqTogene = {'1': 'PB2', '2': 'PB1', '3': 'PA', '4': 'HA',
                 '5': 'NP', '6': 'NA', '7': 'MP', '8': 'NS'}

    for record in seqdict:
        try:
            seqdict[record].gene = seqTogene[seqdict[record].id.split('.')[1]]
        except Exception:
            raise ValueError(
                f"{seqdict[record].gene} has no segemnt number. Check input seqs")
    return(seqdict)

def get_gene(SeqRecord):
    seqTogene = {'1': 'PB2', '2': 'PB1', '3': 'PA', '4': 'HA',
                 '5': 'NP', '6': 'NA', '7': 'MP', '8': 'NS'}
    return(seqTogene[SeqRecord.id.split('.')[1]])

def get_id(SeqRecord):
    return(SeqRecord.id.split('.')[0])

def align(lineage, input_record):
    '''
    align gene to reference
    returns seq_aa and refAA
    '''
    gene = input_record.gene

    refname, ref = get_reference(lineage, gene)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(
        ref=ref, refname=refname, input_gene=gene)
    try:
        seq_aln = codon_align(input_record, refCDS, refAA, 0, cds_end)
        if seq_aln is None:
            raise ValueError("didn't translate properly")

    except Exception:
        raise ValueError("Sequence didn't align, check lineage input")

    seq_aa = safe_translate(seq_aln)

    return(seq_aa, refAA)

def all_variants(sample, ref):
    '''
    Get SNPs from on aligned sequence. Takes one input at a time.

    Parameters : SeqRecord
        sample  SeqRecord for analysis
        ref     SeqRecord sequence used to call SNPS
    Returns : List
        ['T30C', 'G45A', 'A51G', ...] 
    '''

    if len(sample) != len(ref):
        raise ValueError(
            "Sequences are of different lengths. Make sure sequences are aligned")

    snps = []
    for index in range(0, len(ref)):
        nuc = set([ref[index], sample[index]])
        if len(nuc) != 1:
            variant = ref[index] + str(index+1) + sample[index]
            snps.append(variant)

    return(snps)

def specific_variants(sample, gene, lineage):
    '''
    Get aa from aligned sequence.
        NA      H275Y
        PA      S31N
        MP      I38X, 22, 33

    Parameters
        sample : SeqRecord
            aligned amino acid sequence (string)
        gene : str
            name of gene - NA, PA, MP
        lineage: str
            list of h1n1, h3n2, vic, yam

    Return : list containing strings
        'H' or 'Y', 'S' or 'N', 'I' or 'X'
    '''

    if lineage == 'h3n2':
        #gene_pos = {'NA': [273], 'PA': [22, 33, 35, 37, 118, 198], 'MP': [30]}
        gene_pos = {'NA': [273], 'PA': [37], 'MP': [30]}
    else:
        #gene_pos = {'NA': [274], 'PA': [22, 33, 35, 37, 118, 198], 'MP': [30]}
        gene_pos = {'NA': [274], 'PA': [37], 'MP': [30]}
    
    # check input of gene
    if gene.upper() not in gene_pos.keys():
        raise ValueError("Unrecognised Gene:" + gene.upper() + "only NA, PA, MP are allowed")

    nuc = [sample[val] for val in gene_pos[gene]]
    # combine position and aa togther
    #nuc_out = [str(m+1)+n for m, n in zip(gene_pos[gene], nuc)]
    return(nuc)


def get_variants(input_record, lineage):
    '''
    Get variants for any gene segment

    Parameters
        input_record  :  Seq.Record of reference 
        lineage  :  str of lineage
    Return
        list  :  amino acid varants
    '''

    gene = input_record.gene
    if gene in ['PB2', 'PB1', 'NS', 'NP']:
        return()
    try:
        seq_aa, refAA = align(lineage = lineage, input_record = input_record)
    except Exception as e:
        print(f'can not align {e}')
        return()

    if gene == 'HA':
        snps_all = all_variants(seq_aa, refAA)
        return(";".join(snps_all))
    if gene == 'PA':
        pa = specific_variants(seq_aa, input_record.gene, lineage)
        return(";".join(pa))
    if gene == 'NA':
        na = specific_variants(seq_aa, input_record.gene, lineage)
        return(";".join(na))
    if gene == 'MP':
        mp = specific_variants(seq_aa, input_record.gene, lineage)
        return(";".join(mp))

def get_vacc_ref(lineage):
    vacc = {"h1n1" : "A/California/07/2009" , 
    "h3n2" : "A/Beijing/32/1992",
    "bvic" : "B/Hong Kong/02/1993"}
    return(vacc[lineage])