#!/usr/bin/env python3
'''
Variant calling for influenza gene segement. 
'''
from Bio import SeqRecord
from typing import List, Tuple
from .codon_align import get_cds, codon_align, safe_translate
from flukit.utils.utils import get_reference, read_in_mutations

seqTogene = {
    '1': 'PB2', 
    '2': 'PB1', 
    '3': 'PA', 
    '4': 'HA',
    '5': 'NP', 
    '6': 'NA', 
    '7': 'MP', 
    '8': 'NS'}

def set_gene(seqdict: dict) -> dict:
    '''
    Modify dict with SeqRecords adding gene info
    '''
    for record in seqdict:
        try:
            seqdict[record].gene = seqTogene[seqdict[record].id.split('.')[1]]
        except Exception:
            raise ValueError(
                f"{seqdict[record].gene} has no segemnt number. Check input seqs")
    return(seqdict)

def get_gene(SeqRecord: SeqRecord) -> str:
    '''
    converts gene number (.4) to gene abbr (HA) 
    '''
    return(seqTogene[SeqRecord.id.split('.')[1]])

def get_id(SeqRecord: SeqRecord) -> str:
    '''
    get sequence id eg N100500
    '''
    return(SeqRecord.id.split('.')[0])

def align(lineage: str, input_record: SeqRecord) -> Tuple[str, str]:
    '''
    align gene to reference
    returns seq_aa and refAA
    '''

    refname, ref = get_reference(lineage, input_record.gene)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(
        ref=ref, refname=refname, input_gene=input_record.gene)
    try:
        seq_aln = codon_align(input_record, refCDS, refAA, 0, cds_end)
        if seq_aln is None:
            raise ValueError(f"didn't translate properly - {input_record}")

    except Exception:
        raise ValueError("Sequence didn't align, check lineage input")

    seq_aa = safe_translate(seq_aln)

    return(seq_aa, refAA)

def get_ha_snps(sample: SeqRecord, ref: SeqRecord) -> Tuple[str, ...]:
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
            if sample[index] != '-':
                variant = ref[index] + str(index+1) + sample[index]
                snps.append(variant)

    return(tuple(snps))

def get_pa_snps(
    sample: SeqRecord, 
    ref: SeqRecord, 
    gene: str, 
    lineage: str) -> Tuple[str, ...]:
    '''
    get PA mutations of interest
    '''

    gene_pos = read_in_mutations(lineage)
    pa_pos = gene_pos[gene]

    pa_snps = []
    for pos in pa_pos:
        variant = ref[pos] + str(pos+1) + sample[pos]
        pa_snps.append(variant)

    return(tuple(pa_snps))

def specific_variants(sample: SeqRecord, gene: str, lineage: str) -> List[str]:
    '''
    Get aa from aligned sequence.
        NA      H275Y
        MP      S31N
        PA      I38X, 22, 33

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

    gene_pos = read_in_mutations(lineage)
    if gene.upper() not in gene_pos.keys():
        raise ValueError("Unrecognised Gene: " + gene.upper() + ". Only NA, PA, MP are allowed")
    mutations = [sample[val] for val in gene_pos[gene]]

    return(mutations)


def get_variants(input_record: SeqRecord, lineage: str) -> List:
    '''
    Get variants for any gene segment

    Parameters
        input_record  :  SeqRecord of reference 
        lineage       :  str of lineage
    Return
        list          :  amino acid varants
    '''

    gene = input_record.gene
    if gene in ['PB2', 'PB1', 'NS', 'NP']:
        return()
    try:
        seqAA, refAA = align(lineage = lineage, input_record = input_record)
    except Exception as e:
        print(f'can not align {e}')
        return()

    if gene == 'HA':
        snps_all = get_ha_snps(seqAA, refAA)
        return(";".join(snps_all))
    if gene == 'PA':
        pa = get_pa_snps(seqAA, refAA, input_record.gene, lineage)
        return(";".join(pa))
    if gene == 'NA':
        na = specific_variants(seqAA, input_record.gene, lineage)
        return(";".join(na))
    if gene == 'MP':
        mp = specific_variants(seqAA, input_record.gene, lineage)
        return(";".join(mp))

def get_vacc_ref(lineage: str) -> str:
    vacc = {
        "h1n1" : "A/California/07/2009" , 
        "h3n2" : "A/Beijing/32/1992",
        "vic" : "B/Hong Kong/02/1993"
        }
    return(vacc[lineage])
