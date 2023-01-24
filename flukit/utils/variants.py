'''
Variant calling for influenza gene segement. 
'''
import typer
from Bio import SeqRecord
from typing import List
from rich import print
from .align_frames import align
from .utils import read_in_mutations

segements = {
    '1': 'PB2', 
    '2': 'PB1', 
    '3': 'PA', 
    '4': 'HA',
    '5': 'NP', 
    '6': 'NA', 
    '7': 'MP', 
    '8': 'NS'}

def get_ref(lineage: str) -> str:
    vacc = {
        "h1n1" : "A/California/07/2009" , 
        "h3n2" : "A/Beijing/32/1992",
        "vic"  : "B/Hong Kong/02/1993"
        }
    return(vacc[lineage])

def set_gene(seqdict: dict) -> dict:
    '''
    Modify dict with SeqRecords adding gene info
    '''
    for record in seqdict:
        try:
            seqdict[record].gene = segements[seqdict[record].id.split('.')[1]]
        except Exception:
            raise typer.BadParameter(f"Fasta header not formatted properly: {seqdict[record].id}. They must end in a number eg: N1000.4")
    return(seqdict)

def get_gene(SeqRecord: SeqRecord) -> str:
    '''
    converts gene number (.4) to gene abbr (HA) 
    '''
    return(segements[SeqRecord.id.split('.')[1]])

def get_id(SeqRecord: SeqRecord) -> str:
    '''
    get sequence id eg N100500
    '''
    return(SeqRecord.id.split('.')[0])

def get_ha_snps(
    sample: SeqRecord, 
    ref: SeqRecord) -> str:
    '''
    Get SNPs from on aligned sequence. Takes one input at a time.

    Parameters : SeqRecord
        sample  SeqRecord for analysis
        ref     SeqRecord sequence used to call SNPS
    Returns : List
        ['T30C', 'G45A', 'A51G', ...] 
    '''

    if len(sample) != len(ref):
        print(f"[bold yellow]seqAA, refAA of different lengths. Do you have the correct lineage? [/bold yellow]")
        typer.Exit()

    snps = []
    for index in range(0, len(ref)):
        nuc = set([ref[index], sample[index]])
        if len(nuc) != 1:
            if sample[index] != '-':
                variant = ref[index] + str(index+1) + sample[index]
                snps.append(variant)

    return(":".join(snps))

def get_pa_snps(
    sample: SeqRecord, 
    ref: SeqRecord, 
    lineage: str) -> str:
    '''
    get PA mutations of interest
    '''

    gene_pos = read_in_mutations(lineage)
    pa_pos = gene_pos['PA']

    pa_snps = []
    for pos in pa_pos:
        variant = ref[pos] + str(pos+1) + sample[pos]
        pa_snps.append(variant)

    return(";".join(pa_snps))

def get_snps(
    sample: SeqRecord.SeqRecord, 
    gene: str, 
    lineage: str) -> List[str]:
    '''
    Get aa from aligned sequence for
        NA      H275Y
        MP      S31N

    Parameters
        sample : SeqRecord
            aligned amino acid sequence (string)
        gene : str
            name of gene - NA, PA, MP
        lineage: str
            list of h1n1, h3n2, vic, yam

    Return - list[str]
        'H' or 'Y', 'S' or 'N', 'I' or 'X'
    '''

    mutations = read_in_mutations(lineage)
    
    variants = []
    if gene in mutations:
        for pos in mutations[gene]:
            variants.append(sample[pos])
    else:
        variants = ['']
    
    return(";".join(variants))