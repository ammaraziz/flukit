from pandas import DataFrame
from rich import print
from rich.progress import track
from Bio import SeqRecord
from pathlib import Path
from .variants import get_ref, get_snps, get_ha_snps, get_pa_snps
from .align_frames import align
from .utils import write_fasta
from .clades import run_nextclade, update_dataset

# process data

def call_variants(sequences: dict, lineage: str) -> tuple[DataFrame, list]:

    results = DataFrame.from_dict(
        data = {
            "ha_aa":[],
            "NA":[],
            "MP":[],
            "PA":[],
            "vacc_ref":[]
            }, 
            dtype=str
            )

    # List to write out ha sequences
    ha_records = []

    for record in track(sequences, description="Processing..."):
        gene = sequences[record].gene
        
        try:
            seqAA, refAA = align(lineage = lineage, input_record = sequences[record])

            if len(seqAA) == 0:
                print(f"[bold yellow]{record} failed to translate. Out of frame.[/bold yellow]")
                continue            
        except Exception as error:
            print(f"[bold yellow]Error: {error} - {record} [/bold yellow]")
            break

        try:
            if gene in ['PB2', 'PB1', 'NS', 'NP']:
                return('')
            elif gene in ['MP', 'NA']:
                results.at[record, gene] = get_snps(seqAA, gene, lineage)
            elif gene == 'PA':
                results.at[record, gene] = get_pa_snps(seqAA, lineage)
            elif gene == 'HA':
                results.at[record, 'ha_aa'] = get_ha_snps(seqAA, refAA)
                results.at[record, 'vacc_ref'] = get_ref(lineage)
                ha_records.append(sequences[record])
        
        except Exception as error:
            print(f"[bold yellow]Issue with calling variants on sample {sequences[record].id}. Error: {error}[/bold yellow]")
            pass
    # write results, return path to results
    
    return(results, ha_records)

def call_clades(
    ha_records: list, 
    lineage: str,
    output: Path = None,
    update: bool = False):
    '''
    Call clades with nextclade. 
    The list of `sequences` are written to temp file for nextclade. 
    '''

    ha_temp = write_fasta(ha_records)

    if update:
        update_dataset(lineage)
    
    tsv = run_nextclade(ha_temp, lineage, output)
    return(tsv)
    

