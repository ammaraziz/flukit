import typer
from rich import print
from pathlib import Path
from pandas import DataFrame
from rich.progress import track

from .variants import get_ref, get_snps, get_ha_snps, get_pa_snps, get_na_mp_snps, get_na_snps
from .align_frames import align
from .utils import write_temp_fasta, read_meta
from .clades import run_nextclade, update_dataset
from .rename import find_fasta, rename_fasta, write_sequences, write_meta, fuzee_get

def call_variants(sequences: dict, lineage: str) -> tuple[DataFrame, list]:
    results = DataFrame.from_dict(
        data = {"ha_aa":[],"NA":[],"MP":[],"PA":[], "anti-viral":[],"vacc_ref":[]}, 
        dtype='str'
            )

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
            antiviral_str = ""
            if gene in ['PB2', 'PB1', 'NS', 'NP']:
                pass
            elif gene == 'NA':
                results.at[record, gene] = get_na_snps(seqAA, gene, lineage)
                antiviral_str = get_gene_snps(seqAA, refAA, gene, lineage)
            elif gene == 'MP':
                results.at[record, gene] = get_snps(seqAA, gene, lineage)
                antiviral_str =  get_gene_snps(seqAA, refAA, gene, lineage)
            elif gene == 'PA':
                results.at[record, gene] = get_gene_snps(seqAA, refAA, lineage)
                antiviral_str = get_gene_snps(seqAA, refAA, lineage)
            elif gene == 'HA':
                results.at[record, 'ha_aa'] = get_ha_snps(seqAA, refAA)
                results.at[record, 'vacc_ref'] = get_ref(lineage)
                ha_records.append(sequences[record])
            results.at[record, 'anti-viral'] = antiviral_str
        except Exception as error:
            print(f"[bold yellow]Issue with calling variants on sample {sequences[record].id}. Error: {error}[/bold yellow]")
            pass   
    
    return(results, ha_records)

# clade calling
def call_clades(
    ha_records: list, 
    lineage: str,
    output: Path = None,
    update: bool = False) -> DataFrame:
    '''
    Call clades with nextclade. 
    The list of `sequences` are written to temp file for nextclade. 
    '''

    ha_temp = write_temp_fasta(ha_records)

    if update:
        update_dataset(lineage)
    
    tsv = run_nextclade(ha_temp, lineage, output)
    return(tsv)


## work in progress - master function for calling directly
def findrename(
    input_dir: Path, 
    input_meta: Path, 
    output_dir: Path, 
    split_by: str,
    batch_num: str,
    rename: bool,
    ) -> None:
    '''
    Find and rename fasta files
    rename if split_by is gene or multi else do not rename output as is
    '''
    
    if batch_num and input_meta:
        raise typer.BadParameter(f"specify either but not both: batch_num/input_meta")
    if batch_num: # not implemented
        meta = fuzee_get(batch_num)
    else:
        meta = read_meta(input_meta)

    seq_num = list(meta['Seq No'])
    sequences, matched, unmatched = find_fasta(seq_num=seq_num, input_dir=input_dir)

    if rename:
        if split_by == "multi":
            sequences_renamed = rename_fasta(
                sequences=sequences, 
                meta_data=meta, 
                add_gene=True, 
                add_month=True, 
                add_passage=True
            )
        if split_by == "gene":
            sequences_renamed = rename_fasta(
                sequences=sequences, 
                meta_data=meta, 
                add_gene=False, 
                add_month=True, 
                add_passage=True
                )
    else:
        sequences_renamed = sequences

    meta_matched = meta[ meta["Seq No"].isin(matched) ]
    meta_unmatched = meta[ meta["Seq No"].isin(unmatched) ]

    write_sequences(
        sequences=sequences_renamed, 
        output=output_dir, 
        split_by=split_by
        )
    write_meta(
        meta=meta_matched, 
        output=output_dir / "meta_matched.tsv", 
        split_by='multi'
        )
    if not meta_unmatched.empty:
        write_meta(
            meta=meta_unmatched, 
            output=output_dir / "meta_unmatched.tsv", 
            split_by='multi'
            )
    