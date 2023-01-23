import re
import Bio
import rich
import typer
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict

from .utils import read_meta

# meta['Seq No'].values.tolist()

segements_genes = {
    '4' : 'HA',
    '6' : 'NA',
    '1' : 'PB2',
    '2' : 'PB1',
    '3' : 'PA',
    '5' : 'NP',
    '7' : 'MP',
    '8' : 'NS',
    }

def detect_passage(passage: str) -> str:
    '''
    helper function to detect passage, returning abbreviation
    '''
    # cells (siat, mdck) have no abbr
    if re.search('[Ss][Ii][Aa][Tt]', passage):
        return('')
    if re.search('[Mm][Dd][Cc][Kk]', passage):
        return('')
    if re.search('[Oo]riginal', passage):
        return('o')
    if re.search('[Ss]pecimen', passage):
        return('o')
    if re.search('[Ee]\d', passage):
        return('e')
    if re.search('cs', passage):
        return('o')
    else:
        return('')

def sanatise_designations(designation: str) -> str:
    '''
    remove/replace unwanted characters and spaces from designations
    " " -> "_"
    "'" -> ""
    '''
    designation = designation.replace(" ", "_")
    designation = designation.replace("'", "")
    return(designation)

def rename_fasta(
    sequences: Bio.SeqRecord, 
    meta_data: pd.DataFrame, 
    add_gene: bool, 
    add_passage: bool,
    add_month: bool,
    ) -> list[Bio.SeqRecord]:
    '''
    rename fasta
    '''

    meta_data[['id','segment']] = meta_data['Seq No'].str.split('.',expand=True)
    meta_data['gene'] = meta_data['segment'].map(segements_genes)
    meta_data['Month'] = meta_data['Sample Date'].dt.strftime('%b').str.lower()
    meta_data['passage_short'] = meta_data['Passage History'].apply(detect_passage)
    meta_data['new_designation'] = meta_data['Designation'].apply(sanatise_designations)

    if add_passage:
        meta_data['new_designation'] = meta_data['new_designation'] + meta_data['passage_short']
    if add_month:
        # handle missing date/month by replacing nan with ''
        meta_data['new_designation'] = meta_data['new_designation'] + ('_' + meta_data['Month']).fillna('')
    if add_gene:
        meta_data['new_designation'] = meta_data['new_designation'] + '_' + meta_data['gene']
    designations = dict(zip(meta_data['Seq No'], meta_data['new_designation']))

    # rename
    sequences_dict = defaultdict(list)
    for seq in sequences:
        segment = segements_genes[seq.id.split(".")[1]]
        seq.id = designations[seq.id]
        seq.description = ''
        
        sequences_dict[segment].append(seq)
    return(sequences_dict)

def write_meta(meta: pd.DataFrame, output: Path, split_by):
    '''
    write metadata to file
    optionally split by gene
    '''

    if split_by in ['multi', 'ind']:
        meta.to_csv(output, sep='\t', index=False, na_rep='')
    if split_by == 'gene':
        meta_split = [x for _, x in meta.groupby(meta['segment'])]
        for df in meta_split:
            segment = df['segment'].unique()[0]
            gene = segements_genes[segment].lower()
            output_name = output.parent / f"{gene}.tsv"
            df.to_csv(output_name, sep = "\t", index=False)

def write_sequences(
    sequences: dict[Bio.SeqRecord], 
    output: Path, 
    split_by: str,
    ):
    '''
    write sequences to file
    optionally, split output by: ind, gene, multi
    '''
    # individual fasta output
    if split_by == 'ind':
        for seq in sequences:
            SeqIO.write(seq, output / f"{seq.id}.fasta", "fasta")
    # multifasta output
    elif split_by == 'multi':
        with open(output / "multi.fasta", 'w') as handle:
            SeqIO.write(sequences.values(), handle, 'fasta')
    # gene fasta output if not renamed
    elif split_by == 'gene':
        for key, value in sequences.items():
            with open(output / f"{key}.fasta", 'a') as handle:
                SeqIO.write(value, handle, "fasta")
    else:
        print("error hit? search for this haha")

def fuzee_get(
    batch_num: int
    ):
    '''
    get batch data from fuzee api
    '''
    # use the Data > 'GA - Sequencing' page to download all data
    # apply the same parsing settings as utils.read_meta
    raise typer.BadParameter("This feature is not implemented.")

def find_fasta(
    seq_num: list,
    input_dir: Path,
    ) -> tuple[list, set]:
    '''
    find and concat fasta files across multiple dirs from a list of csv inputs 

    seq_no    : list - if not specified then get from fuzee 
    input_dir : Path - specify to search specific dir
    '''
    
    sequence_names =  [num + '.fasta' for num in seq_num] 

    if input_dir:
        fasta_paths = input_dir.glob("*.fasta")
    fasta_names = [p.name for p in fasta_paths]
    #matched = set(fasta_names) & set(sequence_names)
    unmatched = []
    matched = []
    for name in sequence_names:
        if name in fasta_names:
            matched.append(name)
        else:
            rich.print(f"[yellow] WARNING: Missing {input_dir}/{name} [/yellow]")
            unmatched.append(name)
            
    seq_paths = [input_dir / m for m in matched]
    sequences = [SeqIO.read(seq, "fasta") for seq in seq_paths]
    
    return(
        sequences,
        [m.replace(".fasta","") for m in matched], 
        [m.replace(".fasta","") for m in unmatched]
        )