import re
import typer
import pandas as pd
from Bio import SeqIO, SeqRecord
from pathlib import Path
from collections import defaultdict
from typing import List, Set, Tuple, Union, Iterable, Dict
#  meta['Seq No'].values.tolist()

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
    if re.search(r'[Ss][Ii][Aa][Tt]', passage):
        return('')
    if re.search(r'[Mm][Dd][Cc][Kk]', passage):
        return('')
    if re.search(r'[Oo]riginal', passage):
        return('o')
    if re.search(r'[Ss]pecimen', passage):
        return('o')
    if re.search(r'[Ee]\d', passage):
        return('e')
    if re.search(r'cs', passage):
        return('o')
    else:
        return('')

def rename_fasta(
    sequences: List[SeqRecord.SeqRecord],
    meta_data: pd.DataFrame, 
    add_gene: bool, 
    add_passage: bool,
    add_month: bool,
    ) -> Dict[str,SeqRecord.SeqRecord]:
    '''
    rename fasta
    '''

    meta_data[['id','segment']] = meta_data['Seq No'].str.split('.',expand=True)
    meta_data['gene'] = meta_data['segment'].map(segements_genes)
    # what if the date is missing?
    meta_data['Month'] = meta_data['Sample Date'].dt.strftime('%b').str.lower()
    meta_data['passage_short'] = meta_data['Passage History'].apply(detect_passage)
    meta_data['new_designation'] = meta_data['Designation'].replace(" ", "_")

    if add_month:
        meta_data['new_designation'] = meta_data['new_designation'] + meta_data['passage_short']
    if add_passage:
        meta_data['new_designation'] = meta_data['new_designation'] + '_' + meta_data['Month']
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
    return(sequences_dict) # type: ignore
    
def write_meta(meta: pd.DataFrame, output_dir: Path, split_by):
    '''
    write metadata to file
    optionally split by gene
    '''

    if split_by in ['multi', 'individual']:
        meta.to_csv(output_dir / "meta.tsv", sep='\t', index=False, na_rep='')
    if split_by == 'gene':
        meta_split = [x for _, x in meta.groupby(meta['segment'])]
        for df in meta_split:
            segment = df['segment'].unique()[0]
            gene = segements_genes[segment].lower()
            output_name = output_dir / f"{gene}.tsv"
            df.to_csv(output_name, sep = "\t", index=False)

def write_sequences(
    sequences: SeqRecord.SeqRecord, 
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
        for item in sequences:
            with open(output / f"{item}.fasta", 'a') as handle:
                SeqIO.write(sequences[item], handle, "fasta")
    else:
        print("error hit? search for this haha")

def fuzee_get(
    batch_num: str
    ):
    '''
    get batch data from fuzee api
    '''
    # use the Data > 'GA - Sequencing' page to download all data
    # apply the same parsing settings as utils.read_meta
    raise typer.BadParameter("This feature is not implemented.")

def find_fasta(
    seq_num: List[str],
    input_dir: Path,
    ) -> Tuple[List[SeqRecord.SeqRecord], Set]:
    '''
    find and concat fasta files across multiple dirs from a list of csv inputs 
    return - list(SeqIO.SeqRecord), list()

    seq_no    : list - if not specified then get from fuzee 
    input_dir : Path - specify to search specific dir
    '''
    
    sequence_names =  [num + '.fasta' for num in seq_num] 
    fasta_paths = input_dir.glob("*.fasta")
    
    matched = set([p.name for p in fasta_paths]) & set(sequence_names)
    seq_paths = [input_dir / m for m in matched]
    sequences = [SeqIO.read(seq, "fasta") for seq in seq_paths]

    return(sequences, matched)