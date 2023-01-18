import re
import pandas as pd
import Bio
from Bio import SeqIO
from pathlib import Path
from .utils import read_meta

#  meta['Seq No'].values.tolist()

DEFAULT_FASTA_PATHS = [
    Path('S:/Shared/WHOFLU/mol_biol/00-New Sequences/'),
    Path('/mnt/Sdrive/WHOFLU/mol_biol/00-New Sequences/'),
    Path('S:/Shared/WHOFLU/mol_biol/Sequencing-NGS/'),
    Path('/mnt/Sdrive/WHOFLU/mol_biol/Sequencing-NGS/'),
    ]

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

def rename_fasta(
    sequences: Bio.SeqRecord, 
    meta_data: pd.DataFrame, 
    add_gene: bool = True, 
    add_passage: bool = True,
    add_month: bool = True,
    ) -> list[Bio.SeqRecord]:
    '''
    rename fasta
    '''

    meta_data[['id','segment']] = meta_data['Seq No'].str.split('.',expand=True)
    meta_data['gene'] = meta_data['segment'].map(segements_genes)
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

    for seq in sequences:
        seq.id, seq.description  = designations[seq.id], designations[seq.id]
    
    return(sequences)

def write_meta(meta: pd.DataFrame, output_dir: Path, split_by):
    '''
    write metadata to file
    optionally split by gene
    '''

    if split_by in ['multi', 'single']:
        meta.to_csv(output_dir / "meta.tsv", sep='\t', index=False, na_rep='')
    if split_by == 'gene':
        meta_split = [x for _, x in meta.groupby(meta['segment'])]
        for df in meta_split:
            segment = df['segment'].unique()[0]
            gene = segements_genes[segment].lower()
            output_name = output_dir / f"{gene}.tsv"
            df.to_csv(output_name, sep = "\t", index=False)

def write_sequences(
    sequences: list[Bio.SeqRecord], 
    output: Path, 
    split_by: str = None):
    '''
    write sequences to file
    optionally, split output by: single, gene, multi
    '''
    # multifasta output
    if not split_by or split_by == 'multi':
        with open(output / "multi.fasta", 'w') as handle:
            SeqIO.write(sequences, handle, 'fasta')
    # individual fasta output
    elif split_by == 'single':
        for seq in sequences:
            SeqIO.write(seq, output / f"{seq.id}.fasta", "fasta")
    # gene fasta output
    elif split_by == 'gene':
        for seq in sequences:
            segment = seq.id.split(".")[1]
            gene = segements_genes[segment].lower()
            gene_output = output / f"{gene}.fasta"

            with open(gene_output, 'a') as handle:
                SeqIO.write(seq, handle, "fasta")
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
    pass

def find_fasta(
    seq_num: list,
    input_dir: Path = None
    ) -> tuple[list, set]:
    '''
    find and concat fasta files across multiple dirs from a list of csv inputs 
    return - list(SeqIO.SeqRecord), list()

    seq_no    : list - if not specified then get from fuzee 
    input_dir : Path - optional, specify to search specific dir
                     - default to search predefined locations
    '''
    
    sequence_names =  [num + '.fasta' for num in seq_num] 

    if input_dir:
        fasta_paths = input_dir.glob("*.fasta")
    if not input_dir:
        for path in DEFAULT_FASTA_PATHS:
            fasta_paths = path.glob("*.fasta")
    
    matched = set([p.name for p in fasta_paths]) & set(sequence_names)
    seq_paths = [input_dir / m for m in matched]
    sequences = [SeqIO.read(seq, "fasta") for seq in seq_paths]

    return(sequences, matched)