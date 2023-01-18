import re
import pandas as pd

import Bio
from Bio import SeqIO
from datetime import datetime
from pathlib import Path
from utils import read_meta

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

# recode to accept input as SeqRecords instead of reading from path
# record to return renamed SeqRecords
def rename_fasta(
    input_fasta: Path, 
    input_csv: Path, 
    output_fasta: Path,
    add_gene: bool = True, 
    add_passage: bool = True,
    add_month: bool = True,
    ):
    '''
    rename fasta
    '''
    dateparse = lambda x: datetime.strptime(x, '%d/%m/%Y')

    sequences = SeqIO.parse(input_fasta, "fasta")
    # replace with utils.read_meta
    meta = pd.read_csv(
        input_csv, 
        infer_datetime_format = True, 
        parse_dates = ['Sample Date'], 
        date_parser=dateparse, 
        dtype=str,
        na_filter=False,
        )

    meta[['id','segment']] = meta['Seq No'].str.split('.',expand=True)
    meta['gene'] = meta['segment'].map(segements_genes)
    meta['Month'] = meta['Sample Date'].dt.strftime('%b').str.lower()
    meta['passage_short'] = meta['Passage History'].apply(detect_passage)
    meta['new_designation'] = meta['Designation'].replace(" ", "_")

    if add_month:
        meta['new_designation'] = meta['new_designation'] + meta['passage_short']
    if add_passage:
        meta['new_designation'] = meta['new_designation'] + '_' + meta['Month']
    if add_gene:
        meta['new_designation'] = meta['new_designation'] + '_' + meta['gene']

    designations = dict(zip(meta['Seq No'], meta['new_designation']))

    with open(output_fasta, 'w') as handle:
        for seq in sequences:
            seq.id,seq.description  = designations[seq.id], designations[seq.id]
            SeqIO.write(seq, handle, 'fasta')

def write_meta(meta: pd.DataFrame, output: Path, split_by = None):
    '''
    write metadata to file
    optionally split by gene
    '''

    if not split_by:
        meta.to_csv(output, sep=',', index=False, na_rep='')
    if split_by == 'gene':
        meta_split = [x for _, x in meta.groupby(meta['segment'])]
        for df in meta_split:
            segment = df['segment'].unique()[0]
            gene = segements_genes[segment].lower()

            output_name = output / f"{gene}.csv"
            df.to_csv(output_name, sep = ",", index=False)

def write_sequences(
    sequences: list[Bio.SeqRecord], 
    output: Path, 
    split_by: str = None):
    '''
    write sequences to file
    optionally, split output by gene, or single
    '''
    # write out one file
    if not split_by:
        with open(output, 'w') as handle:
            SeqIO.write(sequences, handle, 'fasta')
    # write out individual files
    elif split_by == 'single':
        for seq in sequences:
            print(seq.id)
            print(output / f"{seq.id}.fasta")
            SeqIO.write(seq, output / f"{seq.id}.fasta", "fasta")
    # write out files as genes
    elif split_by == 'gene':
        for seq in sequences:
            segment = seq.id.split(".")[1]
            gene = segements_genes[segment].lower()
            gene_output = output / f"{gene}.fasta"

            with open(gene_output, 'a') as handle:
                SeqIO.write(seq, handle, "fasta")
    else:
        print("error hit")

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
    seq_no: list,
    input_dir: Path = None
    ) -> tuple(list, set):
    '''
    find and concat fasta files across multiple dirs from a list of csv inputs 
    return - list(SeqIO.SeqRecord)

    seq_no    : list - if not specified then get from fuzee 
    input_dir : Path - optional, specify to search specific dir
                     - default to search predefined locations
    '''
    
    # read in meta file
    sequence_names =  [seq + '.fasta' for seq in seq_no] 

    # find fasta files 
    if input_dir:
        fasta_files = input_dir.glob("*.fasta")
    if not input_dir:
        for path in DEFAULT_FASTA_PATHS:
            fasta_files = path.glob("*.fasta")
    
    matched = set([f.name for f in fasta_files]) & set(sequence_names)
    seq_paths = [input_dir / match for match in matched]
    sequences = [SeqIO.read(seq, "fasta") for seq in seq_paths]

    return(sequences, matched)


