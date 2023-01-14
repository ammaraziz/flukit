import re
import pandas as pd
from Bio import SeqIO, SeqRecord
from datetime import datetime
from pathlib import Path

segements_genes = {
    '4': 'HA',
    '6': 'NA',
    '1': 'PB2',
    '2': 'PB1',
    '3': 'PA',
    '5': 'NP',
    '7': 'MP'
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

def rename(
    input_fasta: Path, 
    input_csv: Path, 
    output_fasta: Path,
    ):

    dateparse = lambda x: datetime.strptime(x, '%d/%m/%Y')

    sequences = SeqIO.parse(input_fasta, "fasta")
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
    meta['new_designation'] = meta['Designation'].replace(" ", "_") + meta['passage_short'] + '_' + meta['Month']

    designations = dict(zip(meta['Seq No'], meta['new_designation']))

    #iterate through them to rename and write into new fasta file
    with open(output_fasta, 'w') as handle:
        for seq in sequences:
            seq.id,seq.description  = designations[seq.id], designations[seq.id]
            SeqIO.write(seq, handle, 'fasta')
    
def write_seq(sequences: SeqRecord, output: Path, split_by: str = None):
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

def find(
    input_dir: Path,
    input_csv: Path,
    batch_num: str,
    output_dir: Path):
    '''
    find fasta files in a given directory 
    split_by - str - [gene, all]
    '''

