import os
import re
from datetime import datetime
import pandas as pd
from Bio import SeqIO

def detect_passage(passage: str) -> str:
    '''
    helper function to detect passage
    '''

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
    input_fasta, 
    input_csv, 
    filter_fasta, 
    passage_h, 
    subtype, 
    gene):

    dateparse = lambda x: datetime.strptime(x, '%d/%m/%Y')

    dump_seq = SeqIO.parse(input_fasta, "fasta")
    data = pd.read_csv(
        input_csv, 
        infer_datetime_format = True, 
        parse_dates = ['Sample Date'], 
        date_parser=dateparse, 
        dtype=str,
        na_filter=False,
        )
    
    #param setting
    gene = {
        '4': 'HA',
        '6': 'NA',
        '1': 'PB2',
        '2': 'PB1',
        '3': 'PA',
        '5': 'NP',
        '7': 'MP'
        }

    data[['id','segment']] = data['Seq No'].str.split('.',expand=True)
    data['gene'] = data['segment'].map(gene)
    data['Month'] = data['Sample Date'].dt.strftime('%b').str.lower()
    data['new_designation'] = data['Designation'].replace(" ", "_") + '_' + data['Month']
    data['passage_short'] = data['Passage History'].apply(detect_passage)

    # filter only the row that contains information aboeut Designation, HA clade, Sub Type, and Passage History
    processed_data = data.loc[(data['Designation'].isnull() + data['HA Clade'].isnull() + data['Sub Type'].isnull()+ data['Passage History'].isnull()) ==0 ]

    # filter only the row that match with our query
    processed_data = processed_data.loc[processed_data['Sub Type'] == subtype]
    processed_data = processed_data.loc[processed_data['passage_short'] == passage_h]
    processed_data = processed_data.loc[processed_data['gene'] == gene]

    # use the info from data to subset only sequences you want
    subset_fasta = [] 
    for seq in dump_seq:
        if seq.id in list(processed_data['Seq No']):
            subset_fasta.append(seq)

    #iterate through them to rename and write into new fasta file
    with open(filter_fasta, 'w') as f_fasta:
        for seq in subset_fasta:
            name = list(processed_data.loc[processed_data['Seq No'] == seq.id, 'name'])[0]
            seq.id = name
            seq.description = name
            SeqIO.write(seq, f_fasta, 'fasta')

    #write clean metafile
    colname = ['Seq No','Lab #','Lineage', 'Designation',' Passage History' , 'Sample Date',' Results',	'GISAID #','HA Clade', 'name']
    processed_data[colname].to_csv(filter_csv,index=False)