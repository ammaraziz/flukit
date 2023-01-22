import os
import re
import tempfile
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
from collections import defaultdict
from pandas import read_csv
from pathlib import Path
from importlib_resources import files

config_path = files('flukit').joinpath('config')

def load_features(reference, feature_names=None):
    '''
    Parse genbank file to extract features
    '''

    if not os.path.isfile(reference):
        print("ERROR: reference sequence not found. looking for", reference)
        return None

    features = {}
    for feat in SeqIO.read(reference, 'genbank').features:
        if feat.type == 'CDS':
            if "locus_tag" in feat.qualifiers:
                fname = feat.qualifiers["locus_tag"][0]
                if feature_names is None or fname in feature_names:
                    features[fname] = feat
            elif "gene" in feat.qualifiers:
                fname = feat.qualifiers["gene"][0]
                if feature_names is None or fname in feature_names:
                    features[fname] = feat
        elif feat.type == 'source':  # read 'nuc' as well for annotations - need start/end of whole!
            features['nuc'] = feat

    return features


def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'MX'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    # sequences not mod 3 give messy BiopythonWarning, so avoid by padding.
    if len(sequence) % 3:
        sequence_padded = sequence + "N"*(3-len(sequence) % 3)
    else:
        sequence_padded = sequence
    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq.Seq(sequence_padded).translate(gap='-'))
    except TranslationError:
        translation_exception = True
        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence_padded)
        codons = np.frombuffer(
            str_seq[:len(str_seq) - len(str_seq) % 3].encode(), dtype='S3').astype("U")
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence

def get_reference(input_lineage: str, input_gene: str) -> SeqRecord:
    '''
    Parameters
        input_lineage : str
            'h1n1', 'h3n2', 'vic', 'yam'
        input_gene : str
            'pb2', 'pb1', 'pa', 'ha',' np', 'na', 'mp', 'ns'
    Return
        SeqRecord : containing reference_gene.genbank
    '''

    lineages = ['h1n1', 'h3n2', 'vic', 'yam']
    genes = ['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns']

    if input_lineage.lower() in lineages and input_gene.lower() in genes:
        refname = f"{config_path}/reference_{input_lineage.lower()}_{input_gene.lower()}.gb"
        try:
            ref = SeqIO.read(refname, 'genbank')
        except Exception as e:
            raise ValueError(
                "Error: %s. While reading reference, check file is genbank format and it exists" % e)
    else:
        raise ValueError(
            f"Incorrect lineage or gene entered, check input: {input_lineage}, {input_gene.lower()}")

    return(refname, ref)

def read_in_mutations(lineage: str) -> dict:
    '''
    Reads in tab-seperated file containing subtype specific mutations

    Parameters
    ----------
    clade_file : str
        meta data file

    Returns
    -------
    dict
        {gene : pos, ...}
    '''
    mutation_file = config_path / f'mutations_{lineage}.tsv'
    mutations = defaultdict(list)
    
    df = read_csv(mutation_file, sep='\t',  na_filter=False)
    for index, row in df.iterrows():
        mutations[row.gene].append(row.pos)
    mutations.default_factory = None

    return(mutations)

def write_temp_fasta(sequences: list, output: Path = None) -> Path:
    '''
    Write out SeqIO.dict object to file and return location. 
    If path is None use temp file.

    return - Path   :   location of output file
    '''

    if output is None:
        file = tempfile.NamedTemporaryFile(delete=False)
    else:
        file = output
    SeqIO.write(sequences, file.name, "fasta")

    return(file.name)

def read_meta(meta_data: Path, column: str = None):
    '''
    wrapper function to safely read in fasta files and optionally return a specific column

    Return either pandas dataframe or list
    '''
    import pandas as pd
    from datetime import datetime
    dateparse = lambda x: datetime.strptime(x, '%d/%m/%Y')

    if meta_data.name.split('.')[1] == 'csv':
        sep = ","
    else:
        sep = "\t"

    try:
        meta = pd.read_csv(
            meta_data, 
            infer_datetime_format = True, 
            parse_dates = ['Sample Date'], 
            date_parser=dateparse, 
            dtype=str,
            na_filter=False,
            sep=sep,
            )
            
        if column:
            return(list(meta[column]))
        else:
            return(meta)

    except OSError as error:
        raise OSError(f"File does not exist. Error: {error}")
    except Exception as error:
        raise Exception(f"Uncontrolled error detected: Error: {error}")
 