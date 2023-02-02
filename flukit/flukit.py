#!/usr/bin/env python3

import os
import sys
import typer
import pandas

from pathlib import Path
from rich import print
from Bio import SeqIO
from flukit.utils.variants import set_gene
import flukit.utils.run as run

app = typer.Typer(
    help = "flukit - the influenza surveillance toolkit",
    add_completion=False,
    no_args_is_help=True,
    )

@app.callback()
def main():
    pass

@app.command(no_args_is_help=True, help = "Call variants, mutations and clades on sequences (all genes)")
def variants(
    sequences: Path = typer.Option(
        ...,
        "-s",
        "--sequences",
        help="Path to sequences."), 
    lineage: str = typer.Option(
        ...,
        "-l",
        "--lineage",
        help="lineage options are: h1n1, h3n2, vic."),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="output path"
    ),
    batchNumber: str = typer.Option(
        None,
        "-b",
        "--batchNumber",
        help="prefix used for output files, optional.")
     ):

    # input checks
    if lineage not in ['h1n1', 'h3n2', 'vic']:
        raise typer.BadParameter(f"Linease '{lineage}' is not valid. You must choose from: \n h1n1\n h3n2\n vic")
    if not Path(sequences).resolve():
        raise typer.BadParameter(f"The path is not correct, please check: {sequences}")

    # outputs
    if batchNumber is not None:
        results_out = Path(output) / (batchNumber + "_results.csv")
        clades_out = Path(output) / (batchNumber + "_clades.txt")
    else:
        results_out = Path(output) / "results.csv"
        clades_out = Path(output) / "clades.txt"

    # parse input sequences
    try:
        input_sequences = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
        input_sequences = set_gene(input_sequences)
    except FileNotFoundError:
        raise typer.BadParameter(f"Check input: {sequences} \nFile does not exist.")
    except ValueError as error:
        raise typer.BadParameter(f'''Check input: {sequences} \nError: {error}''')
    except:
        raise typer.BadParameter(
            f"""
            [bold yellow] Error reading in fasta file. Check input: {sequences} \nError: {sys.exc_info()[0]}
            """
        )
    # call variants
    variants, ha_records = run.call_variants(input_sequences, lineage)
    
    # call clades
    clades = run.call_clades(ha_records, lineage)
    
    # combine dataframes
    results = pandas.merge(
        variants, clades, 
        how="left", left_index=True, 
        right_on="seqName", 
        suffixes=(False, False)
        )
    results.insert(0, 'seqName', results.pop('seqName')) # reorder 
    
    # write results
    results.to_csv(results_out, sep=',', index=False)
    print("[bold green]All done![/bold green]")

@app.command(no_args_is_help=True, help = "Find and rename fasta files. Default: find, do not rename, output as a single multi.fasta")
def find(
    input_dir: Path = typer.Option(
        ...,
        "-i",
        "--input-dir",
        help="Input directory containing fasta files"), 
    input_meta: Path = typer.Option(
        ...,
        "-m",
        "--input-meta",
        help="csv/tsv file with the following headers: Seq No, Designation, Sample Date, Passage History"
    ),
    output_dir: Path = typer.Option(
        ...,
        "-o",
        "--output-dir",
        help="Output directory for fasta and meta files"), 
    split_by: str = typer.Option(
        'multi',
        "-sb",
        "--split-by",
        help="Split fasta by: gene, multi"), 
    batch_num: str = typer.Option(
        None,
        "-b",
        "--batch-num",
        help="If specified will retreive meta data from Fuzee via API"), 
    rename: bool = typer.Option(
        False,
        help="Rename fasta"),
        ):

        os.makedirs(output_dir, exist_ok=True)
        # checks
        for dir in [input_dir, input_meta]:
            if not Path(dir).resolve():
                raise typer.BadParameter(f"The path is not correct, please check: {dir}")
        if split_by not in ['gene', 'multi']:
            raise typer.BadParameter(
                f"""Not acceptable value given:  {split_by}
                please choose from:
                [gene]      output file per gene
                [multi]     single multifasta file
                """
                )
        # run``
        run.findrename(
            input_dir=input_dir,
            input_meta = input_meta,
            output_dir=output_dir,
            split_by=split_by,
            batch_num=batch_num,
            rename=rename
        )
        print("[bold green]All done![/bold green]")