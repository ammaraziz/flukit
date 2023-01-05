#!/usr/bin/env python3

import typer
import pandas as pd
from pathlib import Path
from rich.progress import track
from rich import print
from Bio import SeqIO
from .utils.variants import get_variants, set_gene, get_vacc_ref

app = typer.Typer(
    help = "flukit - the influenza surveillance toolkit... kinda",
    add_completion=False)

@app.command(no_args_is_help=True)
def main(
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
        raise typer.BadParameter(
            f"{lineage} is not valid. You must choose from h1n1, h3n2, vic"
            )
    if not Path(sequences).resolve():
        raise typer.BadParameter(
            f"The path is not correct, please check: {sequences}"
        )

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
    except Exception:
        raise typer.BadParameter(
            "Error reading in input sequences. Check input seqs are properly formatted, headers must end in gene number"
            )

    sample_records = pd.DataFrame.from_dict(
        data = {
            "ha_aa":[],
            "na":[],
            "mp":[],
            "pa":[],
            "vacc_ref":[]
            },
        dtype=str)
    
    for record in track(input_sequences, description="Processing..."):
        gene = input_sequences[record].gene
        try:
            if gene == 'MP':
                sample_records.at[record, 'mp'] = get_variants(input_sequences[record], lineage)
            if gene == 'NA':
                sample_records.at[record, 'na'] = get_variants(input_sequences[record], lineage)
            if gene == 'PA':
                sample_records.at[record, 'pa'] = get_variants(input_sequences[record], lineage)
            if gene == 'HA':
                sample_records.at[record, 'ha_aa'] = get_variants(input_sequences[record], lineage)
                sample_records.at[record, 'vacc_ref'] = get_vacc_ref(lineage)
        except Exception as e:
            print(f"Error: {e}")
            pass
    print("[bold green]All done![/bold green]")
    # write out
    sample_records.to_csv(results_out, sep=',', index_label="seqno")
