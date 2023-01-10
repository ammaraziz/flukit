# Flukit - simple variant caller for influenza

Used interally at WHOFLUCC. Not recommended for external usage.

### install

1. clone this repo
2. install using pip

```
cd flukit && python -m pip install .
```
3. Install `nextclade`:

```
conda install -c bioconda nextclade
```

### usage

Run `flukit --help` to see detailed instructions

- input
	- fasta, file with headers ending in `.1`... `.4` which represent the gene. For example: `MySample.4` is the HA gene of MySample.
	- lineage, one of h1n1, h3n2, vic
	- the batch name used as prefix for output files
	- the output path

Example:

```
flukit -s flu.fasta -l h3n2 -b batch150 -o ~/Desktop/
```

- output
	- tsv file in the following format:
	```
	seqno - the fasta header
	ha_aa - HA mutations called against vacc_ref
	na - H275Y mutation
	mp - S31N mutation
	pa - I38X mutation
	vacc_ref - called against ancestral strains
	```

Output is formatted specifically for internal database.
