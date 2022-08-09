# Snakemake pipeline for PopSV

Adapt the configuration file [`config.yaml`](config.yaml).

```sh
snakemake --configfile config.yaml --cores 2
```

## Requirements

- R with the PopSV package installed.
- Python with [snakemake](https://snakemake.readthedocs.io/en/stable/) and [pandas](https://pandas.pydata.org/) installed
    - For example, to install with [pip](https://pip.pypa.io/en/stable/): `pip3 snakemake pandas`

## Test data

This is just to test that the code can run.
It uses a small slice of chr1 and chr2 for one sample, and naively simulate bin counts for others.

```sh
## Download a slice of a BAM from the 1000 Genomes Project
samtools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 1:1-1000000 > samp1.sam
samtools view ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 2:1-1000000 >> samp1.sam
samtools view -b -h samp1.sam > samp1.bam
rm samp1.sam
samtools index samp1.bam

## Simulate bin counts
Rscript sim_dummy_bc.R

## Extract the bin definition from one simulated sample
zcat bc-samp2.tsv.gz | cut -f 1-3 > bin-500bp.tsv

## Snakemake
snakemake --configfile config.yaml -np --cores 2

## Clean up
rm -rf test_*
```

*Remove the `-n` (dry run) to actually run the pipeline.*

## Running the snakemake pipeline in HPC environments

Some relevant links:

- [How to setup a profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and [examples of profiles for some common HPCs](https://github.com/snakemake-profiles/doc).
- [If your HPC uses modules](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html?highlight=hpc#using-environment-modules)
