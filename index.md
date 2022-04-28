* TOC {:toc}

# Using the GRAVI Workflow

## Introduction

Gene Regulatory Analysis using Variable Inputs (GRAVI) is a `snakemake` workflow designed for differential binding analysis and subsequent integration of results across multiple ChIP targets.
These can be distinct ChIP targets within the same biological sample, or the same ChIP target across multiple biological contexts.

- `GRAVI` requires aligned reads as `bam` files as the starting point for the workflow.
- A complete html document will be produced, able to be uploaded to a git repository, along with multiple output files containing peak locations, target genes and integrated results.
- Optional data files include:
    + Results from an RNA-Seq analysis
	+ Key HiC Interactions
	+ Externally derived features, such as those determined from histone marks
	+ Additional coverage tracks for inclusion in all plots

## Quickstart

You will need a snakemake conda environment to begin. 
Please see [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for help setting this up.

1. Create a new `github` repository on your account by going to [the github template repository](https://github.com/steveped/GRAVI/generate)
2. Download the repository to your local server or HPC using `git clone <myrepository>`
3. Place your bam files in the subdirectory `data/bam` or `data/aligned`
    + These should be placed in separate directories for each target, such `data/bam/target1` and `data/bam/target2`
4. Edit `samples.tsv` in the `config` directory. Required columns are: `sample`, `target`, `treat`, `input`
5. Modify any parameters in `config/config.yml`

To run using 16 cores without any queuing system (e.g. on a local machine), enter the following

```
snakemake -p --use-conda --notemp --keep-going --cores 16
```

