# Using the GRAVI Workflow

## Introduction

Gene Regulatory Analysis using Variable Inputs (GRAVI) is a `snakemake` workflow designed for the differential binding analysis and integration of multiple ChIP targets.
These can be distinct ChIP targets within the same biological sample, or the same ChIP target across multiple biological contexts.

- `GRAVI` requires aligned reads as `bam` files as the starting point for the workflow.
- A complete html document will be produced, able to be uploaded to a git repository, along with multiple output files containing peak locations, target genes and integrated results.

## Quickstart

1. Create a new `github` repository using this template
2. Download the repository to your local server or HPC using `git clone <myrepository>`
3. Place your bam files in the subdirectory `data/bam` or `data/aligned`
    + These should be placed in separate directories for each target, such `data/bam/target1` and `data/bam/target2`
4. Edit `samples.tsv` in the `config` directory. Required columns are: `sample`, `target`, `treat`, `input`
5. Modify any parameters in `config/config.yml`

To run using 16 cores without any queuing system (e.g. on a local machine), enter the following

```
snakemake -p --use-conda --notemp --keep-going --cores 16
```

