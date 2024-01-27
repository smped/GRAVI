# GRAVI: Gene Regulatory Analysis using Variable Inputs

This is a `snakemake` workflow for:

1. Performing sample QC
2. Calling ChIP peaks
3. Performing Differential Signal Analysis
4. Comparing pairwise results across ChIP targets

The minimum required input is one ChIP target with two conditions.
Full documentation can be found [here](https://smped.github.io/GRAVI/)

## Snakemake Implementation

The basic workflow is written `snakemake`, requiring at least v7.7, and can be called using the following steps.

Firstly, setup the required conda environments

```
snakemake \
	--use-conda \
	--conda-prefix '/home/stevieped/mambaforge/envs/' \
	--conda-create-envs-only \
	--cores 1
```

It should be noted that the current version requires the R package `extraChIPs` to be updated which requires an active internet connection.
Given internet connectivity is often restricted on HPC clusters, it may be prudent to run this one rule in an interactive session before job submission.

```
snakemake \
	--use-conda \
	--conda-prefix '/home/stevieped/mambaforge/envs/' \
	--allowed-rules update_extrachips \
	--cores 1
```

Finally, the workflow itself can be run using:

```
snakemake \
	-p \
	--use-conda \
	--conda-prefix '/home/stevieped/mambaforge/envs/' \
	--notemp \
	--rerun-triggers mtime \
	--keep-going \
	--cores 16
```

If submitting on an HPC cluster, a snakemake profile will need to be configured first.
Please contact your HPC administrator for assistance with this


Note that this creates common environments able to be called by other workflows and is dependent on the user.
For me, my global conda environments are stored in `/home/stevieped/mambaforge/envs/`.
For other users, this path will need to be modified.

If wishing to tidy the directory after a successful run, you can check which non-essential files can be deleted using `snakemake -n --delete-temp-output --cores 1`.
If the files earmarked for deletion are considered to be non-essential, they can be deleted by removing the `-n` flag from the above code: `snakemake --delete-temp-output --cores 1`.
The bedgraph files produced by `macs2 callpeak` and `macs2 bdgcmp` are typically very large, and once converted to bigwig files, can be safely deleted which will free up considerable disk space.
