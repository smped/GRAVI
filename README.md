# GRAVI: Gene Regulatory Analysis using Variable Inputs

This is a `snakemake` workflow for:

1. Performing sample QC
2. Calling ChIP peaks
3. Performing Differential Binding Analysis
4. Comparing results across ChIP targets
5. Identifying TF motifs (to be added soon)

The minimum required input is one ChIP target with two conditions.

Full documentation can be found [here](https://steveped.github.io/GRAVI/)


## Snakemake Implementation

The basic workflow is written `snakemake` and can be called using the following steps.

Firstly, setup the required conda environments

```
snakemake \
	--use-conda \
	--conda-prefix '/home/steveped/mambaforge/envs/' \
	--conda-create-envs-only \
	--cores 1
```

Secondly, create and inspect the rulegraph

```
snakemake --rulegraph > workflow/rules/rulegraph.dot
dot -Tpdf workflow/rules/rulegraph.dot > workflow/rules/rulegraph.pdf
```

Finally, the workflow itself can be run using:

```
snakemake \
	-p \
	--use-conda \
	--conda-prefix '/home/steveped/mambaforge/envs/' \
	--notemp \
	--keep-going \
	--cores 16
```

Note that this creates common environments able to be called by other workflows and is dependent on the user.
For me, my global conda environments are stored in `/home/steveped/mambaforge/envs/`.
For other users, this path will need to be modified.
