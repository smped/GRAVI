# PDX_HCI-005_AR_ER_GATA3_ChIP

Combines two separate experiments for the PDX line HCI-005

- AR & GATA3 ChIP, as performed by Leila Hosseinzadeh
- ER & H3K27ac ChIP as used for the publication by [Hickey *et al*](https://www.nature.com/articles/s41591-020-01168-7)

## Config Setup

### Samples

This workflow requires a tsv file, nominally `samples.tsv` detailing:

1. Each sample
2. Which ChIP target each samples is associated with
3. Which treatment group it is associated with, and
4. Which Input/Control sample it is associated with

Required columns are: `sample`, `target`, `treat` and `input` (all lower case).
At least one additional column containing replicate information should also be included.
A possible structure is as follows:

```
| sample | target | treat | passage | input |
| ------ | ------ | ----- | ------- | ----- |
| sample1 | H3K27ac | Veh | 1 | input1 |
| sample2 | H3K27ac | E2  | 1 | input1 |
```

It is currently assumed that bam files will be placed in `data/aligned/bam/**` where the value `**` represents each individual ChIP target.
Whilst the root directory can theoretically be changed via `config.yml`, it is not recommended.
Bam files should be named as specified in the `sample` column, with the addition of the `.bam` suffix only.

### Config

The file `config/config.yml` is where users can edit a series of parameters with default values provided.
The example provided should provide a clear guide to the structure.
However, the file can be checked using

```
python scripts/check_yaml.py
```


## Snakemake implementation

The basic workflow is written `snakemake` and can be called using the following steps.

Firstly, setup the required conda environments

```
snakemake \
	--use-conda \
	--conda-prefix '/home/steveped/miniconda3/envs/' \
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
