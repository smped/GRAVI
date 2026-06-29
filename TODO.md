## General

- [ ] Streamline setting species into one locations
- [ ] Separate `regioner_localz_targets.R` so that it only runs one target against the others, then add a new script to merge into the single rds. This should give a 1-2hr speed up in run times
- [ ] Fix bug in labels for Donut plots
- [ ] Improve profile heatmaps
  
## Annotation Setup

  
## Peak Analysis

- [ ] Allow broadPeak as well as narrowPeak
- [ ] Confirm best approach for consensus peaks

## Motif Analysis

- [ ] Add capacity for analysis by cluster
- [ ] Make sure Z-scores are documented as preferred
- [ ] Use both `name` and `altname` columns in modules to ensure all motifs are viable
- [ ] Perform clustering in a standalone Rmd in order to check clustering params in a preliminary run


## Enrichment Testing

- [ ] Enable using custom gene-sets beyond `msigdb`

## Differential Signal

- [ ] Separate setting of normalisation method & analysis model in the config

## Pairwise Analysis

### NFR

## Features (Unlikely)

- [ ] Add 3-way comparisons
- [ ] Additional modules
  - [ ] ROSE?





## Bugs

- [ ] Add `pairwise_comparisons/{t1}_{t2}/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}-de_genes.csv` to output of pairwise_comparisons when RNA-Seq data is provided



