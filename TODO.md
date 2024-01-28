## General

- [ ] Update all Rmd for RELEASE 3.18
    - [x] annotations
    - [ ] macs2
    - [ ] differential signal
    - [ ] pairwise comparisons
- [x] Update r environment for RELEASE 3.18
- [x] Change terminology from differential binding to differential signal
- [ ] Check for compatibility with extraChIPs >= v1.6
- [ ] Redefine all config settings to have a default, with optional overwrite using the target type
- [ ] Add required column to samples.tsv to indicate signal type. This should be restricted to the values ['TF', "H3K27ac'] with the intent of adding 'ATAC' later. Or maybe this is easier just using the above TODO point.
  - This may require removing the 'differential_h3k27ac' module
- [ ] Checks for input file consistency/structure
- [ ] Shift `cowplot` to `patchwork`

## Differential Signal

- [ ] Restrict methods to:
  1. Sliding Windows with SQN (quantro not essential)
  2. Sliding Windows with TMM/RLE (after quantro)
  3. Fixed width using TMM/RLE (after quantro)

This should require writing a DiffSig module for sliding and another for fixed

## Features (Unlikely)

- [ ] Add 3-way comparisons
- [ ] Incorporate motif detection

## Bugs

- [ ] Add `pairwise_comparisons/{t1}_{t2}/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}-de_genes.csv` to output of pairwise_comparisons when RNA-Seq data is provided
- [ ] Find & fix all hacks related to name errors in `plyranges/extraChIPs` coercion & old package versions
- [ ] Cleanup handling of detected genes, particularly in the setup of annotations

## Additional Notes

- `IHW` is currently not installable using conda with Bioc 3.18. Checks for updates may be prudent
