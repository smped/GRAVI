## General

- [ ] Update all Rmd for RELEASE 3.18
    - [x] annotations
    - [x] macs2
    - [ ] differential signal
    - [ ] pairwise comparisons
- [x] Update r environment for RELEASE 3.18
- [x] Change terminology from differential binding to differential signal
- [ ] Check for compatibility with extraChIPs >= v1.7.1
  - [x] macs2
  - [ ] differential signal
  - [ ] pairwise comparisons
- [ ] Checks for input file consistency/structure
- [ ] Shift `cowplot` to `patchwork`
- [ ] Motif Detection
  - [x] Macs2 Summary
  - [x] Shared Peaks
  - [ ] Differential Signal
- [ ] Update index for more complete description of parameters 
  
## Annotation Setup

- [x] Shift all motif URIs to an rds instead of docs/assets. Use a tempdir to create them
  
## Peak Analysis

- [x] Add goseq as well as GREAT

## Differential Signal

- [x] Restrict methods to:
    1. Sliding Windows with SQN (quantro not essential)
    2. Sliding Windows with TMM/RLE (after quantro)
    3. Fixed width using TMM/RLE (after quantro)
- [x] Update main differential signal module
- [x] Update IHW
- [ ] Shift analysis to a separate script & only reporting in the Rmd/HTML
    + This will allow motif/regioner analysis of diff_sig peaks to be included
    + Add params used to the metadata of the results
- [ ] Update RNA-Seq module
    + Wrangle multiple files
- [ ] Rewrite script/rules for creating diff_signal Rmds with updated params
- [ ] Redefine config settings to have a default, with optional overwrite using the target 


## Features (Unlikely)

- [ ] Add 3-way comparisons
- [ ] Additional modules
  - [ ] ROSE
  - [ ] Nucleosome Free Regions (Try HisTrader)

## Bugs

- [ ] Add `pairwise_comparisons/{t1}_{t2}/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}-de_genes.csv` to output of pairwise_comparisons when RNA-Seq data is provided
- [ ] Find & fix all hacks related to name errors in `plyranges/extraChIPs` coercion & old package versions
- [ ] Cleanup handling of detected genes, particularly in the setup of annotations

## Additional Notes

- The package simplifyEnrichment also looks useful for comparing across targets. Examples here: https://jokergoo.github.io/rGREAT_suppl/compare_online_and_local.html

