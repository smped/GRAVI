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
  - [ ] Macs2 Summary
  - [ ] Shared Peaks
  - [ ] Differential Signal
  - [ ] Save Best matches as an rds to keep sequence information

## Differential Signal

- [x] Restrict methods to:
    1. Sliding Windows with SQN (quantro not essential)
    2. Sliding Windows with TMM/RLE (after quantro)
    3. Fixed width using TMM/RLE (after quantro)
- [x] Update main differential signal module
- [x] Update IHW
- [ ] Update RNA-Seq module
- [ ] Rewrite script/rules for creating diff_signal Rmds with updated params
- [ ] Redefine config settings to have a default, with optional overwrite using the target 


## Features (Unlikely)

- [ ] Add 3-way comparisons
- [ ] Additional modules
  - [ ] ROSE
  - [ ] Nucleosome Free Regions

## Bugs

- [ ] Add `pairwise_comparisons/{t1}_{t2}/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}-de_genes.csv` to output of pairwise_comparisons when RNA-Seq data is provided
- [ ] Find & fix all hacks related to name errors in `plyranges/extraChIPs` coercion & old package versions
- [ ] Cleanup handling of detected genes, particularly in the setup of annotations

## Additional Notes

- Add `rGREAT` analysis as the default enrichment for all modules?
- `IHW` is currently not installable using conda with Bioc 3.18. Checks for updates may be prudent
- The package simplifyEnrichment also looks useful for comparing across targets. Examples here: https://jokergoo.github.io/rGREAT_suppl/compare_online_and_local.html

### Motifs & Enrichment Testing

- The Bioconductor package `memes` only implements `AME` out of the tools of interst
- The standalone MEME suite may be able to be added in a modular fashion
    + Centrimo looks highly relevant, as does AME
- HOMER may be simpler and has an R implementation as [`marge`](https://robertamezquita.github.io/marge/index.html)
    + `-nfr ` performed analysis of the nucleosome free region of Histone marks. Not available in `marge`
- Need to figure out `regioneR`
    + Use `regionR` to compare peak sets across targets?
