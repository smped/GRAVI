# To-Do List

## Checks

- [] Add checks for external files in Snakefile
- [] Add checks for YAML input. Maybe put the external file checks in this section

## General Workflow

- [x] Calls to `setup_chunk.Rmd` use `here::here()` which depends on an *.Rproj file or `.here`. This needs to be resolved intelligently, perhaps creating one, or incorporating the point below
- [x] Setup chunk should be written to `analysis/setup_chunk.Rmd` not placed in `workflow/modules`
- [x] Move initiation of `index.Rmd` to be from `workflow/modules`
- [x] Change site classification to use secondary p-values from a point-based H~0~
- [] Add the capacity to download key tables as a csv where viable
- [x] Check bias offsets for goseq (just check results). 
    + May be best to pass the argument `method = c("Wallenius", "Hypergeometric")` via the workflow
    + Perhaps even add Fisher's Non-Central Hypergeometric as an option
- [] Remove github installation for `extraChIPs` & migrate all to R4.2
- [] Check to see if `geom_x/ysidelabel()` has made it on to CRAN yet
- [] Change export of oracle peaks to be bed files, and include `target_treat_` prefix
- [] Add `target_` prefix to consensus peaks. Maybe include score?
- [] Change output paths for `differential_binding` and `pairwise_comparison` modules

## Motifs

- [] Include motif enrichment as part of the standard workflow
- [] Include pathway enrichment

