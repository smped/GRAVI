## Features

- Add 3-way comparisons
- Incorporate motif detection
- Checks for input file consistency/structure
- Add Jon Bischlak's install of R4.2.0
- Update manual to include all recent changes

## Bugs

- Add `pairwise_comparisons/{t1}_{t2}/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}-de_genes.csv` to output of pairwise_comparisons when RNA-Seq data is provided
- Find & fix all hacks related to name errors in `plyranges/extraChIPs` coercion & old package versions
- Cleanup handling of detected genes, particularly in the setup of annotations
