args <- commandArgs(TRUE)
success <- args[[1]]

pkgs <- c(
  "tidyverse", "magrittr", "rtracklayer", "glue", "pander", "scales",
  "plyranges", "yaml", "UpSetR", "BiocParallel", "scico", "extraChIPs",
  "DiagrammeR", "steveped/spBioUtils", "Gviz", "Rsamtools", "GenomicInteractions",
  "cowplot", "GenomeInfoDb", "ngsReports", "csaw", "edgeR", "quantro",
  "qsmooth", "statmod", "IHW", "ggrepel", "rlang", "ggside", "InteractionSet",
  "reactable", "htmltools", "VennDiagram", "Biostrings", "GenomicRanges",
  "steveped/extraChIPs"
)

pks2inst <- setdiff(pkgs, rownames(installed.packages()))
cat("Attempting to install", pks2inst)

BiocManager::install(pks2inst,  update = TRUE, ask = FALSE, force = FALSE)

stopifnot(
  all(basename(pkgs) %in% rownames(installed.packages()))
)

file.create(success)
