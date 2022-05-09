args <- commandArgs(TRUE)
success <- args[[1]]

pkgs <- c(
  "tidyverse", "magrittr", "rtracklayer", "glue", "pander", "scales",
  "plyranges", "yaml", "UpSetR", "BiocParallel", "scico", "DiagrammeR",
  "Gviz", "Rsamtools", "GenomicInteractions", "cowplot",
  "GenomeInfoDb", "ngsReports", "csaw", "edgeR", "quantro", "qsmooth",
  "statmod", "IHW", "ggrepel", "rlang", "ggside", "InteractionSet",
  "reactable", "htmltools", "VennDiagram", "Biostrings", "GenomicRanges"
)

pks2inst <- setdiff(pkgs, rownames(installed.packages()))
cat("Attempting to install", pks2inst)

BiocManager::install(pks2inst,  update = TRUE, ask = FALSE, force = FALSE)

stopifnot(
  all(basename(pkgs) %in% rownames(installed.packages()))
)

remotes::install_github(
  "steveped/extraChIPs",
  ## This is the last commit before shifting to R4.2
  ref = "de3e53111a93d664311a11c2c185ea3286d2a9a7",
  force = FALSE
)

stopifnot("extraChIPs" %in% rownames(installed.packages()))

file.create(success)
