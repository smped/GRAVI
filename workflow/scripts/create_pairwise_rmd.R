#' Create the basic header for the differential signal files.
#'
#' Really just need the wildcards to determine all values
#'
#' Also the output will contain the rmd path
#'
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}
## A function for printing input
cat_list <- function(x, slot = NULL, sep = "\n\t"){
  nm <- setdiff(names(x), "")
  invisible(
    lapply(
      nm,
      \(i) cat("Received", slot, i, sep, paste0( x[[i]], "\n\t"), "\n")
    )
  )
}
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

# all_input <- list(
#   results = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-pairwise_results.rds",
#   dsa1 = "output/differential_signal/AR/AR_E2_E2DHT-differential-signal.rds",
#   dsa2 = "output/differential_signal/ER/ER_E2_E2DHT-differential-signal.rds",
#   localz = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-pairwise_localz.rds",
#   module = "workflow/modules/pairwise_comparison.Rmd",
#   motif_enrichment = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-motif_enrichment.tsv.gz",
#   motif_position = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-motif_position.tsv.gz"
# )
# all_wildcards <- list(
#   tgt1 = "AR",
#   tgt2 = "ER",
#   comp1 = "E2_E2DHT",
#   comp2 = "E2_E2DHT"
# )
# all_output <- list(
#   rmd = "analysis/AR_E2_E2DHT-ER_E2_E2DHT_pairwise_comparison.Rmd"
# )


log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_wildcards <- slot(snakemake, "wildcards")
cat_list(all_input, "input:")
cat_list(all_output, "output", sep = ":")
cat_list(all_wildcards, "wildcards", sep = ":")


## Solidify file paths
# all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

# Set the full comparison
full_comp <- with(
  all_wildcards, paste0(tgt1, "_", comp1, "-", tgt2, "_", comp2)
)
cat_time("full_comp is", full_comp)

cat_time("Loading packages")
library(tidyverse)
library(glue)


cat_time("Writing file header")
hdr <- glue("---\ntitle: 'Pairwise Comparison: {full_comp}'
date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
bibliography: references.bib
link-citations: true
params:
  tgt1: \"{all_wildcards$tgt1}\"
  tgt2: \"{all_wildcards$tgt2}\"
  comp1: \"{all_wildcards$comp1}\"
  comp2: \"{all_wildcards$comp2}\"
  full_comp: \"{full_comp}\"
  dsa1: \"{all_input$dsa1}\"
  dsa2: \"{all_input$dsa2}\"
  localz: \"{all_input$localz}\"
  motif_enrichment: \"{all_input$motif_enrich}\"
  motif_position: \"{all_input$motif_position}\"
  results: \"{all_input$results}\"\n---\n\n"
)
cat(hdr) 
write_lines(hdr, all_output$rmd)


cat_time("Adding main body of file")
## Now add the rest of the module
file.append(all_output$rmd, here::here(all_input$module))
cat_time("Done")
