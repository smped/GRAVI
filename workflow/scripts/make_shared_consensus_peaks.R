## This script takes the complete set of consensus peaks and finds the shared regions
##
## Required inputs are
##
## 1. All target-specific consensus peaks
## 2. The seqinfo object
##
## Output will be
##
## 1. output/peak_analysis/shared/shared_consensus_peaks.bed.gz
##
## First handle any conda weirdness
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

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

## For testing
# all_input <- list(
#   peaks = c(
#     "../GRAVI_testing/output/peak_analysis/AR/AR_consensus_peaks.bed.gz",
#     "../GRAVI_testing/output/peak_analysis/ER/ER_consensus_peaks.bed.gz",
#     "../GRAVI_testing/output/peak_analysis/H3K27ac/H3K27ac_consensus_peaks.bed.gz"
#   ),
#   sq = "../GRAVI_testing/output/annotations/seqinfo.rds"
# )
# all_output <- list(
#   peaks = "../GRAVI_testing/output/peak_analysis/shared/shared_consensus_peaks.bed.gz"
# )
# config <- yaml::read_yaml("../GRAVI_testing/config/config.yml")

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(extraChIPs)
library(plyranges)
library(tidyverse)

cat_time("Loading seqinfo and defining ranges to exclude...\n")
sq <- read_rds(all_input$sq)

cat_time("Loading peaks...\n")
n_targets <- length(all_input$peaks)
cons_peaks <- all_input$peaks %>%
  importPeaks(seqinfo = sq, type = "bed") %>%
  setNames(str_remove(names(.), "_consensus.+"))
## NB: for n = 2, this will be the union peaks
shared_peaks <- cons_peaks %>%
  makeConsensus(p = 1 - 1 / n_targets, method = "coverage") %>%
  filter(n == n_targets)
  
cat_time("Writing", length(shared_peaks), "peaks to", all_output$peaks, "\n")
write_bed(shared_peaks, all_output$peaks)
cat_time("Done\n")
