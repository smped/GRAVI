#' Create the local ZScores enabled by regioneReloaded
#' This is a process requiring large amounts of RAM and multi-threading
#' The recommended number of permutations is 5000
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
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

## Manual lists for testing. Will be overwritten by snakemake objects...
# config <- yaml::read_yaml("config/config.yml")
# all_input <- list(
#   regions = "output/annotations/gene_regions.rds",
#   features = "output/annotations/features.rds",
#   peaks = "output/peak_analysis/AR/AR_consensus_peaks.bed.gz",
#   params = "config/params.yml"
# )
# all_output <- list(rds = "output/peak_analysis/AR/AR_regions_localz.rds")
# all_params <- list(ntimes = 200)
# all_wildcards <- list(target = "AR")

config <- slot(snakemake, "config")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
all_wildcards <- slot(snakemake, "wildcards")
cat_list(all_input, "input")
cat_list(all_output, "output")
cat_list(all_params, "params")
cat_list(all_wildcards, "wildcards")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

## Required input files are the gene-regions and eternal-features, as well as
## a set of peaks to be compared against these regions.
## The enrichment params yaml can also be passed here

## Required params are: adj_p_method
## It will be assumed that 5K permutations will be performed and that window
## sizes are +/-5kb for a 10kb region

## Required output files are:
## file.path(macs2_path, "{target}", "{target}_regions_localz.rds")

cat_time("Loading packages...")
library(regioneReloaded)
library(extraChIPs)
library(plyranges)
library(readr)
library(yaml)
library(rlang)
cat_time("done")

source(here::here("workflow/scripts/custom_functions.R"))
ucsc <- get_ucsc(config$genome$build)

cat_time("Loading all regions...")
regions <- read_rds(all_input$regions)
cat_time("Loading all features...")
features <- read_rds(all_input$features)
cat_time("Forming test_regions...")
test_regions <- c(regions, features)
test_regions <- endoapply(test_regions, granges)
cat_time("Setting genome to be", ucsc$build)
sq <- seqinfo(regions)
genome(sq) <- ucsc$build
seqinfo(test_regions) <- sq
cat_time(" done")

cat_time("Loading peaks from", all_input$peaks)
peaks <- importPeaks(all_input$peaks, seqinfo = sq, type = "bed")
peaks <- unlist(peaks)
cat_time(" done")

cat_time("Loading params from", all_input$params)
regioner_params <- read_yaml(all_input$params)$regioner
cat_time(" done")

threads <- slot(snakemake, "threads")[[1]] - 1
cat_time("Running multiLocalZscore with", threads, "threads...")
mlz_params <- list(
  A = peaks, Blist = test_regions, sampling = FALSE,
  ranFUN = "resampleGenome", evFUN = "numOverlaps",
  max_pv = 1, genome = ucsc$build, mc.cores = threads
)
mlz_params <- c(mlz_params, regioner_params[c("ntimes", "step", "window")])
mlz <- do.call("multiLocalZscore", mlz_params)
cat_time("Done")


cat_time("Writing to", all_output$rds)
write_rds(mlz, all_output$rds, compress = "gz")
cat_time("done")
