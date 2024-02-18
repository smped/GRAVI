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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

## Manual lists for testing. Will be overwritten by snakemake objects...
# config <- yaml::read_yaml("config/config.yml")
# all_input <- list(
#   regions = "output/annotations/gene_regions.rds",
#   features = "output/annotations/features.rds",
#   peaks = "output/macs2/AR/AR_consensus_peaks.bed.gz",
#   params = "config/params.yml"
# )
# all_output <- "output/macs2/AR/AR_regions_localz.rds"

config <- slot(snakemake, "config")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")$rds
cat_list(all_input, "input")
cat("Output will be written to ", all_output, "\n")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- here::here(all_output)

## Required input files are the gene-regions and eternal-features, as well as
## a set of peaks to be compared against these regions.
## The enrichment params yaml can also be passed here

## Required params are: adj_p_method
## It will be assumed that 5K permutations will be performed and that window
## sizes are +/-5kb for a 10kb region

## Required output files are:
## file.path(macs2_path, "{target}", "{target}_regions_localz.rds")

cat("Loading packages...")
library(regioneReloaded)
library(extraChIPs)
library(plyranges)
library(readr)
library(yaml)
library(GenomicRanges)
library(S4Vectors)
cat("done\n")


map <- c(
  hg19 = "hg19", hg38 = "hg38", grch37 = "hg19", grch38 = "hg38",
  mm10 = "mm10", mm39 = "mm39", grcm38 = "mm10", grcm39 = "mm39",
  rn7 = "rn7", mratbn7.2 = "rn7", galgal6 = "galGal6", rhemac10 = "rheMac10",
  canfam5 = "canFam5", susscr11 = "susScr11", pantro6 = "panTro6", dm6 = "dm6"
)
bld <- match.arg(tolower(config$genome$build), names(map))
ucsc_ref <- map[[bld]]

cat("Loading all regions...\n")
regions <- read_rds(all_input$regions)
cat("Loading all features...\n")
features <- read_rds(all_input$features)
cat("Forming test_regions...\n")
test_regions <- c(regions, features)
test_regions <- endoapply(test_regions, granges)
cat("Setting genome to be", ucsc_ref)
sq <- seqinfo(regions)
genome(sq) <- ucsc_ref
seqinfo(test_regions) <- sq
cat(" done\n")

cat("Loading peaks from", all_input$peaks)
peaks <- importPeaks(all_input$peaks, seqinfo = sq, type = "bed")
peaks <- unlist(peaks)
cat(" done\n")

cat("Loading enrichment params from", all_input$params)
enrich_params <- read_yaml(all_input$params)$enrichment
cat(" done\n")

threads <- slot(snakemake, "threads")[[1]] - 1
cat("Starting at", format(Sys.time(), "%H:%M:%S, %d %b %Y"), "\n")
cat("Running multiLocalZscore with", threads, "threads...\n")
mlz <- multiLocalZscore(
  A = peaks, Blist = test_regions,
  ranFUN = "resampleGenome",
  evFUN = "numOverlaps",
  window = 5e3,
  ntimes = 5000,
  adj_pv_method = enrich_params$adj,
  max_pv = 1,
  genome = ucsc_ref,
  mc.cores = threads
)
cat("Finished at", format(Sys.time(), "%H:%M:%S, %d %b %Y"), "\n")

cat("Writing to", all_output, "\n")
write_rds(mlz, all_output, compress = "gz")
cat("done\n")
