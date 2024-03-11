#' This script uses motifTestR to analyse motifs
#'
#' Required input files is:
#' - config/params.yml
#' - output/macs2/{target}/{target}_consensus_peaks.bed.gz
#' - output/annotations/gene_regions.rds
#' - output/annotations/motif_list.rds
#' - output/annotations/exclude_ranges.rds
#' - output/checks/here.chk
#' - output/checks/r-packages.chk
#'
#' Expected output is:
#' - output/macs2/{target}/{target}_motif_enrichment.tsv.gz
#' - output/macs2/{target}/{target}_motif_position.tsv.gz
#'
#' Required config parameters:
#' - genome$build
#'
#' Required enrichment parameters may possibly be
#' - n_iterations
#' - model
#' - peak_width
#' - binwidth
#'
#' Testing will involve
#' - Running `testMotifPos()` setting abs = TRUE. There's no real need to test symmetrically around zero
#' - Running `testMotifEnrich()` using the negbinom model & 100 iterations
#'
#' Can results also be written to a tsv/html during this step, or should that wait?
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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

# all_input <- list(
#   here = "output/checks/here.chk",
#   exclude_ranges = "output/annotations/exclude_ranges.rds",
#   gene_regions = "output/annotations/gene_regions.rds",
#   motifs = "output/annotations/motif_list.rds",
#   packages = "output/checks/r-packages.chk",
#   params = "config/params.yml",
#   peaks = "output/macs2/AR/AR_consensus_peaks.bed.gz"
# )
#
# all_output <- list(
#   enrich = "output/macs2/AR/AR_motif_enrichment.tsv.gz",
#   pos = "output/macs2/AR/AR_motif_position.tsv.gz"
# )
# config <- list(genome = list(build = "GRCh37"))
# threads <- 4

config <- slot(snakemake, "config")
threads <- slot(snakemake, "threads")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
cat_list(all_input, "input -", sep = ":")
cat_list(all_output, "output -", sep = ":")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...")
library(motifTestR)
library(tidyverse)
library(extraChIPs)
library(yaml)
library(scales)
library(universalmotif)

## Contains the function for finding UCSC build info
source(here::here("workflow/scripts/custom_functions.R"))
ucsc <- get_ucsc(config$genome$build)
pkg <- paste(c("BSgenome", ucsc$sp, "UCSC", ucsc$build), collapse = ".")
cat_time("Loading", pkg)
library(pkg, character.only = TRUE)

cat_time("Reading params")
params <- read_yaml(all_input$params)[["motif_analysis"]][["analysis"]]
cat_list(params, "params", ":")

gene_regions <- read_rds(all_input$gene_regions)
sq <- seqinfo(gene_regions)
## Convert to the same as the BSgenome package
genome(sq) <- ucsc$build
seqinfo(gene_regions) <- sq
cat_time("Reading Peaks...")
peaks <- importPeaks(all_input$peaks, type = "bed", seqinfo = sq) |>
  unlist() |>
  granges() |>
  unname() |>
  resize(width = params$peak_width, fix = "center")
peaks$region <- bestOverlap(peaks, gene_regions)
cat_time("Read", comma(length(peaks)), "peaks")

cat_time("Getting peak sequences...")
bs_genome <- get(pkg)
test_seq <- getSeq(bs_genome, peaks)
names(test_seq) <- as.character(peaks)
n_seq <- length(test_seq)
cat_time("done")

motif_list <- read_rds(all_input$motifs)
cat_time("Checking for low frequency matches")
min_matches <- params$ignore_below * n_seq
ignore <- countPwmMatches(motif_list, test_seq, mc.cores = threads) < min_matches
cat_time(
  "Found", sum(ignore), "motifs with matches in <",
  percent(params$ignore_below), "of sequences"
)

cat_time("Started getting best matches")
matches <- getPwmMatches(
  motif_list[!ignore], test_seq, best_only = TRUE, mc.cores = threads
)
cat_time("done")
gc()

cat_time("Testing for Positional Bias")
pos_res <- testMotifPos(
  matches, binwidth = params$binwidth, abs = params$abs, mc.cores = threads
)
gc()
cat_time("done")
cat_time("Writing", all_output$pos)
pos_res |>
  as_tibble(rownames = "altname") |>
  left_join(to_df(motif_list), by = "altname") |>
  dplyr::select(ends_with("name"), cluster, all_of(colnames(pos_res))) |>
  dplyr::select(-contains("consensus")) |>
  write_tsv(all_output$pos)

cat_time("Importing exclude_ranges")
exclude_ranges <- read_rds(all_input$exclude_ranges)

cat_time("Generating RMRanges based on provided gene_regions")
rm_ranges <- makeRMRanges(
  split(peaks, peaks$region), gene_regions, exclude = exclude_ranges,
  n_iter = params$iterations, mc.cores = threads
)
gc()
cat_time("Sampled", comma(length(rm_ranges)), "RMRanges")

cat_time("Extracting RMSeq")
rm_seq <- getSeq(bs_genome, rm_ranges)
mcols(rm_seq) <- mcols(rm_ranges)
cat_time("done")

cat_time("Testing for motif enrichment")
enrich_res <- testMotifEnrich(
  motif_list, test_seq, rm_seq, model = params$model, mc.cores = threads
)
cat_time("Done")
cat_time("Writing", all_output$enrich)
enrich_res |>
  as_tibble(rownames = "altname") |>
  left_join(to_df(motif_list), by = "altname") |>
  dplyr::select(ends_with("name"), cluster, all_of(colnames(enrich_res))) |>
  dplyr::select(-contains("consensus")) |>
  write_tsv(all_output$enrich)
