#' This script uses motifTestR to analyse motifs after DSA
#'
#' Required input files is:
#' - output/differential_signal/{target}/{target}_{ref}_{treat}-differential-signal.rds
#' - output/annotations/motif_list.rds
#' - output/checks/here.chk
#' - output/checks/r-packages.chk
#'
#' Expected output is:
#' - output/differential_signal/{target}/{target}_{ref}_{treat}_motif_enrichment.tsv.gz
#' - output/differential_signal/{target}/{target}_{ref}_{treat}_motif_position.tsv.gz
#'
#' Required config parameters:
#' - genome$build
#'
#' Required enrichment parameters may possibly be
#' - peak_width
#' - binwidth
#'
#' Testing will involve
#' - Running `testMotifPos()` setting abs = TRUE. There's no real need to test symmetrically around zero
#' - Running `testMotifEnrich()` using the hypergeometric distribution and Unchanged sites as baseline
#'
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
#   here = "output/checks/here.chk",
#   motifs = "output/annotations/motif_list.rds",
#   packages = "output/checks/r-packages.chk",
#   results = "output/differential_signal/ER/ER_E2_E2DHT-differential-signal.rds"
# )
#
# all_output <- list(
#   enrich = "output/differential_signal/ER/ER_E2_E2DHT_motif_enrichment.tsv.gz",
#   pos = "output/differential_signal/ER/ER_E2_E2DHT_motif_position.tsv.gz"
# )
# all_params = list(
#   motif_params = list(
#       abs = TRUE,
#       adj = "fdr", # not used here
#       alpha = 0.05, # not used here
#       binwidth = 10,
#       break_ties = "all",
#       peak_width = 400,
#       ignore_below = 0.01,
#       iterations = 100,
#       min_score = "80%",
#       model = "quasipoisson"
#     )
# )
# config <- list(genome = list(build = "GRCh37"))
# threads <- 4

# log <- slot(snakemake, "log")[[1]]
# message("Setting stdout to ", log, "\n")
# sink(log, split = TRUE)
# config <- slot(snakemake, "config")
# threads <- slot(snakemake, "threads")
# all_input <- slot(snakemake, "input")
# all_output <- slot(snakemake, "output")
# all_params <- slot(snakemake, "params")
# all_wildcards <- slot(snakemake, "wildcards")

cat_list(all_input, "input:", sep = "=")
cat_list(all_output, "output:", sep = "=")
cat_list(all_params, "params:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Running with", threads, "threads")
cat_time("Loading packages...")
library(motifTestR)
library(tidyverse)
library(extraChIPs)
library(yaml)
library(scales)
library(universalmotif)
library(plyranges)

motif_params <- all_params$motif_params

## Contains the function for finding UCSC build info
source(here::here("workflow/scripts/custom_functions.R"))
ucsc <- get_ucsc(config$genome$build)
pkg <- paste(c("BSgenome", ucsc$sp, "UCSC", ucsc$build), collapse = ".")
cat_time("Loading", pkg)
library(pkg, character.only = TRUE)

cat_time("Reading Results...")
results <- read_rds(all_input$results)
win_type <- metadata(results)$window_type
genome(results) <- ucsc$build
sq <- seqinfo(results)
rng_col <- c(fixed = "centred_peak", sliding = "keyval_range")[[win_type]]

cat_time("Resizing and recentering ranges")
results <- results |>
  colToRanges(rng_col) |>
  resize(width = motif_params$peak_width, fix = 'center') |>
  set_genome_info(ucsc$build)
results <- splitAsList(results, results$status == "Unchanged")
names(results) <- c("changed", "unchanged")

cat_time("Getting peak sequences...")
bs_genome <- get(pkg)
seq <- results |>
  lapply(\(x) getSeq(bs_genome, x)) |>
  lapply(\(x) x[letterFrequency(x, "N")[,1] == 0]) |>
  lapply(\(x) setNames(x, as.character(x))) |>
  as("DNAStringSetList")
n_seq <- sum(map_int(seq, length))
cat_time("done")

motif_list <- read_rds(all_input$motifs)
cat_time("Checking for low frequency matches")
min_matches <- motif_params$ignore_below * n_seq
counts <- lapply(
  seq, \(x) countPwmMatches(
      motif_list, x, min_score = motif_params$min_score, mc.cores = threads
    )
  )
ignore <- rowSums(do.call("cbind", counts)) < min_matches
cat_time(
  "Found", sum(ignore), "motifs with matches in <",
  percent(motif_params$ignore_below), "of sequences"
)

cat_time("Started getting best matches")
matches <- getPwmMatches(
  motif_list[!ignore], seq$changed, best_only = TRUE,
  min_score = motif_params$min_score, break_ties = motif_params$break_ties,
  mc.cores = threads
)
cat_time("done")
gc()

cat_time("Testing for Positional Bias")
pos_res <- testMotifPos(
  matches, binwidth = motif_params$binwidth, abs = motif_params$abs,
  min_score = motif_params$min_score,  mc.cores = threads
)
cat_time("done\n")
gc()

cat_time("Writing", all_output$pos)
pos_res |>
  as_tibble(rownames = "altname") |>
  left_join(to_df(motif_list), by = "altname") |>
  dplyr::select(ends_with("name"), cluster, all_of(colnames(pos_res))) |>
  dplyr::select(-contains("consensus")) |>
  write_tsv(all_output$pos)

cat_time("Testing for motif enrichment")
enrich_res <- testMotifEnrich(
  motif_list[!ignore], seq$changed, seq$unchanged, model = "hyper",
  min_score = motif_params$min_score, mc.cores = threads
)
cat_time("Done")

cat_time("Writing", all_output$enrich)
enrich_res |>
  as_tibble(rownames = "altname") |>
  left_join(to_df(motif_list), by = "altname") |>
  dplyr::select(ends_with("name"), cluster, all_of(colnames(enrich_res))) |>
  dplyr::select(-contains("consensus")) |>
  write_tsv(all_output$enrich)
cat_time("Done")

