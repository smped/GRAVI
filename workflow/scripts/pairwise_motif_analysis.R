#' Perform motif analysis on all pairwise changed regions
#' Comparisons should be against unchanged-unchanged, using the hypergeometric
#' A minimum size for sites will be required here.
#'
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

# all_wildcards <- list(
#   tgt1 = "AR",
#   comp1 = "E2_E2DHT",
#   tgt2 = "GATA3",
#   comp2 = "E2_E2DHT"
# )
# full_comp <- with(all_wildcards, paste0(tgt1, "_", comp1, "-", tgt2, "_", comp2))
# all_input <- list(
#   rds = file.path(
#     "output", "pairwise_comparisons", full_comp,
#     paste0(full_comp, "-pairwise-results.rds")
#   ),
#   motifs = "output/annotations/motif_list.rds",
#   seqinfo = "output/annotations/seqinfo.rds"
# )
# all_output <- list(
#   enrich_tsv = file.path(
#     "output/pairwise_comparisons", full_comp,
#     paste0(full_comp, "-motif_enrichment.tsv.gz")
#   ),
#   position_tsv = file.path(
#     "output/pairwise_comparisons", full_comp,
#     paste0(full_comp, "-motif_position.tsv.gz")
#   )
# )
# all_params <- list(
#   motif_params = jsonlite::fromJSON("config/json/motif_analysis_param.json")$pairwise
# )
# threads <- 4
# config <- yaml::read_yaml("config/config.yml")
# rm(list = c("full_comp", "bed_groups"))

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_wildcards <- slot(snakemake, "wildcards")
config <- slot(snakemake, "config")
threads <- slot(snakemake, "threads")
all_params <- slot(snakemake, "params")

motif_params <- all_params$motif_params
full_comp <- with(all_wildcards, paste0(tgt1, "_", comp1, "-", tgt2, "_", comp2))

cat_list(all_input, "input:")
cat_list(all_output, "output:")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(motif_params, "motif_params:", "=")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(tidyverse)
library(motifTestR)
library(extraChIPs)
library(plyranges)
library(universalmotif)
library(scales)
library(parallel)

## Contains the function for finding UCSC build info
source(here::here("workflow/scripts/get_ucsc.R"))
ucsc <- get_ucsc(config$genome$build)
pkg <- paste(c("BSgenome", ucsc$sp, "UCSC", ucsc$build), collapse = ".")
cat_time("Loading", pkg)
library(pkg, character.only = TRUE)
bs_genome <- get(pkg)
sq <- read_rds(all_input$seqinfo)
genome(sq) <- ucsc$build

cat_time("Loading motifs")
motif_df <- read_rds(all_input$motifs) %>% to_df() %>% as_tibble()
motif_list <- to_list(motif_df)

cat_time("Loading pairwise results")
pairwise_results <- read_rds(all_input$rds) %>%
  colToRanges("centre") %>%
  resize(width = motif_params$peak_width, fix = 'center') %>%
  splitAsList(.$status) %>%
  endoapply(unique) %>%
  endoapply(granges) %>%
  .[map_int(., length) > 0]
genome(pairwise_results) <- ucsc$build
cat_time("Parsed", length(pairwise_results), "sets of ranges")

cat_time("Getting sequences...")
seq_list <- pairwise_results %>%
  lapply(\(x) setNames(x, as.character(x))) %>%
  mclapply(\(x) getSeq(bs_genome, x), mc.cores = threads) %>%
  mclapply(\(x) x[letterFrequency(x, "N")[,1] == 0], mc.cores = threads) %>%
  as("DNAStringSetList")
n_seq <- map_int(seq_list, length)
cat_time("Extracted", sum(n_seq), "sequences for testing")

cat_time("Checking for low frequency matches")
min_matches <- motif_params$ignore_below * sum(n_seq)
counts <- countPwmMatches(
  motif_list, as(unlist(unlist(seq_list)), "DNAStringSet"),
  min_score = motif_params$min_score, mc.cores = threads
)
ignore <- counts < min_matches
cat_time(
  "Found", sum(ignore), "motifs with matches in <",
  percent(motif_params$ignore_below), "of sequences"
)

cat_time("Started getting best matches")
matches <- seq_list %>%
  lapply(
    \(x) {
      getPwmMatches(
        motif_list[!ignore], x, best_only = TRUE,
        min_score = motif_params$min_score, break_ties = motif_params$break_ties,
        mc.cores = threads
      )
    }
  ) %>% 
  lapply(
    \(x) x[vapply(x, nrow, integer(1)) > 0]
  )
cat_time("done")
gc()


cat_time("Testing for Positional Bias")
pos_res <- matches %>%
  lapply(
    testMotifPos, binwidth = motif_params$binwidth, abs = motif_params$abs,
    min_score = motif_params$min_score, mc.cores = threads
  )
cat_time("done\n")
gc()

cat_time("Writing", all_output$position_tsv)
pos_res %>%
  lapply(as_tibble, rownames = "name") %>%
  lapply(dplyr::select, -contains("consensus")) %>%
  bind_rows(.id = "comparison") %>%
  left_join(motif_df, by = "name") %>%
  dplyr::select(
    comparison, ends_with("name"), cluster, any_of(colnames(pos_res[[1]]))
  ) %>%
  write_tsv(all_output$position_tsv)

#' Enrichment can be against all other sequences, or against unch-unch
#' Currently set to all other groups, making unchanged-unchanged a result
cat_time("Testing for motif enrichment in", length(seq_list), "non-zero sets of sequences")
cat_list(list(groups = names(seq_list)), "non-zero")
enrich_res <-  lapply(
    names(seq_list),
    \(x) {
      cat_time("testing", x)
      bg <- unlist(seq_list[names(seq_list) != x])
      bg <- unique(bg)
      testMotifEnrich(
        motif_list[!ignore], seq_list[[x]], bg,
        model = "hypergeometric", min_score = motif_params$min_score,
        mc.cores = threads
      )
    }
  )
names(enrich_res) <- names(seq_list)
cat_time("Done")

cat_time("Writing", all_output$enrich_tsv)
enrich_res %>%
  lapply(as_tibble, rownames = "name") %>%
  bind_rows(.id = "comparison") %>%
  left_join(motif_df, by = "name") %>%
  dplyr::select(
    comparison, ends_with("name"), cluster, all_of(colnames(enrich_res[[1]]))
  ) %>%
  write_tsv(all_output$enrich_tsv)
cat_time("Done")
