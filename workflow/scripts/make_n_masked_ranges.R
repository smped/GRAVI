#' Finds all Ns in the reference genome & creates a GRanges object for masking 
#' during motif enrichment analysis
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

if ("snakemake" %in% ls()){
  config <- snakemake@config
  all_input <- snakemake@input
  all_output <- snakemake@output
} else{
  os.path.join <- file.path
  config <- yaml::read_yaml(here::here("config.yaml"))
  all_input <- list(
    seqinfo = os.path.join("output", "annotations", "seqinfo.rds"),
    script = os.path.join("workflow", "scripts", "get_ucsc.R")
  )
  all_output <- list(
    rds = os.path.join("output", "annotations", "n_masked_ranges.rds")
  )
}

cat_list(all_input, "input", sep = ":")
cat_list(all_output, "output", sep = ":")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...")
library(tidyverse)
library(GenomicRanges)
library(Biostrings)

## Contains the function for finding UCSC build info
cat_time("Loading UCSC build info function")
source(all_input$script)
ucsc <- get_ucsc(config$genome$build)
pkg <- paste(c("BSgenome", ucsc$sp, "UCSC", ucsc$build), collapse = ".")
cat_time("Loading", pkg)
library(pkg, character.only = TRUE)
bs_genome <- get(pkg)

cat_time("Loading seqinfo")
sq <- read_rds(all_input$seqinfo)

cat_time("Identifying Ns for masking")
n_mask <- vmatchPattern("N", bs_genome) %>%
  GenomicRanges::reduce() |>
  subset(seqnames %in% seqlevels(sq))
seqlevels(n_mask) <- seqlevels(sq)
seqinfo(n_mask) <- sq

cat_time("Saving n_masked_ranges")
write_rds(n_mask, all_output$rds)
cat_time("Done")