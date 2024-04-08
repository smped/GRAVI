#' Handle any conda weirdness
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
  tm <- format(Sys.time(), "%Y-%m-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}


log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(magrittr)
library(glue)
library(extraChIPs)


#### Seqinfo ####
cat_time("Importing seqinfo")
sq <- read_rds(all_input$seqinfo)

#### Check the blacklist for compatibility with sq
cat_time("Checking seqinfo compatability...")
blacklist <- config$externa$blacklist %>% importPeaks(type = "bed")
sq_has_chr <- any(grepl("chr", seqlevels(sq)))
bl_has_chr <- any(grepl("chr", seqlevels(blacklist)))
if (sq_has_chr != bl_has_chr)
  stop("Chromosome identifiers do not match between the blacklist & alignments")
cat_time("done\n")

#### Blacklist ####
## Set all provided files into a single GRanges & export
cat_time("Forming unified blacklist from all files...")
bl <- config$external$blacklist %>%
  unlist() %>%
  here::here() %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  sort()
write_rds(bl, all_output$blacklist)
cat_time("done\n")

