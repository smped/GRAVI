## This script takes the peaks created by merging samples within a target
## and treatment group, then filters down to only the peaks found in 'min_prop'
## of samples that pass QC.
##
## This way, the impact of low quality samples is minimised but having
## individual replicates adds important value to obtaining a high-quality set
## of peaks
##
## Required inputs will be
##
## 1. Paths to all individual replicates
## 2. Paths to the *merged_peaks.narrowPeak files
## 3. Path to the samples_qc.tsv output
## 4. Path to the seqinfo object
## 5. Path to the greylist(s)
##
## Required parameters will be
##
## 1. min_prop
##
## Output will be
##
## 1. output/macs2/{target}/{target}_{treat}_filtered_peaks.narrowPeak
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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

config <- slot(snakemake, "config")
all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")

cat_list(all_wildcards, "wildcards", ":")
cat_list(all_params, "params", ":")
cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(tidyverse)
library(extraChIPs)
library(plyranges)

cat("Loading seqinfo and sampes...\n")
sq <- read_rds(all_input$sq)
samples <- read_tsv(all_input$qc) %>%
    dplyr::filter(treat == all_wildcards$treat, qc == "pass")
rep_paths <- map_chr(samples$sample, \(x) str_subset(all_input$rep, x))
n_rep <- length(rep_paths)
filtered_peaks <- GRanges(seqinfo = sq)
mcnames <- c("score", "signalValue", "pValue", "qValue", "peak")
mcols <- sapply(mcnames, \(x) numeric(), simplify = FALSE)
mcols(filtered_peaks) <- DataFrame(mcols)
if (n_rep > 0) {
    cat("Loading black/grey lists\n")
    bl <- read_rds(all_input$blacklist)
    exclude_ranges <- all_input$greylist %>%
        importPeaks(seqinfo = sq, type = "bed", setNames = FALSE) %>%
        unlist() %>%
        c(bl) %>% 
        GenomicRanges::reduce() 
    cat("Loading merged peaks\n")
    merged_peaks <- all_input$merged %>%
        importPeaks(seqinfo = sq, setNames = FALSE, blacklist = exclude_ranges) %>%
        unlist()
    cat("Loading replicate peaks\n")
    rep_peaks <- importPeaks(rep_paths, seqinfo = sq)
    keep <- countOverlaps(merged_peaks, rep_peaks) > n_rep * all_params$min_prop
    filtered_peaks <- merged_peaks[keep]
}
cat("Exporting", length(filtered_peaks), "filtered_peaks\n")
write_narrowpeaks(filtered_peaks, all_output$peaks)
cat("done\n")

