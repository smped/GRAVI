#' This script runs very simple QC on the results from running macs2 callpeak
#' on individual samples. For each sample group, any sample with > X-fold, or
#' fewer than 1/X-fold peaks when copmared to the median for each sample group
#' is marked for exclusion from the
#' calling of treatment/union peaks, as well as any downstream detection of
#' differential signal. This value 'X' is specified in the main config.yml
#' as the parameter `outlier_threshold`
#'
#' Given that some ChIP targets may yield no peaks under some conditions due to
#' cytoplasmic location, the additional parameter `allow_zero` can be set to
#' `true/false` in config.yml
#'
#' Additionally, the cross-correlation between reads is calculated with a tsv
#' generated for later inclusion in differential_signal and macs2_summary
#' workflows
#'
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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

library(tidyverse)
library(yaml)
library(glue)
library(ngsReports)
library(GenomicRanges)
library(extraChIPs)
library(rtracklayer)
library(Rsamtools)
library(parallel)
library(csaw)
library(plyranges)

target <- slot(snakemake, "wildcards")[["target"]]
threads <- slot(snakemake, "threads")
## Required elements are:
## seqinfo, greylist, peaks, logs, bam, input_bam, gene_regions, colours
all_input <- slot(snakemake, "input")
## Required elements are: qc, cors, treat_peaks_rds, treat_peaks_bed, union_peaks
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
all_wildcards <- slot(snakemake, "wildcards")

## Print values to the log
cat_list(all_wildcards, "wildcards")
cat_list(all_input, "input")
cat_list(all_params, "params")
cat_list(all_output, "output")

## Stabilise paths for input/output
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

outlier_thresh <- as.numeric(all_params$outlier_threshold)
## Might need to wrangle a logical value here
allow_zero <- all_params[["allow_zero"]]

cat("Defining samples...\n")
config <- slot(snakemake, "config")
samples <- here::here(config$samples$file) %>%
  read_tsv() %>%
  dplyr::filter(target == all_wildcards$target)
treat_levels <- unique(samples$treat)
if (!is.null(config$comparisons$contrasts)) {
  ## Ensure levels respect those provided in contrasts
  treat_levels <- config$comparisons$contrasts %>%
    unlist() %>%
    intersect(samples$treat) %>%
    unique()
}
def_cols <- c("sample", "treat", "target", "input")
rep_col <- setdiff(colnames(samples), def_cols)
if (length(rep_col) > 0) {
  samples <- samples %>%
    unite(label, treat, !!sym(rep_col), remove = FALSE) %>%
    mutate(
      "{rep_col}" := as.factor(!!sym(rep_col))
    )
} else {
  samples <- samples %>%
    mutate(label = paste(treat, seq_along(treat), sep = "_"), .by = treat)
}
samples$treat <- factor(samples$treat, levels = treat_levels)

########
## QC ##
########

cat("Defining ranges to exclude...\n")
sq <- read_rds(all_input$seqinfo)
bl <- read_rds(all_input$blacklist)
exclude_ranges <- all_input$greylist %>%
  unlist() %>% 
  importPeaks(seqinfo = sq, type = "bed", setNames = FALSE) %>%
  unlist() %>%
  c(bl) %>% 
  GenomicRanges::reduce() 

cat("Loading peaks\n")
individual_peaks <- all_input$peaks %>%
  importPeaks(seqinfo = sq, blacklist = exclude_ranges, centre = TRUE) %>%
  setNames(str_remove_all(names(.), "_peaks.narrowPeak"))

cat("Loading macs2_logs\n")
macs2_logs <- all_input$logs %>%
  importNgsLogs() %>%
  dplyr::select(
    -contains("file"), -outputs, -n_reads, -alt_fragment_length
  ) %>%
  left_join(samples, by = c("name" = "sample")) %>%
  mutate(
    filtered_peaks = map_int(name, \(x) length(individual_peaks[[x]])),
    prop_passed = filtered_peaks / paired_peaks
  ) %>%
  group_by(treat) %>%
  mutate(
    med = median(filtered_peaks),
    qc = case_when(
      filtered_peaks == 0 & allow_zero ~ "pass",
      abs(log10(filtered_peaks / med)) > log10(outlier_thresh) ~ "fail",
      TRUE ~ "pass"
    ),
    label = case_when(
      qc == "fail" ~ paste(label, "(F)"),
      qc == "pass" ~ label
    )
  ) %>%
  ungroup()
## Now export for use in the merged peak calling
cat("Writing qc_samples\n")
macs2_logs %>%
  dplyr::select(sample = name, any_of(colnames(samples)), qc) %>%
  write_tsv(here::here(all_output$qc))

##################
## Correlations ##
##################
all_bam <- c(all_input$bam, all_input$input_bam)
bfl <- all_bam %>%
  BamFileList() %>%
  setNames(str_remove_all(basename(all_bam), ".bam$"))

## Check if there are any paired end reads
ys <- 1000
cat("Checking for duplicates\n")
anyDups <- mclapply(
  bfl,
  function(x) {
    sbp <- ScanBamParam(
      flag = scanBamFlag(isDuplicate = TRUE),
      which = GRanges(sq)[which.min(seqlengths(sq))],
      what = "qname"
    )
    length(scanBam(x, param = sbp)[[1]]$qname)  > 0
  }, mc.cores = threads
) %>%
  unlist()
cat("Checking for PE reads\n")
anyPE <- mclapply(
  bfl,
  function(x){
    yieldSize(x) <- ys
    open(x)
    flag <- scanBam(x, param = ScanBamParam(what="flag"))[[1]]$flag
    close(x)
    any(bamFlagTest(flag, "isPaired"))
  }, mc.cores = threads
) %>%
  unlist()

fl <- max(macs2_logs$fragment_length)
rp <- readParam(
  pe = ifelse(any(anyPE), "both", "none"),
  dedup = any(anyDups),
  restrict = seqnames(sq)[1:5],
  discard = exclude_ranges,
)
cat("Estimating correlations...\n")
read_corrs <- bfl[seq_along(all_input$bam)] %>%
  path %>%
  mclapply(
    correlateReads, param = rp, max.dist = 5*fl, mc.cores = threads
  ) %>%
  as_tibble() %>%
  mutate(fl = seq_len(nrow(.))) %>%
  pivot_longer(
    cols = all_of(samples$sample),
    names_to = "sample",
    values_to = "correlation"
  )
cat("Exporting correlations\n")
write_tsv(read_corrs, here::here(all_output$cors))
