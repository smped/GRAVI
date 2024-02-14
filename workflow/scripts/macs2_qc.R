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
## seqinfo, blackist, greylist, peaks, logs, bam, input_bam, gene_regions, colours
all_input <- slot(snakemake, "input")
## Required elements are: qc, cors, treat_peaks_rds, treat_peaks_bed, union_peaks
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
outlier_thresh <- as.numeric(all_params$outlier_threshold)
## Might need to wrangle a logical value here
allow_zero <- all_params[["allow_zero"]]
min_prop <- all_params[["min_prop"]] # Not supplied yet!!!

cat("Received all_input\n")
lapply(all_input, cat)

cat("Received all_params\n")
lapply(all_params, cat)

cat("Generating all_output\n")
lapply(all_output, cat)

cat("Defining samples...\n")
config <- slot(snakemake, "config")
samples <- here::here(config$samples$file) %>%
  read_tsv() %>%
  dplyr::filter(target == slot(snakemake, "wildcards")[["target"]])
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
blacklist <- here::here(all_input$blacklist) %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  unname()
greylist <- here::here(all_input$greylist) %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  unname()
exclude_ranges <- c(blacklist, greylist) %>%
  GenomicRanges::reduce()

cat("Loading peaks\n")
individual_peaks <- all_input$peaks %>%
  here::here() %>%
  importPeaks(seqinfo = sq, blacklist = exclude_ranges, centre = TRUE) %>%
  setNames(str_remove_all(names(.), "_peaks.narrowPeak"))

cat("Loading macs2_logs\n")
macs2_logs <- all_input$logs %>%
  here::here() %>%
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
  discard = blacklist,
)
cat("Estimating correlations...\n")
read_corrs <- bfl[seq_along(all_input$bam)] %>%
  path %>%
  mclapply(correlateReads, param = rp, max.dist = 5*fl, mc.cores = threads) %>%
  as_tibble() %>%
  mutate(fl = seq_len(nrow(.))) %>%
  pivot_longer(
    cols = all_of(samples$sample),
    names_to = "sample",
    values_to = "correlation"
  )
write_tsv(read_corrs, here::here(all_output$cors))

######################################
## Add the actual peak merging here ##
######################################

cat("Forming treatment_peaks\n")
n_reps <- summarise(macs2_logs, n = sum(qc == "pass"), .by = treat)
gene_regions <- read_rds(all_input$gene_regions)
regions <- map_chr(gene_regions, \(x) x$region[[1]])
has_features <- !is.null(config$external$features)
treatment_peaks <- treat_levels %>%
  lapply(
    function(x) {
      gr <- macs2_path %>%
        file.path(target, glue("{target}_{x}_merged_peaks.narrowPeak")) %>%
        importPeaks(seqinfo = sq, blacklist = exclude_ranges, centre = TRUE) %>%
        unlist()
      k <- dplyr::filter(n_reps, treat == x)$n * min_prop
      if (k > 0) {
        samp <- dplyr::filter(macs2_logs, treat == x, qc != "fail")$name
        gr$n_reps <- countOverlaps(gr, individual_peaks[samp])
        gr$keep <- gr$n_reps >= k
      } else {
        gr <- GRanges(seqinfo = sq)
      }
      gr
    }
  ) %>%
  setNames(treat_levels) %>%
  GRangesList()
## Exports here are complicated.
## 1. An rds with things how they are
## 2. Individual bed files file.path(macs2_path, target, glue("{target}_{x}_treatment_peaks.bed"))
write_rds(treatment_peaks, all_output$treat_peaks_rds)
for (i in treat_levels) {
  pattern <-  glue("{target}_{i}_treatment_peaks.bed$")
  f <- str_subset(all_output$treat_peaks_bed, pattern)
  export(treatment_peaks[[i]], f)
}
cat("All treatment peaks written to disk\n")


## Now for the union peaks
cat("Forming union_peaks\n")
union_peaks <- treatment_peaks %>%
  endoapply(subset, keep) %>%
  makeConsensus(
    var = c("score", "signalValue", "pValue", "qValue", "centre")
  ) %>%
  mutate(
    score = map_dbl(score, median),
    signalValue = map_dbl(signalValue, median),
    pValue = map_dbl(pValue, median),
    qValue = map_dbl(qValue, median),
    centre = map_dbl(centre, median),
    region = bestOverlap(
      .,
      lapply(gene_regions, select, region) %>%
        GRangesList() %>%
        unlist(),
      var = "region"
    ) %>%
      factor(levels = regions)
  )
export(union_peaks, all_output$union_peaks)
cat("Union peaks exprted\n")
cat("Done\n")
