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

# args <- commandArgs(TRUE)
# target <- args[[1]]
# threads <- args[[2]]
# outlier_thresh <- as.numeric(args[[3]])
# allow_zero <- args[[4]] == "True"

library(tidyverse)
library(yaml)
library(glue)
library(ngsReports)
library(GenomicRanges)
library(extraChIPs)
library(rtracklayer)
library(Rsamtools)
library(BiocParallel)
library(csaw)

target <- snakemake@wildcards[["target"]]
threads <- snakemake@threads
outlier_thresh <- as.numeric(snakemake@params[["outlier_threshold"]])
## Might need to wrangle a logical value here
allow_zero <- snakemake@params[["allow_zero"]]
register(MulticoreParam(workers = threads))

config <- read_yaml(here::here("config", "config.yml"))
samples <- here::here(config$samples$file) %>%
  read_tsv() %>%
  dplyr::filter(target == .GlobalEnv$target)
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

annotation_path <- snakemake@params[["annot_path"]]
macs2_path <- snakemake@params[["macs2_path"]]

sq <- read_rds(snakemake@input[["seqinfo"]])
blacklist <- importPeaks(
  here::here(snakemake@input[["blacklist"]]), type = "bed", seqinfo = sq
)

message("Loading peaks")
individual_peaks <- here::here(
  str_subset(snakemake@input[["indiv_macs2"]], "narrowPeak$")
) %>%
  importPeaks(seqinfo = sq, blacklist = blacklist) %>%
  GRangesList() %>%
  setNames(samples$sample)

message("Loading macs2_logs")
macs2_logs <- here::here(
  str_subset(snakemake@input[["indiv_macs2"]], "callpeak.log$")
) %>% 
  importNgsLogs() %>%
  dplyr::select(
    -contains("file"), -outputs, -n_reads, -alt_fragment_length
  ) %>%
  left_join(samples, by = c("name" = "sample")) %>%
  mutate(
    filtered_peaks = map_int(
      name,
      function(x) {
        length(individual_peaks[[x]])
      }
    ),
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
message("Writing qc_samples")
macs2_logs %>%
  dplyr::select(sample = name, any_of(colnames(samples)), qc) %>%
  write_tsv(
    here::here(snakemake@output[["qc"]])
  )

##################
## Correlations ##
##################
bfl <- bam_path %>%
  file.path(glue("{samples$sample}.bam")) %>%
  c(
    file.path(bam_path, glue("{unique(samples$input)}.bam"))
  ) %>%
  BamFileList() %>%
  setNames(c(samples$sample, unique(samples$input)))

## Check if there are any paired end reads
ys <- 1000
message("Checking for duplicates")
anyDups <- bplapply(
  bfl,
  function(x) {
    sbp <- ScanBamParam(
      flag = scanBamFlag(isDuplicate = TRUE),
      which = GRanges(sq)[which.min(seqlengths(sq))],
      what = "qname"
    )
    length(scanBam(x, param = sbp)[[1]]$qname)  > 0
  }
) %>%
  unlist()
  message("Checking for PE reads")
anyPE <- bplapply(
  bfl,
  function(x){
    yieldSize(x) <- ys
    open(x)
    flag <- scanBam(x, param=ScanBamParam(what="flag"))[[1]]$flag
    close(x)
    any(bamFlagTest(flag, "isPaired"))
  }
) %>%
  unlist()

fl <- max(macs2_logs$fragment_length)
rp <- readParam(
  pe = ifelse(any(anyPE), "both", "none"),
  dedup = any(anyDups),
  restrict = seqnames(sq)[1:5],
  discard = blacklist,
)
message("Estimating correlations")
read_corrs <- bfl[samples$sample] %>%
  path %>%
  # bplapply(correlateReads, param = rp, max.dist = 5*fl) %>%
  lapply(correlateReads, param = rp, max.dist = 5*fl) %>%
  as_tibble() %>%
  mutate(fl = seq_len(nrow(.))) %>%
  pivot_longer(
    cols = all_of(samples$sample),
    names_to = "sample",
    values_to = "correlation"
  )
write_tsv(
  read_corrs,
  here::here(snakemake@output[["cors"]])
)

## Could probably add the actual peak merging here too!!!