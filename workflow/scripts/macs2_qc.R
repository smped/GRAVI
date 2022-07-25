#' This script runs very simple QC on the results from running macs2 callpeak
#' on individual samples. For each sample group, any sample with > X-fold, or
#' fewer than 1/X-fold peaks when copmared to the median for each sample group
#' is marked for exclusion from the
#' calling of oracle/consensus peaks, as well as any downstream detection of
#' differential binding. This value 'X' is specified in the main config.yml
#' as the parameter `outlier_threshold`
#'
#' Given that some ChIP targets may yield no peaks under some conditions due to
#' cytoplasmic location, the additional parameter `allow_zero` can be set to
#' `true/false` in config.yml
#'
#' Additionally, the cross-correlation between reads is calculated with a tsv
#' generated for later inclusion in differential_binding and macs2_summary
#' workflows

if (!"extraChIPs" %in% rownames(installed.packages()))
  BiocManager::install("steveped/extraChIPs", ask = FALSE)

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

args <- commandArgs(TRUE)
target <- args[[1]]
threads <- args[[2]]
register(MulticoreParam(workers = threads))

config <- read_yaml(
  here::here("config", "config.yml")
)
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
rep_col <- setdiff(colnames(samples), c("sample", "treat", "target", "input"))
samples <- samples %>%
  unite(label, treat, !!sym(rep_col), remove = FALSE) %>%
  mutate(
    treat = factor(treat, levels = treat_levels),
    "{rep_col}" := as.factor(!!sym(rep_col))
  )

########
## QC ##
########

outlier_thresh <- config$peaks$qc$outlier_threshold
allow_zero <- as.logical(config$peaks$qc$allow_zero)

annotation_path <- here::here("output", "annotations")
macs2_path <- here::here("output", "macs2", target)
if (!dir.exists(macs2_path)) dir.create(macs2_path)

sq <- read_rds(file.path(annotation_path, "seqinfo.rds"))
blacklist <-  file.path(annotation_path, "blacklist.bed.gz") %>%
  import.bed(seqinfo = sq) %>%
  sort()

individual_peaks <- file.path(
  macs2_path, glue("{samples$sample}_peaks.narrowPeak")
) %>%
  importPeaks(seqinfo = sq, blacklist = blacklist) %>%
  GRangesList() %>%
  setNames(samples$sample)

macs2_logs <- file.path(macs2_path, glue("{samples$sample}_callpeak.log")) %>%
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
macs2_logs %>%
  dplyr::select(sample = name, any_of(colnames(samples)), qc) %>%
  write_tsv(
    file.path(macs2_path, "qc_samples.tsv")
  )

##################
## Correlations ##
##################

bam_path <- here::here(config$paths$bam)
stopifnot(dir.exists(bam_path))
bfl <- bam_path %>%
  file.path(target, glue("{samples$sample}.bam")) %>%
  c(
    file.path(bam_path, "Input", glue("{unique(samples$input)}.bam"))
  ) %>%
  BamFileList() %>%
  setNames(c(samples$sample, unique(samples$input)))

## Check if there are any paired end reads
ys <- 1000
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
read_corrs <- bfl[samples$sample] %>%
  path %>%
  bplapply(correlateReads, param = rp, max.dist = 5*fl) %>%
  as_tibble() %>%
  mutate(fl = seq_len(nrow(.))) %>%
  pivot_longer(
    cols = all_of(samples$sample),
    names_to = "sample",
    values_to = "correlation"
  )
write_tsv(read_corrs, file.path(macs2_path, "cross_correlations.tsv"))
