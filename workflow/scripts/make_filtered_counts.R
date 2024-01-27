#' Takes the section initially in differential signal/signal workflows and
#' forms the `window_counts.rds` and `filtered_counts.rds` objects for easier
#' loading into these workflows. Whilst these windows were initially formed
#' using only the two treatment groups under investigation in the specific
#' workflow, this script creates the file for **all** samples within a
#' target. This will only create any changes where >2 treatments are present,
#' however considerable time will be saved where the same files would previously
#' have been counted several times.
#' The files should be written to `output/differential_*/target` as
#' `{target}_window_counts.rds` and `{target}_filtered_counts.rds`.
#'
#' Files required for setup are:
#'   - config.yml
#'   - samples.tsv
#'   - macs2 summaries (for fragment length)
#'   - Bam Files
#'

library(tidyverse)
library(yaml)
library(csaw)
library(extraChIPs)
library(BiocParallel)
library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(ngsReports)
library(glue)

#' Arguments required to be set during calling from the workflow are
#'   - target
#'   - threads
args <- commandArgs(TRUE)
target <- args[[1]]
threads <- args[[2]]
type <- args[[3]]

#'
#' Setup the key config elements
#'
register(MulticoreParam(workers = threads))
# params <- read_yaml(here::here("config", "params.yml"))
config <- read_yaml(here::here("config", "config.yml"))

#'
#' Setup the paths
#'
macs2_path <- here::here("output", "macs2", target)
annotation_path <- here::here("output", "annotations")
stopifnot(dir.exists(annotation_path))
bam_path <- here::here(config$paths$bam)
stopifnot(dir.exists(bam_path))
out_path <- here::here("output", glue("differential_{type}"), target)
if (!dir.exists(out_path)) dir.create(out_path)
message("Files will be written to:\n", out_path)

#'
#' The sample information. Initially use all treat_levels based on the contrasts
#'
treat_levels <- config$comparisons$contrasts %>%
    unlist() %>%
    unique()
samples <- file.path(macs2_path, "qc_samples.tsv") %>%
    read_tsv() %>%
    dplyr::filter(qc == "pass") %>%
    mutate(treat = fct(treat, levels = treat_levels)) %>%
    droplevels()
treat_levels <- levels(samples$treat)
stopifnot(nrow(samples) > 0)

#'
#' Annotations and black/greylists
#'
sq <- file.path(annotation_path, "seqinfo.rds") %>%
    read_rds()
consensus_peaks <- here::here(macs2_path, "consensus_peaks.bed") %>%
    import.bed(seqinfo = sq)
blacklist <- here::here(annotation_path, "blacklist.bed.gz") %>%
    import.bed(seqinfo = sq) %>%
    sort()
greylist <- file.path(
    annotation_path, glue("{unique(samples$input)}_greylist.bed")
) %>%
    lapply(import.bed, seqinfo = sq) %>%
    GRangesList() %>%
    unlist() %>%
    reduce()

#'
#' Prepare the BamFileList and check for paired reads/duplicates
#'
bfl <- bam_path %>%
    file.path(target, glue("{samples$sample}.bam")) %>%
    c(
        file.path(bam_path, "Input", glue("{unique(samples$input)}.bam"))
    ) %>%
    BamFileList() %>%
    setNames(c(samples$sample, unique(samples$input)))
#'
#' Check for paired reads/duplicates
#'
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
    unlist() %>%
    any()
message(paste(ifelse(anyDups, "", "No"), "Duplicate Reads were found"))
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
    unlist() %>%
    any()
message(paste("Reads were", ifelse(anyPE, "paired", "single"), "end"))


#'
#' Define the parameters for reading in counts
#'
rp <- readParam(
    pe = ifelse(anyPE, "both", "none"),
    dedup = anyDups,
    restrict = seqnames(sq),
    discard = c(granges(blacklist), greylist)
)

#'
#' Set parameters for the windows using the settings in config.yml
#'
macs2_merged_logs <- file.path(
    macs2_path, glue("{treat_levels}_merged_callpeak.log")
) %>%
    importNgsLogs()
fl <- max(macs2_merged_logs$fragment_length)
# win_step <- 10*(1 + fl %/% 60) # Automatically setting things the old way
# win_size <- 3*win_step
win_step <- config$comparisons$windows$step
win_size <- config$comparisons$windows$size

#'
#' Count the windows
#'
window_counts <- windowCounts(
    bam.files = bfl,
    spacing = win_step,
    width = win_size,
    ext = fl,
    filter = length(bfl) - 1,
    param = rp,
    BPPARAM = bpparam()
)
colData(window_counts) <- colData(window_counts) %>%
    as_tibble(rownames = "sample") %>%
    dplyr::select(all_of(c("sample", "bam.files", "totals", "ext", "rlen"))) %>%
    left_join(samples, by = "sample") %>%
    mutate(
        treat = fct_explicit_na(treat, "Input"),
        target = str_replace_na(target, "Input")
    ) %>%
    as.data.frame() %>%
    column_to_rownames("sample") %>%
    DataFrame()
window_counts <- sortSeqlevels(window_counts)
seqinfo(window_counts) <- sq

#'
#' Run the filter
#'
filtered_counts <- dualFilter(
    x = window_counts[, samples$sample],
    bg = window_counts[, samples$input],
    ref = consensus_peaks,
    keep.totals = TRUE,
    q = config$comparisons$filter_q
)
colData(filtered_counts) <- droplevels(colData(filtered_counts))
ip_counts <- window_counts[, unique(samples$input)]
ip_counts <- subsetByOverlaps(ip_counts, filtered_counts)

#'
#' Export the files
#'
write_rds(
    window_counts,
    file.path(out_path, glue("{target}_window_counts.rds")),
    compress = "gz"
)
write_rds(
    filtered_counts,
    file.path(out_path, glue("{target}_filtered_counts.rds")),
    compress = "gz"
)
write_rds(
    ip_counts,
    file.path(out_path, glue("{target}_filtered_input_counts.rds")),
    compress = "gz"
)
