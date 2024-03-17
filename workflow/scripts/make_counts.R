#' Takes the section initially in differential signal workflows and
#' forms the `counts.rds` and `filtered_counts.rds` objects for easier
#' loading into these workflows. Whilst these windows were initially formed
#' using only the two treatment groups under investigation in the specific
#' workflow, this script creates the file for **all** samples within a
#' target. This will only create any changes where >2 treatments are present,
#' however considerable time will be saved where the same files would previously
#' have been counted several times.
#' 
#' Required input files are:
#'   - All bam files for the target
#'   - All bam indexes
#'   - macs2_qc (to determine valid samples)
#'   - macs2 logs (fragment length)
#'   - All consensus peaks
#'   - Checks (packages & here)
#'   - Black & grey lists
#'   - All annotations (mainly seqinfo)
#'
#' Required output files are:
#'   - The counts (put somewhere. Maybe in the peak_analysis folder?)
#' 
#' Required parameters are
#'   - win_type
#'   - win_size
#'   - win_step
#'   - filter_q
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

## A function for printing snakemake input
cat_list <- function(x, slot = NULL, sep = "\n\t"){
    nm <- setdiff(names(x), "")
    invisible(
        lapply(
            nm,
            \(i) cat("Received", slot, i, sep, paste0( x[[i]], "\n\t"), "\n")
        )
    )
}

## Time stamped messages
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log)

all_input <- list(
  bam = file.path(
    "data", "bam", paste0("SRR83151", 74:79, ".bam")
  ),
  bai = file.path(
    "data", "bam", paste0("SRR83151", 74:79, ".bam.bai")
  ),
  peaks = file.path("output", "peak_analysis", "AR_consensus_peaks.bed.gz"),
  macs2_qc = file.path("output", "macs2", "AR", "AR_qc_samples.tsv"),
  macs2_logs = file.path(
    "output", "macs2", "AR", paste0(
      "AR_", c("E2", "E2DHT"), "_merged_callpeak.log"
    )
  )
  blacklist = file.path("output","annotations", "blacklist.rds")
  ## Note this may be a vector of files when being run
  greylist = file.path("output", "annotations", "SRR8315192_greylist.bed.gz"),
  seqinfo = file.path("output", "annotations", "seqinfo.rds")
)
all_output <- list(
  counts = file.path("output", "peak_analysis", "AR", "AR_counts.rds")
)
all_params <- list(
  filter_q = 0,
  win_size =  400, # This should become a global parameter for all
  win_step = 0,
  win_type = "fixed"
)



config <- slot(snakemake, "config")
# all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")
# all_input <- slot(snakemake, "input")
# all_output <- slot(snakemake, "output")

cat_list(all_wildcards, "wildcards", ":")
cat_list(all_params, "params", ":")
cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(yaml)
library(extraChIPs)
library(glue)
library(GenomicRanges)
library(Rsamtools)
library(csaw)
library(ngsReports)
library(BiocParallel)
library(plyranges)

#' Arguments required to be set during calling from the workflow are
#'   - threads
#'   - target
#'   - out_file
#'   - window_type
#'   - window_size
#'   - step (sliding only), use a value of 0 for fixed-with
#'   - filter_q
# args <- c("4", "H3K27ac", "data/counts/H3K27ac_counts.rds","sliding", "240", "80", "0.7")
args <- commandArgs(TRUE)
cat("\nReceived args:", args, "", sep = "\n")
threads <- as.integer(args[[1]])
target <- args[[2]]
counts_file <- here::here(args[[3]])
win_type <- match.arg(args[[4]], c("fixed", "sliding"))
win_size <- as.integer(args[[5]])
win_step <- as.integer(args[[6]])
filter_q <- as.numeric(args[[7]])

#'
#' Setup the key config elements
#'
register(MulticoreParam(workers = threads))
config <- read_yaml(here::here("config", "config.yml"))

#'
#' Setup the paths
#'
macs2_path <- here::here("output", "macs2", target)
annotation_path <- here::here("output", "annotations")
stopifnot(dir.exists(annotation_path))
bam_path <- here::here(config$paths$bam)
stopifnot(dir.exists(bam_path))
out_path <- dirname(counts_file)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
message("Counts will be written to:\n", counts_file)

#'
#' The sample information. Initially use all treat_levels based on the contrasts
#'
treat_levels <- config$comparisons$contrasts %>%
  unlist() %>%
  unique()
samples <- file.path(macs2_path, glue("{target}_qc_samples.tsv")) %>%
  read_tsv() %>%
  dplyr::filter(qc == "pass") %>%
  mutate(treat = fct(treat, levels = treat_levels)) %>%
  droplevels() %>%
  dplyr::select(-all_of("qc"))
treat_levels <- levels(samples$treat)
stopifnot(nrow(samples) > 0)

#'
#' Annotations and black/greylists
#'
sq <- file.path(annotation_path, "seqinfo.rds") %>%
  read_rds()
exclude_gr <- c(
  here::here(config$external$blacklist),
  file.path(
    annotation_path, glue("{unique(samples$input)}_greylist.bed")
  )
) %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  reduce()
peaks <- here::here(macs2_path, paste0(target, "_treatment_peaks.rds")) %>%
  read_rds() %>%
  endoapply(select, centre) %>%
  unlist() %>%
  reduceMC() %>%
  mutate(centre = vapply(centre, \(x) as.integer(mean(x)), integer(1)))

#'
#' Prepare the BamFileList and check for paired reads/duplicates
#'
#'
ids <- c(samples$sample, unique(samples$input))
bfl <- BamFileList(file.path(bam_path, paste0(ids, ".bam")))
names(bfl) <- ids
stopifnot(all(file.exists(path(bfl))))

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
  discard = exclude_gr
)

#'
#' Set parameters for the windows using the settings in config.yml
#'
macs2_merged_logs <- file.path(
  macs2_path, glue("{target}_{treat_levels}_merged_callpeak.log")
) %>%
  importNgsLogs()
fl <- max(macs2_merged_logs$fragment_length)

if (win_type == "fixed") {
  message("Centering & resizing peaks at ", win_size, "bp")
  fw_peaks <- peaks %>%
    mutate(
      range = paste(seqnames, centre, sep = ":"),
      union_peak = granges(.)
    ) %>%
    colToRanges("range") %>%
    resize(width = win_size, fix = "center") %>%
    filter_by_non_overlaps(exclude_gr)
  message("Counting alignments...")
  se <- regionCounts(
    bfl[samples$sample], fw_peaks, ext = fl, param = rp, BPPARAM = bpparam()
  )
  seqlevels(se) <- seqlevels(sq)
  seqinfo(se) <- sq
  message("Updating row/colData...")
  rowRanges(se) <- rowData(se)$union_peak
  rowData(se)$centred_peak <- granges(fw_peaks)
  colData(se) <- colData(se)[c("bam.files", "totals", "ext", "rlen")] %>%
    as_tibble(rownames = "sample") %>%
    left_join(samples, by = "sample") %>%
    as.data.frame() %>%
    DataFrame(row.names = .$sample)
  message("Writing counts to ", counts_file)
  write_rds(se, counts_file, compress = "gz")
}

if (win_type == "sliding") {
  #'
  #' Count the windows
  #'
  message("Counting ", win_size, "bp windows with a step of ", win_step, "bp...")
  window_counts <- windowCounts(
    bam.files = bfl,
    spacing = win_step,
    width = win_size,
    ext = fl,
    filter = length(bfl) - 1,
    param = rp,
    BPPARAM = bpparam()
  )
  message("Updating colData")
  colData(window_counts) <- colData(window_counts) %>%
    as_tibble(rownames = "sample") %>%
    dplyr::select(all_of(c("sample", "bam.files", "totals", "ext", "rlen"))) %>%
    left_join(samples, by = "sample") %>%
    mutate(
      treat = fct_na_value_to_level(treat, "Input"),
      target = str_replace_na(target, "Input")
    ) %>%
    as.data.frame() %>%
    DataFrame(row.names = .$sample)
  window_counts <- sortSeqlevels(window_counts)
  seqinfo(window_counts) <- sq

  #'
  #' Run the filter
  #'
  message("Filtering counts...")
  filtered_counts <- dualFilter(
    x = window_counts[, samples$sample],
    bg = window_counts[, samples$input],
    ref = peaks,
    keep.totals = TRUE,
    q = filter_q
  )
  colData(filtered_counts) <- droplevels(colData(filtered_counts))
  message("Reduced inital ", nrow(window_counts), " windows to ", nrow(filtered_counts))

  message("Writing filtered counts to ", counts_file)
  write_rds(filtered_counts, counts_file, compress = "gz")

}

