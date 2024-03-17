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
sink(log, split = TRUE)

# all_input <- list(
#   bam = file.path(
#     "data", "bam", paste0("SRR83151", 74:79, ".bam")
#   ),
#   bai = file.path(
#     "data", "bam", paste0("SRR83151", 74:79, ".bam.bai")
#   ),
#   peaks = file.path(
#     "output", "peak_analysis", "AR",
#     paste0("AR_", c("E2", "E2DHT") ,"_filtered_peaks.narrowPeak")
#   ),
#   macs2_qc = file.path("output", "macs2", "AR", "AR_qc_samples.tsv"),
#   macs2_logs = file.path(
#     "output", "macs2", "AR", paste0(
#       "AR_", c("E2", "E2DHT"), "_merged_callpeak.log"
#     )
#   ),
#   blacklist = file.path("output","annotations", "blacklist.rds"),
#   ## Note this may be a vector of files when being run
#   greylist = file.path("output", "greylist", "SRR8315192_greylist.bed.gz"),
#   seqinfo = file.path("output", "annotations", "seqinfo.rds")
# )
# all_output <- list(
#   counts = file.path("data", "counts", "AR_counts.rds")
# )
# all_params <- list(
#   filter_q = 0.7,
#   win_size =  180, 
#   win_step = 60,
#   win_type = "sliding",
#   contrasts = list(c("E2", "E2DHT"))
# )
# config <- read_yaml(here::here("config", "config.yml"))
# all_widcards <- list(target = "AR")
# threads <- 1

config <- slot(snakemake, "config")
all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
threads <- slot(snakemake, "threads")
target <- all_wildcards$target

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
library(scales)
cat_time("Configuring for", threads, "threads")
register(MulticoreParam(workers = threads))

cat_time("Checking all parameters")
win_type <- match.arg(all_params$win_type, c("sliding", "fixed"))
if (win_type == "sliding") {
  win_step <- all_params$win_step
  ## Equality will give tiling windows...
  stopifnot(win_step > 0 & win_step <= all_params$win_size)
  stopifnot(all_params$filter_q > 0 & all_params$filter_q <= 1)
}
cat_time("Done")

#'
#' Setup the key config elements
#'
cat_time("Setting treat_levels")
treat_levels <- all_params$contrasts %>% 
  unlist() %>% 
  unique()
samples <- all_input$macs2_qc %>%
  read_tsv() %>%
  dplyr::filter(qc == "pass") %>%
  mutate(treat = fct(treat, levels = treat_levels)) %>%
  droplevels() %>%
  dplyr::select(-all_of("qc"))
treat_levels <- levels(samples$treat)
stopifnot(nrow(samples) > 0)

cat_time("Loading seqinfo")
sq <- read_rds(all_input$seqinfo)
cat_time("Defining black and greylists")
blacklist <- read_rds(all_input$blacklist)
greylist <- all_input$greylist %>% 
  unlist() %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  reduce()
exclude_gr <- reduce(c(blacklist, greylist))
cat_time("Done")

## For the peaks, loading in the treatment-level peaks will provide the
## estimated centres which will be used to recentre & resize the peaks (if reqd)
cat_time("Importing peaks")
peaks <- all_input$peaks %>%
  importPeaks(blacklist = exclude_gr, seqinfo = sq, centre = TRUE) %>%
  unlist() %>%
  select(centre) %>%
  reduceMC() %>%
  mutate(centre = vapply(centre, \(x) as.integer(mean(x)), integer(1)))

#'
#' Prepare the BamFileList and check for paired reads/duplicates
#'
#'
cat_time("Defining BamFileList")
bfl <- BamFileList(c(all_input$bam, all_input$input_bam))
ids <- path(bfl) %>% 
  basename() %>% 
  str_remove_all(".bam$")
names(bfl) <- ids
stopifnot(all(file.exists(path(bfl))))
cat_time("Done")


#'
#' Check for paired reads/duplicates
#'
cat_time("Checking for Paired Reads & duplicates")
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
cat_time(paste(ifelse(anyDups, "", "No"), "Duplicate Reads were found"))
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
cat_time(paste("Reads were", ifelse(anyPE, "paired", "single"), "end"))


#'
#' Define the parameters for reading in counts
#'
cat_time("Defining readParam")
rp <- readParam(
  pe = ifelse(anyPE, "both", "none"),
  dedup = anyDups,
  restrict = seqnames(sq),
  discard = exclude_gr
)

#'
#' Set parameters for the windows using the settings in config.yml
#'
cat_time("Estimating fragment length from macs2 logs")
macs2_merged_logs <- all_input$macs2_logs %>%
  importNgsLogs()
fl <- max(macs2_merged_logs$fragment_length)

if (win_type == "fixed") {

  cat_time("Centering & resizing peaks at ", all_params$win_size, "bp")
  fw_peaks <- peaks %>%
    mutate(
      range = paste(seqnames, centre, sep = ":"),
      union_peak = granges(.)
    ) %>%
    colToRanges("range") %>%
    resize(width = all_params$win_size, fix = "center") %>%
    filter_by_non_overlaps(exclude_gr)

  cat_time("Counting alignments...")
  se <- regionCounts(
    bfl[samples$sample], fw_peaks, ext = fl, param = rp, BPPARAM = bpparam()
  )
  seqlevels(se) <- seqlevels(sq)
  seqinfo(se) <- sq

  cat_time("Updating row/colData...")
  rowRanges(se) <- rowData(se)$union_peak
  rowData(se)$centred_peak <- granges(fw_peaks)
  colData(se) <- colData(se)[c("bam.files", "totals", "ext", "rlen")] %>%
    as_tibble(rownames = "sample") %>%
    left_join(samples, by = "sample") %>%
    as.data.frame() %>%
    DataFrame(row.names = .$sample)

  cat_time("Writing counts to ", all_output$rds)
  write_rds(se, all_output$rds, compress = "gz")
}

if (win_type == "sliding") {
  #'
  #' Count the windows
  #'
  cat_time("Counting ", all_params$win_size, "bp windows with a step of ", all_params$win_step, "bp...")
  window_counts <- windowCounts(
    bam.files = bfl,
    spacing = all_params$win_step,
    width = all_params$win_size,
    ext = fl,
    filter = length(bfl) - 1,
    param = rp,
    BPPARAM = bpparam()
  )
  cat_time("Counted", comma(nrow(window_counts)), "windows")

  cat_time("Updating colData")
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
  cat_time("Filtering counts...")
  filtered_counts <- dualFilter(
    x = window_counts[, samples$sample],
    bg = window_counts[, samples$input],
    ref = peaks,
    keep.totals = TRUE,
    q = all_params$filter_q
  )
  cat_time("Updating metadata")
  colData(filtered_counts) <- droplevels(colData(filtered_counts))
  cat_time(
    "Reduced inital ", comma(nrow(window_counts)), " windows to ", 
    comma(nrow(filtered_counts)), "filtered windows"
  )

  cat_time("Writing filtered counts to ", all_output$rds)
  write_rds(filtered_counts, all_output$rds, compress = "gz")

}
cat_time("Done counting")

