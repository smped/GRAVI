#' Create the local ZScores enabled by regioneReloaded comparing targets
#' to other targets
#'
#' This is a process requiring large amounts of RAM and multi-threading
#' The recommended number of permutations is 5000
#'
#' Required input files will be:
#'
#' 1. All union peaks
#' 2. The seqinfo object
#'
#' Required output files will be:
#'
#' 1. "output/macs2/shared/all_consensus_localz.rds")
#'
#' The enrichment params yaml can also be passed here
#'
#' It will be assumed that 5K permutations will be performed and that window
#' sizes are +/-5kb for a 10kb region
#'
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
sink(log, split = TRUE)

## Manual lists for testing. Will be overwritten by snakemake objects...
# config <- yaml::read_yaml("config/config.yml")
# all_input <- list(
#   peaks = c(
#             "output/macs2/AR/AR_consensus_peaks.bed.gz",
#             "output/macs2/ER/ER_consensus_peaks.bed.gz",
#             "output/macs2/H3K27ac/H3K27ac_consensus_peaks.bed.gz"
#           ),
#   sq = "output/annotations/seqinfo.rds"
# )
# all_output <- list(rds = "output/macs2/shared/all_consensus_localz.rds")
# all_params <- list(ntimes = 200)

config <- slot(snakemake, "config")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
cat_list(all_input, "input")
cat_list(all_output, "output")
cat_list(all_params, "params")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)


cat("Loading packages...")
library(regioneReloaded)
library(extraChIPs)
library(plyranges)
library(readr)
library(yaml)
cat("done\n")

cat("Determining the UCSC compatible reference.. ")
map <- c(
  hg19 = "hg19", hg38 = "hg38", grch37 = "hg19", grch38 = "hg38",
  mm10 = "mm10", mm39 = "mm39", grcm38 = "mm10", grcm39 = "mm39",
  rn7 = "rn7", mratbn7.2 = "rn7", galgal6 = "galGal6", rhemac10 = "rheMac10",
  canfam5 = "canFam5", susscr11 = "susScr11", pantro6 = "panTro6", dm6 = "dm6"
)
bld <- match.arg(tolower(config$genome$build), names(map))
ucsc_ref <- map[[bld]]
cat("done\n")

cat("Setting genome to be", ucsc_ref)
sq <- read_rds(all_input$sq)
genome(sq) <- ucsc_ref
cat(" done\n")

cat("Loading peaks from", paste("\n\t", all_input$peaks))
peaks <- importPeaks(all_input$peaks, seqinfo = sq, type = "bed")
names(peaks) <- gsub("_consensus.+", "", names(peaks))
cat(" done\n")

threads <- slot(snakemake, "threads")[[1]] - 1
cat("Starting at", format(Sys.time(), "%H:%M:%S, %d %b %Y"), "\n")
cat("Running multiLocalZscore with", threads, "threads...\n")
mlz_list <- names(peaks) %>%
  lapply(
    \(i) {
      multiLocalZscore(
        A = peaks[[i]],
        Blist = peaks[setdiff(names(peaks), i)],
        ranFUN = "resampleGenome",
        evFUN = "numOverlaps",
        window = 5e3,
        ntimes = all_params$ntimes,
        max_pv = 1,
        genome = ucsc_ref,
        mc.cores = threads
      )
    }
  ) %>%
  setNames(names(peaks))
cat("Finished at", format(Sys.time(), "%H:%M:%S, %d %b %Y"), "\n")


cat("Writing to", all_output$rds, "\n")
write_rds(mlz_list, all_output$rds, compress = "gz")
cat("done\n")
