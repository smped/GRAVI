#' This simply merges peaks within a given range
#' No centres are able to be retained so this will simply produce a bed file
#'
#' Required inputs: A narrowPeak file
#' Required outputs: A bed file
#' Required params: Distance to merge within
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
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

## Setup for testing
# all_input <- list(
#   peaks = "output/peak_analysis/H3K27ac/H3K27ac_E2_filtered_peaks.narrowPeak"
# )
# all_output <- list(
#   bed = "output/peak_analysis/H3K27ac/H3K27ac_E2_merged_peaks.bed"
# )
# all_params <- list(
#   within = 300
# )

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
cat_list(all_input, "input:", sep = "=")
cat_list(all_output, "output:", sep = "=")
cat_list(all_output, "params:", sep = "=")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(extraChIPs)
library(plyranges)

cat_time("Loading peaks")
peaks <- all_input$peaks |>
  importPeaks() |>
  unlist()

cat_time("Merging peaks with", all_params$within, "bp")
merged <- peaks |>
  plyranges::select(score) |>
  reduceMC(min.gapwidth = all_params$within) |>
  plyranges::mutate(score = vapply(score, max, numeric(1)))
cat_time("Reduced", length(peaks), "peaks to", length(merged), "peaks")

cat_time("Experting to", all_output$bed)
write_bed(merged, all_output$bed)
cat_time("Done")
