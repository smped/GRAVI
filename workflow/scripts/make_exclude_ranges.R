#' Uses AnnotationHub to find any defined regions for excluding from sequence
#' data. These are N-rich regions which can cause troubles when working with
#' PWMs. Mostly, these represent centromeres, telomeres & heterochromatin
#'
#' This needs to be run as a local rule
#'
#' The only required input is
#' - output/checks/r-packages.chk
#' - output/annotations/seqinfo.rds
#'
#' The output will be placed in
#' - output/annotations/exclude_ranges.rds
#'
#' Maybe this script should be extended later to automate downloading of the
#' blacklist using the same strategy, although this will require checking for
#' an externally provided one.
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
message("Setting stdout to ", log, "\n")
sink(log)

# all_input <- list(
#   packages = "output/checks/r-packages.chk",
#   seqinfo = "output/annotations/seqinfo.rds"
# )
#
# all_output <- list(
#   rds = "output/annotations/exclude_ranges.rds"
# )
# config <- list(genome = list(build = "GRCh37"))
# threads <- 4

config <- slot(snakemake, "config")
threads <- slot(snakemake, "threads")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
cat_list(all_input, "input -", sep = ":")
cat_list(all_output, "output -", sep = ":")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...")

library(AnnotationHub)
library(GenomicRanges)
library(readr)

cat_time("Getting UCSC build information")
source(here::here("workflow/scripts/custom_functions.R"))
ucsc <- get_ucsc(config$genome$build)
cat_time("done")
sq <- read_rds(all_input$seqinfo)
ucsc_sq <- sq
genome(ucsc_sq) <- ucsc$build

cat_time("Loading the AnnotationHub")
ah <- AnnotationHub()
cat_time("Subsetting AnnotationHub")
hubs <- ah |>
  subset(preparerclass == "excluderanges") |>
  subset(genome == ucsc$build) |>
  subset(grepl("telomere|centromere|heterochromatin", description)) |>
  names()
exclude_ranges <- GRanges(seqinfo = ucsc_sq)
if (length(hubs) == 0) {
  cat_time("No relevant excluderanges objects found")
} else {
  cat_time("Downloading", length(hubs), "sets of ranges")
  exclude_ranges <- lapply(hubs, \(x) ah[[x]]) |>
    lapply(trim) |>
    lapply(\(x) {
      seqlevels(x) <- seqnames(ucsc_sq)
      seqinfo(x) <- ucsc_sq
      mcols(x) <- mcols(x)["type"]
      x
    }) |>
    GRangesList() |>
    unlist() |>
    sort() |>
    reduce()
}
cat_time("Exporting ranges")
write_rds(exclude_ranges, all_output$rds, compress = "gz")
