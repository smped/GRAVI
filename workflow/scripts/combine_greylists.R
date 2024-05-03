#' Combine all grey lists into a single object
#'
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

# samples <- readr::read_tsv("config/samples.tsv")
# all_input <- list(
#   gl = file.path(
#     "output", "greylist", "{unique(samples$input)}_greylist.bed.gz"
#   ) |> glue::glue(),
#   script = file.path("workflow", "scripts", "combine_greylists.R"),
#   sq = file.path("output", "annotations", "seqinfo.rds")
# )
# rm(samples)
# all_output <- list(
#   rds = file.path("output", "greylist", "greylists.rds")
# )

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(extraChIPs)

cat_time("Importing all grey-lists")
sq <- read_rds(all_input$sq)
gl <- all_input$gl %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  setNames(str_remove_all(names(.), "_greylist.bed.gz$")) %>%
  endoapply(granges)

cat_time("Exporting to", all_output$rds)
write_rds(gl, all_output$rds, compress = "gz")


