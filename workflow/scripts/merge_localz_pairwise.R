#' This is a simple script which takes all group-wise results, discards 
#' the nulls, and returnes the remainder as a list
#' 
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

## Get all the snakemake vals
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

cat_time("Loading packages")
library(readr)
library(regioneReloaded)

cat_time("Reading all RDS files")
all_rds <- lapply(all_input$rds, read_rds)
nulls <- vapply(all_rds, is.null, logical(1))
cat_time("Retaining", sum(!nulls), "of", length(nulls), "with results")
all_rds <- all_rds[!nulls]

cat_time("Writing to", all_output$rds)
write_rds(all_rds, all_output$rds, compress = "gz")
cat_time("Done")

