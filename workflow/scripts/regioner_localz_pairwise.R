#' Create the local ZScores enabled by regioneReloaded
#' This is a process requiring large amounts of RAM and multi-threading
#' The recommended number of permutations is 5000
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

# ## Manual lists for testing. Will be overwritten by snakemake objects...
# config <- yaml::read_yaml("config/config.yml")
# all_input <- list(
#   regions = "output/annotations/gene_regions.rds",
#   features = "output/annotations/features.rds",
#   bed = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-increased_unchanged.bed.gz"
# )
# all_output <- list(
#   rds = "output/pairwise_comparisons/AR_E2_E2DHT-ER_E2_E2DHT/AR_E2_E2DHT-ER_E2_E2DHT-increased_unchanged_localz.rds"
# )
# all_params <- yaml::read_yaml("config/params.yml")
# threads <- 4

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
config <- slot(snakemake, "config")
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
threads <- slot(snakemake, "threads")[[1]] - 1

## Print all input
regioner_params <- all_params$regioner
cat_list(all_input, "input:", "=")
cat_list(all_output, "output:", "=")
cat_list(regioner_params, "regioner params:", "=")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...")
library(regioneReloaded)
library(extraChIPs)
library(plyranges)
library(readr)
cat_time("done")

source(here::here("workflow/scripts/custom_functions.R"))
ucsc <- get_ucsc(config$genome$build)

cat_time("Loading regions...")
regions <- read_rds(all_input$regions)
cat_time("Loading all features...")
features <- read_rds(all_input$features) %>%
  endoapply(select, feature) %>%
  unlist() %>%
  names_to_column("source") %>%
  mutate(feature = paste(source, feature, sep = ": ")) %>%
  splitAsList(.$feature)

cat_time("Forming test_regions...")
test_regions <- c(regions, features)
test_regions <- endoapply(test_regions, granges)

cat_time("Setting genome to be", ucsc$build)
sq <- seqinfo(regions)
genome(sq) <- ucsc$build
seqinfo(test_regions) <- sq
cat_time(" done")

cat_time("Importing region set")
peaks <- all_input$bed %>% 
  importPeaks(type = 'bed', seqinfo = sq) %>% 
  unlist() %>% 
  unname() %>% 
  granges() %>% 
  unique()

if (length(peaks) < regioner_params$min_regions){

  cat_time("Too few regions found for testing. An empty object will be output")
  mlz <- NULL

} else{

  cat_time("Running multiLocalZscore with", threads, "threads...")
  mlz_params <- list(
    A = peaks, Blist = test_regions, sampling = FALSE,
    ranFUN = "resampleGenome", evFUN = "numOverlaps",
    max_pv = 1, genome = ucsc$build, mc.cores = threads
  ) %>%
    c(
      regioner_params[c("ntimes", "step", "window")]
    )
  mlz <- do.call("multiLocalZscore", mlz_params)
  cat_time("Done")

}


cat_time("Writing to", all_output$rds)
write_rds(mlz, all_output$rds, compress = "gz")
cat_time("done")
