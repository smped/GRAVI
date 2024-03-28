## This script takes the complete set of consensus peaks and finds the shared regions
##
## Required inputs are
##
## 1. All target-specific consensus peaks
## 2. The seqinfo object
##
## Output will be
##
## 1. output/peak_analysis/shared/shared_consensus_peaks.bed.gz
##
## First handle any conda weirdness
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

## For testing
# all_input <- list(
#   features = "output/annotations/features.rds",
#   gtf_gene = "output/annotations/gtf_gene.rds",
#   hic = "output/annotations/hic.rds",
#   regions = "output/annotations/gene_regions.rds",
#   peaks = c(
#     "../GRAVI_testing/output/peak_analysis/AR/AR_consensus_peaks.bed.gz",
#     "../GRAVI_testing/output/peak_analysis/ER/ER_consensus_peaks.bed.gz",
#     "../GRAVI_testing/output/peak_analysis/H3K27ac/H3K27ac_consensus_peaks.bed.gz"
#   ),
#   sq = "../GRAVI_testing/output/annotations/seqinfo.rds",
#   yaml = "../GRAVI_testing/config/params.yml"
# )
# all_output <- list(
#   bed = "../GRAVI_testing/output/peak_analysis/shared/shared_consensus_peaks.bed.gz",
#   rds = "../GRAVI_testing/output/peak_analysis/shared/shared_consensus_peaks.rds"
# )
# config <- yaml::read_yaml("../GRAVI_testing/config/config.yml")

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(extraChIPs)
library(plyranges)
library(tidyverse)
library(yaml)

cat_time("Loading seqinfo and defining ranges to exclude...\n")
sq <- read_rds(all_input$sq)

cat_time("Loading peaks...\n")
n_targets <- length(all_input$peaks)
shared_peaks <- all_input$peaks %>%
  importPeaks(seqinfo = sq, type = "bed") %>%
  setNames(str_remove(names(.), "_consensus.+"))
## NB: for n = 2, this will be the union peaks
shared_peaks <- shared_peaks %>%
  makeConsensus(p = 1 - 1 / n_targets, method = "coverage") %>%
  filter(n == n_targets)

cat_time("Writing", length(shared_peaks), "peaks to", all_output$bed, "\n")
write_bed(shared_peaks, all_output$bed)
cat_time("Done\n")

## Map to genes, feature & regions
cat_time("Loading all annotations")
gtf_gene <- read_rds(all_input$gtf_gene)
gene_regions <- read_rds(all_input$regions)
region_levels <- map_chr(gene_regions, \(x) x$region[1]) %>%
  setNames(names(gene_regions))
features <- read_rds(all_input$features)
hic <- read_rds(all_input$hic)
mapping_params <- all_input$yaml %>%
  read_yaml() %>%
  pluck("mapping")

## Find if there are any regions in the features which can be matched
## to promoters or enhancers
cat_time("Checking for promoters/enhancers in the features")
which_prom <- grepl("prom", str_to_lower(names(features)))
feat_prom <- features[which_prom] %>%
  unlist() %>%
  GenomicRanges::reduce()

which_enh <- grepl("enhanc", str_to_lower(names(features)))
feat_enh <- features[which_enh] %>%
  unlist() %>%
  GenomicRanges::reduce()

cat_time("Mapping peaks to regions and features")
shared_peaks$region <- bestOverlap(shared_peaks, gene_regions)
shared_peaks$region <- factor(region_levels[shared_peaks$region], unname(region_levels))
if (length(features)) 
  shared_peaks$feature <- bestOverlap(shared_peaks, features, missing = "no_feature")

save.image(here::here("output/envs/shared_peaks.RData"))

cat_time("Mapping peaks to genes")
prom <- GenomicRanges::reduce(c(feat_prom, granges(gene_regions$promoter)))
shared_peaks <- mapByFeature(
  shared_peaks, gtf_gene,
  prom = prom,
  enh = feat_enh,
  gi = hic,
  gr2gene = mapping_params$gr2gene,
  prom2gene = mapping_params$prom2gene,
  enh2gene = mapping_params$enh2gene,
  gi2gene = mapping_params$gi2gene
)
cat_time("Writing mapped peaks to", all_output$rds)
write_rds(shared_peaks, all_output$rds, compress = "gz")
cat_time("Done")



