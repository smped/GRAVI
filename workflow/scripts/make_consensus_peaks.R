## This script takes the peaks retained after filtering merged peaks and forms
## a set of consensus peaks, as the union of all ranges covered by a peak
##
## Required inputs are
##
## 1. output/peak_analysis/{target}/{target}_{treat}_filtered_peaks.narrowPeak
## 2. The greylist
## 3. The seqinfo object
##
##
## Output will be
##
## 1. output/mapeak_analysiscs2/{target}/{target}_consensus_peaks.bed
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
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

## For testing
# all_input <- list(
#     peaks = c(
#         "../GRAVI_testing/output/peak_analysis/AR/AR_E2DHT_filtered_peaks.narrowPeak",
#         "../GRAVI_testing/output/peak_analysis/AR/AR_E2_filtered_peaks.narrowPeak"
#     ),
#     sq = "../GRAVI_testing/output/annotations/seqinfo.rds",
#     blacklist = "../GRAVI_testing/output/annotations/blacklist.rds",
#     features = "../GRAVI_testing/output/annotations/features.rds",
#     greylist = "../GRAVI_testing/output/greylist/SRR8315192_greylist.bed.gz",
#     gtf_gene = "../GRAVI_testing/output/annotations/gtf_gene.rds",
#     hic = "../GRAVI_testing/output/annotations/hic.rds",
#     regions = "../GRAVI_testing/output/annotations/gene_regions.rds",
#     yaml = "../GRAVI_testing/config/params.yml"
# )
# all_output <- list(
#     peaks = "../GRAVI_testing/output/peak_analysis/AR/AR_consensus_peaks.bed.gz",
#     rds = "../GRAVI_testing/output/peak_analysis/AR/AR_consensus_peaks.rds"
# )
# all_wildcards <- list(target = "AR")
# config <- yaml::read_yaml("../GRAVI_testing/config/config.yml")

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")
all_wildcards <- slot(snakemake, "wildcards")

cat_list(all_input, "input")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(extraChIPs)
library(plyranges)
library(yaml)

cat_time("Loading seqinfo and defining ranges to exclude...\n")
sq <- read_rds(all_input$sq)
bl <- read_rds(all_input$blacklist)
exclude_ranges <- all_input$greylist %>%
    unlist() %>%
    importPeaks(seqinfo = sq, type = "bed", setNames = FALSE) %>%
    unlist() %>%
    c(bl) %>%
    GenomicRanges::reduce()

cat_time("Checking peak type")
peak_type <- "narrow"
var <- c("score", "centre")
if (any(str_detect(all_input$peaks, "(bed|bed.gz)$"))) {
  peak_type <- "bed"
  vars <- "score"
}

cat_time("Loading peaks/ranges using type =", peak_type)
filtered_peaks <- all_input$peaks %>%
  importPeaks(
    type = peak_type, seqinfo = sq, blacklist = exclude_ranges,
    nameRanges = FALSE, centre = TRUE
  )
cons_peaks <- filtered_peaks %>% makeConsensus(var = vars)

if ("score" %in% vars) {
  cat_time("Taking the maximum score for each peak")
  cons_peaks$score <- map_dbl(cons_peaks$score, max)
}
if ("centre" %in% vars) {
  cat_time("Taking the median centre for each peak")
  cons_peaks$centre <- floor(map_int(cons_peaks$centre, median))
}
cons_peaks <- plyranges::select(cons_peaks, any_of(vars))

cat_time("Writing", length(cons_peaks), "ranges to", all_output$peaks, "\n")
cons_peaks %>%
  plyranges::select(any_of("score")) %>%
  write_bed(all_output$peaks)
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
cons_peaks$region <- bestOverlap(cons_peaks, gene_regions)
cons_peaks$region <- factor(region_levels[cons_peaks$region], unname(region_levels))
if (length(features))
  cons_peaks$feature <- bestOverlap(cons_peaks, features, missing = "no_feature")

cat_time("Mapping peaks to genes")
prom <- GenomicRanges::reduce(c(feat_prom, granges(gene_regions$promoter)))
cons_peaks <- mapByFeature(
  cons_peaks, gtf_gene,
  prom = prom,
  enh = feat_enh,
  gi = hic,
  gr2gene = mapping_params$gr2gene,
  prom2gene = mapping_params$prom2gene,
  enh2gene = mapping_params$enh2gene,
  gi2gene = mapping_params$gi2gene
)
cat_time("Writing mapped peaks to", all_output$rds)
write_rds(cons_peaks, all_output$rds)
cat_time("Done")

