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

## For testing
# all_input <- list(
#     peaks = c(
#         "../GRAVI_testing/output/nfr/H3K27ac/H3K27ac_E2.nfr.bed.gz",
#         "../GRAVI_testing/output/nfr/H3K27ac/H3K27ac_E2DHT.nfr.bed.gz"
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
#     bed = "../GRAVI_testing/output/nfr/H3K27ac/H3K27ac_consensus_nfr.bed.gz",
#     rds = "../GRAVI_testing/output/nfr/H3K27ac/H3K27ac_consensus_nfr.rds"
# )
# all_wildcards <- list(target = "H3K27ac")
# all_params <- list(
#   method = 'coverage',
#   min_width = 75,
#   p = 1,
#   min_gapwidth = 52
# )
# config <- yaml::read_yaml("../GRAVI_testing/config/config.yml")

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")
all_wildcards <- slot(snakemake, "wildcards")
all_params <- slot(snakemake, "params")

cat_list(all_input, "input")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(all_params, "params:", "=")
cat_list(all_output, "output:", "=")

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
vars <- c("score", "centre")
if (any(str_detect(all_input$peaks, "(bed|bed.gz)$"))) peak_type <- "bed"


cat_time("Loading peaks/ranges using type =", peak_type)
filtered_peaks <- all_input$peaks %>%
  importPeaks(
    type = peak_type, seqinfo = sq, blacklist = exclude_ranges,
    nameRanges = FALSE, centre = TRUE
  )
vars <- intersect(vars, colnames(mcols(filtered_peaks[[1]])))
if (length(vars) == 0) vars <- NULL

cat_time("Checking params")
# From  names(Rdpack::S4formals("reduce", c(x = "GenomicRanges")))
valid_args <- list(formals(makeConsensus), formals(reduceMC)) %>%
  lapply(names) %>%
  unlist() %>%
  unique() %>%
  setdiff("...")
cons_params <- list(
  ## These all need to be set on a cluster, but not when running interactively
  ## Don't know why...
  x = filtered_peaks, var = vars, simplify = FALSE, ignore.strand = TRUE,
  p = 0, method = 'union'
) %>%
  .[!names(.) %in% names(all_params)] %>%
  c(all_params) %>%
  .[names(.) %in% valid_args]
cat_time("Forming consensus peaks")
cons_peaks <- do.call("makeConsensus", cons_params)

if ("score" %in% vars) {
  cat_time("Taking the maximum score for each peak")
  cons_peaks$score <- map_dbl(cons_peaks$score, max)
}
if ("centre" %in% vars) {
  cat_time("Taking the median centre for each peak")
  cons_peaks$centre <- floor(map_dbl(cons_peaks$centre, median))
}
cons_peaks <- plyranges::select(cons_peaks, any_of(vars))

cat_time("Writing", length(cons_peaks), "ranges to", all_output$bed, "\n")
if ("score" %in% colnames(mcols(cons_peaks))) {
  cons_peaks %>%
    plyranges::select(any_of("score")) %>%
    write_bed(all_output$bed)
} else {
    write_bed(granges(cons_peaks), all_output$bed)
}
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
feat_prom <- features %>%
  endoapply(subset, grepl("prom", str_to_lower(feature))) %>%
  unlist() %>%
  GenomicRanges::reduce()
feat_enh <- features %>%
  endoapply(subset, grepl("enh", str_to_lower(feature))) %>%
  unlist() %>%
  GenomicRanges::reduce()

cat_time("Mapping peaks to regions")
cons_peaks$region <- bestOverlap(cons_peaks, gene_regions)
cons_peaks$region <- factor(
  region_levels[cons_peaks$region], unname(region_levels)
)

if (length(features)) {
  cat_time("Mapping peaks to features")
  feat_df <- features %>%
    lapply(
      \(x) bestOverlap(cons_peaks, x, var = "feature")
    ) %>%
    as_tibble() %>%
    mutate(range = as.character(cons_peaks))
  if (length(features) > 1) {
    cat_time("Merging features across source files")
    feat_df <- feat_df %>%
      nest(data = all_of(names(features))) %>%
      mutate(
        feature = lapply(
          data, \(x) {
            x <- unlist(x)
            x[!is.na(x)]
          }
        )
      ) %>%
      unnest(data, keep_empty = TRUE) %>%
      dplyr::select(feature, all_of(names(features))) %>%
      as.data.frame()
    mcols(cons_peaks) <- cbind(mcols(cons_peaks), feat_df)
    cons_peaks$feature <- CharacterList(cons_peaks$feature)
  } else {
    cons_peaks$feature <- str_replace_na(
      feat_df[[names(features)]], "no_feature"
    )
  }
}

cat_time("Mapping peaks to genes")
prom <- GenomicRanges::reduce(c(feat_prom, granges(gene_regions$promoter)))
mapping_params <- c(
  mapping_params,
  list(
    gr = cons_peaks, genes = gtf_gene, prom = prom, enh = feat_enh, gi = hic
  )
)
cons_peaks <- do.call("mapByFeature", mapping_params)

cat_time("Writing mapped peaks to", all_output$rds)
write_rds(cons_peaks, all_output$rds)
cat_time("Done")

