#' This script takes two sets of differential signal results and makes all
#' required pairwise outputs. These should be
#'
#' 1. Pairwise results as an rds with all initial ranges
#' 2. All pairwise changed ranges as bed files for motif testing
#'     + Increased-Increased
#'     + Increased-Decreased
#'     + Increased-Unchanged
#'     + Decreased-Increased
#'     + Decreased-Decreased
#'     + Decreased-Unchanged
#'     + Unchanged-Increased
#'     + Unchanged-Decreased
#'     + Unchanged-Unchanged (as a control set during testing)
#'
#' Mappings to genes, regions and features should replicate the methods from
#' all earlier steps
#'
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
# all_wildcards <- list(
#   tgt1 = "AR",
#   tgt2 = "ER",
#   comp1 = "E2_E2DHT",
#   comp2 = "E2_E2DHT"
# )
# full_comp <- with(all_wildcards, paste0(tgt1, "_", comp1, "-", tgt2, "_", comp2))
# bed_groups <- list(
#   c("increased", "decreased", "unchanged"), c("increased", "decreased", "unchanged")
# ) |>
#   expand.grid() |>
#   as.matrix() |>
#   apply(1, \(x) paste(x[2], x[1], sep = "_")) |>
#   paste0(".bed.gz")
# all_input <- list(
#   sq = "output/annotations/seqinfo.rds",
#   blacklist = "output/annotations/blacklist.rds",
#   features = "output/annotations/features.rds",
#   greylist = "output/greylist/greylists.rds",
#   gtf_gene = "output/annotations/gtf_gene.rds",
#   hic = "output/annotations/hic.rds",
#   regions = "output/annotations/gene_regions.rds",
#   results1 = file.path(
#     "output/differential_signal", all_wildcards$tgt1,
#     paste0(all_wildcards$tgt1, "_", all_wildcards$comp1, "-differential-signal.rds")
#   ),
#   results2 = file.path(
#     "output/differential_signal", all_wildcards$tgt2,
#     paste0(all_wildcards$tgt2, "_", all_wildcards$comp2, "-differential-signal.rds")
#   ),
#   yaml = "config/params.yml"
# )
# all_output <- list(
#   rds = file.path(
#     "output/pairwise_comparisons", full_comp, paste0(full_comp, "-pairwise_results.rds")
#   ),
#   bed = file.path(
#     "output/pairwise_comparisons", full_comp, paste(full_comp, bed_groups, sep = "-")
#   )
# )
# all_params <- list(
#   pairwise_params = list(
#     adj = "fdr",
#     alpha = 0.05
#   )
# )
# pw_params <- all_params$pairwise_params
# rm(list = c("bed_groups", "full_comp"))
# config <- yaml::read_yaml("config/config.yml")
# threads <- 2

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_wildcards <- slot(snakemake, "wildcards")
config <- slot(snakemake, "config")
threads <- slot(snakemake, "threads")
all_params <- slot(snakemake, "params")
pw_params <- all_params$pairwise_params

cat_list(all_input, "input:", "-")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(pw_params, "pairwise_params:", "=")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(extraChIPs)
library(tidyverse)
library(glue)
library(yaml)
library(scales)
library(rlang)
library(plyranges)
library(parallel)

cat_time("Setting comparison names")
both_comps <- all_wildcards[c("tgt1", "comp1", "tgt2", "comp2")] |>
  unlist() |>
  matrix(nrow = 2, byrow = TRUE) |>
  apply(1, paste, collapse = "_") |>
  (\(x) split(x, x))() |>
  setNames(paste0("comp", 1:2))
full_comp <- both_comps |>
  unlist() |>
  paste(collapse = "-")

cat_time("Loading samples")
treat_levels <- all_wildcards[c("comp1", "comp2")] %>%
  lapply(str_split_fixed, pattern = "_", n = 2) %>%
  lapply(as.character)
samples <- read_tsv(here::here(config$samples$file)) %>%
  dplyr::filter(
    (target == all_wildcards$tgt1 & treat %in% treat_levels$comp1) |
      (target == all_wildcards$tgt2 & treat %in% treat_levels$comp2)
  )


cat_time("Loading results")
grl_cols <- c("logCPM", "logFC", "PValue", "FDR", "status")
dsa_results <- list2(
  "{both_comps$comp1}" := all_input$results1,
  "{both_comps$comp2}" := all_input$results2,
) %>%
  lapply(read_rds) %>%
  mclapply(
    \(x) {
      ## Set FDR as the common column, now IHW is sorted out
      fdr_col <- ifelse(metadata(x)$ihw == "none", "FDR", "fdr_ihw")
      ## Also switch any p_mu0 column to be PValue for declaring categories
      p_col <- ifelse(metadata(x)$fc == 1, "PValue", "p_mu0")
      x %>%
        mutate(FDR = !!sym(fdr_col), PValue = !!sym(p_col)) %>%
        select(all_of(grl_cols))
    },
    mc.cores = 2
  ) %>%
  GRangesList()

cat_time("Merging results")
combined_results <- dsa_results %>%
  endoapply(mutate, centre = start + 0.5 * width) %>%
  mapGrlCols(var = c("centre", grl_cols))

cat_time("Removing any black/grey-listed regions")
bl <- read_rds(all_input$blacklist)
gl <- read_rds(all_input$greylist)[unique(samples$input)] %>% unlist()
exclude_ranges <- c(bl, gl)
combined_results <- combined_results[!overlapsAny(combined_results, exclude_ranges)]

cat_time("Updating mcols")
## Updated the status to undetected where apropriate
mc <- mcols(combined_results)
mc[str_ends(names(mc), "status")] <- mc[str_ends(names(mc), "status")] %>%
  mclapply(fct_na_value_to_level, "Undetected", mc.cores = 2)
## Distance between gr centres
mc$d <- mc[str_ends(names(mc), "centre")] %>%
  as.matrix() %>%
  rowDiffs() %>%
  abs() %>%
  as.numeric()

## Reclassify the status
lambda <- dsa_results %>%
  lapply(metadata) %>%
  lapply(pluck, "fc") %>%
  lapply(log2) %>%
  lapply(abs)
either_sig <- mc[str_ends(names(mc), "status")] %>%
  mclapply(str_detect, "(In|De)creased", mc.cores = threads) %>%
  as.data.frame() %>%
  as.matrix() %>%
  rowAnys()
cat_time("Updating status for", both_comps$comp1)
stat1_col <- paste0(both_comps$comp1, "_status")
p1_col <- paste0(both_comps$comp1, "_PValue")
lfc1_col <- paste0(both_comps$comp1, "_logFC")
mc[[stat1_col]] <- case_when(

  ## No change if significant nowhere
  !either_sig ~ mc[[stat1_col]],

  ## No change if the adjusted p (mu0) is > alpha
  p.adjust(mc[[p1_col]], pw_params$adj) >= pw_params$alpha ~ mc[[stat1_col]],

  ## The remaining sites will be significant somewhere & have a significant mu0
  ## p-value. Make significant if |lfc| > |lambda|
  mc[[lfc1_col]] > lambda[[both_comps$comp1]] ~ "Increased",
  mc[[lfc1_col]] < -lambda[[both_comps$comp1]] ~ "Decreased",

  ## Anything else is as it was
  TRUE ~ mc[[stat1_col]]

) %>%
  factor(levels = levels(dsa_results[[1]]$status))

cat_time("Updating status for", both_comps$comp2)
stat2_col <- paste0(both_comps$comp2, "_status")
p2_col <- paste0(both_comps$comp2, "_PValue")
lfc2_col <- paste0(both_comps$comp2, "_logFC")
mc[[stat2_col]] <- case_when(

  ## No change if significant nowhere
  !either_sig ~ mc[[stat2_col]],

  ## No change if the adjusted p (mu0) is > alpha
  p.adjust(mc[[p2_col]], pw_params$adj) >= pw_params$alpha ~ mc[[stat2_col]],

  ## The remaining sites will be significant somewhere & have a significant mu0
  ## p-value. Make significant if |lfc| > |lambda|
  mc[[lfc2_col]] > lambda[[both_comps$comp2]] ~ "Increased",
  mc[[lfc2_col]] < -lambda[[both_comps$comp2]] ~ "Decreased",

  ## Anything else is as it was
  TRUE ~ mc[[stat2_col]]

) %>%
  factor(levels = levels(dsa_results[[2]]$status))
mc$status <- fct_cross(
  mc[[stat1_col]], mc[[stat2_col]], sep = " - ", keep_empty = TRUE
)
mcols(combined_results) <- mc[!str_detect(names(mc), "centre")]

## There may be duplicated ranges at this point
dup_ranges <- combined_results[duplicated(combined_results)] %>%
  granges() %>%
  GenomicRanges::reduce()
cat_time("Found", length(dup_ranges), "duplicated ranges")
## Choose the most significant pairing will be chosen, with ties split by which are closest
merged_dup_ranges <- combined_results %>%
  subsetByOverlaps(dup_ranges) %>%
  arrange(1 / str_count(status, "(In|De)creased"), d) %>%
  as_tibble() %>%
  distinct(range, .keep_all = TRUE) %>%
  colToRanges("range", seqinfo = seqinfo(combined_results)) %>%
  sort()
cat_time("Removing duplicated ranges")
combined_results <- combined_results %>%
  filter_by_non_overlaps(merged_dup_ranges) %>%
  c(merged_dup_ranges) %>%
  sort()

cat_time("Re-mapping to regions")
regions <- read_rds(all_input$regions)
region_levels <- map_chr(regions, \(x) x$region[1]) %>% setNames(names(regions))
combined_results$region <- combined_results %>%
  bestOverlap(unlist(regions), var = "region") %>%
  factor(levels = region_levels)

cat_time("Checking for features")
features <- read_rds(all_input$features)
if (length(features)) {
  cat_time("Mapping peaks to features")
  feat_df <- features %>%
    lapply(
      \(x) bestOverlap(combined_results, x, var = "feature")
    ) %>%
    as_tibble() %>%
    mutate(range = as.character(combined_results))
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
    mcols(combined_results) <- cbind(mcols(combined_results), feat_df)
    combined_results$feature <- CharacterList(combined_results$feature)
  } else {
    combined_results$feature <- str_replace_na(
      feat_df[[names(features)]], "no_feature"
    )
  }
}

## Map to genes, feature & regions
cat_time("Loading all annotations")
gtf_gene <- read_rds(all_input$gtf_gene)
hic <- read_rds(all_input$hic)
mapping_params <- all_input$yaml %>%
  read_yaml() %>%
  pluck("mapping")

## Find if there are any regions in the features which can be matched
## to promoters or enhancers
cat_time("Checking for promoters/enhancers in provided features")
feat_prom <- features %>%
  endoapply(subset, grepl("(prom|tssa$)", str_to_lower(feature))) %>%
  unlist() %>%
  GenomicRanges::reduce()
cat_time("Found", length(feat_prom), "promoters in provided features")
prom <- GenomicRanges::reduce(c(feat_prom, granges(regions$promoter)))
feat_enh <- features %>%
  ## Exclude any 'weak enhancers'
  endoapply(subset, grepl("enh[^w]", str_to_lower(feature))) %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  filter_by_non_overlaps(prom)
cat_time("Found", length(feat_enh), "enhancers in provided features")

cat_time("Mapping peaks to genes")
mapping_params <- c(
  mapping_params,
  list(
    gr = combined_results,
    genes = gtf_gene, prom = prom, enh = feat_enh, gi = hic
  )
)
combined_results <- do.call("mapByFeature", mapping_params)

cat_time("Exporting complete set of results to", all_output$rds)
if (!dir.exists(dirname(all_output$rds)))
  dir.create(dirname(all_output$rds), recursive = TRUE)
write_rds(combined_results, all_output$rds, compress = "gz")

cat_time("Exporting all joint status bed files")
split_ranges <- combined_results %>%
  splitAsList(.$status) %>%
  endoapply(granges)
names(split_ranges) %>%
  str_subset("Undetected", negate = TRUE) %>%
  lapply(
    \(x) {
      tag <- str_to_lower(x) %>% str_replace_all(" - ", "_")
      bed <- str_subset(all_output$bed, tag)
      write_bed(split_ranges[[x]], bed, )
    }
  )
cat_time("done")
