#' Handle any conda weirdness
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
  tm <- format(Sys.time(), "%Y-%m-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}


# ## Required input for this workflow is:
# all_input <- list(
#   gene_regions = "output/annotations/gene_regions.rds",
#   ## NFRs in here somehow?
#   seqinfo = "output/annotations/seqinfo.rds"
# )
# all_output <- list(rds = "output/annotations/features.rds")
# config <- here::here("config", "config.yml") |>
#   yaml::read_yaml()

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input:")
cat_list(all_output, "output:")
cat_list(config$external, "external:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(tidyverse)
library(rtracklayer)
library(extraChIPs)
library(plyranges)

#' When parsing multiple files, separate into GTF & bed format then parse
#' separately. Bed files will use the bed file name as the feature, unless
#' named in the config yaml.
cat_time("Loading existing annotations")
sq <- read_rds(all_input$seqinfo)
gene_regions <- read_rds(all_input$gene_regions)

## File existence is already handled by the main snakemake file
cat_time("Checking files can be identified as BED or GTF")
all_files <- unlist(config$external$features, recursive = TRUE)
if (is.null(names(all_files))) {
  names(all_files) <- str_remove_all(basename(all_files), "\\.gz$")
}
all_bed <- all_files[grepl("\\.bed", all_files)]
all_gtf <- all_files[grepl("\\.gtf", all_files)]
chk <- c(
  all_files %in% c(all_bed, all_gtf), length(intersect(all_bed, all_gtf)) == 0
)
if (!all(chk)) {
  cat("Unable to identify all files as BED or GTF format")
  stop()
}

#' Checks need to be:
#' 1. [x] non-empty names
#' 2. [x] matching chromosome identifiers
#' 3. [ ] names which do not include gene regions
#' 4. [ ] Regions which don't lie outside of the chromosome boundaries as might
#' occur when a file from the wrong build is provided

bed_features <- GRangesList()
seqinfo(bed_features) <- sq
if (length(all_bed)) {
  cat_time("Loading features provided as BED files")
  bed_features <- all_bed %>%
    importPeaks(type = 'bed', nameRanges = FALSE, seqinfo = sq) %>%
    lapply(select, any_of("name"))
  names(bed_features) <- names(all_bed)
  for (i in names(bed_features)) {
    if ("name" %in% colnames(mcols(bed_features[[i]]))) {
      bed_features[[i]]$feature <- bed_features[[i]]$name
      bed_features[[i]] <- select(bed_features[[i]], feature)
    } else {
      cat("Couldn't find feature names. Setting the file name as feature for", i)
      bed_features[[i]] <- granges(bed_features[[i]])
      bed_features[[i]]$feature <- i
    }
  }
  bed_features <- GRangesList(bed_features)
}

gtf_features <- GRangesList()
seqinfo(gtf_features) <- sq
if (length(all_gtf)) {
  cat_time("Loading features provided as GTF files")
  gtf_features <- all_gtf %>%
    lapply(import.gff, which = GRanges(sq)) %>%
    lapply(select, any_of(c("name", "feature")))
  for (i in names(gtf_features)) {
    if ("feature" %in% names(mcols(gtf_features[[i]]))) {
      gtf_features[[i]] <- dplyr::select(gtf_features[[i]], feature)
    } else {
      if ("type" %in% names(mcols(gtf_features[[i]]))) {
        cat("Setting the feature to be 'type' for ", i)
        gtf_features[[i]] <- dplyr::select(gtf_features[[i]], feature = type)
      } else {
        cat("Couldn't find feature names. Setting the file name as feature for", i)
        gtf_features[[i]] <- granges(gtf_features[[i]])
        bed_fgtf_featureseatures[[i]]$feature <- i
      }
    }
  }
  gtf_features <- GRangesList(gtf_features)
  seqinfo(gtf_features) <- sq
}
all_features <- unlist(c(bed_features, gtf_features))

##########################
## Should NFRs go here? ##
##########################

if (length(all_features)) {
  cat_time("Finding overlaps with gene_regions")
  ol_df <- gene_regions %>%
    lapply(\(x) propOverlap(all_features, x)) %>%
    DataFrame()
  mcols(all_features) <- cbind(mcols(all_features), ol_df)
  all_features <- all_features %>%
    splitAsList(names(.)) %>%
    endoapply(unname)
}

if (!length(all_features)) cat("No features provided. An empty object will be written")
cat_time("Exporting features to", all_output$rds)
write_rds(all_features, all_output$rds, compress = "gz")

