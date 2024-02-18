#' This script runs through all the external files and checks they
#'
#' 1. Exist
#' 2. Have been supplied in a matching format
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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

config <- slot(snakemake, "config")
all_external <- config$external
all_external <- lapply(all_external, unlist)
all_external <- all_external[vapply(all_external, length, integer(1)) > 0]
all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")[[1]]
cat_list(all_input, "input")
cat_list(all_external, "external files")
cat("Output will be written to ", all_output, "\n")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)
all_external <- lapply(all_external, here::here)

cat("Loading packages...\n")
library(Rsamtools)
library(tidyverse)
library(extraChIPs)
library(rtracklayer)
library(plyranges)

#### Seqinfo won't have been created yet ####
cat("Parsing genome info...")
sq <- all_input$bam %>%
  BamFileList() %>%
  seqinfo() %>%
  sortSeqlevels() %>%
  as.data.frame() %>%
  .[rownames(.) %in% paste0("chr", c(1:22, "X", "Y")),] %>%  # This covers mouse & rat
  mutate(
    isCircular = FALSE, genome = config$genome$build
  ) %>%
  as("Seqinfo")
cat("done\n")

#### Check the blacklist for compatibility
cat("Checking blacklist...")
blacklist <- all_external$blacklist %>%
  importPeaks(type = "bed")
sq_has_chr <- any(grepl("chr", seqlevels(sq)))
bl_has_chr <- any(grepl("chr", seqlevels(blacklist)))
if (sq_has_chr != bl_has_chr)
  stop("Chromosome identifiers do not match between the blacklist & alignments")
cat("done\n")

#### GTF ####
cat("Checking GTF annotations...\n")
stopifnot(length(all_external$gtf) == 1)
stopifnot(file.exists(all_external$gtf))
reqd_cols <- c(
  "type", "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
cat("Importing ", all_external$gtf, "\n")
gtf <- all_external$gtf %>%
  import.gff(
    which = GRanges(sq)[1], # THIS WILL FAIL if incompatible
  ) %>%
  select(all_of(reqd_cols)) %>% # WILL ABORT if any are missing
  mutate(gene_id = str_remove_all(gene_id, "\\..+$")) %>%
  sort() %>%
  subset(seqnames %in% seqlevels(sq))
cat("GTF imported successfully...\n")
if (length(gtf) == 0)
  stop(
    "No valid ranges found in the provided GTF.\nPlease check for compatible ",
    "chromosome identifiers"
  )
cat("done\n")

#### Optional RNA-Seq ####
if (!is.null(all_external$rnaseq)) {
  cat("Checking RNA-Seq data...")
  stopifnot(length(all_external$rnaseq) == 1)
  stopifnot(file.exists(all_external$rnaseq))
  if (str_detect(all_external$rnaseq, "tsv$"))
    rnaseq <- read_tsv(all_external$rnaseq)
  if (str_detect(all_external$rnaseq, "csv$"))
    rnaseq <- read_csv(all_external$rnaseq)
  if (!"gene_id" %in% colnames(rnaseq)) stop("Supplied RNA-Seq data must contain the column 'gene_id'")
  cat("RNA-Seq imported and contains gene_ids...\n")
  shared_ids <- intersect(rnaseq$gene_id, gtf$gene_id)
  if (length(shared_ids) == 0)
    stop("RNA-Seq gene ids do not mach those in the GTF")
  cat("done\n")
} else {
  cat("No RNA-Seq data specified\n")
}

## FEATURES ##
if (!is.null(all_external$features)) {
  cat("Checking features...")
  feat_exists <- file.exists(all_external$features)
  if (any(!feat_exists))
    stop("Couldn't find specified features as ", all_external$features[!feat_exists], "\n")
  cat("Found feature files:\n\t", paste0(all_external$features, "\n\t"))
  ## This will fail if chromosome ids are incompatible
  feat_gtf <- lapply(
    all_external$features, import.gff, which = GRanges(sq)[1:5]
  )
  has_feat_col <- map_lgl(feat_gtf, \(x) "feature" %in% colnames(mcols(x)))
  if (any(!has_feat_col))
    stop("The required column 'feature' is missing from ", all_external$features[!has_feat_col])

  region_cols <- c(
    "promoter", "upstream_promoter", "intron", "exon", "proximal_intergenic",
    "distal_intergenic", "gene_body", "intergenic"
  )
  not_permitted <- intersect(region_cols, colnames(mcols(feat_gtf)))
  if (length(not_permitted) > 0)
        stop("Disallowed feature names:", paste0("\n\t", not_permitted))

  ## But to make sure
  empty_gtf <- map_lgl(feat_gtf, \(x) length(x) == 0)
  if (any(empty_gtf))
    stop(
      "No features were found in\n\t",
      paste(all_external$features[empty_gtf], "\n\t")
    )
    cat("done\n")
} else {
  cat("No features provided\n")
}

## HIC ##
## Write checks for a bedpe file...

## COVERAGE ##
## Write checks for bw files

cat("All checks passed. Writing ", all_output[[1]], "\n")
checked <- unlist(all_external, recursive = TRUE)
writeLines(checked, all_output[[1]])
cat("Done")
