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
all_output <- slot(snakemake, "output")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(Rsamtools)
library(tidyverse)
library(extraChIPs)
library(rtracklayer)

#### Seqinfo won't have been created yet ####
sq <- all_input$bam %>%
  BamFileList() %>%
  seqinfo() %>%
  sortSeqlevels() %>%
  as.data.frame() %>%
  .[rownames(.) %in% paste0("chr", c(1:22, "X", "Y")),] %>%  # This covers mouse & rat
  mutate(
    isCircular = FALSE,
    genome = config$genome$build
  ) %>%
  as("Seqinfo")

#### Check the blacklist for compatibility
blacklist <- config$external$blacklist %>%
  unlist() %>%
  here::here() %>%
  importPeaks(type = "bed")
sq_has_chr <- any(grepl("chr", seqlevels(sq)))
bl_has_chr <- any(grepl("chr", seqlevels(blacklist)))
if (sq_has_chr != bl_has_chr)
  stop("Chromosome identifiers do not match between the blacklist & alignments")

#### GTF ####
gtf <- here::here(config$external$gtf)
stopifnot(file.exists(gtf))
reqd_cols <- c(
  "type", "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
cat("Importing ", gtf, "\n")
gtf_gene <- gtf %>%
  import.gff(
    which = GRanges(sq)[1], # THIS WILL FAIL if incompatible
  ) %>%
  select(all_of(reqd_cols)) %>% # WILL ABORT if any are missing
  mutate(
    gene_id = str_remove_all(gene_id, "\\..+$")
  ) %>%
  sort() %>%
  subset(seqnames %in% seqlevels(sq))
cat("GTF imported successfully...\n")
if (length(gtf_gene) == 0)
  stop(
    "No valid ranges found in the provided GTF.\nPlease check for compatible ",
    "chromosome identifiers"
  )

#### Optional RNA-Seq ####
rna_path <- here::here(config$external$rnaseq[[1]])
rnaseq <- tibble(gene_id = character())
if (length(rna_path) > 0) {
  stopifnot(file.exists(rna_path))
  if (str_detect(rna_path, "tsv$")) rnaseq <- read_tsv(rna_path)
  if (str_detect(rna_path, "csv$")) rnaseq <- read_csv(rna_path)
  if (!"gene_id" %in% colnames(rnaseq)) stop("Supplied RNA-Seq data must contain the column 'gene_id'")
  cat("RNA-Seq imported and contains gene_ids...\n")
  shared_ids <- intersect(rnaseq$gene_id, all_gtf$gene$gene_id)
  if (length(shared) == 0) stop("RNA-Seq gene ids do not mach those in the GTF")
}

## FEATURES ##
cat("Checking features...\n")
feat_files <- here::here(config$external$features) %>% unlist()
any_features <- length(feat_files > 0)
if (any_features) {
  feat_exists <- file.exists(feat_files)
  if (any(!feat_exists))
    stop("Couldn't find specified features as ", feat_files[!feat_exists], "\n")
  cat("Found feature files:\n\t", paste0(feat_files, "\n\t"))
  ## This will fail if chromosome ids are incompatible
  feat_gtf <- lapply(feat_files, import.gff, which = GRanges(sq)[1:5])
  ## But to make sure
  empty_gtf <- map_lgl(feat_gf, \(x) length(x) == 0)
  if (any(empty_gtf))
    stop("No features were found in\n\t", paste(feat_files[empty_gtf], "\n\t"))
}

## HIC ##
## Write checks for a bedpe file...

## COVERAGE ##
## Write checks for bw files

checked <- c(
  unlist(config$external$blacklist),
  unlist(config$external$gtf),
  unlist(config$external$rnaseq),
  unlist(config$external$features)
)
writeLines(checked, all_output[[1]])
