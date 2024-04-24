#' This script defines
#'
#' - Seqinfo (sq)
#' - Chromosome sizes (chrom_sizes for bedGraphToBigWig)
#' - Transcript Models for plotting with Gviz
#' - GRanges for genes, transcripts & exons, taken directly from the gtf
#' - TSS
#' - Unique gene-centric regions
#'
#' Running this as a stand-alone script removes any dependency on config.yml
#' which reduces the number of times it is re-run by snakemake
#'
#'
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


log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input:")
cat_list(all_output, "output:")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(glue)
library(plyranges)
library(yaml)
library(Rsamtools)
library(extraChIPs)
params <- read_yaml(all_input$yaml)

#### Seqinfo ####
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
write_rds(sq, all_output$seqinfo)
cat_time("Seqinfo exported...\n")

#### chrom_sizes ####
## For bedGraphToBigWig
sq %>%
  as_tibble() %>%
  dplyr::select(seqnames, seqlengths) %>%
  write_tsv(all_output$chrom_sizes, col_names = FALSE)
cat_time("chrom_sizes exported...\n")

#### GTF ####
## Perhaps set to this to also take a named list of bed files?
gtf <- here::here(config$external$gtf)[[1]]
stopifnot(file.exists(gtf))
reqd_cols <- c(
  "type", "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
cat_time("Importing ", gtf, "\n")
all_gtf <- gtf %>%
  import.gff(
    which = GRanges(sq), # THIS WILL FAIL if incompatible
    feature.type = c("gene", "transcript", "exon")
  ) %>%
  select(all_of(reqd_cols)) %>% # WILL ABORT if any are missing
  mutate(
    gene_id = str_remove_all(gene_id, "\\..+$"),
    transcript_id = str_remove_all(transcript_id, "\\..+$"),
    exon_id = str_remove_all(exon_id, "\\..+$"),
  ) %>%
  sort() %>%
  subset(seqnames %in% seqlevels(sq)) %>%
  splitAsList(f = .$type)
gtf_lens <- map_int(all_gtf[c("gene", "transcript", "exon")], length)
if (any(gtf_lens == 0)) {
	cat("Zero length categories found for", names(all_gtf)[gtf_lens == 0])
	stop()
}
cat_time("GTF imported successfully...\n")
seqlevels(all_gtf) <- seqlevels(sq)
seqinfo(all_gtf) <- sq

cat_time("Exporting gene, transcript and exon-level objects\n")
write_rds(all_gtf$gene, all_output$gtf_gene, compress = "gz")
write_rds(all_gtf$transcript, all_output$gtf_transcript, compress = "gz")
write_rds(all_gtf$exon, all_output$gtf_exon, compress = "gz")
cat_time("All gtf_*.rds objects written successfully...\n")

#### Transcript Models (Gviz) ####
trans_models <- all_gtf$exon %>%
  select(
    type, gene = gene_id, exon = exon_id, transcript = transcript_id,
    symbol = gene_name
  )
write_rds(trans_models, all_output$trans_models, compress = "gz")
cat_time("trans_models.rds written successfully...\n")

#### TSS ####
tss <- all_gtf$transcript %>%
  resize(width = 1) %>%
  reduceMC() %>%
  mutate(region = "TSS") %>%
  select(region, everything(), -type)
write_rds(tss, all_output$tss, compress = "gz")
cat_time("TSS regions exported...\n")

#### Promoters ####
cat_time("Defining gene_regions...\n")
region_params <- params$gene_regions %>%
  lapply(unlist) %>%
  lapply(unname) %>%
  c(
    all_gtf[c("gene", "transcript", "exon")] %>%
      as.list() %>%
      setNames(paste0(names(.), "s"))
  )
gene_regions <- do.call("defineRegions", region_params)

cat_time("Exporting gene_regions...\n")
write_rds(gene_regions, all_output$regions, compress = "gz")


