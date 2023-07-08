#' This script defines
#'
#' - Seqinfo (sq)
#' - Chromosome sizes (chrom_sizes for bedGraphToBigWig)
#' - Transcript Models for plotting with Gviz
#' - GRanges for genes, transcripts & exons, taken directly from the gtf
#' - TSS
#' - Unique gene-centric regions
#'
#' The trimming of upstream promoter regions takes about 90min with 16 cores
#'
#' Running this as a stand-alone script removes any dependency on config.yml
#' which reduces the number of times it is re-run by snakemake
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(glue)
library(plyranges)
library(yaml)
library(Rsamtools)
library(extraChIPs)

args <- commandArgs(TRUE)
config <- read_yaml(here::here("config", "config.yml"))
params <- read_yaml(here::here("config", "params.yml"))
samples <- here::here(config$samples$file) %>%
  read_tsv()

#### Paths ####
annotation_path <- here::here(args[[1]])
if (!dir.exists(annotation_path)) dir.create(annotation_path, recursive = TRUE)
all_out <- list(
  chrom_sizes = file.path(annotation_path, "chrom.sizes"),
  gene_regions = file.path(annotation_path, "gene_regions.rds"),
  gtf_gene = file.path(annotation_path, "gtf_gene.rds"),
  gtf_trans = file.path(annotation_path, "gtf_transcript.rds"),
  gtf_exon = file.path(annotation_path, "gtf_exon.rds"),
  seqinfo = file.path(annotation_path, "seqinfo.rds"),
  transcript_models = file.path(annotation_path, "trans_models.rds"),
  tss = file.path(annotation_path, "tss.rds")
)

#### Seqinfo ####
sq <- samples %>%
  mutate(
    path = here::here(config$paths$bam, glue("{sample}.bam"))
  ) %>%
  .[["path"]] %>%
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
write_rds(sq, all_out$seqinfo)
cat("Seqinfo exported...\n")

#### Check the blacklist for compatibility
blacklist <- here::here(config$external$blacklist) %>%
  import.bed()
sq_has_chr <- any(grepl("chr", seqlevels(sq)))
bl_has_chr <- any(grepl("chr", seqlevels(blacklist)))
if (sq_has_chr != bl_has_chr)
  stop("Chromosome identifiers do not match between the blacklist & alignments")


#### chrom_sizes ####
## For bedGraphToBigWig
sq %>%
  as_tibble() %>%
  dplyr::select(seqnames, seqlengths) %>%
  write_tsv(all_out$chrom_sizes, col_names = FALSE)
cat("chrom_sizes exported...\n")

#### GTF ####
gtf <- here::here(config$external$gtf)
stopifnot(file.exists(gtf))
reqd_cols <- c(
  "type",
  "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
cat("Importing ", gtf, "\n")
all_gtf <- gtf %>%
  import.gff(
    which = GRanges(sq), # Should fail if incompatible
    feature.type = c("gene", "transcript", "exon")
  ) %>%
  select(all_of(reqd_cols)) %>% # Will abort if any are missing
  mutate(
    gene_id = str_remove_all(gene_id, "\\..+$"),
    transcript_id = str_remove_all(transcript_id, "\\..+$"),
    exon_id = str_remove_all(exon_id, "\\..+$"),
  ) %>%
  sort() %>%
  subset(seqnames %in% seqlevels(sq)) %>%
  keepSeqlevels(seqlevels(sq)) %>%
  splitAsList(f = .$type)
cat("GTF imported successfully...\n")

if (all(vapply(all_gtf, length, integer(1)) == 0))
  stop(
    "No valid ranges found in the provided GTF.\nPlease check for compatible ",
    "chromosome identifiers"
    )
seqinfo(all_gtf) <- sq

cat("Exporting gene, transcript and exon-level objects\n")
write_rds(all_gtf$gene, all_out$gtf_gene, compress = "gz")
write_rds(all_gtf$transcript, all_out$gtf_trans, compress = "gz")
write_rds(all_gtf$exon, all_out$gtf_exon, compress = "gz")
cat("All gtf_*.rds objects written successfully...\n")

#### Transcript Models (Gviz) ####
trans_models <- all_gtf$exon %>%
  select(
    type, gene = gene_id, exon = exon_id, transcript = transcript_id, 
    symbol = gene_name
  )
write_rds(trans_models, all_out$transcript_models, compress = "gz")
cat("trans_models.rds written successfully...\n")

#### Optional RNA-Seq ####
# rna_path <- here::here(config$external$rnaseq)
# rnaseq <- tibble(gene_id = character())
# if (length(rna_path) > 0) {
#   stopifnot(file.exists(rna_path))
#   if (str_detect(rna_path, "tsv$")) rnaseq <- read_tsv(rna_path)
#   if (str_detect(rna_path, "csv$")) rnaseq <- read_csv(rna_path)
#   if (!"gene_id" %in% colnames(rnaseq)) stop("Supplied RNA-Seq data must contain the column 'gene_id'")
#   cat("RNA-Seq imported. detected will be an informative column...\n")
# }
# tx_col <- intersect(c("tx_id", "transcript_id"), colnames(rnaseq))
# rna_gr_col <- ifelse(length(tx_col) > 0, "transcript_id", "gene_id")
# rna_col <- c(tx_col, "gene_id")[[1]]

#### TSS ####
tss <- all_gtf$transcript %>%
  resize(width = 1) %>%
  reduceMC() %>%
  mutate(region = "TSS") %>%
  select(region, everything(), -type)
write_rds(tss, all_out$tss, compress = "gz")
cat("TSS regions exported...\n")

#### Promoters ####
cat("Defining gene_regions...\n")
gr_params <- params$gene_regions
gene_regions <- defineRegions(
  genes = all_gtf$gene, transcripts = all_gtf$transcript, exons = all_gtf$exon,
  promoter = unlist(gr_params$promoter), upstream = gr_params$upstream,
  proximal = gr_params$intergenic
)


write_rds(gene_regions, all_out$gene_regions, compress = "gz")
cat("Exporting gene_regions...\n")
all_exist <- vapply(all_out, file.exists, logical(1))
if (!all(all_exist)) {
  nm <- names(all_exist)[!all_exist]
  stop("\nFailed to create:\n\t", paste(nm, collapse = "\n\t"))
}
cat("Data export completed at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

