#' This script defines
#'
#' - Seqinfo (sq)
#' - Chromosome sizes (chrom_sizes for bedGraphToBigWig)
#' - Transcript Models for plotting with Gviz
#' - GRanges for genes, transcripts & exons, taken directly from the gtf
#' - TSS
#' - Unique gene-centric regions
#' - Motifs
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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

cat_list(all_input, "input")
cat_list(all_output, "output")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat("Loading packages...\n")
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(glue)
library(plyranges)
library(yaml)
library(Rsamtools)
library(extraChIPs)
library(MotifDb)
library(universalmotif)
params <- read_yaml(all_input$yaml)
samples <- here::here(config$samples$file) %>%
  read_tsv()

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
write_rds(sq, all_output$sq)
cat("Seqinfo exported...\n")

### Blacklist ###
## Set all provided files into a single GRanges & export
cat("Forming unified blacklist from all files...")
bl <- config$external$blacklist %>%
  unlist() %>%
  here::here() %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  sort()
write_rds(bl, all_output$blacklist)
cat("done\n")

#### chrom_sizes ####
## For bedGraphToBigWig
sq %>%
  as_tibble() %>%
  dplyr::select(seqnames, seqlengths) %>%
  write_tsv(all_output$chrom_sizes, col_names = FALSE)
cat("chrom_sizes exported...\n")

#### GTF ####
gtf <- here::here(config$external$gtf)[[1]]
stopifnot(file.exists(gtf))
reqd_cols <- c(
  "type", "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
cat("Importing ", gtf, "\n")
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
cat("GTF imported successfully...\n")
seqlevels(all_gtf) <- seqlevels(sq)
seqinfo(all_gtf) <- sq

cat("Exporting gene, transcript and exon-level objects\n")
write_rds(all_gtf$gene, all_output$genes, compress = "gz")
write_rds(all_gtf$transcript, all_output$transcripts, compress = "gz")
write_rds(all_gtf$exon, all_output$exons, compress = "gz")
cat("All gtf_*.rds objects written successfully...\n")

#### Transcript Models (Gviz) ####
trans_models <- all_gtf$exon %>%
  select(
    type, gene = gene_id, exon = exon_id, transcript = transcript_id,
    symbol = gene_name
  )
write_rds(trans_models, all_output$trans_models, compress = "gz")
cat("trans_models.rds written successfully...\n")

#### TSS ####
tss <- all_gtf$transcript %>%
  resize(width = 1) %>%
  reduceMC() %>%
  mutate(region = "TSS") %>%
  select(region, everything(), -type)
write_rds(tss, all_output$tss, compress = "gz")
cat("TSS regions exported...\n")

#### Promoters ####
cat("Defining gene_regions...\n")
gr_params <- params$gene_regions
gene_regions <- defineRegions(
  genes = all_gtf$gene, transcripts = all_gtf$transcript, exons = all_gtf$exon,
  promoter = unlist(gr_params$promoter), upstream = gr_params$upstream,
  proximal = gr_params$intergenic
)

cat("Exporting gene_regions...\n")
write_rds(gene_regions, all_output$regions, compress = "gz")

#### Features ####
feat <- GRangesList()
seqinfo(feat) <- sq
if (!is.null(config$external$features)) {
  fl <- unlist(config$external$features)
  cat("Parsing features from", fl, "\n")
  feat <- lapply(fl, import.gff, which = GRanges(sq))
  feat <- lapply(feat, select, feature)
  feat <- unlist(GRangesList(feat))
  seqlevels(feat) <- seqlevels(sq)
  seqinfo(feat) <- sq
  cat("Finding overlap with gene regions...")
  ol <- lapply(gene_regions, \(x) propOverlap(feat, x))
  mcols(feat) <- cbind(mcols(feat), DataFrame(ol))
  cat("done\n")
  feat <- splitAsList(feat, feat$feature)
  cat("Features have split into a GRangesList of length", length(feat), "\n")
  cat("Features provided appear to be:", paste0("\n\t", names(feat)), "\n")
  cat("Writing to", all_output$features, "...")
} else {
  cat(
    "No features provided. Writing an empty object to",
    all_output$features, "..."
  )
}
write_rds(feat, all_output$features, compress = "gz")
cat("done\n")

#### Motifs ####
motif_params <- params$motif_analysis$motifdb
cat("Converting to Universal Motif format\n")
db <- convert_motifs(MotifDb) |> to_df()
if (is.null(motif_params$data_source))
  stop("No data source provided for transription factors")
db <- subset(db, dataSource %in% motif_params$data_source)
if (!is.null(motif_params$organism))
  db <- subset(db, organism %in% motif_params$organism)
cat("Database has been subset to", nrow(db), "motifs\n")

cat("Calculating correlations between motifs...")
cormat <- db |>
  to_list() |>
  compare_motifs(method = "PCC")
cormat[cormat < motif_params$cluster_above] <- 0
cormat[is.na(cormat)] <- 0
cat("done\n")
cat("Clustering motifs...\n")
cl <- as.dist(1 - cormat) %>%
  hclust() %>%
  cutree(h = 0.9)
cat(max(cl), "clusters formed\n")

cat("Exporting to:", all_output$motifs, "\n")
db |>
  mutate(cluster = cl[name]) |>
  to_list() |>
  write_rds(all_output$motifs, compress = "gz")
cat("done\n")


## Now exit confirming everything
all_exist <- map_lgl(all_output, file.exists)
if (!all(all_exist)) {
  nm <- names(all_exist)[!all_exist]
  stop("\nFailed to create:\n\t", paste(nm, collapse = "\n\t"), )
}
cat("Data export completed at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))
