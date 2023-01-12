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
library(BiocParallel)
library(extraChIPs)

args <- commandArgs(TRUE)
gtf <- args[[1]]
# threads <- args[[2]]
# register(MulticoreParam(threads))

config <- read_yaml(here::here("config", "config.yml"))
params <- read_yaml(here::here("config", "params.yml"))
samples <- here::here(config$samples$file) %>%
  read_tsv()

#### Paths ####
annotation_path <- here::here("output", "annotations")
if (!dir.exists(annotation_path)) dir.create(annotation_path, recursive = TRUE)
all_out <- list(
  chrom_sizes = file.path(annotation_path, "chrom.sizes"),
  gene_regions = file.path(annotation_path, "gene_regions.rds"),
  gtf  = file.path(annotation_path, "all_gr.rds"),
  seqinfo = file.path(annotation_path, "seqinfo.rds"),
  transcript_models = file.path(annotation_path, "trans_models.rds"),
  tss = file.path(annotation_path, "tss.rds")
)

#### Seqinfo ####
sq <- samples %>%
  mutate(
    path = here::here(config$paths$bam, target, glue("{sample}.bam"))
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

#### chrom_sizes ####
## For bedGraphToBigWig
sq %>%
  as_tibble() %>%
  dplyr::select(seqnames, seqlengths) %>%
  write_tsv(all_out$chrom_sizes, col_names = FALSE)
cat("chrom_sizes exported...\n")

#### GTF ####
reqd_cols <- c(
  "type",
  "gene_id", "gene_type", "gene_name",
  "transcript_id", "transcript_type", "transcript_name",
  "exon_id"
)
all_gr <- gtf %>%
  import.gff(
    which = GRanges(sq),
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
seqinfo(all_gr) <- sq
write_rds(all_gr, all_out$gtf, compress = "gz")
cat("all_gr.rds written successfully...\n")

#### Transcript Models (Gviz) ####
trans_models <- all_gr$exon %>%
  as_tibble(rangeAsChar = FALSE) %>%
  group_by(transcript_id) %>%
  mutate(exon = paste0(transcript_id, "_", seq_along(transcript_id))) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = sq) %>%
  mutate(feature = as.character(type)) %>%
  select(type, gene = gene_id, exon, transcript = transcript_id, symbol = gene_name) %>%
  sort() %>%
  setNames(.$transcript)
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
tss <- all_gr$transcript %>%
  resize(width = 1) %>%
  reduceMC() %>%
  mutate(region = "TSS") %>%
  select(region, everything(), -type)
write_rds(tss, all_out$tss, compress = "gz")
cat("TSS regions exported...\n")

#### Promoters ####
prom_params <- params$gene_regions$promoters
gene_regions <- list(
  promoters = all_gr$transcript %>%
    promoters(
      upstream = prom_params$upstream,
      downstream = prom_params$downstream
    ) %>%
    select(-ends_with("type"), -starts_with("exon")) %>%
    reduceMC(ignore.strand = TRUE) %>%
    mutate(
      region = glue(
        "Promoter (-{prom_params$upstream}/+{prom_params$downstream}bp)"
      )
    ) %>%
    select(region, everything())
)
cat("Promoters defined...\n")

#### Upstream Promoters ####
cat("Began defining upstream promoters at ", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))
## Remove any which are shorter than 1% of the requested distance as these are likely
## chopped by an exon
min_width <- 0.01 * (params$gene_regions$upstream - params$gene_regions$promoters$upstream)
gene_regions$upstream <- all_gr$transcript %>%
  promoters(upstream = params$gene_regions$upstream, downstream = 0)  %>%
  select(-ends_with("type"), -starts_with("exon")) %>%
  setdiffMC(gene_regions$promoters, ignore.strand = TRUE)%>%
  # Ignoring any exon overlaps ensures an upstream promoter is annotated as
  # being more important than an exon, given that exons refer to RNA structure
  # setdiffMC(all_gr$exon, ignore.strand = TRUE) %>%
  reduceMC() %>%
  subset(width > min_width) %>%
  mutate(
    region = glue("Upstream Promoter (<{params$gene_regions$upstream/1e3}kb)")
  ) %>%
  select(region, everything())
  ## Some of these will possibly extend into other genes.
  ## The only real solution is to cut any sections from the upstream ranges
  ## which overlap other genes, whilst retaining those ranges which are internal
  ## to the gene-of-origin. The only viable way to do that is to manually exclude
  ## the gene-of-origin then take the setdiff. This will take an hour or two
  # split(f = seq_along(.)) %>%
  # bplapply(
  #   function(x) {
  #     gr <- subset(all_gr$gene, !gene_id %in% unlist(x$gene_id))
  #     setdiffMC(x, gr, ignore.strand = TRUE)
  #   },
  #   BPPARAM = bpparam()
  # ) %>%
  # GRangesList() %>%
  # unlist() %>%
  # setNames(c())
cat("Finished defining upstream promoters at ", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

#### Exons ####
gene_regions$exons <- all_gr$exon %>%
  unstrand() %>%
  select(-ends_with("type")) %>%
  reduceMC() %>%
  setdiffMC(
    ## Exons overlapping promoter regions are assigned as promoters
    lapply(gene_regions[c("promoters", "upstream")], granges) %>%
      GRangesList() %>%
      unlist
  ) %>%
  mutate(region = "Exon") %>%
  select(region, everything())
cat("Exons defined at ", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

#### Introns ####
gene_regions$introns <- all_gr$gene %>%
  unstrand() %>%
  select(-ends_with("type"), -starts_with("trans"), -starts_with("exon")) %>%
  setdiffMC(
    lapply(gene_regions[c("promoters", "upstream", "exons")], granges) %>%
      GRangesList() %>%
      unlist()
  ) %>%
  reduceMC() %>%
  mutate(region = "Intron") %>%
  select(region, everything())
cat("Introns defined at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

#### Intergenic Distal ####
suppressWarnings(
  gene_regions$distal <- sq %>%
    GRanges() %>%
    setdiff(
      all_gr$gene %>%
        resize(
          width  = width(.) + 2 * params$gene_regions$intergenic,
          fix = 'center'
        ) %>%
        trim() %>%
        unstrand()
    ) %>%
    mutate(
      region = glue("Intergenic (>{params$gene_regions$intergenic/1e3}kb)")
    )
)
cat("Distal Intergenic Regions defined at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

#### Intergenic Proximal ####
gene_regions$proximal <- setdiff(
  GRanges(sq),
  lapply(
    gene_regions[c("promoters", "upstream", "exons", "introns", "distal")],
    granges
  ) %>%
    GRangesList() %>%
    unlist()
) %>%
  join_nearest(all_gr$gene) %>%
  mutate(
    region = glue("Intergenic (<{params$gene_regions$intergenic/1e3}kb)"),
  ) %>%
  select(region, gene_id, gene_name)
cat("Proximal Intergenic Regions defined at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

#### Final Tidy Up ####
##Remove any columns which are just NA
gene_regions <- gene_regions %>%
  lapply(
    function(x) {
      keep <- vapply(
        mcols(x),
        function(y) {
          if (!is(y, "list_OR_List")) {
            if (all(is.na(y))) return(FALSE)
          }
          return(TRUE)
        },
        logical(1)
      )
      mcols(x) <- mcols(x)[keep]
      x
    }
  )
## Place in the correct order
o <- c(
  "promoters", "upstream", "exons", "introns", "proximal",  "distal"
)
gene_regions <- gene_regions[o]
write_rds(gene_regions, all_out$gene_regions, compress = "gz")
cat("Data export completed at", format(Sys.time(), "%H:%M:%S, %d %b %Y\n"))

