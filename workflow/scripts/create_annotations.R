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
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%m-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}


log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
all_params <- slot(snakemake, "params")
config <- slot(snakemake, "config")

cat_list(all_input, "input:")
cat_list(all_output, "output:")
cat_list(all_params, "params:", "-")

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
library(MotifDb)
library(universalmotif)
library(msigdbr)
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
write_rds(sq, all_output$seqinfo)
cat_time("Seqinfo exported...\n")

#### Blacklist ####
## Set all provided files into a single GRanges & export
cat_time("Forming unified blacklist from all files...")
bl <- config$external$blacklist %>%
  unlist() %>%
  here::here() %>%
  importPeaks(type = "bed", seqinfo = sq) %>%
  unlist() %>%
  GenomicRanges::reduce() %>%
  sort()
write_rds(bl, all_output$blacklist)
cat_time("done\n")

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
cat_time("GTF imported successfully...\n")
seqlevels(all_gtf) <- seqlevels(sq)
seqinfo(all_gtf) <- sq

cat_time("Exporting gene, transcript and exon-level objects\n")
write_rds(all_gtf$gene, all_output$genes, compress = "gz")
write_rds(all_gtf$transcript, all_output$transcripts, compress = "gz")
write_rds(all_gtf$exon, all_output$exons, compress = "gz")
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
gr_params <- params$gene_regions
gene_regions <- defineRegions(
  genes = all_gtf$gene, transcripts = all_gtf$transcript, exons = all_gtf$exon,
  promoter = unlist(gr_params$promoter), upstream = gr_params$upstream,
  proximal = gr_params$intergenic
)

cat_time("Exporting gene_regions...\n")
write_rds(gene_regions, all_output$regions, compress = "gz")

#### Features ####
feat <- GRangesList()
seqinfo(feat) <- sq
if (!is.null(config$external$features)) {
  fl <- unlist(config$external$features)
  cat_time("Parsing features from", fl, "\n")
  feat <- lapply(fl, import.gff, which = GRanges(sq))
  feat <- lapply(feat, select, feature)
  feat <- unlist(GRangesList(feat))
  seqlevels(feat) <- seqlevels(sq)
  seqinfo(feat) <- sq
  cat_time("Finding overlap with gene regions...")
  ol <- lapply(gene_regions, \(x) propOverlap(feat, x))
  mcols(feat) <- cbind(mcols(feat), DataFrame(ol))
  cat_time("done\n")
  feat <- splitAsList(feat, feat$feature)
  cat_time("Features have split into a GRangesList of length", length(feat), "\n")
  cat_time("Features provided appear to be:", paste0("\n\t", names(feat)), "\n")
  cat_time("Writing to", all_output$features, "...")
} else {
  cat_time(
    "No features provided. Writing an empty object to",
    all_output$features, "..."
  )
}
write_rds(feat, all_output$features, compress = "gz")
cat_time("done\n")

#### Motifs ####
motif_params <- params$motifdb
cat_time("Converting to Universal Motif format\n")
db <- convert_motifs(MotifDb) |> to_df()
if (is.null(motif_params$data_source))
  stop("No data source provided for transription factors")
db <- subset(db, dataSource %in% motif_params$data_source)
if (!is.null(motif_params$organism))
  db <- subset(db, organism %in% motif_params$organism)
cat_time("Database has been subset to", nrow(db), "motifs\n")

cat_time("Calculating correlations between motifs...")
cormat <- db |>
  to_list() |>
  compare_motifs(method = "PCC", use.type = "ICM")
cormat[cormat < motif_params$cluster_above] <- 0
cormat[is.na(cormat)] <- 0
cat_time("done\n")
cat_time("Clustering motifs...\n")
cl <- as.dist(1 - cormat) %>%
  hclust() %>%
  cutree(h = 0.9)
cat_time(max(cl), "clusters formed\n")

cat_time("Exporting to:", all_output$motifs, "\n")
db |>
  mutate(cluster = cl[name]) |>
  to_list() |>
  write_rds(all_output$motifs, compress = "gz")
cat_time("done\n")

cat_time("Creating IC Matrix Thumbnails as uri strings")
img_path <- tempdir()
cat_time("Writing motifs to", img_path)
motif_uri <- db |>
  to_list() |>
  lapply(
    \(x) {
      w <- 30 * ncol(x)
      nm <- slot(x, "altname")
      png_out <- file.path(img_path, paste0(nm, ".png"))
      png(png_out, height = 150, width = w)
      p <- view_motifs(x) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
        )
      print(p)
      dev.off()
      knitr::image_uri(png_out)
    }
  ) |>
  setNames(db$altname)
cat_time("Removing", img_path)
unlink(img_path, recursive = TRUE)
cat_time("Done")

cat_time("Writing", all_output$motif_uri)
write_rds(motif_uri, all_output$motif_uri, compress = "gz")
cat_time("Done")

#### MSigDB ####
cat_time("Preparaing MSigDB using msigdbr...")
msigdb_params <- params$msigdb
msigdb <- msigdbr(msigdb_params$species) %>%
  dplyr::filter(
    gs_cat %in% msigdb_params$gs_cat | gs_subcat %in% msigdb_params$gs_subcat,
    ensembl_gene %in% all_gtf$gene$gene_id
  ) %>%
  dplyr::filter(
    dplyr::n() >= min(msigdb_params$size),
    dplyr::n() <= max(msigdb_params$size),
    .by = gs_name
  )
cat_time("Loaded", length(unique(msigdb$gs_name)), "gene-sets")

cat_time("Updating Gene-Set URLs...")
gs_url <- msigdb %>%
  distinct(gs_cat, gs_subcat, gs_name, gs_url, gs_exact_source) %>%
  mutate(
    gs_url = case_when(
      gs_subcat == "CP:REACTOME" ~ str_remove_all(gs_url, "\\|.+"),
      gs_subcat == "CP:KEGG" ~ paste0("https://www.genome.jp/pathway/", gs_exact_source),
      gs_subcat == "CP:WIKIPATHWAYS" ~ paste0(
        "https://www.wikipathways.org/pathways/", gs_exact_source, ".html"
      ),
      gs_url == "" ~ "http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp",
      TRUE ~ gs_url
    ) %>%
      setNames(gs_name)
  ) %>%
  pull(gs_url)
msigdb$gs_url <- gs_url[msigdb$gs_name]

cat_time("Exporting to", all_output$msigdb)
write_rds(msigdb, all_output$msigdb, compress = "gz")
cat_time("Done")

cat_time("Data export completed")

#### Prepare the RMD ####
cat_time("Forming annotation_description.Rmd")
all_output <- lapply(
  all_output,
  str_remove_all,
  pattern = paste0(here::here(), .Platform$file.sep)
)
ln <- glue(
  "
	---
	title: 'Description of Annotations'
	date: \"`r format(Sys.Date(), '%d %B, %Y')`\"
	bibliography: references.bib
	link-citations: true
	params:
	  chrom_sizes: \"{{all_output$chrom_sizes}}\"
	  colours: \"{{all_params$colours}}\"
	  features: \"{{all_output$features}}\"
	  gene_regions: \"{{all_output$regions}}\"
	  gtf_exon: \"{{all_output$exons}}\"
	  gtf_gene: \"{{all_output$genes}}\"
	  gtf_transcript: \"{{all_output$transcripts}}\"
	  motif_list: \"{{all_output$motifs}}\"
	  motif_uri: \"{{all_output$motif_uri}}\"
	  seqinfo: \"{{all_output$seqinfo}}\"
	  trans_models: \"{{all_output$trans_models}}\"
	  tss: \"{{all_output$tss}}\"
	---

	",
  .open = "{{",
  .close = "}}"
)
readr::write_lines(ln, all_output$rmd)

cat_time("Written YAML header; Appending Module")
file.append(all_output$rmd, all_input$module)

cat_time("Done")

