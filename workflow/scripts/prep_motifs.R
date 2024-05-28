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
library(glue)
library(yaml)
library(MotifDb)
library(universalmotif)
params <- read_yaml(all_input$yaml)


#### Motifs ####
## If provided in external, ignore all settings in params.yml & use that
## Only MEME and JASPAR format are supported for manually provided files
db <- NULL
if (!is.null(config$external$motifdb)) {

  f <- config$external$motifdb[[1]] # Ignore all but the first
  ln <- read_lines(f, n_max = 1)

  ## Check MEME format by looking for MEME in the first line
  if (grepl("MEME", ln))  {
    db <- read_meme(f) |> to_df()
    cat_time("Imported", nrow(db), "motifs in MEME format")
  }

  ## Check for JASPAR format by looking for `>` at the start of every 5th line
  if (is.null(db) & grepl("^>", ln))  {
    db <- read_jaspar(f) |> to_df()
    cat_time("Imported", nrow(db), "motifs in JASPAR format")
  }

  ## Check that the `altname` column has values
  na_alts <- is.na(db$altname)
  db$altname[na_alts] <- db$name[na_alts]
  ## As a test, also try extracting the first word from the name column
  db$name[na_alts] <- str_extract(db$name, "^[A-Z0-9]+")

  if (is.null(db))
    cat_time("Unable to determine motif format. Using settings from params")

}

motif_params <- params$motifdb

if (is.null(db)) {

  cat_time("Converting to Universal Motif format from MotifDb\n")
  db <- convert_motifs(MotifDb) |> to_df()
  if (is.null(motif_params$data_source))
    stop("No data source provided for transcription factors")
  db <- subset(db, dataSource %in% motif_params$data_source)
  if (!is.null(motif_params$organism))
    db <- subset(db, organism %in% motif_params$organism)
  cat_time("Database has been subset to", nrow(db), "motifs\n")

}

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


cat_time("Data export completed")


