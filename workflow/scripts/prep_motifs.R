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

if ("snakemake" %in% ls()) {
  log <- slot(snakemake, "log")[[1]]
  message("Setting stdout to ", log, "\n")
  sink(log, split = TRUE)
  all_input <- slot(snakemake, "input")
  all_output <- slot(snakemake, "output")
  config <- slot(snakemake, "config")
} else {
  os.path.join <- here::here
  annotation_path = os.path.join("output", "annotations")
  all_input <- list(
    script = os.path.join("workflow", "scripts", "prep_motifs.R"),
    yaml = os.path.join("config", "params.yml")
  )
  all_output <- list(
    motifs = os.path.join(annotation_path, "motif_list.rds"),
    motif_uri = os.path.join(annotation_path, "motif_uri.rds")
  )
  config <- here::here("config", "config.yml") |> yaml::read_yaml()
}


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
library(motifTestR)
params <- read_yaml(all_input$yaml)

#' Add a function for clustering which handles any NA values in the distance
#' matrix. This has been updated in later versions of motifTestR
clusterMotifs <- function(
        motifs, type = c("PPM", "ICM"),
        method = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN", "ALLR_LL", "WEUCL", "WPCC"),
        power = 1, agglom = "complete", thresh = 0.2, return_d = FALSE,
        plot = FALSE, labels = FALSE, cex = 1, main = "Motif Cluster Dendrogram",
        sub = NULL, xlab = NULL, ylab = "Height", linecol = "red", ...
){
    # Convert to universal motif, if a list is passed
    if (all(vapply(motifs, is, logical(1), class2 = "matrix"))) {
        motifs <- lapply(
            names(motifs), \(x) create_motif(motifs[[x]], name = x)
        )
    }
    stopifnot(all(vapply(motifs, is, logical(1), class2 = "universalmotif")))
    method <- match.arg(method)
    type <- match.arg(type)
    if (type == "ICM" & method %in% c("ALLR", "ALLR_LL"))
        stop("Cannot use ICM with ALLR or ALLR_LL")
    is_dist <- method %in% c("EUCL", "KL", "HELL", "SEUCL", "MAN", "WEUCL")
    mat <- compare_motifs(motifs, use.type = type, method = method, ...)
    ## This is really only useful for correlations
    if (power != 1 & method %in% c("PCC", "WPCC")) mat <- mat^power
    ## Make a distance/dissimilarity matrix
    mat <- abs(mat) / max(abs(mat), na.rm = TRUE) # Scale to be in [0,1]
    mat[is.na(mat)] <- 1 # Set NA to be the maximum distance
    if (!is_dist) mat <- 1 - abs(mat)
    d <- as.dist(mat)
    cl <- hclust(d, method = agglom)
    if (plot) {
        plot(
            cl, labels = labels, cex = cex, main = main, sub = sub, xlab = xlab,
            ylab = ylab
        )
        abline(a = thresh, b = 0, col = linecol)
    }
    cl <- cutree(cl, h = thresh)
    if (return_d) {
        nm <- lapply(split(cl, cl), names)
        d_split <- lapply(nm, \(x) as.matrix(d)[x, x])
        return(list(cl = cl, d = d_split))
    }
    cl
}

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

## A better thing to do here would be to run an md5sum on the motif matrices 
## and remove duplicates based on that, but for now just remove duplicates 
## based on the name column
cat_time("Ensuring no duplicated names")
db <- distinct(db, name, .keep_all = TRUE)
cat_time("Database has been subset to", nrow(db), "motifs\n")

cat_time("Clustering motifs\n")
cluster_params <- list(
  motifs = to_list(db), plot = FALSE, return_d = TRUE, nthreads = 1
) %>% 
  c(
    params$motif_clustering %>% 
      setNames(str_replace_all(names(.), "_", "."))
  ) 
cluster_ids <- do.call("clusterMotifs", cluster_params)
db$cluster <- cluster_ids$cl
cat_time("Motifs clustered into", length(unique(cluster_ids$cl)), "clusters\n")

cat_time("Exporting to:", all_output$motifs, "\n")
db |>
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


