#' This checks that all provided arguments are valid across all yaml files
#' Those passed in the main config file can be taken directly from the
#' snakemake object
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
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

# all_input <- list(
#   params = "config/params.yml",
#   colours = "config/colours.yml"
# )
# all_output <- list("output/checks/args.chk")
# config <- yaml::read_yaml("config/config.yml")

all_input <- slot(snakemake, "input")
all_output <- slot(snakemake, "output")
config <- slot(snakemake, "config")

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

cat_list(all_input, "input:", "=")

all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages")
library(yaml)
library(MotifDb)

params <- read_yaml(all_input$params)
colours <- read_yaml(all_input$colours)
err <- c()

cat_time("Checking params for gene regions")
if (!is.logical(params$gene_regions$intron)) err <- c(err, "gene_regions$intron must be logical")
if (!all(
  vapply(lapply(params$gene_regions[c("promoters", "upstream", "proximal")], unlist), is.numeric, logical(1))
  ))
  err <- c(err, "Region parameters malformed. These must match the arguments for defineRegions")

cat_time("Checking params for mapping to genes")
if (!all(unlist(params$mapping) >= 0)) err <- c(err, "Mapping parameters must be all >= 0")

cat_time("Checking MotifDb Params")
if (!any(params$motifdb$data_source %in% mcols(MotifDb)$dataSource))
  err <- c(err, "No specified Motif data sources match those in the data base")
if (!any(params$motifdb$organism %in% mcols(MotifDb)$organism))
  err <- c(err, "No specified Motif organisms match those in the data base")

cat_time("Checking enrichment params")
if (!params$enrichment$method %in% c("great", "gene_id"))
  err <- c(err, "Invalid enrichment method. Must be 'great' or 'gene_id'")


cat_time("Checking motif parameters")
motif_mods <- config$motif_analysis |>
  lapply(\(x) x$model) |>
  unlist() |>
  unique()
if (!all(motif_mods %in% eval(formals(motifTestR::testMotifEnrich)$model)))
  err <- c(err, "Invalid models for enrichment testing")

cat_time("Checking all significance-related params")
adj <- params |>
  lapply(\(x) x$adj) |>
  unlist() |>
  c(unlist(lapply(config$motif_analysis, \(x) x$adj))) |>
  unique()
if (!all(adj %in% p.adjust.methods)) err <- c(err, "Invalid p-value adjustment method")

cat_time("Checking Differential Signal params")
ds_win_type <- config$differential_signal |>
  lapply(\(x) x$window_type) |>
  unlist()
if (!all(ds_win_type %in% c("fixed", "sliding")))
  err <- c(err, paste(
    "Invalid window type for",
    paste(names(ds_win_type)[!ds_win_type %in% c("sliding", "fixed")], collapse = ", ")
  ))

ds_norm <- config$differential_signal |>
  lapply(\(x) x$norm) |>
  unlist()
valid_norm <- c(eval(formals(edgeR::calcNormFactors.default)$method), "sq")
if (!all(ds_norm %in% valid_norm))
  err <- c(err, paste(
    "Invalid normalisation for",
    paste(names(ds_norm)[!ds_norm %in% valid_norm], collapse = ", ")
  ))

ds_methods <- config$differential_signal |>
  lapply(\(x) x$method) |>
  unlist()
valid_methods <- c("qlf", "lt", "wald")
if (!all(ds_methods %in% valid_methods))
  err <- c(err, paste(
    "Invalid analytic method for",
    paste(names(ds_methods)[!ds_methods %in% valid_methods], collapse = ", ")
  ))

ds_ihw <- config$differential_signal |>
  lapply(\(x) x$ihw) |>
  unlist()
valid_ihw <- c("none", "targets", "features", "regions")
if (!all(ds_ihw %in% valid_ihw))
  err <- c(err, paste(
    "Invalid IHW method for",
    paste(names(ds_ihw)[!ds_ihw %in% valid_ihw], collapse = ", ")
  ))

cat_time("Checking profile heatmaps")
bw_type <- config$profile_heatmaps |>
  lapply(\(x) x$bw_type) |>
  unlist()
valid_types <- c("FE", "coverage")
if (!all(bw_type %in% valid_types))
  err <- c(err, paste(
    "Invalid BigWig type for",
    paste(names(bw_type)[!bw_type %in% valid_types], collapse = ", ")
  ))

prof_cols <- config$profile_heatmaps |>
  lapply(\(x) col2rgb(x$gradient)) |>
  vapply(is.matrix, logical(1))
if (!all(prof_cols))
  err <- c(err, paste(
    "Invalid profile heatmaps colour specification for",
    paste(names(which(!prof_cols)), collapse = ", ")
  ))

cat_time("Checking colours")
valid_cols <- colours |>
  lapply(unlist) |>
  lapply(col2rgb) |>
  vapply(is.matrix, logical(1))
## This will automatically error...

#################################
## Have a think about the rest ##
#################################

if (length(err)) {
  cat_time("Errors found")
  cat(err, sep = "\n")
  stop()
}

file.create(all_output[[1]])
