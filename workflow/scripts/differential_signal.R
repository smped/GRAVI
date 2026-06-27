#' This script runs the differential signal analysis and provides all output
#' Providing all peaks/bed files as output enables downstream results to be
#' included in the reporting section of the Rmd
#'
#' Key outputs will be:
#' 1. DiffSig Results
#' 2. IHW
#' 3. Increased Regions
#' 4. Decreased Regions
#'
#' The results can have key parameters added to the metadata
#'
#' Key Inputs:
#' 1. counts
#' 2. gtf_gene
#' 3. gene_regions
#' 4. features
#' 5. hic
#' 6. sq
#' 7. All Consensus peaks!!!
#'
#' Key Params
#' 1. diff_sig_params
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
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}



if ("snakemake" %in% ls()) {

  log <- slot(snakemake, "log")[[1]]
  message("Setting stdout to ", log, "\n")
  sink(log, split = TRUE)
  all_input <- slot(snakemake, "input")
  all_output <- slot(snakemake, "output")
  config <- slot(snakemake, "config")
  all_wildcards <- slot(snakemake, "wildcards")
  all_params <- slot(snakemake, "params")
  threads <- slot(snakemake, "threads")

} else {
  
  ## For testing
  target <- "H3K27ac"
  all_input <- list(
    counts = "output/differential_signal/{target}/{target}_counts.rds",
    gtf = "output/annotations/gtf.rds",
    hic = "output/annotations/hic.rds",
    features = "output/annotations/features.rds",
    regions = "output/annotations/gene_regions.rds",
    sq = "output/annotations/seqinfo.rds",
    yaml = "config/params.yml"
  )
  all_input <- lapply(all_input, glue::glue)
  all_input$peaks <- vapply(
    c("AR", "ER", "GATA3", "H3K27ac"),
    \(x) file.path(
      "output", "peak_analysis", x, paste0(x, "_consensus_peaks.bed.gz")
    ),
    character(1)
  )
  all_output <- list(
    changed = "output/differential_signal/{target}/{target}_E2_E2DHT-changed.bed.gz",
    decreased = "output/differential_signal/{target}/{target}_E2_E2DHT-decreased.bed.gz",
    increased = "output/differential_signal/{target}/{target}_E2_E2DHT-increased.bed.gz",
    ihw = "output/differential_signal/{target}/{target}_E2_E2DHT-ihw.rds",
    rds = "output/differential_signal/{target}/{target}_E2_E2DHT-differential-signal.rds"
  ) |>
    lapply(glue::glue)
  all_params <- list(
    diff_sig_params = jsonlite::fromJSON("config/json/differential_signal_param.json")[[target]],
    peak_calling_params = jsonlite::fromJSON("config/json/peak_calling_param.json")[[target]]
  )
  all_wildcards <- list(target = target, ref = "E2", treat = "E2DHT")
  rm(target)
  config <- here::here("config/config.yml") |> yaml::read_yaml()
  threads <- 4

}

diff_sig_params <- all_params$diff_sig_params
peak_params <- all_params$peak_calling_params

cat_list(all_input, "input")
cat_list(all_output, "output")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(diff_sig_params, "DiffSig Params:", "=")
cat_list(peak_params, "Peak Params:", "=")

## Solidify file paths
all_input <- lapply(all_input, here::here)
all_output <- lapply(all_output, here::here)

cat_time("Loading packages...\n")
library(tidyverse)
library(extraChIPs)
library(edgeR)
library(yaml)
library(quantro)
library(qsmooth)
library(doParallel)
library(BiocParallel)
library(IHW)
library(glue)
library(plyranges)
library(GenomicInteractions)
cat_time("Setting to run using", threads, "threads")
register(MulticoreParam(workers = threads))

cat_time("Checking parameters")
win_type <- match.arg(diff_sig_params$window_type, c("sliding", "fixed"))
method <- match.arg(diff_sig_params$method, c("qlf", "lt", "wald"))
norm <- match.arg(
  diff_sig_params$norm, c("TMM", "TMMwsp", "RLE", "upperquartile", "none", "sq")
)
if (norm == "sq") {
  if (!method == "lt" & !win_type == "sliding") {
    cat("SQ Normalisation only enabled for limma-trend using sliding windows")
    stop()
  }
}
ihw_method <- match.arg(
  diff_sig_params$ihw, c("none", "targets", "regions", "features")
)
fdr_alpha <- diff_sig_params$alpha

cat_time("Loading annotations")
gtf <- read_rds(all_input$gtf)
sq <- read_rds(all_input$sq)
regions <- read_rds(all_input$regions)
features <- read_rds(all_input$features)
has_features <- length(features) > 0
hic <- read_rds(all_input$hic)
mapping_params <- all_input$yaml %>%
  read_yaml() %>%
  pluck("mapping")

cat_time("Loading counts")
counts <- read_rds(all_input$counts)
counts <- counts[,counts$treat %in% unlist(all_wildcards[c("ref", "treat")])]
colData(counts) <- droplevels(colData(counts))

# cat_time("Adding logCPM assay")
assay_name <- ifelse(!method == "lt", "counts", "logCPM")
## Sliding windows (i.e. sq-lt) will already have a logCPM assay which avoids
## any concerns about double normalisation here
if (!"logCPM" %in% assayNames(counts)) {
  dge <- calcNormFactors(counts, method = norm)
  dge$samples$lib.size <- counts$totals
  lcpm <- cpm(dge, log = TRUE)
  rownames(lcpm) <- NULL
  assay(counts, "logCPM") <- lcpm
}

cat_time("Checking logCPM distributions using quantro")
quantro_p <- NULL
if (norm != "none") {
  qtest <- metadata(counts)$quantro
  quantro_p <- c(
    perm = quantroPvalPerm(qtest),
    anova = anova(qtest)[["Pr(>F)"]][[1]]
  )
  cat(sprintf("Lowest p-value is %.3e \n", min(quantro_p)))
} else {
  cat_time("Q-test not required (norm = 'none')")
}
if (any(quantro_p < diff_sig_params$quantro_alpha) & norm != "sq") {
  cat_time("Quantro-test rejected H0. Setting normalisation to none")
  norm <- "none"
}

qs <- NULL
if (norm == "sq") {
  cat_time("Performing Smooth Quantile Normalisation")
  assay_name <- "qsmooth"
  qs <- qsmooth(assay(counts, "logCPM"), group_factor = counts$treat)
  assay(counts, assay_name) <- qsmoothData(qs)
}

cat_time("Defining model parameters")
nesting <- NULL
if (!is.null(diff_sig_params$nesting)) {
  nesting <- match.arg(diff_sig_params$nesting, colnames(colData(counts)))
}
fm <- as.formula(
  ifelse(is.null(nesting), "~treat", paste("~", nesting, "+treat"))
)
cat_time("Model formula set as", as.character(fm))
X <- model.matrix(fm, data = colData(counts))
colnames(X) <- str_remove_all(colnames(X), "treat")
colData(counts)$design <- X
paired_cors <- block <- txt <- NULL
if (!is.null(nesting) & method == "lt") {
  ## These will be passed to fitAssayDiff. This in turn passes these to
  ## lmFit, although when method is qlf they will be passed to glmQLFit.
  ## As they are not parameters for that modelling approach, they will be
  ## ignored
  cat_time("Calculating correlations")
  block <- colData(counts)[[nesting]]
  n_max <- min(1e4, nrow(counts))
  set.seed(1e6)
  ind <- sample.int(nrow(counts), n_max, replace = FALSE)
  paired_cors <- duplicateCorrelation(
    object = assay(counts, assay_name)[ind, ],
    design = X,
    block = block
  )$consensus.correlation
}
cat_time("Fitting model")
fit <- fitAssayDiff(
  counts, assay = assay_name, design = X, coef = all_wildcards$treat,
  method = method, norm = ifelse(norm == "sq", "none", norm),
  fc = diff_sig_params$fc, block = block, correlation = paired_cors
)
pcols <- c("PValue", "p_mu0")

if (win_type == "sliding") {

  cat_time("Merging windows")
  #' This is set to merge within 2 window steps, with minimum window set via params.
  #' Perhaps a more rational solution would be to use the 'merge_within' value
  #' from peak calling. This should've been set for broader signal targets,
  #' whilst will often be zero for narrow-type targets. Given there is a
  #' parallel between peak calling & detection of windows, this may be a simple
  #' strategy that won't require setting of any additional parameters. It will
  #' however, require passing the peak-callaing parameters to the snakemake rule
  #'
  #' Hard-wired to return the 'adaptive' region that's changed as keyval_range
  results <- mergeByHMP(
    fit, pval = pcols,
    merge_within = peak_params$merge_within,
    hm_pre = "", keyval = "merged", min_win = diff_sig_params$min_win
  ) %>%
    select(
      starts_with("n_"), keyval_range, starts_with("log"), any_of(pcols),
      FDR = PValue_fdr
    ) %>%
    addDiffStatus(alpha = fdr_alpha)

  ## Map genes, features & regions, which are otherwise propagated through
  cat_time("Mapping merged windows to regions")
  # results$region <- bestOverlap(results$keyval_range, unlist(regions), var = "region")
  results$region <- bestOverlap(results, unlist(regions), var = "region")
  results$region <- factor(
    results$region, levels = map_chr(regions, \(x) x$region[1])
  )
  if (has_features) {
    cat_time("Mapping merged windows to features")
    feat_df <- features %>%
      lapply(
        # \(x) bestOverlap(results$keyval_range, x, var = "feature")
        \(x) bestOverlap(results, x, var = "feature")
      ) %>%
      as_tibble() %>%
      mutate(range = as.character(results))
    if (length(features) > 1) {
      cat_time("Merging features across source files")
      feat_df <- feat_df %>%
        nest(data = all_of(names(features))) %>%
        mutate(
          feature = lapply(
            data, \(x) {
              x <- unlist(x)
              x[!is.na(x)]
            }
          )
        ) %>%
        unnest(data, keep_empty = TRUE) %>%
        dplyr::select(feature, all_of(names(features))) %>%
        as.data.frame()
      existing <- setdiff(colnames(mcols(results)), colnames(feat_df))
      mcols(results) <- cbind(mcols(results)[existing], feat_df)
      results$feature <- CharacterList(results$feature)
    } else {
      results$feature <- str_replace_na(
        feat_df[[names(features)]], "no_feature"
      )
    }
  }

  cat_time("Defining promoters & enhancers")
  feat_prom <- features %>%
    endoapply(subset, grepl("(prom|tssa$)", str_to_lower(feature))) %>%
    unlist() %>%
    GenomicRanges::reduce()
  prom <- GenomicRanges::reduce(c(feat_prom, granges(regions$promoter)))
  cat_time("Found", length(feat_prom), "promoters in provided features")
  feat_enh <- features %>%
    ## Exclude any 'weak enhancers'
    endoapply(subset, grepl("enh[^w]", str_to_lower(feature))) %>%
    unlist() %>%
    GenomicRanges::reduce() %>%
    filter_by_non_overlaps(prom)
  cat_time("Found", length(feat_enh), "enhancers in provided features")

  cat_time("Mapping to genes")

  map_by_feature_args <- formals(mapByFeature)
  map_by_feature_args[names(mapping_params)] <- mapping_params
  map_by_feature_args <- map_by_feature_args[map_lgl(map_by_feature_args, is, "atomicVector")]
  map_by_feature_args$prom <- prom
  map_by_feature_args$enh <- feat_enh
  map_by_feature_args$gi <- hic
  map_by_feature_args$gr <- results
  map_by_feature_args$genes <- gtf$gene

  results <- do.call("mapByFeature", map_by_feature_args)


} else {
  results <- rowRanges(fit) %>% addDiffStatus(alpha = fdr_alpha)
}
fdr_column <- "FDR"
cat_time("Done")

if (ihw_method != "none") {

  cat_time("Running Independent hypothesis weighting")
  if (ihw_method == "regions") {
    ## Regions are unique, so we can use existing mappings
    covariate <- results$region
  }

  if (ihw_method == "targets") {
    cat_time("Importing consensus peaks")
    ihw_gr <- all_input$peaks %>%
      importPeaks(
        type = "bed", seqinfo = sq, glueNames = "{basename(dirname(x))}"
      ) %>%
      .[names(.) != all_wildcards$target]
    ## Hits may occur to multiple targets
    cat_time("Mapping to other targets")
    hits <- vapply(
      ihw_gr, \(x) overlapsAny(results, x), logical(length(results))
    )
    covariate <- apply(
      hits, MARGIN = 1, \(x) paste(colnames(hits)[x], collapse = " + ")
    )
    covariate[covariate == ""] <- "None"

  }

  if (ihw_method == "features") {
    if (has_features) {
      covariate <- mcols(results)[names(features)] %>%
        lapply(str_replace_na, "no_feature") %>%
        lapply(fct_infreq) %>%
        lapply(fct_lump_min, min = 1e3) %>%
        as_tibble() %>%
        unite(feature, all_of(names(features)), sep = "; ") %>%
        pull("feature")
    } else {
      covariate <- rep_len("no feature", length(results))
    }
  }

  cat_time("Grouping ranges by", ihw_method)
  results$ihw_covariate <- covariate %>%
    fct_infreq() %>%
    fct_lump_min(min = 1e3)

  ## Check the merging by fct_lump_min has left all groups > 1e3
  cat_time("Setting final IHW groups")
  if (any(fct_count(results$ihw_covariate)$n < 1e3)) {
    lv_to_drop <- fct_count(results$ihw_covariate) %>%
      dplyr::filter(f != "Other") %>%
      dplyr::filter(n == min(n)) %>%
      pull("f") %>%
      as.character() %>%
      c("Other")
    results$ihw_covariate <- results$ihw_covariate %>%
      fct_other(drop = lv_to_drop, other_level = "Other")
  }

  ihw_proceed <- length(levels(results$ihw_covariate)) > 1
  if (ihw_proceed) fdr_column <- "fdr_ihw"
  ihw <- NULL
  if (ihw_proceed) {
    cat_time("Running IHW")
    ihw <- ihw(
      pvalues = results$PValue,
      covariates = results$ihw_covariate,
      alpha <- fdr_alpha,
      covariate_type = "nominal"
    )
    cat_time("Updating status column")
    results <- mutate(results, fdr_ihw = adj_pvalues(ihw))
    results <- addDiffStatus(results, sig_col = fdr_column, alpha = fdr_alpha)
  } else {
    cat_time("No viable groupings. Not performing IHW")
    ihw_method <- "none"
    results$ihw_covariate <- NULL
  }

}

cat_time("Updating metadata")
vals <- c(
  "alpha", "fc", "filter_q", "method", "nesting", "window_type",
  "window_size","window_step"
)
metadata(results) <- c(
  diff_sig_params[vals],
  list(
    norm = norm, ref = all_wildcards$ref, treat = all_wildcards$treat,
    ihw = ihw_method
  )
) %>%
  .[sort(names(.))]
metadata(results)$quantro_p <- quantro_p
metadata(results)$description <- glue(
    "Differential Signal for {all_wildcards$target} was assessed using {win_type} windows of {diff_sig_params$window_size}bp ",

    ifelse(
        win_type == "sliding",
        "with a step-size of {diff_sig_params$window_step}. Windows within {peak_params$merge_within}bp were merged after testing using the
        harmonic-mean p-value [@Wilson2019-ln] to obtain a representative p-value for merged regions. ",
        "centred at the estimated peak-centres returned by `macs2 callpeak` [@Zhang18798982]. "
    ),

    ifelse(
        !is.null(quantro_p),
        glue(
            "Distributions of counts between treatment groups were first checked using `quantro` [@HicksQuantro2015] and ",
            ifelse(
                any(quantro_p < 0.05),
                "counts were found to be from different distributions ",
                "no difference in the underlying distributions of counts was found "
            ),
            "(p~perm~ = {round(quantro_p[['perm']], 3)}; p~anova~ = ",
            "{round(quantro_p[['anova']], 3)}). "
        ),
        ""
    ),

    ifelse(norm == "none", "No normalisation was applied. ", "{str_to_upper(norm)}-normalisation was applied "),
    case_when(
        norm == "sq" ~ "[@HicksSQN2017]. ",
        norm == "RLE" ~ "[@Anders2010-sd]. ",
        str_detect(norm, "TMM") ~ "[@Robinson2010-qp]. ",
        TRUE ~ ""
    ),

    "Read totals across the complete genome were always taken as the representative library size for each sample. ",

    ifelse(is.null(nesting), "", "Samples were nested within {nesting}. "),
    "\n\nStatistical analysis was performed using ",
    case_when(
        method == "qlf" ~ "Quasi-Likelihood fits [@LunSmythGLMQL2017] on counts ",
        method == "lt" ~ "*limma-trend* [@LawVoom2014] on normalised logCPM values ",
        method == "wald" ~ "the negative binomial Wald Test on counts [@Love2014Wald] "
    ),

    ifelse(
        diff_sig_params$fc > 0,
        "and a range-based H~0~, setting changed signal within the range [-{round(log2(diff_sig_params$fc), 3)}, {round(log2(diff_sig_params$fc), 3)}] as not being of interest [@McCarthyTreat2009]. ",
        "and a conventional H~0~, testing whether any change in signal is zero or non-zero. "
    ),

    "The analysis tested the treatment {all_wildcards$treat} against the baseline condition of {all_wildcards$ref}. ",

    ifelse(
        ihw_method == "none", "",
        sprintf(
            "P-values after all testing were then weighted using IHW [@IgnatiadisIHW2016], setting overlap with %s as the covariate. ",
            case_when(
                ihw_method == "regions" ~ "genomic regions",
                ihw_method == "targets" ~ "consensus peaks from alternative targets",
                ihw_method == "features" ~ "supplied features"
            )
        )
    ),
    sprintf(
        "Significant differential signal was determined using %s-adjusted p-values < {fdr_alpha} along with direction of change.",
        ifelse(ihw_method == "none", "FDR", "FDR~IHW~")
    )
)


## Think a bit more carefully here. Maybe just update in the Rmd
# cat_time("Exporting counts with updated assays")
# write_rds(counts, all_input$counts, compress = "gz")
cat_time("Exporting results")
write_rds(results, all_output$rds, compress = "gz")
cat_time("Exporting changed sites")
results %>%
  filter(grepl("(In|De)creased", status)) %>%
  write_bed(all_output$changed)
cat_time("Exporting increased sites")
results %>%
  filter(status == "Increased") %>%
  write_bed(all_output$increased)
cat_time("Exporting decreased sites")
results %>%
  filter(status == "Decreased") %>%
  write_bed(all_output$decreased)
cat_time("Exporting IHW results")
write_rds(ihw, all_output$ihw)


cat_time("Done")
