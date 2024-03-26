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

log <- slot(snakemake, "log")[[1]]
message("Setting stdout to ", log, "\n")
sink(log, split = TRUE)

## For testing
all_input <- list(
  counts = "output/counts/ER_counts.rds",
  gtf_gene = "output/annotations/gtf_gene.rds",
  hic = "output/annotations/hic.rds",
  features = "output/annotations/features.rds",
  peaks = vapply(
    c("AR", "ER", "H3K27ac"),
    \(x) file.path(
      "output", "peak_analysis", x, paste0(x, "_consensus_peaks.bed.gz")
    ), character(1)
  ),
  regions = "output/annotations/gene_regions.rds",
  sq = "output/annotations/seqinfo.rds",
  yaml = "config/params.yml"
)
all_output <- list(
  decreased = "output/differential_signal/ER/ER_E2_E2DHT-decreased.bed.gz",
  increased = "output/differential_signal/ER/ER_E2_E2DHT-increased.bed.gz",
  ihw = "output/differential_signal/ER/ER_E2_E2DHT-ihw.rds",
  rds = "output/differential_signal/ER/ER_E2_E2DHT-differential-signal.rds"
)
all_params <- list(
  alpha = 0.05,
  contrasts = structure(c("E2", "E2DHT"), dim = 1:2),
  fc = 1.2,
  filter_q = NULL,
  ihw = "targets",
  method = "qlf",
  norm = "TMM",
  pair_column = NULL,
  rna_toptable = "data/external/ZR75_DHT_StrippedSerum_RNASeq_topTable.tsv",
  window_type = "fixed",
  window_size = 400L,
  window_step = NULL,
  target = "ER"
)
all_wildcards <- list(target = "ER", ref = "E2", treat = "E2DHT")
config <- yaml::read_yaml("../GRAVI_testing/config/config.yml")
threads <- 4

# all_input <- slot(snakemake, "input")
# all_output <- slot(snakemake, "output")
# config <- slot(snakemake, "config")
# all_wildcards <- slot(snakemake, "wildcards")
# all_params <- slot(snakemake, "params")$diff_sig_params
# threads <- slot(snakemake, "threads")

cat_list(all_input, "input")
cat_list(all_output, "output")
cat_list(all_wildcards, "wildcards:", "=")
cat_list(all_params, "DiffSig Params:", "=")

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
cat_time("Setting to run using", threads, "threads")
register(MulticoreParam(workers = threads))

cat_time("Checking parameters")
win_type <- match.arg(all_params$window_type, c("sliding", "fixed"))
method <- match.arg(all_params$method, c("qlf", "lt"))
norm <- match.arg(
  all_params$norm, c("TMM", "TMMwsp", "RLE", "upperquartile", "none", "sq")
)
if (norm == "sq") {
  if (method == "qlf" | win_type == "fixed") {
    cat("SQ Normalisation only enabled for limma-trend using sliding windows")
    stop()
  }
}
ihw_method <- match.arg(
  all_params$ihw, c("none", "targets", "regions", "features")
)
fdr_alpha <- all_params$alpha

cat_time("Loading annotations")
sq <- read_rds(all_input$sq)
regions <- read_rds(all_input$regions)
features <- read_rds(all_input$features)
hic <- read_rds(all_input$hic)
mapping_params <- all_input$yaml %>%
  read_yaml() %>%
  pluck("mapping")

cat_time("Loading counts")
counts <- read_rds(all_input$counts)

cat_time("Adding logCPM assay")
assay_name <- ifelse(method == "qlf", "counts", "logCPM")
## Sliding windows (i.e. sq-lt) will already have a logCPM assay
if (!"logCPM" %in% assayNames(counts)) {
  dge <- calcNormFactors(counts, method = norm)
  dge$samples$lib.size <- counts$totals
  lcpm <- cpm(dge, log = TRUE)
  rownames(lcpm) <- NULL
  assay(counts, "logCPM") <- lcpm
}

cat_time("Checking count distributions using quantro")
quantro_p <- NULL
if (norm != "none") {
  registerDoParallel(threads)
  qtest <- quantro(assay(counts, "counts"), counts$treat, B = 1e3)
  quantro_p <- c(
    perm = quantroPvalPerm(qtest),
    anova = anova(qtest)[["Pr(>F)"]][[1]]
  )
  cat(sprintf("Lowest p-value is %.3e", min(quantro_p)))
} else {
  cat_time("Q-test not required (norm = 'none')")
}
if (any(quantro_p < 0.05) & norm != "sq") {
  cat_time("Q-test failed. Setting normalisation to none")
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
pair_col <- NULL
if (!is.null(all_params$pair_column)) {
  pair_col <- match.arg(all_params$pair_column, colnames(colData(counts)))
}
fm <- as.formula(
  ifelse(is.null(pair_col), "~treat", paste("~", pair_col, "+treat"))
)
cat_time("Model formula set as", as.character(fm))
X <- model.matrix(fm, data = colData(counts))
colnames(X) <- str_remove_all(colnames(X), "treat")
colData(counts)$design <- X
paired_cors <- block <- txt <- NULL
if (!is.null(pair_col) & method == "lt") {
  ## These will be passed to fitAssayDiff. This in turn passes these to
  ## lmFit, although when method is qlf they will be passed to glmQLFit.
  ## As they are not parameters for that modelling approach, they will be
  ## ignored
  cat_time("Calculating correlations")
  block <- colData(counts)[[pair_col]]
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
  fc = all_params$fc, block = block, correlation = paired_cors
)
pcols <- c("PValue", "p_mu0")

if (win_type == "sliding") {
  cat_time("Merging windows")
  results <- mergeByHMP(
    fit, pval = pcols,
    merge_within = floor(1 + 2 * all_params$window_size / 3),
    hm_pre = ""
  ) %>%
    plyranges::select(
      starts_with("n_"), keyval_range, starts_with("log"), any_of(pcols),
      FDR = PValue_fdr
    ) %>%
    addDiffStatus(alpha = fdr_alpha)

  ## Map genes, features & regions, which are otherwise propagated through
  cat_time("Mapping merged windows to regions")
  results$region <- bestOverlap(results, unlist(regions), var = "region")
  results$region <- factor(
    results$region, levels = map_chr(regions, \(x) x$region[1])
  )
  if (has_features) {
    cat_time("Mapping merged windows to features")
    results$feature <- bestOverlap(results, features, missing = "no_feature")
  }

  cat_time("Defining promoters & enhancers")
  which_prom <- grepl("prom", str_to_lower(names(features)))
  feat_prom <- features[which_prom] %>%
    unlist() %>%
    GenomicRanges::reduce()
  which_enh <- grepl("enhanc", str_to_lower(names(features)))
  feat_enh <- features[which_enh] %>%
    unlist() %>%
    GenomicRanges::reduce()

  cat_time("Mapping to genes")
  prom <- GenomicRanges::reduce(c(feat_prom, granges(regions$promoter)))
  results <- mapByFeature(
    results, gtf_gene,
    prom = prom,
    enh = feat_enh,
    gi = hic,
    gr2gene = mapping_params$gr2gene,
    prom2gene = mapping_params$prom2gene,
    enh2gene = mapping_params$enh2gene,
    gi2gene = mapping_params$gi2gene
  )
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
    ## Hits may occur to multiple features
    hits <- vapply(
      features, \(x) overlapsAny(results, x), logical(length(results))
    )
    covariate <- apply(
      hits, MARGIN = 1, \(x) paste(colnames(hits)[x], collapse = " + ")
    )
    covariate[covariate == ""] <- "None"
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
  }

}
cat_time("Updating metadata")
vals <- c(
  "alpha", "fc", "filter_q", "method", "pair_column", "window_type",
  "window_size","window_step"
)
metadata(results) <- c(
  all_params[vals],
  list(
    norm = norm, ref = all_wildcards$ref, treat = all_wildcards$treat,
    ihw = ihw_method
  )
) %>%
  .[sort(names(.))]
metadata(results)$description <- glue(
    "Differential Signal for {all_wildcards$target} was assessed using {win_type} windows of {all_params$window_size}bp ",
    ifelse(
        win_type == "sliding",
        "with a step-size of {all_params$window_step}. Windows were merged after testing using the
        harmonic-mean p-value [@Wilson2019-ln] to obtain a representative p-value for merged regions. ",
        "centred at the estimated peak-centres returned by `macs2 callpeak` [@Zhang18798982]. "
    ),
    ifelse(
        !is.null(quantro_p),
        glue(
            "Distributions of counts between treatment groups were checked using quantro [@HicksQuantro2015] and ",
            ifelse(
                any(quantro_p < 0.05),
                "counts were found to be from different distributions. ",
                "no difference in the underlying distributions of counts was found. "
            )
        ),
        ""
    ),
    ifelse(norm == "none", "No normalisation was applied. ", "{norm}-normalisation was applied "),
    case_when(
        norm == "sq" ~ "[@HicksSQN2017]. ",
        norm == "RLE" ~ "[@Anders2010-sd]. ",
        str_detect(norm, "TMM") ~ "[@Robinson2010-qp]. ",
        TRUE ~ ""
    ),
    "Read totals across the complete genome were always taken as the representative library size for each sample. ",
    ifelse(is.null(pair_col), "", "Samples were nested within {pair_col}. "),
    "\n\nStatistical analysis was performed using ",
    case_when(
        method == "qlf" ~ "Quasi-Likelihood fits [@LunSmythGLMQL2017] on counts ",
        method == "lt" ~ "Limma-Trend [@LawVoom2014] on normalised logCPM values ",
    ),
    ifelse(
        all_params$fc > 0,
        "and a range-based H~0~, setting changed signal within the range [-{round(log2(all_params$fc), 3)}, {round(log2(all_params$fc), 3)}] as not being of interest [@McCarthyTreat2009]. ",
        "and a conventional H~0~, testing whether any change in signal is zero or non-zero. "
    ),
    "The analysis tested the treatment {all_wildcards$treat} against the baseline condition of {all_wildcards$ref}. ",
    ifelse(
        ihw_method == "none", "",
        sprintf(
            "P-values after all testing were then weighted using IHW [@IgnatiadisIHW2016] setting overlap with %s as the covariate. ",
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



