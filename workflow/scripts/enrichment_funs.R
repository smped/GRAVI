#' Prepare a wrapper for running enrichment testing either using rGREAT or
#' a Hypergeometric test. The hypergeometric should generally be performed
#' with accounting for bias, but in the case of checking for overall hits to
#' genes, gene_length bias can be included by setting the method to "Wallenius"
#'
#' The simplest strategy may be to ensure that all ranges passed have gene_ids
#' attached. These can be ignored by rGREAT but will be used by goseq
#'
#' @param test_ranges The peaks/ranges to be tested. Must have a gene_id column
#' @param db The msigdb database
#' @param gtf Gene-based GTF
#' @param bg A set of background ranges if not testing against the genome.
#' Must have a gene_id column
#' @param genome The UCSC genome
#' @param min_sig Discard any significant hits below this value
#' @param adj_method One p.adjust.methods
#' @param threads The number of cores
#' @param method The goseq method to use
.geneid_enrich <- function(
    test_ranges, db, gtf, bg = NULL, genome, min_sig, threads,
    adj_method = p.adjust.methods, method = c("Hypergeometric", "Wallenius"),
    ...
){
  method <- match.arg(method)
  adj_method <- match.arg(adj_method)
  stopifnot("gene_id" %in% colnames(mcols(test_ranges)))

  # Make the pwf
  ## The background must include the test set
  bg_ids <- unique(gtf$gene_id)
  test_ids <- unlist(test_ranges$gene_id)
  if (!is.null(bg)) {
    bg_ids <- unique(unlist(bg$gene_id))
    test_ids <- intersect(test_ids, bg_ids)
  }

  ## Setup the pathways
  db <- dplyr::filter(db, ensembl_gene %in% bg_ids)
  db <- droplevels(db)

  ## Form an empty object for no sig genes
  res <- tibble::tibble(
    gs_name = character(), observed_gene_hits = integer(),
    gene_set_size = integer(),  p = numeric(), adj_p = numeric(),
    gene_id = list(), genes_with_hits = character()
  )
  if (length(test_ids) == 0 | all(bg_ids %in% test_ids)) return(res)

  if (method == "Wallenius") {
    library(goseq)
    ## This is classic goseq
    gs_list <- split(db$gs_name, db$ensembl_gene)
    tbl <- tibble::tibble(gene_id = gtf$gene_id, w = width(gtf))
    tbl <- tbl[tbl$gene_id %in% bg_ids,]
    tbl <- dplyr::arrange(tbl, desc(w))
    tbl <- tbl[!duplicated(tbl$gene_id),]
    de <- tbl$gene_id %in% test_ids
    names(de) <- tbl$gene_id
    pwf <- nullp(
      de, genome = genome, bias.data = log10(tbl$w), plot.fit = FALSE
    )
    if (sum(pwf$DEgenes) == 0) return(res)
    res <- goseq(
      pwf = pwf, genome = genome, gene2cat = gs_list, method = method
    )
    res <- as_tibble(res)
    res <- dplyr::select(
      res, gs_name = category, observed_gene_hits = numDEInCat,
      gene_set_size = numInCat, p = over_represented_pvalue
    )

  } else {

    ## This is a slightly more direct hypergeometric test based on goseq
    n_de <- sum(bg_ids %in% test_ids)
    n_genes <- length(bg_ids)
    if (n_de == 0) return(res)
    gs_list <- split(db$ensembl_gene, db$gs_name)
    list_res <- mclapply(
      gs_list,
      function(x){
        observed_gene_hits <- sum(x %in% test_ids)
        gene_set_size <- sum(x %in% bg_ids)
        p = dhyper(
          observed_gene_hits, gene_set_size, n_genes - gene_set_size, n_de
        ) + phyper(
          observed_gene_hits, gene_set_size, n_genes - gene_set_size, n_de,
          lower.tail = FALSE
        )
        tibble(observed_gene_hits, gene_set_size, p)
      }, mc.cores = threads
    )
    res <- dplyr::bind_rows(list_res, .id = "gs_name")

  }

  test_db <- dplyr::filter(db, ensembl_gene %in% test_ids)
  genes_by_gs <- split(test_db$gene_symbol, test_db$gs_name)
  res <- dplyr::arrange(res, p)
  res$adj_p <- p.adjust(res$p, adj_method)
  res <- res[res$observed_gene_hits >= min_sig,]
  res$gene_id <- genes_by_gs[res$gs_name]
  res$genes_with_hits <- vapply(
    res$gene_id, \(x) paste(sort(unique(x)), collapse = "; "), character(1)
  )
  tibble::as_tibble(res)

}


#' @param test_ranges The peaks/ranges to be tested
#' @param db The msigdb database
#' @param gtf Gene-based GTF
#' @param bg A set of background ranges if not testing against the genome
#' @param genome The UCSC genome
#' @param min_sig Discard any significant hits below this value
#' @param adj_method One p.adjust.methods
#' @param threads The number of cores
#' @param params Upstream & downstream distances for promoters. Taken from the yaml
#' @param min_size Discard any gene_sets smaller than this before testing
#'
#'
.great_enrich <- function(
    test_ranges, db, gtf, bg = NULL, genome, min_sig, threads, params,
    adj_method = p.adjust.methods, ...
) {

  res <- tibble::tibble(
    gs_name = character(), genome_fraction = numeric(),
    observed_region_hits = integer(), mean_tss_dist = numeric(),
    gene_set_size = integer(), observed_gene_hits = integer(),
    fold_enrichment = numeric(), gene_with_hits = character(), gene_id = list(),
    p = numeric(), adj_p = numeric()
  )
  if (length(test_ranges) == 0) return(res)

  library(rGREAT)
  ## Update any genome info
  genome(test_ranges) <- genome
  if (!is.null(bg)) {
    genome(bg) <- genome
    ## The bg set must include all of the test set
    stopifnot(any(overlapsAny(test_ranges, bg)))
  }

  ## Define the TSS
  id2gene <- setNames(gtf$gene_name, gtf$gene_id)
  gtf <- plyranges::select(gtf, gene_id)
  genome(gtf) <- genome
  tss <- extendTSS(
    gtf, genome = genome, gene_id_type = 'ENSEMBL',
    basal_upstream = params$upstream, basal_downstream = params$downstream
  )
  valid_ids <- subsetByOverlaps(tss, test_ranges)$gene_id

  ## Setup the rest of the objects
  adj_method <- match.arg(adj_method)
  gs_list <- split(db$ensembl_gene, db$gs_name)
  hit_db <- dplyr::filter(db, ensembl_gene %in% valid_ids)
  hit_id_list <- split(hit_db$ensembl_gene, hit_db$gs_name)

  ## Run rGREAT
  great_res <- great(
    gr = test_ranges, gene_sets = gs_list, extended_tss = tss,
   cores = threads, background = bg
  )

  rm <- c("fold_enrichment", "p_adjust", "p_value")
  tbl <- getEnrichmentTable(great_res)
  tbl <- tibble::as_tibble(tbl[!names(tbl) %in% rm])
  tbl <- tbl[tbl$observed_gene_hits >= min_sig,]
  tbl <- dplyr::arrange(tbl, p_value_hyper, desc(fold_enrichment_hyper))
  tbl$adj_p <- p.adjust(tbl$p_value_hyper, adj_method)
  tbl$gene_id <- hit_id_list[tbl$id]
  tbl$genes_with_hits <- vapply(
    tbl$gene_id,
    \(x) paste(sort(unique(id2gene[x])), collapse = "; "),
    character(1)
  )

  ## Return with the columns needed
  names(tbl) <- gsub("_hyper", "", names(tbl))
  dplyr::select(
    tbl, gs_name = id, genome_fraction, observed_region_hits, mean_tss_dist,
    gene_set_size, observed_gene_hits, fold_enrichment, genes_with_hits,
    gene_id, p = p_value, adj_p
  )

}

#' Setup the columns for enrichment tables
enrich_cols <- list(
  gs_name = colDef(
    "GeneSet", minWidth = 180,
    cell = function(value) htmltools::tags$a(
      href = gs_url[[value]],
      target = "_blank",
      str_replace_all(value, "_", " ")
    ),
    html = TRUE
  ),
  genome_fraction = colDef(show = FALSE),
  observed_region_hits = colDef(show = FALSE),
  mean_tss_dist = colDef(
    "Mean TSS Distance (kb)", cell = \(value) round(value / 1e3, 2),
    filterMethod = js_greater
  ),
  gene_set_size = colDef(
    "Gene Set Size", maxWidth = 100, filterMethod = js_greater
  ),
  observed_gene_hits = colDef(
    "Gene Hits", maxWidth = 100, filterMethod = js_greater
  ),
  fold_enrichment = colDef(
    name = "Fold Enrichment", format = colFormat(digits = 3),
    minWidth = 110, filterMethod = js_greater
  ),
  genes_with_hits = colDef(
    name = "Genes With Associated Peaks", minWidth = 200,
    cell = \(value) with_tooltip(value, width = 60)
  ),
  gene_id = colDef(show = FALSE),
  p = colDef(show = FALSE),
  adj_p = colDef(
    name = glue("P<sub>adj</sub>"),
    html = TRUE, cell = \(value) sprint_pval(value),
    maxWidth = 110, filterMethod = js_less
  )
)


#' Produce the plain text describing the enrichment strategy...
desc_enrichment <- function(method = c("great", "gene_id")){
  method <- match.arg(method)

  list(
    great = "The package `rGREAT` [@Gu2023-rp] was used to test for enrichment of pathways
  and genesets, as an R implementation of the GREAT methodology [@McLean2010-gw].
  This approach takes a set of 'background' genomic ranges and assesses relative
  enrichment of matches to a given pathway amongst the provided peaks.

  All `rGREAT` results were taken from the hypergeometric tests returned by the
  algorithm.
  ",
    gene_id = "Enrichment testing was performed using hypergeometric approaches using gene-ids
  mapped to peaks, instead of the ranges themselves. This allowed for custom
  approaches when assigning genes to peaks as likely regulatory targets.
  A set of background gene-ids is used as the reference set to determine the
  significance of any enrichment.
  "
  )[[method]]

}
