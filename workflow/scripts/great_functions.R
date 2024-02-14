#' Map standard msigdb categories to the format required for rGREAT
#'
#' @param x vector of categories to convert
#' @param pre Prefix to add. See the vignette
#'
#' @return Character vector
#'
map_great_cats <- function(x = NULL, pre = "msigdb:") {
  ## The gs_cat entries are complete compatible. Only the subcats need to be mapped
  df <- msigdbr::msigdbr_collections()
  cats <- unique(dplyr::filter(df, gs_cat %in% x)$gs_cat)
  sub_cats <- dplyr::filter(df, !gs_cat %in% cats, gs_subcat %in% x)
  sub_cats <- paste(sub_cats$gs_cat, sub_cats$gs_subcat, sep = ":")
  paste0(pre, c(cats, sub_cats))
}
#' Map common references to the options set by rGREAT
#' @param x Genome. For example "GRCh37"
#' @return The UCSC equivalent
map_great_refs <- function(x = NULL){
  map <- c(
    hg19 = "hg19", hg38 = "hg38", grch37 = "hg19", grch38 = "hg38",
    mm10 = "mm10", mm39 = "mm39", grcm38 = "mm10", grcm39 = "mm39",
    rn7 = "rn7", mratbn7.2 = "rn7", galgal6 = "galGal6", rhemac10 = "rheMac10",
    canfam5 = "canFam5", susscr11 = "susScr11", pantro6 = "panTro6", dm6 = "dm6"
  )
  x <- match.arg(tolower(x), names(map))
  map[[x]]
}
#' The parameters are the same as rGREAT::great with the following additions
#' @param tss Must be prepared using extendTSS!!!
#' @param max_gene_set_size Remove gene sets which are too large and likely
#' @param adj Allows flexibility in the p-values adjustment method
#' uninformative
#' @param min_hits Remove gene sets with fewer than this number of hits. Gene
#' sets are removed after adjusting p-values
#' @param great_cores Passed to each gene-set level run of `great()`
multi_set_great <- function(
    gr, gene_sets, tss, min_gene_set_size = 5, max_gene_set_size = Inf,
    basal_upstream = 5000, basal_downstream = 1000, extension = 1e+06,
    background = NULL, exclude = "gap", great_cores = 1,
    adj = p.adjust.methods, min_hits = 1, verbose = TRUE
){
  adj <- match.arg(adj)
  gl <- lapply(
    gene_sets,
    \(x) {
      great(
        gr, x, extended_tss = tss, min_gene_set_size = min_gene_set_size,
        mode = "basalPlusExt", basal_upstream = basal_upstream,
        basal_downstream = basal_downstream, extension = extension,
        background = background, exclude = exclude, cores = great_cores,
        verbose = verbose
      )
    }
  )

  ## Base the output on the first element, then add the new objects
  slot_names <- slotNames(gl[[1]])
  new_obj <- lapply(slot_names, \(i) slot(gl[[1]], i))
  names(new_obj) <- slot_names

  ## Tidy up the table
  tbl <- bind_rows(lapply(gl, slot, "table"))
  tbl <- dplyr::filter(tbl, gene_set_size < max_gene_set_size)
  tbl$p_adjust <- p.adjust(tbl$p_value, adj)
  tbl$p_adjust_hyper <- p.adjust(tbl$p_value_hyper, adj)
  tbl <- dplyr::filter(tbl, observed_gene_hits >= min_hits)
  tbl <- arrange(tbl, p_value)

  ## Form the remaining new objects
  sets <- do.call(c, unname(lapply(gl, slot, "gene_sets")))
  sets <- sets[names(sets) %in% tbl$id]
  set_names <- as.character(unlist(gene_sets))

  # And return the merged object
  new_obj$table <- as.data.frame(tbl)
  new_obj$gene_sets <- sets
  new_obj$gene_sets_name <- set_names
  do.call("GreatObject", new_obj)
}
