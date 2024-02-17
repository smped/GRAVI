#' Combines two steps in separating strings and capitalising
#'
#' @param x a xharaxter vector
#' @param pattern the separator to replace
#' @param replacement The replacement pattern
#' @param exclude Do not operate on these coplete strings
#'
#' @return a character vector
str_sep_to_title <- function(
  x, pattern = "[^[:alnum:]]+", replacement = " ",
  exclude = c("DNA", "RNA")
) {
  out <- x
  out[!x %in% exclude] <- out[!x %in% exclude] %>%
    stringr::str_replace_all(pattern, replacement) %>%
    stringr::str_to_title()
  out
}

#' Plot a simple colour scheme as a simple bar
#'
#' @param x A vector of **named** colours
#' @param xlab,ylab Default names for axis labels
#'
#' @return a ggplot2 object
.plotScheme <- function(x, xlab = "", ylab = "") {
  x %>%
    as_tibble() %>%
    pivot_longer(cols = everything()) %>%
    mutate(name = factor(name, levels = names(x))) %>%
    ggplot(
      aes(name, 1, fill = name)
    ) +
    geom_raster() +
    scale_fill_manual(values = unlist(x)) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_y_continuous(expand = expansion(c(0, 0))) +
    labs(x = xlab, y = ylab) +
    guides(fill = "none") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
}

#' Create a tbl_graph object from goseq results
#'
#' @details
#' Distances between genesets are created by using the overlap-coefficient
#'
#' @param res data.frame with goseq results
#' @param gs List of gene-sets with gene ids with each element. Names are gene-sets
#' @param p_col Column name with (adjusted) p-values
#' @param min_dist remove any edges with distances < min_dist
#' @param alpha Threshold for filtering goseq results
#' @param max_gs The maximum number of genesets to draw
#'
make_tbl_graph <- function(
    res, gs, p_col = "adj_p", gs_col = "gs_name",
    hit_col = "numDEInCat", size_col = "numInCat",
    alpha = enrich_alpha,
    min_dist = max_network_dist,
    max_gs = max_network_size
) {
  p_col <- match.arg(p_col, colnames(res))
  gs_col <- match.arg(gs_col, colnames(res))
  res <- dplyr::filter(res, !!sym(p_col) < alpha)
  res <- arrange(res, !!sym(p_col))
  ## If the gene sets are identical, retain the first (most-significant) only
  gs <- gs[res[[gs_col]]]
  nm <- gs %>%
    map_chr(paste, collapse = "; ") %>%
    .[!duplicated(.)] %>%
    names() %>%
    .[seq_len(min(c(length(.), max_gs)))]
  if (length(nm) < 2) return(tbl_graph())
  combs <- combn(nm, 2)
  d <- lapply(
    seq_len(ncol(combs)),
    function(i) {
      pair <- combs[,i]
      all <- unlist(gs[pair])
      n <- min(vapply(gs[pair], length, integer(1)))
      ## Use the overlap coefficient
      ol <- sum(table(all) == 2) / n
      1 - ol
    }
  )
  edges <- combs %>%
    t() %>%
    set_colnames(c("from", "to")) %>%
    as_tibble() %>%
    mutate(d = unlist(d), oc = 1 - d) %>%
    dplyr::filter(d <= min_dist)
  ## Subsume any with d == 0 into the more highly ranked
  ## This needs to happen sequentially
  omit <- c()
  # continue <- any(edges$d == 0)
  # while (continue) {
  #   omit <- c(omit, edges$to[which(edges$d == 0)[1]])
  #   edges <- dplyr::filter(edges, !from %in% omit, !to %in% omit)
  #   continue <- any(edges$d == 0)
  # }
  nodes <- tibble(label = setdiff(nm, omit)) %>%
    left_join(res, by = c("label" = gs_col)) %>%
    mutate(prop = !!sym(hit_col) / !!sym(size_col))
  node_ids <- setNames(seq_along(nodes$label), nodes$label)
  edges$from <- node_ids[edges$from]
  edges$to <- node_ids[edges$to]
  tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
}


#' Run an enrichment analysis using the Hypergeometric Distribution
#'
#' @details
#' Stripped down version of \link[goseq](goseq) with only the Hypergeomtric
#' test implemented. The original output column over_represented_pvalue has
#' been renamed to simply pval. The output is also a tibble as opposed to
#' a standard `data.frame`
#'
#' @param pwf a Simple data.frame as would be produced by \link[goseq](nullp).
#' Rownames must be gene ids and the column `DEgenes` is mandatory.
#' @param gene2cat List of gene_ids with pathway mappings included for each element
#' @param cores The number of threads to use when running in parallel
#'
#' @importFrom parallel mclapply
#'
goseq_hyper <- function (pwf, gene2cat = NULL, cores = threads) {
  ## Stripped down version of goseq, fixed to the Hypergeometric Distribution
  ## The pwf should be matrix/df with column "DEgenes" and rownames as gene_ids
  ## The columns bias.data & pwf are ignored for this stripped down version
  ##
  ## Should really be rewritten to match values using slightly more readable code
  ## All we really need is a vector of genes for each geneset
  stopifnot("DEgenes" %in% colnames(pwf))
  mapped_genes <- intersect(rownames(pwf), names(gene2cat))
  gene2cat <- gene2cat[mapped_genes]
  cat2gene <- goseq:::reversemapping(gene2cat)
  cat2gene <- cat2gene[vapply(cat2gene, length, integer(1)) > 0]
  cat2gene <- lapply(cat2gene, unique)
  de <- intersect(rownames(subset(pwf, DEgenes == 1)), mapped_genes)
  n_de <- length(de)
  n_genes <- length(mapped_genes)

  f <- lapply
  if (cores > 1){
    f <- parallel::mclapply
    formals(f)$mc.cores <- cores
  }

  list_res <- f(
    cat2gene,
    function(x){
      df <- tibble(
        numDEInCat = sum(x %in% de),
        numInCat =  length(x),
        pval = dhyper(numDEInCat, numInCat, n_genes - numInCat, n_de) +
          phyper(numDEInCat, numInCat, n_genes - numInCat, n_de, lower.tail = FALSE),
        DEgenes = list(intersect(x, de))
      )
      df
    }
  )
  df <- bind_rows(list_res, .id = "gs_name")
  df <- arrange(df, pval)
  dplyr::select(df, gs_name, pval, everything())

}
#' Convert any supplied references to the UCSC format needed by rGREAT
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
