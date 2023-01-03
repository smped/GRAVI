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
#' @param min_dist remove any edges with distances < min_dist
#' @param alpha Threshold for filtering goseq results
#'
make_tbl_graph <- function(res, gs, min_dist = 0.9, alpha = enrich_alpha) {
  res <- dplyr::filter(res, adj_p < alpha)
  if (nrow(res) < 2) return(tbl_graph())
  nm <- res$gs_name
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
    mutate(d = unlist(d)) %>%
    dplyr::filter(d < min_dist)
  ## If there is a distance of zero, take the most highly ranked
  omit <- dplyr::filter(edges, d == 0)$to
  nodes <- tibble(id = seq_along(nm), label = nm) %>%
    dplyr::filter(!id %in% omit) %>%
    left_join(res, by = c("label" = "gs_name")) %>%
    mutate(prop = numDEInCat / numInCat)
  edges <- dplyr::filter(edges, d > 0)
  node_ids <- setNames(nodes$id, nodes$label)
  edges$from <- node_ids[edges$from]
  edges$to <- node_ids[edges$to]
  tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
}
