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
