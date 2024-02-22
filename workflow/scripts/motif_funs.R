#' Functions for working with PWMs and sequences
#'
#' Find all PWM matches within an XStringSet
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param min_score The minimum score to return a match
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param with_score Return scores for each match
#' @param ... Passed to \link[Biostrings]{matchPWM}
#'
#' @return A DataFrame with columns: `seq`, `start`, `end`, `width`, `match`,
#' `score` and `direction`
#' The first four columns indicate the position of the match within the query
#' stringset, whilst `match` returns the matching sequence as an `XStringSet`.
#' The `direction` column denotes `F` and `R` as a forward or reverse match
#' respectively
#'
library(Biostrings)
getMatches <- function(
    pwm, stringset, rc = TRUE, min_score = "80%", with_score = TRUE,  ...
){

  ## Checks
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
  stopifnot(is(stringset, "XStringSet"))

  # Form a map from the original sequences to a views object
  map <- data.frame(end = cumsum(width(stringset)), width = width(stringset))
  map$start <- map$end - map$width + 1
  map$names <- names(stringset)

  # Form the entire XStringSetList into a Views object
  views <- Views(
    unlist(stringset), start = map$start, width = map$width, names = map$names
  )
  hits <- matchPWM(
    pwm, views, min.score = min_score, with.score = with_score, ...
  )
  mcols(hits)$direction <- rep_len("F", length(hits))
  if (rc) {
    hits_rev <- matchPWM(
      reverseComplement(pwm), views, min.score = min_score,
      with.score = with_score, ...
    )
    mcols(hits_rev)$direction <- rep_len("R", length(hits_rev))
    hits <- c(hits, hits_rev)
  }
  if (length(hits) == 0) {
    cols <- c(
      "seq", "start", "end", "width", "match", "score", "direction",
      "from_centre"
    )
    out <- setNames(lapply(cols, \(x) integer()), cols)
    out$seq <- character()
    out$match <- DNAStringSet()
    out$direction <- factor(levels = c("F", "R"))
    out[c("score", "from_centre")] <- lapply(
      out[c("score", "from_centre")], \(x) numeric()
    )
    return(DataFrame(out))
  }
  ## Map back to the original Views
  hits_to_map <- findInterval(start(hits), map$start)
  w <- width(hits)
  ## Setup the output
  hits_df <- list()
  hits_df$seq <- hits_to_map
  if (!is.null(names(stringset))) {
    hits_df$seq <- factor(map$names[hits_to_map], map$names)
    hits_df$seq <- droplevels(hits_df$seq)
  }
  hits_df$start <- as.integer(start(hits) - c(0, map$end)[hits_to_map])
  hits_df$end <- as.integer(hits_df$start + w - 1)
  hits_df$width <- w
  hits_df$match <- as(hits, "XStringSet")
  to_rev <- mcols(hits)$direction == "R"
  hits_df$match[to_rev] <- reverseComplement(hits_df$match[to_rev])
  if (with_score) hits_df$score <- mcols(hits)$score
  hits_df$direction <- factor(mcols(hits)$direction, levels = c("F", "R"))
  o <- order(hits_df$seq, hits_df$start)
  if (is.factor(hits_df$seq)) hits_df$seq <- as.character(hits_df$seq)

  ## Add the distance from the centre
  seq_widths <- width(stringset[hits_df$seq])
  hits_df$from_centre <- 0.5 * (hits_df$start + hits_df$end) - seq_widths / 2

  ## The final object
  DataFrame(hits_df)[o,]
}

#' Get the best matching PWF position with a sequence
#'
#' @details
#' The best match to the PWF within each sequence is returned, using the highest
#' score as the criteria determining the best match.
#' In the case of matches with tied scores, the method of choosing the match
#' can be selected, with 'random' being the most highly recommended.
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param ties Choose how to resolve matches with tied scores
#' @param ... Passed to \link[Biostrings]{matchPWM}
getBestMatch <- function(
    pwm, stringset, rc = TRUE, min_score = "80%",
    ties = c("random", "first", "last", "central"),  ...
) {
  ## Get the matches
  matches <- getMatches(pwm, stringset, rc, min_score, with_score = TRUE, ...)
  if (nrow(matches) == 0) return(matches)

  ## Decide how to break ties
  ties <- match.arg(ties)
  if (ties == "first") pos <- matches$start
  if (ties == "last") pos <- width(stringset[matches$seq]) - matches$end
  if (ties == "central") {
    pos <- 0.5 * (matches$start + matches$end) - 0.5 * width(stringset[matches$seq])
  }
  if (ties == "random") pos <- sample.int(nrow(matches), nrow(matches))
  coerce_fact <- is.character(matches$seq)
  if (coerce_fact)
    matches$seq <- factor(matches$seq, levels = unique(matches$seq))

  ## Pick the best match using the values we have
  o <- order(matches$seq, 1 / matches$score, pos)
  matches <- matches[o,]
  matches <- matches[!duplicated(matches$seq), ]
  if (coerce_fact) matches$seq <- as.character(matches$seq)

  # ## Add the distance from the centre, taking the mean of start/end
  # seq_widths <- width(stringset[matches$seq])
  # matches$from_centre <- 0.5 * (matches$start + matches$end) - seq_widths / 2
  matches

}
#' Test for a Uniform Distribution across a set of best matches
#'
#' @details
#' This function tests for an evenness of motif matches across a set of
#' sequences, using the assumption that if there is no positional bias, then
#' matches should be evenly spread across all positions within a set of sequences.
#' Conversely, if there is positional bias, typically but not neccessarily near
#' the centre of a range, this function intends to detect this signal.
#'
#' Firstly, the best match within a sequence is found.
#' The positions of all best matches within a set of sequences are then grouped
#' into bins of size specified, and a binomial test performed on each bin.
#' This tests for a number of successes matching the expected probability.
#' By default, the hypothesis test is set to detect bins with greater than
#' expected 'success'.
#' The p-values obtain across all bins are then merged into a single p-value
#' using the Harmonic mean p-value to obtain a representative p-value across
#' the set of sequences.
#' The bin with the lowest raw p-value is returned as the most enriched position.
#'
#' @return
#' A data.frame with columns `start`, `end`, `total_matches`, `matches_in_bin`,
#' `expected`, `enrichment`, and `p`
#' The total matches represent the total number of matches within the set of
#' sequences, whilst the number observed in the specific bin are also given.
#'
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#' @param ties Choose how to resolve matches with tied scores
#' @param binwidth Width of bins across the range to group data into
#' @param abs Use absolute positions around zero to find symmetrical enrichment
#' @param alt Alternative hypothesis for the binomial test
#' @param ... Passed to \link[Biostrings]{matchPWM}
testUniform <- function(
    pwm, stringset, rc = TRUE, min_score = "80%",
    ties = c("random", "first", "last", "central"), binwidth = 10, abs = FALSE,
    alt = c("greater", "two.sided", "less"), ...
) {
  w <- unique(width(stringset))
  if (length(w) > 1) stop("All queried sequences must be the same width")
  bm <- getBestMatch(pwm, stringset, rc, min_score, ties)
  n_matches <- nrow(bm)
  if (n_matches == 0) return(DataFrame())

  all_breaks <- c(
    seq(-w/2, w/2, by = binwidth) - binwidth / 2, w / 2 + binwidth / 2
  )
  if (abs) {
    ## Set all positions as absolute, i.e. distance from centre
    bm$from_centre <- abs(bm$from_centre)
    all_breaks <- seq(0, w/2, by = binwidth / 2)
  }
  bm$bin <- cut(bm$from_centre, breaks = all_breaks, include.lowest = TRUE)
  ## Run the analysis using the range as a single range
  df <- summarise(
    as.data.frame(bm[,c("from_centre", "bin")]),
    matches_in_bin = dplyr::n(),
    start = min(from_centre), end = max(from_centre),
    .by = bin
  )
  df$width <- df$end - df$start + 1
  n_bins <- nrow(df)
  alt <- match.arg(alt)
  df$expected <- n_matches * df$width / sum(df$width)
  df$p <- map_dbl(
    seq_len(n_bins),
    \(i) {
      x <- as.list(df[i,])
      binom.test(x$matches_in_bin, n_matches, x$width / sum(df$width), alt)$p.value
    }
  )
  n <- nrow(df)
  hmp <- extraChIPs:::.ec_HMP(df$p, rep_along(df$p, 1 / n))
  df <- subset(df, p < hmp | p == min(p))
  df <- summarise(
    df,
    start = all_breaks[min(as.integer(bin))],
    end = all_breaks[max(as.integer(bin)) + 1], n_bins = dplyr::n(),
    matches_in_bin = sum(matches_in_bin), expected = sum(expected)
  )
  df$total_matches <- n_matches
  df$enrichment <- df$matches_in_bin / df$expected
  df$width <- df$end - df$start
  df$centre <- df$start + 0.5 * df$width
  cols <- c(
    "start", "end", "width", "centre",
    "n_bins", "matches_in_bin", "total_matches", "expected", "enrichment", "hmp"
  )
  df$hmp <- hmp
  df[, cols]

}

## The following tests show potential.
## The next strategy is to use `nullranges` to select multiple sets of ranges
## which match the overlapping features, permute 1000-500 times for a viable
## Null distribution & Z-scores for general motif enrichment.

all_unif_test <- merged_filtered_db %>%
  mclapply(
    \(x) {
      pwm <- x@motif
      df <- testUniform(pwm, seq, binwidth = 40)
      if (length(df) == 0) df <- list()
      df$name <- x@name
      df$altname <- x@altname
      as.data.frame(df)
    },
    mc.cores = 4
  ) %>%
  bind_rows() %>%
  as_tibble() %>%
  dplyr::select(ends_with("name"), everything()) %>%
  dplyr::filter(!is.na(hmp)) %>%
  mutate(
    adj_p = p.adjust(hmp, "fdr")
  ) %>%
  arrange(hmp)

all_abs_test <- merged_filtered_db %>%
  mclapply(
    \(x) {
      pwm <- x@motif
      df <- testUniform(pwm, seq, binwidth = 40, abs = TRUE)
      if (length(df) == 0) df <- list()
      df$name <- x@name
      df$altname <- x@altname
      as.data.frame(df)
    },
    mc.cores = 4
  ) %>%
  bind_rows() %>%
  as_tibble() %>%
  dplyr::select(ends_with("name"), everything()) %>%
  dplyr::filter(!is.na(hmp))%>%
  mutate(adj_p = p.adjust(hmp, "fdr")) %>%
  arrange(hmp)



