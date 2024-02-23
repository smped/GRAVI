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
  ## Add the distance from the centre
  seq_widths <- width(stringset[hits_df$seq])
  hits_df$from_centre <- 0.5 * (hits_df$start + hits_df$end) - seq_widths / 2
  ## The match itself
  hits_df$match <- as(hits, "XStringSet")
  to_rev <- mcols(hits)$direction == "R"
  hits_df$match[to_rev] <- reverseComplement(hits_df$match[to_rev])
  hits_df <- c(as.list(mcols(hits)), hits_df)
  hits_df$direction <- factor(mcols(hits)$direction, levels = c("F", "R"))
  o <- order(hits_df$seq, hits_df$start)
  if (is.factor(hits_df$seq)) hits_df$seq <- as.character(hits_df$seq)

  ## The final object
  DataFrame(hits_df)[o,c("seq", setdiff(names(hits_df), "seq"))]
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
#' using the Harmonic mean p-value (HMP) to obtain a representative p-value
#' across the set of sequences.
#' All bins with p-values below the HMP are merged in the final result, with
#' start, end, width and centre indicating the range of the merged bins.
#' The column `n_bins` also indicates how many bins were merged to obtain the
#' final values.
#'
#' @return
#' A data.frame with columns `start`, `end`, `centre`, `consensus_motif`,
#' `n_bins`, `matches_in_bin`, `total_matches`, `expected`, `bin_enrichment`
#' and `hmp`
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
  df <- dplyr::summarise(
    as.data.frame(bm[,c("from_centre", "bin")]),
    matches_in_bin = dplyr::n(),
    start = min(from_centre), end = max(from_centre),
    .by = bin
  )
  df$width <- df$end - df$start + 1
  n_bins <- length(levels(bm$bin))
  alt <- match.arg(alt)
  df$expected <- n_matches * df$width / sum(df$width)
  df$p <- map_dbl(
    seq_along(df$bin),
    \(i) {
      x <- as.list(df[i,])
      binom.test(x$matches_in_bin, n_matches, x$width / sum(df$width), alt)$p.value
    }
  )
  hmp <- harmonicmeanp::p.hmp(df$p, L = n_bins)
  # hmp <- extraChIPs:::.ec_HMP(df$p, rep_along(df$p, 1 / n_bins))
  df <- subset(df, p < hmp | p == min(p))
  df <- dplyr::summarise(
    df,
    start = all_breaks[min(as.integer(bin))],
    end = all_breaks[max(as.integer(bin)) + 1], n_bins = dplyr::n(),
    matches_in_bin = sum(matches_in_bin), expected = sum(expected)
  )
  df$total_matches <- n_matches
  df$bin_enrichment <- df$matches_in_bin / df$expected
  df$width <- df$end - df$start
  df$centre <- df$start + 0.5 * df$width
  motif_cols <- strsplit(consensusString(bm$match), "")[[1]]
  consensus <- consensusMatrix(bm$match)[rownames(pwm),]
  colnames(consensus) <- motif_cols
  df$consensus_motif <- list(consensus)
  cols <- c(
    "start", "end", "width", "centre", "n_bins", "matches_in_bin",
    "total_matches", "expected", "bin_enrichment", "hmp", "consensus_motif"
  )
  df$hmp <- hmp
  df[, cols]

}
#' Prepare a set of matching ranges for bootstrapping or permuting
#'
#' Prepare a set of ranges from y which exactly match those in x
#'
#' @param x GRanges with ranges to be matched
#' @param y GRanges with ranges to select random matching ranges from
#' @param exclude GRanges of ranges to omit from testing
#' @param n_iter The number of times to repeat the random selection process
#' @param replace Sample with our without replacement from the random ranges
#' @param ... Not used
setGeneric(
  "matchRandomRanges",
  function(x, y, ...) standardGeneric("matchRandomRanges")
)
setMethod(
  "matchRandomRanges", signature = signature(x = "GRanges", y = "GRanges"),
  function(
    x, y, exclude = GRanges(), n_total = NULL, n_iter = NULL, replace = TRUE, ...
  ) {

    # Make sure any ranges don't fall off the end of chromosomes
    stopifnot(is(exclude, "GRanges"))
    sq_gr <- GRanges(seqinfo(x))
    w <- width(x)
    chr_ends <- c(
      resize(sq_gr, width = max(w) / 2, fix = "start"),
      resize(sq_gr, width = max(w) / 2, fix = "end")
    )
    exclude <- c(granges(exclude), chr_ends)
    y <- setdiff(y, exclude)

    ## Set the iterations
    if (!is.null(n_iter)) {
      iter_size <- length(x)
      if (!is.null(n_total)) message("Recalculating totals for iterations")
      n_total <- as.integer(n_iter * iter_size)
    }

    ## Now the random sampling
    gpos <- GPos(y)
    i <- sample(seq_along(gpos), n_total, replace = replace)
    rand_w <- sample(w, n_total, replace = TRUE)
    gr <- GRanges(gpos[i])
    gr <- resize(gr, width = w, fix = "center")
    if (!is.null(n_iter)) gr$iteration <- rep(seq_len(n_iter), times = iter_size)
    sort(gr)

  }
)
setMethod(
  "matchRandomRanges",
  signature = signature(x = "GRangesList", y = "GRangesList"),
  function(
    x, y, exclude = GRanges(), n_iter = 1, replace = TRUE, mc.cores = 1, ...
  ){

    stopifnot(length(x) <= length(y))
    nm_x <- names(x)
    nm_y <- names(y)
    if (!is.null(nm_x) & !is.null(nm_y)) {
      stopifnot(all(nm_x %in% nm_y))
      y <- y[nm_x]
    }
    i <- seq_along(x)
    out <- parallel::mclapply(
      i,
      \(i) matchRandomRanges(
        x[[i]], y[[i]], exclude = exclude,  n_iter = n_iter, replace = replace,
        ...
      ), mc.cores = mc.cores
    )
    if (!is.null(nm_x)) names(out) <- nm_x
    GRangesList(out)
  }
)
#' Count the matches to a PWM within an XStringSet
#'
#' Count the matches to a PWM within an XStringSet
#'
#' @details
#' Will simply count the matches within an XStringSet and return an integer
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param rc logical(1) Also find matches using the reverse complement of pwm
#' @param min_score The minimum score to return a match
#'
countPwmMatches <- function(pwm, stringset, rc = TRUE, min_score = "80%", ...) {

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
  width <- ncol(pwm)
  n_matches <- countPWM(pwm, views)
  if (rc) n_matches <- c(n_matches, countPWM(reverseComplement(pwm), views))
  as.integer(sum(n_matches))

}
#' @title Test Single Motif Enrichment Using an Iterative approach
#' @description
#' Iterate through background sequences to form a null distribution and test
#' enrichment
#'
#' @details
#' This is a computationally intensive process which steps through a set of
#' background sequences, which, ideally, are well matched to the test sequences.
#' This counts matches within each set of of iterations and builds up estimates
#' of the mean and sd for the expected counts. These are then used to calculate
#' a Z-score, under the Central Limit Theorem and a p-value derived from this
#' Z-score.
#'
#' Preliminary testing has shown that background counts follow an approximate
#' Poisson distribution and as such, the function [poisMotifEnrich()] may be
#' less computationally demanding.
#'
#' @return
#' A data.frame with columns: `matches`, `n`, `mean_bg`, `sd_bg`, `n_iter`, `Z`
#' and `p`
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param bg An XStringSet
#' @param by A column in the mcols element of bg, usually denoting an iteration
#' number
#' @param mc.cores Passed to \link[parallel]{mclapply}
#'
iterMotifEnrich <- function(
    pwm, stringset, bg, by = "iteration", mc.cores = 1, ...
) {

  ## Checks
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
  stopifnot(is(stringset, "XStringSet"))
  stopifnot(is(bg, "XStringSet"))
  stopifnot(by %in% colnames(mcols(bg)))

  ## Find the matches
  matches <- countPwmMatches(pwm, stringset, ...)
  bglist <- split(bg, mcols(bg)[[by]])
  bg_matches <- parallel::mclapply(
    bglist, \(x) countPwmMatches(pwm, x, ...), mc.cores = mc.cores
  )
  bg_matches <- as.integer(unlist(bg_matches))
  mean_bg <- mean(bg_matches)
  sd_bg <- sd(bg_matches)
  n_iter <- length(bg_matches)
  n <- length(stringset)
  Z <- (matches - mean_bg) / sd_bg# / sqrt(n))
  # perm_p <- (sum(bg_matches > matches) + 1) / n_iter
  p <- 1 - pchisq(Z^2, 1)
  data.frame(matches, n, mean_bg, sd_bg, n_iter, Z, p)#, perm_p)

}
# pwm <- merged_filtered_db %>%
#   to_df() %>%
#   subset(name == "ESR1") %>%
#   to_list() %>%
#   .[[1]] %>%
#   slot("motif")
# iterMotifEnrich(pwm, seq, bg_seq, by = "iteration", mc.cores = 4)
# # Hmm
# enrich_df <- merged_filtered_db %>%
#   lapply(
#     \(x) {
#       nm <- x@name
#       alt <- x@altname
#       df <- iterMotifEnrich(x@motif, seq, bg_seq, by = "iteration", mc.cores = 4)
#       mutate(df, name = nm, altname = alt)
#     }
#   ) %>%
#   bind_rows() %>%
#   as_tibble() %>%
#   dplyr::select(ends_with("name"), everything()) %>%
#   mutate(adj_p = p.adjust(p, "fdr"))
# ## Took about an hour for all 639

#' @title Test Single Motif Enrichment Using an Poisson Model
#' @description
#' Test enrichment after estimating parameters of a (null) Poisson Model
#'
#' @details
#' This assume that motif matches within the set of background sequences follow
#' a Poisson Distribution, and with enough background sequences, the rate
#' parameter can be accurately estimated. For confidence, the number of
#' background sequences should be far larger than the set of tested sequences.
#'
#' The observed counts of motif matches are then tested against a Poisson Model
#' using the estimated rate and a two-sided test.
#' Whilst being orders of magnitude faster than [interMotifEnrich()], large
#' background sets can still be computationally demanding
#'
#' @return
#' A data.frame with columns: `matches`, `sequences`, `expected`, `enrichment`,
#' `est_rate`, `approxZ` and `p`.
#' The enrichment column is simply the ratio between observed and expected,
#' whilst the Zscore uses the estimated Poisson rate to estimate a Zscore.
#' P-values are obtained from [poisson.test()]
#'
#' @param pwm A Position Weight Matrix
#' @param stringset An XStringSet
#' @param bg An XStringSet
#'
poisMotifEnrich <- function(pwm, stringset, bg, ...){

  ## Checks
  colsums <- round(colSums(pwm), 4) # Precision error
  if (!all(colsums == 1))
    stop("PWMs must be a matrix with columns that sum to 1")
  stopifnot(rownames(pwm) == c("A", "C", "G", "T"))
  stopifnot(is(stringset, "XStringSet"))
  stopifnot(is(bg, "XStringSet"))
  stopifnot(by %in% colnames(mcols(bg)))

  # Get the counts
  obs <- countPwmMatches(pwm, stringset, ...)
  bg_rate <- countPwmMatches(pwm, bg, ...) / length(bg)
  n <- length(stringset)
  p <- poisson.test(obs, n, bg_rate)$p.value
  expected <- bg_rate * n
  Z  <- (obs - expected) / sqrt(expected) ## Assumes a strict poisson
  data.frame(
    matches = obs, sequences = n, expected, enrichment = obs / expected,
    est_rate = bg_rate, approxZ = Z, p = p
  )
}
# poisMotifEnrich(pwm, seq, bg_seq)
# enrich_pois_df <- merged_filtered_db %>%
#   parallel::mclapply(
#     \(x) {
#       nm <- x@name
#       alt <- x@altname
#       df <- poisMotifEnrich(x@motif, seq, bg_seq,)
#       mutate(df, name = nm, altname = alt)
#     }, mc.cores = 4
#   ) %>%
#   bind_rows() %>%
#   as_tibble() %>%
#   dplyr::select(ends_with("name"), everything()) %>%
#   mutate(adj_p = p.adjust(p, "fdr"))

## The following tests show potential.
## The next strategy is to use `nullranges` to select multiple sets of ranges
## which match the overlapping features, permute 1000-500 times for a viable
## Null distribution & Z-scores for general motif enrichment.
#
# all_unif_test <- merged_filtered_db %>%
#   parallel::mclapply(
#     \(x) {
#       pwm <- x@motif
#       df <- testUniform(pwm, seq, binwidth = 40)
#       if (length(df) == 0) df <- list()
#       df$name <- x@name
#       df$altname <- x@altname
#       as.data.frame(df)
#     },
#     mc.cores = 4
#   ) %>%
#   bind_rows() %>%
#   as_tibble() %>%
#   dplyr::select(ends_with("name"), everything()) %>%
#   dplyr::filter(!is.na(hmp)) %>%
#   mutate(
#     adj_p = p.adjust(hmp, "fdr")
#   ) %>%
#   arrange(hmp)
#
# all_abs_test <- merged_filtered_db %>%
#   parallel::mclapply(
#     \(x) {
#       pwm <- x@motif
#       df <- testUniform(pwm, seq, binwidth = 40, abs = TRUE)
#       if (length(df) == 0) df <- list()
#       df$name <- x@name
#       df$altname <- x@altname
#       as.data.frame(df)
#     },
#     mc.cores = 4
#   ) %>%
#   bind_rows() %>%
#   as_tibble() %>%
#   dplyr::select(ends_with("name"), everything()) %>%
#   dplyr::filter(!is.na(hmp))%>%
#   mutate(adj_p = p.adjust(hmp, "fdr")) %>%
#   arrange(hmp)



