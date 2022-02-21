#' @param x, y GenomicRanges objects
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param compress If any columns are to be returned as List columns, return
#' these as a CompressedList. If set to FALSE a SimpleList will be returned for
#' these columns
#' @importFrom GenomicRanges setdiff findOverlaps
setdiffWithMcols <- function(x, y, ignore.strand = FALSE, compress = TRUE, ...) {
  #####################################################################
  ## This should be setup as an S4 method for GRanges, GRangesList & ##
  ## maybe GInteractions objects?                                    ##
  #####################################################################
  if (ncol(mcols(x)) == 0) return(setdiff(x, y, ignore.strand))

  gr <- setdiff(x, y, ignore.strand)
  hits <- findOverlaps(gr, x, ignore.strand = ignore.strand)
  i <- queryHits(hits)
  j <- subjectHits(hits)

  ## If mcols only has one column, a vector will be returned so this handles
  ## the multi-column and single-column case
  DF <- DataFrame(mcols(x)[j,])
  DF <- setNames(DF, names(mcols(x)))

  ## Return columns as Compressed/SimpleList objects if the ranges map back
  ## to multiple ranges in the original object
  needs_list <- any(duplicated(i))
  if (needs_list) {
    ## This will work for mcols which are not lists. Need to figure an approach
    ## out which handles columns which are already lists!!!
    DF <- lapply(DF, splitAsList, f = as.factor(i))
    DF <- lapply(DF, setNames, nm = c())
    if (compress) DF <- lapply(DF, as, "CompressedList")
    ## Check for any columns which can be unlisted
    cols_to_unlist <- vapply(
      DF,
      function(x) {
        l <- vapply(x, function(y) length(unique(y)), integer(1))
        all(l == 1)
      },
      logical(1)
    )
    if (any(cols_to_unlist)) DF[cols_to_unlist] <- lapply(
      DF[cols_to_unlist], unlist
    )
  }
  mcols(gr) <- DF
  gr

}
all_gr$gene[1:40] %>%
  select(gene_id, gene_name) %>%
  setdiffWithMcols(all_gr$exon)

