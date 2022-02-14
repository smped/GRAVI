#' @title Auto-Detect File-Type and Read
#'
#' @description Automatically detects a file-type and imports
#'
#' @details
#' This will use the file extension to the automatically detect a file-type
#' and will load accordingly.
#' The currently implemented file types are tsv/txt, csv, rds, ged, gff
#'
#' @param file Path to a file
#' @param ... Passed to all import functions. Unless all files are the same format, this may cause issues
#'
#' @import readr
#' @importFrom rlang call2 !!!
autoload <- function(file, ...) {
	stopifnot(file.exists(file))
	n <- length(file)
	out <- lapply(file, .autoload_single, ...)
	names(out) <- file
  if (n == 1) out <- out[[1]]
	out
}

.autoload_single <- function(file, ...){

  stopifnot(length(file) == 1)
  fname <- basename(file)

  fun <- case_when(
    str_detect(fname, "tsv|txt") ~ "read_tsv",
    str_detect(fname, "csv") ~ "read_csv",
    str_detect(fname, "rds") ~ "read_rds",
    str_detect(fname, "bed") ~ "import.bed",
    str_detect(fname, "gtf") ~ "import.gff",
    TRUE ~ NA_character_
  )
  if (is.na(fun)) stop("File type not detected")
  if (str_detect(fun, "read_")) {
    args <- c(
        list(file = file),  list(...)
      )
  }
  if (str_detect(fun, "import")) {
      args <- c(
      list(con = file),  list(...)
    )
  }

  f <- rlang::call2(fun, !!!args)
  eval(f)

}
