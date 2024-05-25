#' Get UCSC Genome IDs from any build
#' Every internal function here shold be from base
get_ucsc <- function(build) {
  map <- c(
    hg19 = "hg19", hg38 = "hg38", grch37 = "hg19", grch38 = "hg38",
    mm10 = "mm10", mm39 = "mm39", grcm38 = "mm10", grcm39 = "mm39",
    rn7 = "rn7", mratbn7.2 = "rn7", galgal6 = "galGal6",
    # rhemac10 = "rheMac10", canfam5 = "canFam5", # Not available as a BSgenome
    susscr11 = "susScr11", pantro6 = "panTro6", dm6 = "dm6"
  )
  map_sp <- c(
    hg19 = "Hsapiens", hg38 = "Hsapiens", mm10 = "Mmusculus", mm39 = "Mmusculus",
    rn7 = "Rnorvegicus", galGal6 = "Ggallus", susScr11 = "Sscrofa",
    panTro6 = "Ptroglodytes", dm6 = "Dmelanogaster"
  )
  bld <- match.arg(tolower(build), names(map))
  ucsc_build <- map[[bld]]
  sp <- map_sp[[ucsc_build]]
  list(build = ucsc_build, sp = sp)
}
