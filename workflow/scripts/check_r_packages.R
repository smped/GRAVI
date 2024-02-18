# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log)

cat(.libPaths(),"\n")

## All installed packages
all_inst <- rownames(installed.packages())
stopifnot(length(system2('which', 'egrep', stdout = TRUE)) > 0)

## All required packages from modules
wd <- here::here()
cmd <- paste0("'library\\(.+\\)' ", file.path(wd, "workflow/modules/*Rmd"))
reqd_mod <- system2("egrep", cmd, stdout = TRUE)

## All required packages from scripts
cmd <- paste0("'library\\(.+\\)' ", file.path(wd, "workflow/scripts/*R"))
reqd_scripts <- system2("egrep", cmd, stdout = TRUE)

## Both sets of packages
all_reqd <- sort(
  unique(
    gsub(".+library\\((.+)\\)", "\\1", c(reqd_mod, reqd_scripts))
  )
)

## Check the BSgenome.<sp>.UCSC.<build> package, required for regioneR
## This can also be used for spitting out motifs as .fa files for MEME
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
config <- slot(snakemake, "config")
bld <- match.arg(tolower(config$genome$build), names(map))
ucsc_build <- map[[bld]]
sp <- map_sp[[ucsc_build]]
pkg <- paste(c("BSgenome", sp, "UCSC", ucsc_build), collapse = ".")
all_reqd <- c(all_reqd, pkg)

## Any missing packages
not_installed <- setdiff(all_reqd, all_inst)
cat(length(all_reqd), " packages are required\n")
cat(length(not_installed), " package(s) not found\n")
if (length(not_installed)) {
	cat(
    "Attempting to install:\n",
    paste(not_installed, collapse = "\n"), "\n"
  )
	BiocManager::install(not_installed, update = FALSE, force = FALSE)
}

## Check the updated list
all_inst <- rownames(installed.packages())
not_installable <- setdiff(all_reqd, all_inst)
if (length(not_installable))
	stop("unable to install:\n", paste(not_installable, collapse = "\n"))
cat("All required packages have been installed\n")

## Set the minimum version for extraChIPs
updateEC <- packageVersion("extraChIPs") < "1.7.1"
if (updateEC) {
  cat("Updating extraChIPs to a suitable version\n")
  BiocManager::install(
    "smped/extraChIPs", ref = "devel", update = FALSE, force = FALSE
  )
  stopifnot(packageVersion("extraChIPs") >= "1.7.1")
}

cat("Writing check file\n")
outfile <- slot(snakemake, "output")[[1]]
writeLines(all_reqd, outfile)
