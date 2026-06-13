# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}
cat_time <- function(...){
  tm <- format(Sys.time(), "%Y-%b-%d %H:%M:%S\t")
  cat(tm, ..., "\n")
}

log <- slot(snakemake, "log")[[1]]
cat("Setting stdout to ", log, "\n")
sink(log)

cat("libPaths set to\n", .libPaths(),"\n\n")
config <- slot(snakemake, "config")

all_params <- slot(snakemake, "params")

## All installed packages
all_inst <- rownames(installed.packages())
stopifnot(length(system2('which', 'egrep', stdout = TRUE)) > 0)

## All required packages from modules
wd <- here::here()
cmd <- paste0("'library\\(.+\\)' ", file.path(wd, "workflow/modules/*Rmd"))
reqd_mod <- system2("egrep", cmd, stdout = TRUE)

## All required packages from scripts
cmd <- paste0("'library\\(.+\\)' ", file.path(wd, "workflow/scripts/*R"))
cmd2 <- paste0("'library\\(.+\\)' ", file.path(wd, "scripts/*R"))
reqd_scripts <- system2("egrep", c(cmd, cmd2), stdout = TRUE)

## Both sets of packages
all_reqd <- sort(
  unique(
    gsub(".+library\\((.+)\\)", "\\1", c(reqd_mod, reqd_scripts))
  )
)
all_reqd <- all_reqd[!grepl("character.only", all_reqd)]

## Check the BSgenome.<sp>.UCSC.<build> package, required for regioneR
source(here::here("workflow/scripts/get_ucsc.R"))
ucsc <- get_ucsc(config$genome$build)
bs_pkg <- paste(c("BSgenome", ucsc$sp, "UCSC", ucsc$build), collapse = ".")
all_reqd <- c(all_reqd, bs_pkg)

## SimpleUpset is required for this version of extraChIPs & it's only on CRAN
all_reqd <- c(all_reqd, "SimpleUpset")

## Any missing packages
not_installed <- setdiff(all_reqd, all_inst)
cat_time(length(all_reqd), " packages are required")
cat_time(length(not_installed), " package(s) not found")
if (length(not_installed)) {
	cat_time(
    "Attempting to install:\n", paste(not_installed, collapse = "\n")
  )
	BiocManager::install(not_installed, update = FALSE, force = FALSE)
}

## Check the updated list
all_inst <- rownames(installed.packages())
not_installable <- setdiff(all_reqd, all_inst)
if (length(not_installable))
	stop("unable to install:\n", paste(not_installable, collapse = "\n"))
cat_time("All required packages have been installed")



cat_time("Writing check file")
outfile <- slot(snakemake, "output")[[1]]
writeLines(all_reqd, outfile)
