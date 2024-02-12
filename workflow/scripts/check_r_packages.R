# Handle any conda weirdness
conda_pre <- system2("echo", "$CONDA_PREFIX", stdout = TRUE)
if (conda_pre != "") {
  conda_lib_path <- file.path(conda_pre, "lib", "R", "library")
  if (!dir.exists(conda_lib_path)) conda_lib_path <- NULL
  prev_paths <- .libPaths()
  paths_to_set <- unique(c(conda_lib_path, prev_paths))
  .libPaths(paths_to_set)
}

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
all_reqd <- unique(gsub(".+library\\((.+)\\)", "\\1", c(reqd_mod, reqd_scripts)))


## Any missing packages
not_installed <- setdiff(all_reqd, all_inst)
message(length(all_reqd), " packages are required")
message(length(missing), " package(s) not found")
if (length(not_installed)) {
	message("Attempting to install:\n", paste(not_installed, collapse = "\n"))
	BiocManager::install(not_installed, update = FALSE, force = FALSE)
}

## Check the updated list
all_inst <- rownames(installed.packages())
not_installable <- setdiff(all_reqd, all_inst)
if (length(not_installable))
	stop("unable to install:\n", paste(not_installable, collapse = "\n"))
message("All required packages have been installed")

## Set the minimum version for extraChIPs
updateEC <- packageVersion("extraChIPs") >= "1.7.1"
if (updateEC) {
  message("Updating extraChIPs to a suitable version")
  BiocManager::install(
    "smped/extraChIPs", ref = "devel", update = FALSE, force = FALSE
  )
  stopifnot(packageVersion("extraChIPs") >= "1.7.1")
}

message("Writing check file")
d <- here::here("output/checks")
if (!dir.exists(d)) dir.create(d, recursive = TRUE)
writeLines(all_reqd, file.path(d, "r-packages.chk"))
