args <- commandArgs(TRUE)
f <- args[[1]]
reqd_version <- args[[2]]

cat("Checking for extraChIPs >= v", reqd_version, "\n")
inst_version <- packageVersion("extraChIPs")
cat("Installed version is ", as.character(inst_version), "\n")

needs_update <- inst_version < reqd_version
cat(ifelse(needs_update, "Updating...", "No update required"))
if (needs_update) BiocManager::install("extraChIPs", ask = FALSE)

success <- packageVersion("extraChIPs") >= reqd_version
stopifnot(success)

file.create(f)
