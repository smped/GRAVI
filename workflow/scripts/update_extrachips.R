args <- commandArgs(TRUE)
f <- args[[1]]
version <- args[[2]]

cat("Checking for extraChIPs >= v", version)
needs_update <- packageVersion("extraChIPs") >= version
cat(ifelse(needs_update, "Updating...", "No update required"))
if (needs_update) BiocManager::install("extraChIPs", ask = FALSE)

success <- packageVersion("extraChIPs") >= version
stopifnot(success)

file.create(f)
