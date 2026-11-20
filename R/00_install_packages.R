# OPTIONAL, ONE-TIME SETUP. Run only if the required packages are missing.
# This may download and install software. It does not run the RNA analysis.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c("muscData", "edgeR"), update = FALSE, ask = FALSE)
cat("Setup finished. The analysis entry point is R/03_rerun_checked.R.
")
