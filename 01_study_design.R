# PROJECT 2, PART 1: understand the paired donor-condition study design.
# Source: Kang et al. (2018), GSE96583, as documented in muscData/muscat.
# This small table reconstructs the DOCUMENTED DESIGN, not expression data.
# It must later be checked against the downloaded cell/sample metadata.
# Run from PBMC_IFNB_project. No packages or RNA-data download are needed here.
if (!file.exists("PBMC_IFNB_project.Rproj")) {
  stop("Open PBMC_IFNB_project.Rproj before running this script.")
}

# These public donor IDs appear in the official muscat sample table.
donor_ids <- c("101", "1015", "1016", "1039", "107", "1244", "1256", "1488")

# Each donor contributes one control and one stimulated sample.
study_design <- data.frame(
  Donor_ID = rep(donor_ids, each = 2),
  Condition = rep(c("ctrl", "stim"), times = length(donor_ids)),
  stringsAsFactors = FALSE
)

print(study_design)
number_of_donors <- length(unique(study_design$Donor_ID))
number_of_sample_rows <- nrow(study_design)
cat("Number of donors:", number_of_donors, "\n")
cat("Number of donor-condition rows:", number_of_sample_rows, "\n")

# Counts here are sample rows, not cell counts or RNA measurements.
donor_condition_table <- table(study_design$Donor_ID, study_design$Condition)
print(donor_condition_table)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(study_design, "results/tables/documented_study_design.csv", row.names = FALSE)

# First plot: both conditions contain the SAME eight donors.
samples_per_condition <- table(study_design$Condition)
barplot(samples_per_condition,
  col = c("grey70", "steelblue"), ylim = c(0, 10),
  xlab = "Condition", ylab = "Donor-condition sample rows",
  main = "Documented paired study design",
  sub = "The same donors contribute to both conditions")

# Save a separate PNG copy with the same measurements and settings.
png("results/figures/documented_study_design.png", width = 1400, height = 1000, res = 160)
barplot(samples_per_condition,
  col = c("grey70", "steelblue"), ylim = c(0, 10),
  xlab = "Condition", ylab = "Donor-condition sample rows",
  main = "Documented paired study design",
  sub = "The same donors contribute to both conditions")
dev.off()

cat("Saved the design table and plot. This design-only script does not analyze RNA.\n")
