# PBMC IFN-beta project: checked, self-contained rerun
# Prepared from Seema Sharma's uploaded 02_load_data.R on 20 November 2026.
# AI-assisted code review and reconstruction; the numerical model is unchanged.
#
# HOW TO RUN
# 1. Open PBMC_IFNB_project.Rproj in RStudio.
# 2. Save your open scripts. Use Session > Restart R for a fresh session.
# 3. In the Console run source(file.choose()) and select this file.
#
# Packages must already be installed. This script does not install packages.
# Kang18_8vs8() normally reuses the ExperimentHub cache; missing cache data
# may require an internet connection and a download.
# New outputs go to a unique results/rerun_* folder. Existing results and
# source scripts are not overwritten. Analysis objects live inside a function
# and do not depend on objects restored from .RData or the Global Environment.
#
# CHANGES FROM THE UPLOADED SCRIPT
# - Moved retained-profile review after creation of qc_cd14_qc.
# - Removed repeated plotting, subset, summary and >=300-model blocks.
# - Removed the package-install command from the analysis workflow.
# - Created every required output directory before writing figures/tables.
# - Added alignment, aggregation, paired-design and historical-result checks.
# - Saved software versions, session information, methods and findings.
#
# VALIDATION STATUS BEFORE DELIVERY
# Source reviewed and R syntax checked with a tree-sitter parser.
# No R interpreter or dataset was used to execute this script here.
# The clean-session rerun on your computer remains to be verified.
#
# Scientific scope: supplied CD14 labels, raw-count pseudobulk by donor and
# condition, TMM, edgeR paired model ~ Donor_ID + Condition, modern QL fit.
# Results are relative RNA condition associations in this supplied dataset.
# Two condition-specific pools mean condition and shared pool effects cannot
# be separated by donor blocking. QC sensitivity does not resolve that limit.

local({
run_checked_analysis <- function() {
  if (!file.exists("PBMC_IFNB_project.Rproj")) {
    stop("Open PBMC_IFNB_project.Rproj, then run this file from the project root.")
  }
  required <- c("muscData", "SummarizedExperiment", "Matrix", "edgeR")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) {
    stop("Required packages are missing: ", paste(required[!available], collapse = ", "),
         ". Stop here and report this message; no packages were installed.")
  }
  if (!"normLibSizes" %in% getNamespaceExports("edgeR")) {
    stop("This workflow needs the edgeR normLibSizes() API; report your edgeR version.")
  }
  dir.create("results", showWarnings = FALSE)
  run_dir <- tempfile(
    pattern = paste0("rerun_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = "results"
  )
  if (!dir.create(run_dir, recursive = FALSE)) stop("Could not create a new results folder.")
  for (folder in c("tables", "figures", "logs", "docs")) {
    if (!dir.create(file.path(run_dir, folder))) stop("Could not create output subfolders.")
  }

  log_file <- file.path(run_dir, "logs/rerun_console.txt")
  initial_sinks <- sink.number()
  initial_devices <- dev.list()
  sink(log_file, split = TRUE)
  on.exit({
    devices_now <- dev.list()
    for (device in setdiff(unname(devices_now), unname(initial_devices))) {
      try(grDevices::dev.off(device), silent = TRUE)
    }
    while (sink.number() > initial_sinks) sink()
  }, add = TRUE)

  cat("New output folder:", normalizePath(run_dir, winslash = "/"), "\n")
  writeLines("Run started; completion has not yet been verified.",
             file.path(run_dir, "logs/run_status.txt"))
  packages_recorded <- c("muscData", "edgeR", "limma", "SingleCellExperiment",
                         "SummarizedExperiment", "Matrix", "ExperimentHub")
  software_versions <- data.frame(
    Software = c("R", packages_recorded),
    Version = c(as.character(getRversion()), vapply(packages_recorded,
      function(pkg) as.character(utils::packageVersion(pkg)), character(1)))
  )
  write.csv(software_versions, file.path(run_dir, "tables/software_versions.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(run_dir, "logs/session_info_at_start.txt"))

  tryCatch({
# Project 2: load the public RNA dataset
library(muscData)

sce_ifnb <- Kang18_8vs8()
stopifnot(
  "counts" %in% SummarizedExperiment::assayNames(sce_ifnb),
  all(c("ind", "stim", "cell", "multiplets") %in%
      names(SummarizedExperiment::colData(sce_ifnb))),
  !anyNA(sce_ifnb$ind), !anyNA(sce_ifnb$stim),
  !anyNA(sce_ifnb$multiplets),
  setequal(as.character(unique(sce_ifnb$stim)), c("ctrl", "stim"))
)

# Inspect the object and its dimensions
print(sce_ifnb)
print(dim(sce_ifnb))

# Save a small text record of the loaded object
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  capture.output(sce_ifnb),
  file.path(run_dir, "logs/loaded_data_summary.txt")
)
print(dim(sce_ifnb))

# Inspect the downloaded cell metadata
cell_metadata <- as.data.frame(
  SummarizedExperiment::colData(sce_ifnb)
)

print(head(cell_metadata, n = 6))

# Count cell profiles for each donor and condition
donor_condition_cells <- table(
  cell_metadata$ind,
  cell_metadata$stim,
  useNA = "ifany"
)

print(donor_condition_cells)

# Save this metadata check
dir.create(file.path(run_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

write.csv(
  as.data.frame(donor_condition_cells),
  file.path(run_dir, "tables/donor_condition_cell_counts_before_qc.csv"),
  row.names = FALSE
)

# Inspect metadata fields and supplied flags
print(names(cell_metadata))

print(head(cell_metadata, n = 6))

print(table(
  cell_metadata$multiplets,
  useNA = "ifany"
))

# Keep profiles classified as singlets
keep_singlet <- sce_ifnb$multiplets == "singlet"

sce_singlets <- sce_ifnb[, keep_singlet]

print(dim(sce_singlets))

# Check donor-condition coverage after this filter
singlet_donor_condition <- table(
  sce_singlets$ind,
  sce_singlets$stim,
  useNA = "ifany"
)

print(singlet_donor_condition)

# Save the original flag counts
write.csv(
  as.data.frame(table(cell_metadata$multiplets, useNA = "ifany")),
  file.path(run_dir, "tables/multiplet_flags_before_qc.csv"),
  row.names = FALSE
)

# Percentage retained after selecting singlet-labelled profiles
retained_percent <- 100 *
  singlet_donor_condition / donor_condition_cells

print(round(retained_percent, digits = 1))


write.csv(
  as.data.frame(retained_percent),
  file.path(run_dir, "tables/singlet_retention_percent.csv"),
  row.names = FALSE
)

png(
  file.path(run_dir, "figures/singlet_retention_percent.png"),
  width = 1700,
  height = 1100,
  res = 160
)

par(mar = c(5, 5.5, 4, 1.5))

barplot(
  t(retained_percent),
  beside = TRUE,
  col = c("grey70", "steelblue"),
  ylim = c(0, 100),
  xlab = "Donor ID",
  ylab = "Profiles retained (%)",
  main = "Retention after the singlet-label filter",
  legend.text = c("ctrl", "stim"),
  args.legend = list(x = "top", horiz = TRUE, bty = "n")
)

dev.off()

# Retrieve the raw gene-by-cell count matrix
rna_counts <- SummarizedExperiment::assay(
  sce_singlets,
  "counts"
)

# Create a metadata table for the retained profiles
qc_metadata <- as.data.frame(
  SummarizedExperiment::colData(sce_singlets)
)

# Calculate two QC measurements per profile
qc_metadata$total_counts <- Matrix::colSums(rna_counts)

qc_metadata$detected_genes <- Matrix::colSums(rna_counts > 0)

# Inspect their distributions
print(summary(qc_metadata$total_counts))
print(summary(qc_metadata$detected_genes))

# Save the measurements together with donor and condition labels
write.csv(
  qc_metadata,
  file.path(run_dir, "tables/cell_qc_metrics_before_thresholds.csv"),
  row.names = TRUE
)

# Find profiles with the smallest number of detected genes
lowest_gene_profiles <- qc_metadata[
  qc_metadata$detected_genes == min(qc_metadata$detected_genes),
  c("ind", "stim", "cell", "total_counts", "detected_genes")
]

print(lowest_gene_profiles)

write.csv(
  lowest_gene_profiles,
  file.path(run_dir, "tables/lowest_gene_profiles.csv"),
  row.names = TRUE
)

# Save a plot comparing the two measurements
png(
  file.path(run_dir, "figures/counts_vs_detected_genes.png"),
  width = 1700,
  height = 1100,
  res = 160
)

par(mar = c(5, 5.5, 4, 1.5))

plot(
  qc_metadata$total_counts,
  qc_metadata$detected_genes,
  pch = 16,
  cex = 0.3,
  col = "grey60",
  xlab = "Total UMI counts per profile",
  ylab = "Detected genes per profile",
  main = "RNA count depth and detected genes"
)

points(
  lowest_gene_profiles$total_counts,
  lowest_gene_profiles$detected_genes,
  pch = 16,
  col = "red",
  cex = 1
)

dev.off()


# Find this profile's position in our QC table
profile_index <- which(
  rownames(qc_metadata) == "GCCACTACAGGTTC-1"
)

# Extract its counts across all genes
if (length(profile_index) != 1L) {
  stop("The historical profile barcode was not found uniquely; review the dataset version.")
}
profile_counts <- rna_counts[, profile_index]

# Keep nonzero counts and order them from largest to smallest
detected_profile_counts <- sort(
  profile_counts[profile_counts > 0],
  decreasing = TRUE
)

print(detected_profile_counts)
print(sum(detected_profile_counts))


# Identify mitochondrial genes by their MT- prefix
mito_genes <- grepl("^MT-", rownames(rna_counts))
print(rownames(rna_counts)[mito_genes])

# Stop if no mitochondrial genes were identified
stopifnot(sum(mito_genes) > 0)
# No mitochondrial threshold is applied: the supplied rows were all zero
# in the historical run. Their current totals are checked below and recorded.

# Calculate mitochondrial UMI percentage for each profile
mito_counts <- Matrix::colSums(
  rna_counts[mito_genes, , drop = FALSE]
)

qc_metadata$mito_percent <- 100 *
  mito_counts / qc_metadata$total_counts

print(summary(qc_metadata$mito_percent))

# Inspect this metric for our reviewed profile
print(qc_metadata[
  profile_index,
  c("total_counts", "detected_genes", "mito_percent")
])

# Update the saved QC table with the new measurement
write.csv(
  qc_metadata,
  file.path(run_dir, "tables/cell_qc_metrics_before_thresholds.csv"),
  row.names = TRUE
)

# Total counts for each mitochondrial gene across retained profiles
singlet_mito_totals <- Matrix::rowSums(
  rna_counts[mito_genes, , drop = FALSE]
)

print(singlet_mito_totals)

# Check the original downloaded object, before our singlet filter
original_counts <- SummarizedExperiment::assay(
  sce_ifnb,
  "counts"
)

original_mito_genes <- grepl("^MT-", rownames(original_counts))

original_mito_totals <- Matrix::rowSums(
  original_counts[original_mito_genes, , drop = FALSE]
)

print(original_mito_totals)

# Save evidence for the mitochondrial QC limitation
mito_diagnostic <- data.frame(
  Gene = names(original_mito_totals),
  Original_counts = as.numeric(original_mito_totals),
  Singlet_counts = as.numeric(
    singlet_mito_totals[names(original_mito_totals)]
  )
)

write.csv(
  mito_diagnostic,
  file.path(run_dir, "tables/mitochondrial_count_diagnostic.csv"),
  row.names = FALSE
)


# Select the supplied CD14+ Monocytes annotation
keep_cd14 <- qc_metadata$cell %in% "CD14+ Monocytes"

# Keep matching RNA columns and metadata rows
sce_cd14 <- sce_singlets[, keep_cd14]

qc_cd14 <- qc_metadata[keep_cd14, , drop = FALSE]

# Check dimensions and donor-condition coverage
print(dim(sce_cd14))

cd14_donor_condition <- table(
  qc_cd14$ind,
  qc_cd14$stim,
  useNA = "ifany"
)

print(cd14_donor_condition)

# Inspect measurements within our target population
print(summary(qc_cd14$total_counts))
print(summary(qc_cd14$detected_genes))

# Create labels such as 101_ctrl and 101_stim
qc_cd14$qc_group <- interaction(
  qc_cd14$ind,
  qc_cd14$stim,
  sep = "_",
  lex.order = TRUE,
  drop = TRUE
)

# Save the donor-condition profile counts
write.csv(
  as.data.frame(cd14_donor_condition),
  file.path(run_dir, "tables/cd14_counts_before_thresholds.csv"),
  row.names = FALSE
)

# Save two QC panels
png(
  file.path(run_dir, "figures/cd14_qc_by_donor_condition.png"),
  width = 2000,
  height = 1600,
  res = 160
)

par(mfrow = c(2, 1), mar = c(7, 5.5, 3, 1))

boxplot(
  total_counts ~ qc_group,
  data = qc_cd14,
  col = c("grey75", "steelblue"),
  log = "y",
  las = 2,
  cex.axis = 0.8,
  xlab = "",
  ylab = "Total UMI count (log scale)",
  main = "CD14-labelled profiles: count depth"
)

boxplot(
  detected_genes ~ qc_group,
  data = qc_cd14,
  col = c("grey75", "steelblue"),
  log = "y",
  las = 2,
  cex.axis = 0.8,
  xlab = "",
  ylab = "Detected genes (log scale)",
  main = "CD14-labelled profiles: detected genes"
)

dev.off()

# Compare candidate gene-detection cutoffs by donor and condition
candidate_qc <- data.frame(
  Group = levels(qc_cd14$qc_group),
  
  Profiles = as.integer(table(qc_cd14$qc_group)),
  
  Below200 = as.integer(tapply(
    qc_cd14$detected_genes < 200,
    qc_cd14$qc_group,
    sum
  )),
  
  Below300 = as.integer(tapply(
    qc_cd14$detected_genes < 300,
    qc_cd14$qc_group,
    sum
  ))
)

# Calculate the percentage flagged within each group
candidate_qc$Pct200 <- round(
  100 * candidate_qc$Below200 / candidate_qc$Profiles,
  digits = 1
)

candidate_qc$Pct300 <- round(
  100 * candidate_qc$Below300 / candidate_qc$Profiles,
  digits = 1
)

print(candidate_qc)

write.csv(
  candidate_qc,
  file.path(run_dir, "tables/candidate_gene_qc_cutoffs.csv"),
  row.names = FALSE
)


# Retrieve counts for the CD14-labelled population
cd14_counts <- SummarizedExperiment::assay(sce_cd14, "counts")
stopifnot(
  identical(colnames(cd14_counts), rownames(qc_cd14)),
  all(c("HBB", "HBA1", "HBA2") %in% rownames(cd14_counts)),
  all(qc_cd14$total_counts > 0)
)

# Calculate the fraction contributed by three hemoglobin genes
globin3_counts <- Matrix::colSums(
  cd14_counts[c("HBB", "HBA1", "HBA2"), , drop = FALSE]
)

qc_cd14$globin3_percent <- 100 *
  globin3_counts / qc_cd14$total_counts

# Review profiles below 200 detected genes
low_gene_review <- qc_cd14[
  qc_cd14$detected_genes < 200,
  c("ind", "stim", "total_counts",
    "detected_genes", "globin3_percent")
]

print(low_gene_review)



write.csv(
  low_gene_review,
  file.path(run_dir, "tables/cd14_low_gene_review.csv"),
  row.names = TRUE
)

# Record the baseline rule for every CD14-labelled profile
qc_cd14$keep_ge200 <- qc_cd14$detected_genes >= 200

# Apply the same selection to RNA columns and metadata rows
sce_cd14_qc <- sce_cd14[, qc_cd14$keep_ge200]

qc_cd14_qc <- qc_cd14[
  qc_cd14$keep_ge200,
  ,
  drop = FALSE
]

print(dim(sce_cd14_qc))

# Check donor-condition coverage after filtering
cd14_counts_after_qc <- table(
  qc_cd14_qc$ind,
  qc_cd14_qc$stim,
  useNA = "ifany"
)

print(cd14_counts_after_qc)

# Inspect hemoglobin contribution among retained profiles
print(summary(qc_cd14_qc$globin3_percent))

# Save decisions for all profiles, including excluded ones
write.csv(
  qc_cd14,
  file.path(run_dir, "tables/cd14_gene_qc_decisions.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(cd14_counts_after_qc),
  file.path(run_dir, "tables/cd14_counts_after_gene_qc.csv"),
  row.names = FALSE
)

# Order retained profiles by decreasing hemoglobin contribution
globin_order <- order(
  qc_cd14_qc$globin3_percent,
  decreasing = TRUE
)

# Inspect the ten highest values
globin_review <- head(
  qc_cd14_qc[
    globin_order,
    c("ind", "stim", "total_counts",
      "detected_genes", "globin3_percent")
  ],
  n = 10
)

print(globin_review)

write.csv(
  globin_review,
  file.path(run_dir, "tables/cd14_retained_globin_review.csv"),
  row.names = TRUE
)

# Define a small marker panel
review_marker_genes <- c("CD14", "LYZ", "S100A9")

# Check that these genes are available
marker_presence <- data.frame(
  Gene = review_marker_genes,
  Present = review_marker_genes %in% rownames(sce_cd14_qc)
)

print(marker_presence)

present_markers <- review_marker_genes[marker_presence$Present]

# Select the same ten profiles, in the same order as globin_review
review_indices <- head(globin_order, n = 10)

review_marker_counts <- as.matrix(
  SummarizedExperiment::assay(sce_cd14_qc, "counts")[
    present_markers,
    review_indices,
    drop = FALSE
  ]
)

# Add marker counts beside each profile's QC information
marker_review <- cbind(
  globin_review,
  as.data.frame(t(review_marker_counts))
)

print(marker_review)

write.csv(
  marker_review,
  file.path(run_dir, "tables/cd14_globin_marker_review.csv"),
  row.names = TRUE
)

# Record a descriptive review flag: three globin genes exceed 50%
qc_cd14_qc$globin_majority_flag <-
  qc_cd14_qc$globin3_percent > 50

print(table(qc_cd14_qc$globin_majority_flag))

write.csv(
  qc_cd14_qc,
  file.path(run_dir, "tables/cd14_retained_qc_metadata.csv"),
  row.names = TRUE
)

# Retrieve raw counts for the baseline population
cd14_qc_counts <- SummarizedExperiment::assay(
  sce_cd14_qc,
  "counts"
)

sample_group_ids <- levels(qc_cd14_qc$qc_group)

# Create an empty gene-by-sample table
pseudobulk_counts <- matrix(
  0,
  nrow = nrow(cd14_qc_counts),
  ncol = length(sample_group_ids),
  dimnames = list(
    rownames(cd14_qc_counts),
    sample_group_ids
  )
)

# Sum each gene across profiles belonging to each sample group
for (group_name in sample_group_ids) {
  profiles_in_group <- qc_cd14_qc$qc_group == group_name
  
  pseudobulk_counts[, group_name] <- Matrix::rowSums(
    cd14_qc_counts[, profiles_in_group, drop = FALSE]
  )
}

# Inspect dimensions and check that all counts were preserved
print(dim(pseudobulk_counts))
stopifnot(sum(pseudobulk_counts) == sum(cd14_qc_counts))

write.csv(
  pseudobulk_counts,
  file.path(run_dir, "tables/cd14_pseudobulk_counts_baseline.csv"),
  row.names = TRUE
)

# Extract one metadata row per donor-condition group
sample_metadata <- unique(
  qc_cd14_qc[, c("qc_group", "ind", "stim")]
)

colnames(sample_metadata) <- c(
  "Sample_ID", "Donor_ID", "Condition"
)

# Match metadata order to the pseudobulk columns
sample_metadata <- sample_metadata[
  match(colnames(pseudobulk_counts), sample_metadata$Sample_ID),
  ,
  drop = FALSE
]

rownames(sample_metadata) <- as.character(sample_metadata$Sample_ID)

# Prepare categorical variables for the later paired model
sample_metadata$Donor_ID <- factor(sample_metadata$Donor_ID)

sample_metadata$Condition <- factor(
  sample_metadata$Condition,
  levels = c("ctrl", "stim")
)

# Add profile numbers and total UMI counts
sample_metadata$Profiles <- as.integer(
  table(qc_cd14_qc$qc_group)[
    as.character(sample_metadata$Sample_ID)
  ]
)

sample_metadata$Total_UMIs <- as.numeric(
  colSums(pseudobulk_counts)
)

print(sample_metadata)

# Confirm exact alignment
print(identical(
  as.character(sample_metadata$Sample_ID),
  colnames(pseudobulk_counts)
))

write.csv(
  sample_metadata,
  file.path(run_dir, "tables/pseudobulk_sample_metadata.csv"),
  row.names = FALSE
)


library(edgeR)

# Explicitly identify the categories and reference condition
sample_metadata$Donor_ID <- factor(sample_metadata$Donor_ID)

sample_metadata$Condition <- factor(
  sample_metadata$Condition,
  levels = c("ctrl", "stim")
)

# Create the paired model's design matrix
design <- model.matrix(
  ~ Donor_ID + Condition,
  data = sample_metadata
)

stopifnot(
  identical(colnames(pseudobulk_counts), rownames(design)),
  all(sample_metadata$Profiles > 0),
  all(table(sample_metadata$Donor_ID, sample_metadata$Condition) == 1L)
)
print(dim(design))
print(colnames(design))

# Check the number of independent model columns
print(qr(design)$rank)
stopifnot(qr(design)$rank == ncol(design))

write.csv(
  design,
  file.path(run_dir, "tables/paired_design_matrix.csv"),
  row.names = TRUE
)

# Store raw pseudobulk counts in an edgeR object
dge <- DGEList(
  counts = pseudobulk_counts,
  samples = sample_metadata,
  group = sample_metadata$Condition
)

# Identify genes with sufficient expression across samples
gene_filter <- filterByExpr(
  dge,
  group = sample_metadata$Condition,
  min.count = 10,
  min.total.count = 15
)

gene_filter_record <- data.frame(
  Gene = rownames(dge),
  Retained = gene_filter
)

print(table(gene_filter))

# Select gene rows and recalculate library totals
dge <- dge[gene_filter, , keep.lib.sizes = FALSE]

# Calculate TMM library normalization factors
dge <- normLibSizes(dge, method = "TMM")

print(dim(dge))

# Record the normalization information
normalization_record <- data.frame(
  Sample_ID = colnames(dge),
  Library_UMIs = dge$samples$lib.size,
  TMM_factor = dge$samples$norm.factors,
  Effective_library_size =
    dge$samples$lib.size * dge$samples$norm.factors
)

print(normalization_record)

write.csv(
  gene_filter_record,
  file.path(run_dir, "tables/pseudobulk_gene_filter.csv"),
  row.names = FALSE
)

write.csv(
  normalization_record,
  file.path(run_dir, "tables/pseudobulk_normalization.csv"),
  row.names = FALSE
)

# Display the gene-filtering results
print(table(gene_filter))
print(dim(dge))

# Assign colours to conditions
sample_colours <- ifelse(
  sample_metadata$Condition == "ctrl",
  "grey40",
  "blue"
)

# Save the sample-level MDS plot
png(
  file.path(run_dir, "figures/cd14_pseudobulk_mds.png"),
  width = 1800,
  height = 1200,
  res = 160
)

par(mar = c(5, 5.5, 4, 2))

mds_result <- plotMDS(
  dge,
  labels = as.character(sample_metadata$Sample_ID),
  col = sample_colours,
  top = min(500, nrow(dge)),
  method = "logFC",
  gene.selection = "pairwise",
  cex = 0.8,
  main = "CD14-labelled pseudobulk samples",
  sub = paste(
    nrow(dge), "retained genes;",
    ncol(dge), "donor-condition samples"
  )
)

legend(
  "topright",
  legend = c("Control", "Stimulated"),
  text.col = c("grey40", "blue"),
  bty = "n"
)

dev.off()

# Fit the paired quasi-likelihood count model
fit <- glmQLFit(
  dge,
  design = design,
  robust = TRUE,
  legacy = FALSE
)

# Test the stimulated-versus-control coefficient
ql_test <- glmQLFTest(
  fit,
  coef = "Conditionstim"
)

# Retrieve all tested genes with BH-adjusted p-values
de_results <- topTags(
  ql_test,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue"
)$table

de_results <- data.frame(
  Gene = rownames(de_results),
  de_results,
  row.names = NULL
)

print(head(de_results, n = 10))

cat("Genes tested:", nrow(de_results), "\n")
cat(
  "Genes with BH FDR < 0.05:",
  sum(de_results$FDR < 0.05, na.rm = TRUE),
  "\n"
)

write.csv(
  de_results,
  file.path(run_dir, "tables/paired_expression_results_baseline.csv"),
  row.names = FALSE
)

# Save a model-variability diagnostic
png(
  file.path(run_dir, "figures/paired_model_variability.png"),
  width = 1700,
  height = 1100,
  res = 160
)

par(mar = c(5, 5.5, 4, 1.5))
plotQLDisp(fit, main = "Paired model: variability estimates")
dev.off()

result_summary <- data.frame(
  Genes_tested = nrow(de_results),
  
  FDR_below_005 = sum(
    de_results$FDR < 0.05,
    na.rm = TRUE
  ),
  
  Higher_in_stim = sum(
    de_results$FDR < 0.05 & de_results$logFC > 0,
    na.rm = TRUE
  ),
  
  Lower_in_stim = sum(
    de_results$FDR < 0.05 & de_results$logFC < 0,
    na.rm = TRUE
  )
)

print(result_summary)

write.csv(
  result_summary,
  file.path(run_dir, "tables/baseline_result_summary.csv"),
  row.names = FALSE
)


# Examine NEXN across the eight donor pairs
# Uses the existing dge and sample_metadata objects

stopifnot(
  "NEXN" %in% rownames(dge),
  identical(
    colnames(dge),
    as.character(sample_metadata$Sample_ID)
  )
)

# TMM-normalized log2 counts per million for display
log_cpm <- edgeR::cpm(
  dge,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

nexn_plot_data <- data.frame(
  Sample_ID = as.character(sample_metadata$Sample_ID),
  Donor_ID = as.character(sample_metadata$Donor_ID),
  Condition = as.character(sample_metadata$Condition),
  Log2_CPM = as.numeric(log_cpm["NEXN", ])
)

write.csv(
  nexn_plot_data,
  file.path(run_dir, "tables/nexn_paired_plot_data.csv"),
  row.names = FALSE
)

# Define the plot so it can be saved and displayed
draw_nexn_plot <- function() {
  par(mar = c(5.5, 6.5, 4.5, 1.5))
  plot(
    NA,
    xlim = c(0.7, 2.3),
    ylim = range(nexn_plot_data$Log2_CPM) + c(-0.5, 0.5),
    xaxt = "n",
    xlab = "Condition",
    ylab = "TMM-normalized log2 CPM",
    main = "NEXN: eight donor pairs"
  )
  
  axis(1, at = c(1, 2), labels = c("Control", "Stimulated"))
  
  for (donor in unique(nexn_plot_data$Donor_ID)) {
    donor_data <- nexn_plot_data[
      nexn_plot_data$Donor_ID == donor, ,
      drop = FALSE
    ]
    
    donor_data <- donor_data[
      match(c("ctrl", "stim"), donor_data$Condition), ,
      drop = FALSE
    ]
    
    lines(
      c(1, 2), donor_data$Log2_CPM,
      col = "grey65", lwd = 1.5
    )
    
    points(
      c(1, 2), donor_data$Log2_CPM,
      pch = 19, col = c("grey35", "steelblue"), cex = 1.2
    )
  }
  
  mtext(
    "Each line connects samples from the same donor",
    side = 3, line = 0.3, cex = 0.8
  )
}

# Save the figure
png(
  file.path(run_dir, "figures/nexn_paired_expression.png"),
  width = 1800, height = 1300, res = 180
)
draw_nexn_plot()
dev.off()

{
  # Prepare counts and filtering flag
  counts_for_sensitivity <- SummarizedExperiment::assay(
    sce_cd14, "counts"
  )
  
  stopifnot(identical(
    colnames(counts_for_sensitivity),
    rownames(qc_cd14)
  ))
  
  keep_ge300 <- qc_cd14$detected_genes >= 300
  metadata_ge300 <- sample_metadata
  
  # Aggregate counts for each donor-condition sample
  pseudobulk_ge300 <- matrix(
    0,
    nrow = nrow(counts_for_sensitivity),
    ncol = nrow(metadata_ge300),
    dimnames = list(
      rownames(counts_for_sensitivity),
      as.character(metadata_ge300$Sample_ID)
    )
  )
  
  for (sample_id in colnames(pseudobulk_ge300)) {
    selected_profiles <- keep_ge300 &
      as.character(qc_cd14$qc_group) == sample_id
    
    pseudobulk_ge300[, sample_id] <- Matrix::rowSums(
      counts_for_sensitivity[, selected_profiles, drop = FALSE]
    )
  }
  
  metadata_ge300$Profiles <- as.integer(
    table(qc_cd14$qc_group[keep_ge300])[
      as.character(metadata_ge300$Sample_ID)
    ]
  )
  
  metadata_ge300$Total_UMIs <- as.numeric(
    colSums(pseudobulk_ge300)
  )
  
  stopifnot(
    all(metadata_ge300$Profiles > 0),
    identical(colnames(pseudobulk_ge300), rownames(design))
  )
  
  # Filter genes and normalize
  dge_ge300 <- edgeR::DGEList(
    counts = pseudobulk_ge300,
    samples = metadata_ge300,
    group = metadata_ge300$Condition
  )
  
  filter_ge300 <- edgeR::filterByExpr(
    dge_ge300,
    group = metadata_ge300$Condition,
    min.count = 10,
    min.total.count = 15
  )
  
  dge_ge300 <- dge_ge300[
    filter_ge300, , keep.lib.sizes = FALSE
  ]
  
  dge_ge300 <- edgeR::normLibSizes(
    dge_ge300, method = "TMM"
  )
  
  # Fit and test the paired model
  fit_ge300 <- edgeR::glmQLFit(
    dge_ge300,
    design = design,
    robust = TRUE,
    legacy = FALSE
  )
  
  test_ge300 <- edgeR::glmQLFTest(
    fit_ge300,
    coef = "Conditionstim"
  )
  
  results_ge300 <- edgeR::topTags(
    test_ge300,
    n = Inf,
    adjust.method = "BH",
    sort.by = "PValue"
  )$table
  
  results_ge300 <- data.frame(
    Gene = rownames(results_ge300),
    results_ge300,
    row.names = NULL
  )
  
  # Summarize and save
  summary_ge300 <- data.frame(
    Analysis = "At least 300 detected genes",
    Profiles = sum(keep_ge300),
    Genes_tested = nrow(results_ge300),
    FDR_below_005 = sum(results_ge300$FDR < 0.05),
    Higher_in_stim = sum(
      results_ge300$FDR < 0.05 & results_ge300$logFC > 0
    ),
    Lower_in_stim = sum(
      results_ge300$FDR < 0.05 & results_ge300$logFC < 0
    )
  )
  
  write.csv(
    metadata_ge300,
    file.path(run_dir, "tables/sample_metadata_ge300.csv"),
    row.names = FALSE
  )
  
  write.csv(
    results_ge300,
    file.path(run_dir, "tables/paired_expression_results_ge300.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary_ge300,
    file.path(run_dir, "tables/result_summary_ge300.csv"),
    row.names = FALSE
  )
  
  print(summary_ge300)
}

{
  # Sensitivity check 2: exclude globin-majority profiles
  # Start from the baseline population with >=200 detected genes
  
  counts_globin_check <- SummarizedExperiment::assay(
    sce_cd14_qc, "counts"
  )
  
  stopifnot(identical(
    colnames(counts_globin_check),
    rownames(qc_cd14_qc)
  ))
  
  keep_globin_check <- qc_cd14_qc$globin3_percent <= 50
  stopifnot(!anyNA(keep_globin_check))
  
  metadata_globin <- sample_metadata
  
  # Aggregate raw counts for the same 16 samples
  pseudobulk_globin <- matrix(
    0,
    nrow = nrow(counts_globin_check),
    ncol = nrow(metadata_globin),
    dimnames = list(
      rownames(counts_globin_check),
      as.character(metadata_globin$Sample_ID)
    )
  )
  
  for (sample_id in colnames(pseudobulk_globin)) {
    selected_profiles <- keep_globin_check &
      as.character(qc_cd14_qc$qc_group) == sample_id
    
    pseudobulk_globin[, sample_id] <- Matrix::rowSums(
      counts_globin_check[, selected_profiles, drop = FALSE]
    )
  }
  
  metadata_globin$Profiles <- as.integer(
    table(qc_cd14_qc$qc_group[keep_globin_check])[
      as.character(metadata_globin$Sample_ID)
    ]
  )
  
  metadata_globin$Total_UMIs <- as.numeric(
    colSums(pseudobulk_globin)
  )
  
  stopifnot(
    all(metadata_globin$Profiles > 0),
    identical(colnames(pseudobulk_globin), rownames(design))
  )
  
  # Repeat the same gene filtering and normalization
  dge_globin <- edgeR::DGEList(
    counts = pseudobulk_globin,
    samples = metadata_globin,
    group = metadata_globin$Condition
  )
  
  filter_globin <- edgeR::filterByExpr(
    dge_globin,
    group = metadata_globin$Condition,
    min.count = 10,
    min.total.count = 15
  )
  
  dge_globin <- dge_globin[
    filter_globin, , keep.lib.sizes = FALSE
  ]
  
  dge_globin <- edgeR::normLibSizes(
    dge_globin, method = "TMM"
  )
  
  # Fit and test the same paired model
  fit_globin <- edgeR::glmQLFit(
    dge_globin,
    design = design,
    robust = TRUE,
    legacy = FALSE
  )
  
  test_globin <- edgeR::glmQLFTest(
    fit_globin,
    coef = "Conditionstim"
  )
  
  results_globin <- edgeR::topTags(
    test_globin,
    n = Inf,
    adjust.method = "BH",
    sort.by = "PValue"
  )$table
  
  results_globin <- data.frame(
    Gene = rownames(results_globin),
    results_globin,
    row.names = NULL
  )
  
  # Summarize and save
  summary_globin <- data.frame(
    Analysis = "Exclude globin-majority profiles",
    Profiles = sum(keep_globin_check),
    Genes_tested = nrow(results_globin),
    FDR_below_005 = sum(results_globin$FDR < 0.05),
    Higher_in_stim = sum(
      results_globin$FDR < 0.05 & results_globin$logFC > 0
    ),
    Lower_in_stim = sum(
      results_globin$FDR < 0.05 & results_globin$logFC < 0
    )
  )
  
  write.csv(
    metadata_globin,
    file.path(run_dir, "tables/sample_metadata_globin_sensitivity.csv"),
    row.names = FALSE
  )
  
  write.csv(
    results_globin,
    file.path(run_dir, "tables/paired_expression_results_globin_sensitivity.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary_globin,
    file.path(run_dir, "tables/result_summary_globin_sensitivity.csv"),
    row.names = FALSE
  )
  
  print(summary_globin)
}

{
  # Compare individual genes with the baseline
  compare_to_baseline <- function(alternative) {
    
    stopifnot(
      !anyDuplicated(de_results$Gene),
      !anyDuplicated(alternative$Gene)
    )
    
    shared <- merge(
      de_results[, c("Gene", "logFC", "FDR")],
      alternative[, c("Gene", "logFC", "FDR")],
      by = "Gene",
      suffixes = c("_baseline", "_alternative")
    )
    
    baseline_sig <- shared$FDR_baseline < 0.05
    alternative_sig <- shared$FDR_alternative < 0.05
    
    same_direction <- sign(shared$logFC_baseline) ==
      sign(shared$logFC_alternative)
    
    metrics <- c(
      Shared_genes = nrow(shared),
      
      Tested_only_in_baseline = length(setdiff(
        de_results$Gene, alternative$Gene
      )),
      
      Tested_only_in_alternative = length(setdiff(
        alternative$Gene, de_results$Gene
      )),
      
      LogFC_rank_correlation = round(cor(
        shared$logFC_baseline,
        shared$logFC_alternative,
        method = "spearman"
      ), 4),
      
      Same_direction_percent = round(
        100 * mean(same_direction), 2
      ),
      
      FDR_below_005_in_both = sum(
        baseline_sig & alternative_sig
      ),
      
      FDR_below_005_baseline_only = sum(
        baseline_sig & !alternative_sig
      ),
      
      FDR_below_005_alternative_only = sum(
        !baseline_sig & alternative_sig
      ),
      
      Significant_in_both_opposite_direction = sum(
        baseline_sig & alternative_sig & !same_direction
      )
    )
    
    list(metrics = metrics, genes = shared)
  }
  
  comparison_ge300 <- compare_to_baseline(results_ge300)
  comparison_globin <- compare_to_baseline(results_globin)
  
  sensitivity_agreement <- data.frame(
    Metric = names(comparison_ge300$metrics),
    GE300 = unname(comparison_ge300$metrics),
    Globin_excluded = unname(comparison_globin$metrics)
  )
  
  write.csv(
    comparison_ge300$genes,
    file.path(run_dir, "tables/baseline_vs_ge300_by_gene.csv"),
    row.names = FALSE
  )
  
  write.csv(
    comparison_globin$genes,
    file.path(run_dir, "tables/baseline_vs_globin_by_gene.csv"),
    row.names = FALSE
  )
  
  write.csv(
    sensitivity_agreement,
    file.path(run_dir, "tables/sensitivity_gene_agreement.csv"),
    row.names = FALSE
  )
  
  print(sensitivity_agreement, row.names = FALSE)
}


# Summarize all three analyses in one table.
analysis_summary <- rbind(
  data.frame(Analysis = "Baseline: >=200 detected genes",
             Profiles = ncol(sce_cd14_qc), result_summary),
  summary_ge300, summary_globin
)
write.csv(analysis_summary, file.path(run_dir, "tables/all_analysis_summaries.csv"), row.names = FALSE)

# Compare a rerun with the numerical observations already supplied in chat.
# These are historical checks, not evidence that statistical assumptions hold.
observed <- c(
  Features_loaded = nrow(sce_ifnb), Profiles_loaded = ncol(sce_ifnb),
  Singlets = ncol(sce_singlets), CD14_before_gene_QC = ncol(sce_cd14),
  Baseline_profiles = ncol(sce_cd14_qc), Baseline_tested = nrow(de_results),
  Baseline_significant = result_summary$FDR_below_005,
  Baseline_higher = result_summary$Higher_in_stim,
  Baseline_lower = result_summary$Lower_in_stim,
  GE300_profiles = sum(keep_ge300), GE300_tested = nrow(results_ge300),
  GE300_significant = summary_ge300$FDR_below_005,
  GE300_higher = summary_ge300$Higher_in_stim, GE300_lower = summary_ge300$Lower_in_stim,
  Globin_profiles = sum(keep_globin_check), Globin_tested = nrow(results_globin),
  Globin_significant = summary_globin$FDR_below_005,
  Globin_higher = summary_globin$Higher_in_stim, Globin_lower = summary_globin$Lower_in_stim,
  Original_mitochondrial_total = sum(original_mito_totals),
  Flagged_globin_majority = sum(qc_cd14_qc$globin_majority_flag),
  GE300_rank_correlation = unname(comparison_ge300$metrics["LogFC_rank_correlation"]),
  Globin_rank_correlation = unname(comparison_globin$metrics["LogFC_rank_correlation"]),
  GE300_same_direction_percent = unname(comparison_ge300$metrics["Same_direction_percent"]),
  Globin_same_direction_percent = unname(comparison_globin$metrics["Same_direction_percent"]),
  GE300_significant_both = unname(comparison_ge300$metrics["FDR_below_005_in_both"]),
  Globin_significant_both = unname(comparison_globin$metrics["FDR_below_005_in_both"]),
  GE300_opposite_significant = unname(comparison_ge300$metrics["Significant_in_both_opposite_direction"]),
  Globin_opposite_significant = unname(comparison_globin$metrics["Significant_in_both_opposite_direction"])
)
expected <- c(35635,29065,24679,5385,5355,3697,2506,1252,1254,
              5254,3687,2497,1246,1251,5348,3696,2525,1260,1265,
              0,7,0.9999,0.9997,99.73,99.78,2484,2502,0,0)
stopifnot(length(observed) == length(expected))
historical_checks <- data.frame(
  Metric = names(observed), Expected = expected, Observed = unname(observed),
  Matches = abs(unname(observed) - expected) < 1e-10
)
write.csv(historical_checks, file.path(run_dir, "tables/historical_checks.csv"), row.names = FALSE)

# Check full result tables against existing CSVs when they are available.
# Read references from the original results/tables, never from the new run.
reference_files <- c(
  Baseline = "paired_expression_results_baseline.csv",
  GE300 = "paired_expression_results_ge300.csv",
  Globin_excluded = "paired_expression_results_globin_sensitivity.csv"
)
current_results <- list(Baseline = de_results, GE300 = results_ge300,
                        Globin_excluded = results_globin)
reference_checks <- do.call(rbind, lapply(names(reference_files), function(label) {
  path <- file.path("results/tables", reference_files[[label]])
  if (!file.exists(path)) {
    return(data.frame(Analysis = label, Status = "REFERENCE_NOT_FOUND",
                      Same_gene_set = NA, Numeric_match = NA))
  }
  old <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL)
  columns <- c("Gene", "logFC", "logCPM", "F", "PValue", "FDR")
  if (is.null(old) || !all(columns %in% names(old)) || anyDuplicated(old$Gene)) {
    return(data.frame(Analysis = label, Status = "INVALID_REFERENCE",
                      Same_gene_set = NA, Numeric_match = NA))
  }
  fresh <- current_results[[label]]
  same_genes <- setequal(old$Gene, fresh$Gene)
  same_values <- FALSE
  if (same_genes) {
    old <- old[match(fresh$Gene, old$Gene), columns, drop = FALSE]
    numeric_columns <- setdiff(columns, "Gene")
    same_values <- all(vapply(numeric_columns, function(column) {
      isTRUE(all.equal(as.numeric(old[[column]]), as.numeric(fresh[[column]]),
                       tolerance = 1e-8, check.attributes = FALSE))
    }, logical(1)))
  }
  data.frame(Analysis = label,
             Status = if (same_genes && same_values) "MATCH" else "DIFFERENCE",
             Same_gene_set = same_genes, Numeric_match = same_values)
}))
write.csv(reference_checks, file.path(run_dir, "tables/saved_result_comparison.csv"), row.names = FALSE)

# Record what was reviewed and what remains uncertain.
methods <- c(
  "# Methods and interpretation",
  "", "## Question and dataset",
  "How do relative RNA profiles of supplied CD14-labelled monocytes differ between control and IFN-beta-stimulated samples across eight donors?",
  "Public Kang18_8vs8 dataset from muscData; supplied donor, condition, cell-type and multiplet metadata were used.",
  "", "## Analysis",
  "Kept supplied singlet-labelled profiles and the supplied CD14+ Monocytes annotation.",
  "Baseline: at least 200 detected genes. No additional mitochondrial threshold was applied.",
  "Summed raw counts by donor-condition group, producing 16 pseudobulk samples.",
  "Filtered genes with filterByExpr(group=Condition, min.count=10, min.total.count=15).",
  "Applied TMM and the paired design ~ Donor_ID + Condition.",
  "Used glmQLFit(robust=TRUE, legacy=FALSE), glmQLFTest(coef='Conditionstim'), and BH FDR <0.05.",
  "Positive logFC denotes higher relative RNA in stimulated samples.",
  "Two separate QC sensitivity analyses used >=300 detected genes, or excluded baseline profiles with HBB+HBA1+HBA2 above 50% of counts.",
  "Per-gene sensitivity comparisons include only genes tested in both analyses.",
  "", "## Observations and decisions",
  sprintf("Baseline retained %d profiles and tested %d genes; %d met FDR <0.05 (%d higher, %d lower in stimulated).",
          ncol(sce_cd14_qc), nrow(de_results), result_summary$FDR_below_005,
          result_summary$Higher_in_stim, result_summary$Lower_in_stim),
  sprintf("GE300: %d profiles, %d tested genes and %d at FDR <0.05.",
          sum(keep_ge300), nrow(results_ge300), summary_ge300$FDR_below_005),
  sprintf("Globin sensitivity: %d profiles, %d tested genes and %d at FDR <0.05.",
          sum(keep_globin_check), nrow(results_globin), summary_globin$FDR_below_005),
  "See sensitivity_gene_agreement.csv for observed direction agreement, rank correlations and significance changes.",
  "Retain the >=200 analysis as the baseline, report both alternatives, and flag any mismatch with historical checks for review.",
  "", "## NEXN example",
  "NEXN was selected after reviewing the baseline results; it is an illustration, not independent validation.",
  "Each plotted point is a donor-condition pseudobulk sample, linked within donor. Display uses TMM-normalized log2 CPM with prior.count=2.",
  "", "## Limitations",
  "The supplied 13 MT- gene rows were all zero in the historical object, including before singlet filtering. The upstream cause is unresolved; zero values do not establish healthy cells.",
  "Cell-type labels were supplied, not independently re-annotated. The small marker review does not validate all identities.",
  "Genetic multiplet flags do not identify every same-donor doublet. No ambient-RNA correction was performed.",
  "The three-globin fraction is a descriptive count fraction, not an estimated contamination fraction or validated threshold.",
  "Donor aliquots were pooled by condition into one control and one stimulated pool. Donor blocking cannot separate condition from shared pool effects.",
  "These are study-specific relative RNA associations, not protein measurements, absolute RNA per cell, or an unconfounded causal effect.",
  "No new biological discovery or independent validation is claimed. Coding and troubleshooting were AI-assisted.",
  "", "## Sources",
  "https://bioconductor.org/packages/release/data/experiment/html/muscData.html",
  "https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf",
  "https://stacks.cdc.gov/view/cdc/79371/cdc_79371_DS1.pdf"
)
writeLines(methods, file.path(run_dir, "docs/methods_and_findings.md"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "logs/session_info.txt"))
write.csv(data.frame(File = list.files(run_dir, recursive = TRUE)),
          file.path(run_dir, "logs/output_inventory.csv"), row.names = FALSE)

cat("\nALL ANALYSIS SUMMARIES\n")
print(analysis_summary, row.names = FALSE)
cat("\nHISTORICAL CHECKS\n")
print(historical_checks, row.names = FALSE)
cat("\nCOMPARISON WITH EXISTING RESULT CSV FILES\n")
print(reference_checks, row.names = FALSE)
matches <- all(historical_checks$Matches) && all(reference_checks$Status == "MATCH")
status <- if (matches) {
  "COMPLETED: historical checks and all three saved-result comparisons match."
} else {
  "COMPLETED WITH REVIEW NEEDED: inspect historical_checks.csv and saved_result_comparison.csv."
}
writeLines(status, file.path(run_dir, "logs/run_status.txt"))
cat("\n", status, "\n", sep = "")
cat("NEW output folder: ", normalizePath(run_dir, winslash = "/"), "\n", sep = "")
cat("Existing results were not overwritten. No packages were installed.\n")
invisible(run_dir)
  }, error = function(e) {
    error_text <- paste("STOPPED:", conditionMessage(e))
    writeLines(error_text, file.path(run_dir, "logs/run_status.txt"))
    writeLines(capture.output(sessionInfo()), file.path(run_dir, "logs/session_info.txt"))
    cat("\n", error_text, "\nOutput folder: ", normalizePath(run_dir, winslash = "/"), "\n", sep = "")
    stop(conditionMessage(e), call. = FALSE)
  })
}
run_checked_analysis()
})
