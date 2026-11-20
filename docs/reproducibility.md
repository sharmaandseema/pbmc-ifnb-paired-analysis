# Reproducibility and provenance

## Recorded execution

The learner executed the checked analysis on 20 November 2026 after the instructed R restart. The saved output identifies Windows 11, R 4.6.1, muscData 1.26.0, edgeR 4.10.7 and limma 3.68.5. ExperimentHub snapshotDate was 2026-04-21; the data loaded from cache. Full records are under results/logs/ and software versions are under results/tables/.

The rerun output folder was named rerun_20261120_112453_1d083a8797. The submitted archive contains the verified outputs, the initial study-design script, exploratory analysis script and earlier notes. The checked script had been downloaded separately and was not inside that archive; this package includes the previously delivered checked script, unchanged.

## What matched

- All 29 comparisons with previously observed dimensions, counts and rounded sensitivity metrics matched.
- All three full gene-result CSV comparisons recorded MATCH, including identical gene sets and numerical comparisons of logFC, logCPM, F, PValue and FDR using all.equal with tolerance 1e-8.
- The assistant's archive review additionally checked sample-column alignment, pseudobulk totals against sample metadata, nonnegative integer counts, paired design rank, QC subset membership, BH adjustments and the sensitivity metrics from the CSVs.

No R interpreter was used for the assistant's full-data execution. The full R analysis was run on the learner's computer. The archived matrix is pseudobulk, so raw-cell aggregation was not independently recalculated by the assistant.

## Code review

The exploratory script used qc_cd14_qc before creating it, repeated several blocks, embedded an edgeR installation command, and relied on earlier creation of some output directories. The checked script reorders the dependency, removes repeated execution, separates setup, creates output directories, preserves existing results and writes failure/completion records. The numerical model, baseline rule and two sensitivity definitions were retained.

The optional study-design script has only its final console wording updated to clarify that it is design-only. The unused exploratory script is retained in the learner's original project rather than presented as an active analysis entry point.

## Packaged output choices

Verified output files were copied from the submitted rerun folder to results/tables, results/figures and results/logs. Statistical CSVs were not edited. The new top-level README, notebook, methods commentary and explanation documents summarize the completed work rather than the earlier unfilled templates.

For display, singlet_retention_percent.png uses the learner's separately supplied Rplot01.png export of the same retention plot because the archived rerun export partly clips its x-axis title. The original rerun export is preserved as singlet_retention_percent_rerun.png. documented_study_design.png uses the supplied Rplot.png. The matching design CSV reconstructs the documented 16-row table and agrees with observed donor-condition coverage. All other PNGs are unchanged archived rerun exports. A few stimulated-sample labels overlap in the MDS image; the sample metadata table provides their exact identities.

FILE_MANIFEST_SHA256.csv records packaged-file hashes, excluding the manifest itself. It checks file integrity; it is not a numerical rerun test. The old output_inventory.csv is the original rerun inventory and does not enumerate newly authored documentation or the curated folder.

## Limits

This is consistency in one recorded software environment, not independent biological replication, new-data validation or proof of model assumptions. There is no dependency lockfile. Installing newer package versions can change results. Three package-build warnings were reported for generics, dbplyr and matrixStats; the recorded run still completed and matched.

The original input/cache, full cell-level count matrix and RDS workspaces are not included. The small derived pseudobulk count matrix is included. A new computer needs internet access and the required packages to retrieve the original data through muscData.

## Rerun procedure

Open the packaged Rproj, restart R, and source R/03_rerun_checked.R. Inspect the NEW output folder printed at completion. Its logs/run_status.txt should report completion, tables/historical_checks.csv should match the historical checks, and tables/saved_result_comparison.csv should match the three packaged reference result tables. Investigate differences rather than silently updating the reference values.
