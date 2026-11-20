# Table dictionary

All paths below are under results/tables/. Raw counts are distinct from normalized displays and model estimates.

| File or group | Row unit | Use |
|---|---|---|
| documented_study_design.csv | Documented donor-condition sample | Design-only learning table, checked against observed coverage. |
| donor_condition_cell_counts_before_qc.csv | Donor-condition group | Loaded profile counts. |
| multiplet_flags_before_qc.csv | Supplied multiplet-label category | Counts of singlet, doublet and ambiguous profiles. |
| singlet_retention_percent.csv | Donor-condition group | Percentage retained by supplied singlet labels. |
| cell_qc_metrics_before_thresholds.csv | Supplied singlet profile | Metadata, UMI total, detected genes and uninformative mitochondrial percentage. |
| lowest_gene_profiles.csv | Reviewed minimum-gene profile | Evidence for the initial unusual profile. |
| mitochondrial_count_diagnostic.csv | MT- gene | Counts in original and singlet objects. |
| cd14_counts_before_thresholds.csv; cd14_counts_after_gene_qc.csv | Donor-condition group | Coverage before and after >=200 detected-gene filtering. |
| candidate_gene_qc_cutoffs.csv | Donor-condition group | Numbers and percentages below 200 or 300 detected genes. |
| cd14_gene_qc_decisions.csv | Supplied CD14-labelled singlet profile | QC measurements and baseline keep decision for all 5,385 profiles. |
| cd14_retained_qc_metadata.csv | Baseline retained profile | Metadata and descriptive globin flag for 5,355 profiles. |
| cd14_low_gene_review.csv; cd14_retained_globin_review.csv; cd14_globin_marker_review.csv | Reviewed profile | Selected diagnostic QC/marker values; not a validation of all annotations. |
| cd14_pseudobulk_counts_baseline.csv | Gene/feature | Raw summed counts in 16 donor-condition columns, including untested genes. |
| pseudobulk_sample_metadata.csv; sample_metadata_*.csv | Pseudobulk sample | Donor, condition, retained profile number and pre-gene-filter UMI total. |
| paired_design_matrix.csv | Pseudobulk sample | Paired model matrix; includes donor terms and Conditionstim. |
| pseudobulk_gene_filter.csv | Feature | Baseline expression-filter decision. |
| pseudobulk_normalization.csv | Pseudobulk sample | Post-gene-filter library total, TMM factor and effective library size. Totals differ from pre-gene-filter totals. |
| paired_expression_results_*.csv | Tested gene | Gene, logFC, logCPM, QL F statistic, raw p-value and BH FDR. Positive logFC means higher in stimulated samples. |
| baseline_result_summary.csv; result_summary_*.csv; all_analysis_summaries.csv | Analysis | Tested/significant/up/down counts and, where included, retained profiles. |
| nexn_paired_plot_data.csv | Donor-condition sample | Prior-adjusted normalized log2 CPM for the NEXN display. |
| baseline_vs_*_by_gene.csv | Shared tested gene | Baseline and alternative logFC/FDR values, joined by gene name. |
| sensitivity_gene_agreement.csv | Agreement metric | Shared-gene denominator, direction and rank agreement, and threshold crossings. |
| historical_checks.csv | Historical observation | Expected/observed comparison for 29 recorded quantities. |
| saved_result_comparison.csv | Analysis | Match of full rerun and earlier result tables within the specified tolerance. |
| software_versions.csv | Software package | Recorded version; full session details are under results/logs/. |

Missing genes in a filtered DE table were not tested in that analysis. Their absence is not evidence of no expression or no biological effect. FDR threshold crossings can occur with small changes in estimates, variability or the tested gene set.
