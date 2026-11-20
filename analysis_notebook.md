# Completed analysis notebook

Seema Sharma | Public Kang18_8vs8 RNA dataset | Work recorded 19-20 November 2026

This record was completed from the learner's outputs and the submitted verified-rerun archive. Observations are separated from interpretation and decisions. Coding and documentation were AI-assisted.

## 1. Study design and data units

**Question:** What are the biological replicates and paired measurements?

**Observed:** 35,635 features by 29,065 profiles, with metadata ind, stim, cluster, cell and multiplets. All eight documented donors have both ctrl and stim profiles. There are 16 donor-condition combinations.

**Interpretation:** A profile is a cell-associated barcode, a pseudobulk sample is one donor-condition aggregate, and a donor is a person. Thousands of profiles do not provide thousands of independent donor replicates.

**Limit:** The original study pooled samples by condition. Condition and shared pool effects cannot be separated by this model.

**Decision / next check:** Keep pairing in the design; inspect supplied flags and group coverage before aggregation.

**Evidence:** donor_condition_cell_counts_before_qc.csv and documented_study_design.csv under results/tables/.

## 2. Supplied singlet labels

**Question:** Which supplied demultiplexing labels should define the analysis population?

**Observed:** 24,679 singlet, 3,169 doublet and 1,217 ambiguous profiles. Keeping singlet-labelled profiles excludes 4,386 profiles. Both conditions remain represented for all donors; retention ranges from about 79% to 89% across donor-condition groups.

**Interpretation:** The dataset supplies usable genetic multiplet labels. Retention percentages are a property of this filtering step, not proof of cell health or biological response.

**Limit:** Genetic labels do not detect every same-donor doublet. No new doublet algorithm was run.

**Decision / next check:** Keep supplied singlets and inspect per-profile RNA count depth and detected genes.

**Evidence:** multiplet_flags_before_qc.csv, singlet_retention_percent.csv and the matching figure.

## 3. Unusual profiles and unavailable mitochondrial QC

**Question:** Can count totals or mitochondrial percentages alone establish quality?

**Observed:** Across singlets, total counts range from 562 to 11,819 and detected genes from 12 to 2,760. Barcode GCCACTACAGGTTC-1 has 763 counts but only 12 detected genes; HBB, HBA1 and HBA2 contribute 750 counts (about 98.3%). It carries a supplied Megakaryocytes label. All 13 MT- rows have zero totals in the original object and after singlet filtering.

**Interpretation:** A large count total can be dominated by a few genes. Mitochondrial percentage cannot distinguish profiles in this supplied matrix.

**Limit:** This review cannot establish whether the globin-rich profile is an erythroid cell, a mixed droplet or affected by ambient RNA. The upstream reason for mitochondrial zeros is unresolved.

**Decision / next check:** Do not apply a mitochondrial cutoff or interpret zero as good quality. Examine the selected cell population by donor and condition.

**Evidence:** lowest_gene_profiles.csv, mitochondrial_count_diagnostic.csv and counts_vs_detected_genes.png.

## 4. Target population and baseline gene QC

**Question:** Which supplied CD14-labelled profiles have sufficient gene detection for the baseline?

**Observed:** 5,385 singlet profiles carry the supplied CD14+ Monocytes label. Thirty have <200 detected genes: 29 from donor 1039 and one from donor 1488. The 29 donor-1039 profiles have three-globin fractions of roughly 76%-97%; the donor-1488 profile has 193 detected genes and 0% three-globin fraction. The >=200 rule retains 5,355 profiles, with 112-779 profiles in each donor-condition group.

**Interpretation:** The baseline removes a small, investigated low-gene tail while preserving all donor pairs. The threshold is a dataset-specific choice, not a universal guarantee of quality.

**Limit:** Supplied cell labels remain provisional. Three genes were inspected only in the ten retained profiles with the highest globin fractions.

**Decision / next check:** Use >=200 detected genes as baseline and plan >=300 as one sensitivity analysis.

**Evidence:** candidate_gene_qc_cutoffs.csv, cd14_low_gene_review.csv, cd14_gene_qc_decisions.csv and cd14_counts_after_gene_qc.csv.

## 5. Retained globin-rich profiles

**Question:** Do some retained profiles remain dominated by globin counts?

**Observed:** Seven baseline profiles have HBB+HBA1+HBA2 above 50% of counts. Nine of the ten highest-globin retained profiles express at least one of CD14, LYZ or S100A9; one expresses none of these three genes.

**Interpretation:** The supplied CD14 annotation has partial RNA marker support in this small reviewed group, but mixed or ambient contributions remain possible.

**Limit:** Marker zeros do not prove absence of a cell identity. The three-globin percentage is not an ambient contamination estimate.

**Decision / next check:** Retain all 5,355 baseline profiles, flag the seven descriptively, and exclude them only in a separate planned sensitivity analysis.

**Evidence:** cd14_globin_marker_review.csv and cd14_retained_qc_metadata.csv.

## 6. Pseudobulk and model setup

**Question:** Are the measurements and sample metadata suitable for a donor-aware comparison?

**Observed:** Summing raw counts produces a 35,635-by-16 table. The total equals the selected cell-count total, metadata align with columns, and every donor contributes one control and one stimulated sample. The paired design is 16-by-9 with rank 9. Gene filtering retains 3,697 genes. TMM factors range from about 0.794 to 1.251.

**Interpretation:** The table and design encode the intended paired comparison. Filtering removes low-information genes before testing; normalization accounts for library composition through effective library scaling.

**Limit:** Full rank does not separate shared pool effects from condition. Sixteen sample columns are still eight paired donors.

**Decision / next check:** Fit the paired modern edgeR QL model and inspect sample relationships and model variability.

**Evidence:** cd14_pseudobulk_counts_baseline.csv, pseudobulk_sample_metadata.csv, paired_design_matrix.csv, pseudobulk_gene_filter.csv and pseudobulk_normalization.csv.

## 7. Baseline results and visual evidence

**Question:** What condition differences does the fitted model estimate?

**Observed:** 2,506 of 3,697 tested genes have BH FDR <0.05, including 1,252 higher and 1,254 lower in stimulated samples. MDS separates the conditions mainly along its first dimension (76% in that display), with donor 1039 higher on dimension two in both conditions. The variability plot shows raw and moderated estimates around an expression-dependent trend. NEXN is higher in all eight donor pairs; its model logFC is approximately 6.6575.

**Interpretation:** There is a large study-specific condition-associated expression pattern. NEXN illustrates a consistent direction across the inspected donor pairs.

**Limit:** MDS is not a statistical test; 76% is not the fraction of responding genes. The variability diagnostic does not validate all model assumptions. NEXN was selected after inspecting results, and its plotting prior means displayed differences are not the model coefficient. No independent validation or fold-change confidence intervals were produced.

**Decision / next check:** Preserve baseline results and assess the two planned QC alternatives without choosing an analysis merely to maximize discoveries.

**Evidence:** paired_expression_results_baseline.csv, nexn_paired_plot_data.csv and the MDS, variability and NEXN figures.

## 8. Sensitivity findings

**Question:** Do the main conclusions depend strongly on the two chosen QC alternatives?

**Observed:** The >=300 analysis retains 5,254 profiles, tests 3,687 genes and finds 2,497 at FDR <0.05. Excluding seven globin-majority profiles retains 5,348 profiles, tests 3,696 genes and finds 2,525 at FDR <0.05. Among shared tested genes, logFC rank correlations with baseline are 0.9999 and 0.9997; direction agreement is 99.73% and 99.78%. No gene significant in both analyses has opposite directions in either comparison.

**Interpretation:** The overall pattern is stable under these two specific QC choices. Some significance classifications change near the testing threshold.

**Limit:** Shared-gene comparisons omit genes filtered out in an alternative. Ten baseline genes are untested under >=300, including six baseline-significant genes. One baseline gene is untested under globin exclusion and it was not baseline-significant. Both alternatives share the original annotations, normalization and model assumptions.

**Decision / next check:** Retain the >=200 baseline and report both alternatives. Proceed to a clean-session rerun and output comparison.

**Evidence:** all_analysis_summaries.csv, sensitivity_gene_agreement.csv and both baseline_vs_*_by_gene.csv files.

## 9. Rerun and final record

**Question:** Does the checked script reproduce the saved analysis without relying on session objects?

**Observed:** The learner ran the checked script after the instructed restart. The output loaded data from cache, passed all 29 historical checks and recorded MATCH for each of the three full saved result tables within the 1e-8 comparison tolerance. The archive contains these tables, figures, software versions and session information.

**Interpretation:** Rerun consistency is supported for the recorded environment. Review of the CSVs additionally confirmed key alignments, library sums, BH values and sensitivity metrics.

**Limit:** The assistant did not rerun the full dataset in R. Cross-machine portability, future versions, independent biological replication and raw-cell reaggregation by the assistant were not tested. Package-build warnings for generics, dbplyr and matrixStats were reported, although this run completed and matched.

**Decision:** Preserve the checked script and verified outputs; prepare a curated GitHub copy. Publication remains a separate manual step.

## Template for future projects

Question -> observed output -> interpretation -> limitation -> decision -> next check.

Record the script, input object, data unit, parameters, before/after counts and output path. State whether a value was expected or actually observed. Save the first error and the correction. Saving a script with Ctrl+S does not execute it; Ctrl+Enter runs selected code.
