# Methods and interpretation

## Question and dataset
How do relative RNA profiles of supplied CD14-labelled monocytes differ between control and IFN-beta-stimulated samples across eight donors?
Public Kang18_8vs8 dataset from muscData; supplied donor, condition, cell-type and multiplet metadata were used.

## Analysis
Kept supplied singlet-labelled profiles and the supplied CD14+ Monocytes annotation.
Baseline: at least 200 detected genes. No additional mitochondrial threshold was applied.
Summed raw counts by donor-condition group, producing 16 pseudobulk samples.
Filtered genes with filterByExpr(group=Condition, min.count=10, min.total.count=15).
Applied TMM and the paired design ~ Donor_ID + Condition.
Used glmQLFit(robust=TRUE, legacy=FALSE), glmQLFTest(coef='Conditionstim'), and BH FDR <0.05.
Positive logFC denotes higher relative RNA in stimulated samples.
Two separate QC sensitivity analyses used >=300 detected genes, or excluded baseline profiles with HBB+HBA1+HBA2 above 50% of counts.
Per-gene sensitivity comparisons include only genes tested in both analyses.

## Observations and decisions
Baseline retained 5355 profiles and tested 3697 genes; 2506 met FDR <0.05 (1252 higher, 1254 lower in stimulated).
GE300: 5254 profiles, 3687 tested genes and 2497 at FDR <0.05.
Globin sensitivity: 5348 profiles, 3696 tested genes and 2525 at FDR <0.05.
See sensitivity_gene_agreement.csv for observed direction agreement, rank correlations and significance changes.
Retain the >=200 analysis as the baseline, report both alternatives, and flag any mismatch with historical checks for review.

## NEXN example
NEXN was selected after reviewing the baseline results; it is an illustration, not independent validation.
Each plotted point is a donor-condition pseudobulk sample, linked within donor. Display uses TMM-normalized log2 CPM with prior.count=2.

## Limitations
The supplied 13 MT- gene rows were all zero in the historical object, including before singlet filtering. The upstream cause is unresolved; zero values do not establish healthy cells.
Cell-type labels were supplied, not independently re-annotated. The small marker review does not validate all identities.
Genetic multiplet flags do not identify every same-donor doublet. No ambient-RNA correction was performed.
The three-globin fraction is a descriptive count fraction, not an estimated contamination fraction or validated threshold.
Donor aliquots were pooled by condition into one control and one stimulated pool. Donor blocking cannot separate condition from shared pool effects.
These are study-specific relative RNA associations, not protein measurements, absolute RNA per cell, or an unconfounded causal effect.
No new biological discovery or independent validation is claimed. Coding and troubleshooting were AI-assisted.

## Sources
https://bioconductor.org/packages/release/data/experiment/html/muscData.html
https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
https://stacks.cdc.gov/view/cdc/79371/cdc_79371_DS1.pdf


## Verified rerun and saved-file review

On 20 November 2026, the learner reported completion after the instructed R restart. The archived output contains 29 matching historical checks and three MATCH rows for the full saved gene-result tables. The analysis itself was run on the learner's Windows computer, not in the assistant's environment.

Review of the submitted archive additionally checked metadata-column alignment, pseudobulk library totals, paired design rank, QC subset membership, BH adjustments from the saved p-values, sensitivity correlations and all eight NEXN directions. These CSV checks do not reproduce the raw-cell aggregation independently because the full gene-by-cell matrix is not in the archive.

## Statistical interpretation

The model includes donor effects and one common condition coefficient. It does not fit eight donor-specific condition coefficients or test response heterogeneity. The NEXN figure is the directly inspected donor-level example; the DE list alone does not establish agreement in every donor for every gene.

BH FDR <0.05 describes a multiple-testing procedure under its assumptions. It is not a 5% probability that an individual result is false and does not imply that exactly 5% of this observed list is false. No minimum fold-change threshold was imposed in the reported count of discoveries.

TMM changes effective library scaling rather than the original raw counts. All alternatives use the same normalization and count-model framework; their agreement is limited evidence about the selected QC choices.
