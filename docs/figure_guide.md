# How to read the saved figures

| Figure under results/figures | What it shows | Supported observation | What it cannot establish |
|---|---|---|---|
| documented_study_design.png | Bars count documented donor-condition rows, not cells. | Eight rows per condition, representing the same eight people. | Gene expression or a biological effect. |
| singlet_retention_percent.png | Percent of original profiles retained by the supplied singlet label, by donor and condition. | Retention varies by group; both conditions remain represented for all donors. | Cell health or a stimulation response. |
| counts_vs_detected_genes.png | One point per supplied singlet; raw UMI total on x, detected genes on y. The reviewed minimum-gene profile is red. | A low-gene tail exists; count depth alone does not ensure broad gene detection. | The cause or identity of every unusual profile. |
| cd14_qc_by_donor_condition.png | Count-depth and detected-gene boxplots in supplied CD14-labelled singlets before baseline gene filtering; both y axes are logarithmic. | Donor 1039 has low-gene tails in both conditions; donor 1488 control has a lower central detected-gene distribution. | Automatic exclusion rules or validated cell identities. |
| cd14_pseudobulk_mds.png | Sixteen donor-condition pseudobulk samples, distances based on leading log-fold changes; grey control, blue stimulated. | Conditions separate chiefly on dimension 1; donor 1039 is high on dimension 2 in both conditions. | Causal response, validated normalization or a fraction of responding genes. The 76% label describes the displayed MDS dimension. |
| paired_model_variability.png | Each gene has raw black and moderated red variability estimates; blue is the expression-dependent trend. | Moderation reduces noisy scatter while retaining some highly variable genes. | A pass/fail certificate for every model assumption. |
| nexn_paired_expression.png | One point per donor-condition pseudobulk sample; lines pair donors. Y is TMM-normalized log2 CPM with prior.count=2. | NEXN is higher in stimulated samples for every inspected donor pair. | Independent validation, absolute RNA per cell, protein abundance or the plotted line difference being identical to the model logFC. |

Some points or labels overlap. Use the source tables for exact values. The original rerun retention export is also preserved as singlet_retention_percent_rerun.png; its x-axis title is partly clipped. The display export has a readable title and the same underlying retention percentages.
