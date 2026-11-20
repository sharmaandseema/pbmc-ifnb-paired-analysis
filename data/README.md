# Data source and storage

Original public data: Kang et al. (2018), GEO GSE96583, provided as muscData::Kang18_8vs8() through ExperimentHub EH2259. The recorded object contains a counts assay with 35,635 features and 29,065 profiles. Metadata include ind, stim, cluster, cell and multiplets.

The successful rerun used cached data with ExperimentHub snapshotDate 2026-04-21. The accessor can download the object on a new computer and normally reuses its cache afterward. The cache is outside the project unless separately configured.

Full original counts, cache files and RDS workspaces are not bundled. The small derived raw-count pseudobulk matrix is available under results/tables/cd14_pseudobulk_counts_baseline.csv. Its values sum counts over selected profiles and are not log-normalized expression.

Public donor IDs are dataset labels. This repository did not collect patient samples or generate the original sequencing data. See the original data source for applicable usage terms.
