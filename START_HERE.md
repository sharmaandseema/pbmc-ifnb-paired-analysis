# Start here

This folder is a curated copy of the completed PBMC IFN-beta learning project. It does not change your original PBMC_IFNB_project folder.

## Read and keep the work

1. Extract the downloaded ZIP using Windows **Extract All**.
2. Keep the extracted PBMC_IFNB_GitHub folder in Documents or another chosen location.
3. Open docs/PBMC_IFNB_Observation_Guide.pdf for the project explanation.
4. Open README.md on GitHub or in RStudio for the overview; docs/analysis_notebook.md contains the detailed observation record.

The six analysis parts and the matched rerun are complete. GitHub publication remains a separate step. There is no need to repeat the full analysis merely to read or upload this package.

## Publish to a new GitHub repository

Suggested repository name: pbmc-ifnb-paired-analysis

Suggested description: Paired-donor pseudobulk RNA analysis of public Kang18 data in R and edgeR, with QC sensitivity checks and a verified rerun.

1. Sign in to the sharmaandseema account and create a new repository with the suggested name or your preferred name.
2. Select Public if you want it visible as a portfolio project. Leave the automatic README, .gitignore and license options unchecked because the package supplies its own README and .gitignore; no reuse license has been selected here.
3. Use **uploading an existing file**, or **Add file > Upload files** in an existing repository.
4. Open the extracted PBMC_IFNB_GitHub folder. Drag its CONTENTS (R, docs, data, results and the root files) onto the upload page. Do not upload the ZIP or the outer folder itself as a single file.
5. Check that README.md and PBMC_IFNB_project.Rproj are at the repository root and that .gitignore is included. Commit with a message such as: Add paired IFN-beta RNA analysis and reproducibility records.
6. Review the rendered README, open its NEXN figure and check the docs links. Pin the repository on your profile if you want it featured.

Nothing in this download has been uploaded to GitHub automatically.

## Rerun only when needed

Open this folder's PBMC_IFNB_project.Rproj. If the packages are already installed, restart R and run source("R/03_rerun_checked.R"). A fresh results/rerun_* folder is created and compared with the packaged reference tables. The README explains optional installation and reproducibility limits.

Keep your original downloaded-data/cache backup separately if you want to preserve it. The GitHub copy excludes the full original cell-count matrix, cache and RDS workspaces.
