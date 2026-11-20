# Explain this project

## A one-minute explanation

I completed a guided R analysis of public single-cell RNA data from eight donors with lupus, measured under control and IFN-beta-stimulated conditions. I focused on profiles carrying a supplied CD14 monocyte annotation. After reviewing quality measurements, I retained 5,355 profiles and summed raw counts into 16 donor-condition samples. I used edgeR with a paired-donor model. Among 3,697 tested genes, 2,506 met BH FDR below 0.05. Two planned QC alternatives produced very similar fold-change rankings and directions. I also reran the checked workflow and verified agreement with the saved results. This is a learning and reproducibility project with AI-assisted coding. The main limitation is that condition and shared pool effects cannot be separated, and cell identities were not independently validated.

## Questions you may be asked

**Is this proteomics?** No. It measures RNA counts. I did not measure protein abundance or validate protein-level responses.

**Why aggregate cells?** The biological replication comes from donors. Summing counts within each donor and condition gives pseudobulk samples suitable for a donor-aware count model; treating every cell as an independent donor would overstate replication.

**Why include donor in the model?** Each donor contributes both conditions. The donor term accounts for baseline differences between people while estimating a common stimulated-versus-control coefficient.

**What does a positive logFC mean?** Higher estimated relative RNA abundance in stimulated samples under the model. A log2 fold change of 1 corresponds to a twofold ratio; a large estimate still needs scrutiny of expression and assumptions.

**What does FDR mean?** BH is a procedure to control the expected false-discovery proportion under its assumptions. FDR <0.05 does not mean a particular gene has a 5% probability of being false.

**Why was mitochondrial filtering not used?** The MT- rows were present but all counts were zero, including in the original downloaded object. The metric had no discriminatory information. I documented this instead of interpreting zero as good cell quality.

**Why were globin-rich profiles retained in baseline?** I flagged seven descriptively and tested their exclusion separately. Globin fraction alone does not identify the cause or validate an annotation. The main pattern remained very similar under their exclusion.

**Did every gene respond in every donor?** That was not established. The model estimates a common condition coefficient; I directly inspected NEXN as one post hoc paired example. All eight NEXN pairs rose in the normalized display.

**What did the reproducibility check establish?** The checked workflow regenerated results matching earlier saved CSVs in the recorded environment. It did not establish independent biological validation, causal identification or identical results with all software versions.

**What did you do yourself, and where was help used?** I worked through the data units, QC decisions, output interpretation, saved records and rerun in RStudio with guidance. Coding, troubleshooting and documentation were AI-assisted. This project documents developing skills rather than independent mastery.

## Evidence to show during a discussion

1. README results table and the study design.
2. NEXN paired figure, explaining points, lines, scale and post hoc selection.
3. Sensitivity agreement table, explaining its shared-gene denominator.
4. The checked R script and saved reproducibility checks.
5. A notebook entry that separates observation, interpretation and limits.

## A concise portfolio description

Guided paired-donor pseudobulk analysis of public Kang18_8vs8 scRNA-seq data in R and edgeR. Reviewed QC and supplied annotations, aggregated raw counts for eight paired donors, fitted a donor-adjusted QL model, examined two QC sensitivities, and verified agreement on rerun. AI-assisted coding; study-design and annotation limitations documented.
