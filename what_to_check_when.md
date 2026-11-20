# What to check, when and why

Use these checks to decide whether the next step is supported. They describe the completed six-part project; they are not universal QC thresholds.

| Part | Check before proceeding | Evidence in this project | Decision or limit |
|---|---|---|---|
| 1. Design | What is a profile, sample and donor? Are measurements paired? | Eight donors, both conditions, 16 donor-condition combinations. | Preserve donor pairing. Shared condition-specific pooling remains confounded. |
| 2. Quality | Are supplied flags usable? Do metrics reveal unusual profiles? | Singlet labels, detected genes, count depth, globin review and all-zero MT rows. | Keep supplied singlets and CD14 labels; baseline >=200 genes; no mitochondrial threshold; annotations remain provisional. |
| 3. Aggregation | Are raw counts being summed and metadata aligned? | 35,635-by-16 pseudobulk table, preserved count total, positive group sizes and exact column alignment. | Proceed with donor-condition aggregates, not cell-level pseudo-replication. |
| 4. Model | Is donor encoded? Is the design full rank? What remains confounded? | 16-by-9 design, rank 9; group-based expression filtering, TMM, MDS and QL diagnostic. | Use paired model; full rank cannot separate pool from condition. Diagnostics do not prove all assumptions. |
| 5. Results | What is the direction and magnitude? Are findings sensitive to choices? | Full logFC/FDR tables, paired NEXN display and two QC alternatives. | Report estimates and limits; no confidence intervals or independent validation were generated. |
| 6. Reproduction | Can saved code recreate results and evidence? | 29 historical checks and three matching result tables; software and session records. | Rerun consistency verified in the recorded environment; cross-machine/version portability was not tested. |

## Before any new command

Ask: what question does this step answer; what data unit does it use; which object holds it; what result would change the next decision?

## If an error appears

Save the first exact error. Check that the creating command ran before commands that print, summarize or save its object. Select a complete dependent block. Ctrl+S saves the script; Ctrl+Enter runs selected code. Package installation is separate from loading a package in a fresh session.

## After a figure appears

State the axes, units, point/bar meaning, colors and transformations. Record one directly visible observation before proposing an explanation. A display can illustrate a finding without independently validating it.

## Carry-forward limits

RNA is not protein. Supplied labels are not validated identities. More cells cannot replace donors. Reproducing output cannot resolve confounding. Agreement under two QC choices does not validate assumptions shared by both analyses.
