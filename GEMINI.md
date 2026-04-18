# Collaboration Log: RNA-seq Pipeline Enhancement

This file tracks the step-by-step collaborative work between Antigravity and the User to improve the RNA-seq pipeline based on the mentor's scripts.

## Protocol
- **Step-by-Step Execution**: No actions will be taken beyond the current agreed-upon step.
- **Verification**: Each step must be verified before moving to the next.
- **Execution Rule**: Always ask the user to activate the `mamba` environment (`ma cancer_rnaseq`) first, then only provide the `Rscript` command for the user to run. Do not execute visualization scripts automatically.

## RNA-Seq Pipeline System Instructions

The following strict principles govern the maintenance, synchronization, and modification of scripts within this repository. 

### 1. Absolute Identicality (Line-by-Line)
Test scripts **MUST** strictly mirror their corresponding production scripts (`scripts_upstream/`, `scripts_downstream/`) line-by-line. 
- Do not alter the structure of `echo` statements, resource checks, conditional blocks, or native loops.
- Do not append `[TEST]`, `(TEST)`, or any other custom tags to output messages.
- The visual and structural output of a Git diff between the Test script and the Production script must be as minimal as possible.

### 2. Explicit Deviation Tracking
Only when dealing with **Test vs. Production script synchronization**, any fundamental divergence MUST be explicitly flagged.
- The modified line must be immediately preceded by a comment using the exact format: 
  `# --- DEVIATION: [Describe the deviation rationale]`
- If writing an R script, use the standard `#` delimiter.

### 3. No Uncontrolled Automation 
Do not use blanket `sed`, unstructured string replacements, or automated multi-file patching.
- Any synchronization from Prod to Test must preserve the structure of the file seamlessly.

### 4. Preservation of Integrity Checks
Do not strip or relax native integrity checks simply because it is a "test environment". 
- Validations such as `gzip -t` verifications and memory constraints (`JAVA_MEM`) must remain fully intact.

### 5. Bioinformatician Perspective: Subset/Mock Parameter Appropriateness
When evaluating the pipeline through the lens of a **Subset Mock**, production parameters may be mathematically incompatible. The expert role must assess and apply logical overrides:
- **Genome Indexing**: STAR indexing for small genomes MUST use `--genomeSAindexNbases` calculated correctly.
- **Statistical Power**: Significance thresholds (p-values) MUST be adjusted in test scripts to ensure downstream blocks execute for verification.

### 6. HPC Operations (HPC-01)
Operating on the high-performance computing environment (HPC-01) requires extreme caution.
- **MANDATORY REVIEW**: All HPC scripts MUST be scrutinized before execution.
- **RESTRICTED ACCESS**: Strictly adhere to the authorized directories (`/mnt/22T/phongdinh`, `/mnt/22T/vuphong`, `/home/phongdinh`).
- **Job Submission**: Use `sbatch` for substantial loads; negligible loads can use direct bash.
- **Package Management**: Do not automate package installations. Notify the user if a dependency is missing.

---

## Status Log
1. **[2026-04-16]** Initialized `GEMINI.md`. Awaiting first instruction.
2. **[2026-04-16]** **Analysis of Mentor's `featurecounts` folder:**
    - **Structure:** Identical folders for `human` and `dog`.
    - **File Formats:** Standard `featureCounts` outputs (`.tsv`) and metadata files (`.tsv`).
    - **Count Matrices:** Found `human_featurecounts_counts.tsv` (~79k genes, 24 samples). Column 1: `Geneid`, Col 2-6: Annotations, Col 7+: raw counts.
    - **Metadata:** Found `human_sample_metadata.tsv` (26 rows). Maps `sample_id` to conditions (`cell_line`, `treatment`, `drug`, `replicate`).
    - **Gap Discovery:** Mentor's data is already summarized. Your current `01_data_prep.R` expects individual group files in `../_data/counts/`. To use his data, we can create a "bridge" script that reformats his consolidated matrices into your `.RData` pipeline format.
3. **[2026-04-16]** **Strategy Comparison (Species vs. Group):**
    - **User (Group-based):** Highly modular and automated, but siloes data. Misses out on statistical power from cross-cell-line comparisons.
    - **Mentor (Species-based):** Higher statistical power using multifactorial designs (`~ cell_line + treatment`). Enables discovery of *conserved* drug mechanisms. Necessary for Human-Dog comparisons.
    - **Decision:** Pivot toward a **Species-level** modular approach. We will treat "Human" and "Dog" as the main units of analysis to leverage the mentor's multifactorial power while maintaining the user's automated script style.
4. **[2026-04-16]** **Data Bridging Complete:**
    - **Created:** `scripts_downstream/00_bridge.R` to adapt mentor's pre-consolidated matrices.
    - **Execution:** Successfully processed `human` (24 samples) and `dog` (24 samples) data.
    - **Output:** Generated `.RData` objects in `scripts_downstream/.RData/human/` and `scripts_downstream/.RData/dog/` to align with the pipeline's execution context.
5. **[2026-04-16]** **Multifactorial DGE Success:**
    - **Upgraded:** `scripts_downstream/02_deseq2_dge.R` to use `~ cell_line + condition`.
    - **Integrated:** High-performance GTF-to-symbol mapping (using Ensembl 113).
    - **Results:** Generated multifactorial DGE tables for Human and Dog with readable gene symbols (e.g., RCC2, TNFRSF8).
6. **[2026-04-16]** **PCA Alignment Complete:**
    - **Upgraded:** `scripts_downstream/04_01_pca.R` for species-level visualization.
    - **Logic:** Now uses both `color = condition` and `shape = cell_line` + per-cell-line sub-PCAs.
    - **Results:** Confirmed clear separation of cell-lines and consistency of drug-induced shifts.
7. **[2026-04-16]** **Volcano Plot Upgrade:**
    - **Upgraded:** `scripts_downstream/04_02_volcano.R` to reuse results from the DGE step.
    - **Logic:** Ported mentor's `ggplot2` + `ggrepel` labeling heuristics (top genes by significance).
    - **Visuals:** Maintained user's preferred color palette (Red/Blue) while standardizing high-resolution informative output.
