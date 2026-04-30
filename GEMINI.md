# Collaboration Log: RNA-seq Pipeline Enhancement

This file tracks the step-by-step collaborative work between Antigravity and the User to improve the RNA-seq pipeline based on the mentor's scripts.

## Protocol
- **Step-by-Step Execution**: No actions will be taken beyond the current agreed-upon step.
- **Verification**: Each step must be verified before moving to the next.
- **User-Led Progression**: Never proactively move to the next phase/pillar or ask to move on. Wait for the user's explicit command to proceed.
- **Research Constraints**: Do not run terminal commands for independent data verification or research unless explicitly requested as a specific step.
- **Execution Rule**: Always ask the user to activate the `mamba` environment (`ma cancer_rnaseq`) first, then only provide the `Rscript` command for the user to run. Do not execute visualization scripts automatically.

## RNA-Seq Pipeline System Instructions

The following strict principles govern the maintenance, synchronization, and modification of scripts within this repository. 

### 1. Absolute Identicality (Line-by-Line)
Test scripts **MUST** strictly mirror their corresponding production scripts (located in `01_core_alignment/`, `02_dge_analysis/`, and `03_as_analysis/`) line-by-line. 
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

### 7. Documentation-First Integrity
Whenever a script is being modified or created, you **MUST** refer to the corresponding documentation (located in the `docs/` subdirectory of each module) to maintain technical and logical integrity. 
- Ensure the script implementation aligns exactly with the documented purpose, mathematical parameters, and workflow sequence.
- If a modification creates a discrepancy with existing documentation, the documentation MUST be updated simultaneously to reflect the new truth.
