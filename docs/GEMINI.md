# RNA-Seq Pipeline System Instructions

The following strict principles govern the maintenance, synchronization, and modification of scripts within this repository. 

These rules MUST be adhered to by any automated agent or collaborating developer working on the test suite (`test_upstream/`, `test_downstream/`).

## 1. Absolute Identicallity (Line-by-Line)
Test scripts **MUST** strictly mirror their corresponding production scripts (`scripts_upstream/`, `scripts_downstream/`) line-by-line. 
- Do not alter the structure of `echo` statements, resource checks, conditional blocks, or native loops.
- Do not append `[TEST]`, `(TEST)`, or any other custom tags to output messages (`echo`, `message()`, etc.).
- The visual and structural output of a Git diff between the Test script and the Production script must be as incredibly minimal as mathematically possible.

## 2. Explicit Deviation Tracking
Only when dealing with **Test vs. Production script synchronization** (e.g., matching `test_upstream/` to `scripts_upstream/`), any fundamental divergence MUST be explicitly flagged.
- This rule does NOT apply globally to utility or standalone scripts.
- The modified line must be immediately preceded by a comment using the exact format: 
  `# --- DEVIATION: [Describe the deviation rationale]`
- If writing an R script, use the standard `#` delimiter. If shell, use `#`. Do not leave unmarked logic changes.

## 3. No Uncontrolled Automation 
Do not use blanket `sed`, unstructured string replacements, or automated multi-file patching to blindly force changes or "save time."
- Any synchronization from Prod to Test must preserve the structure of the file seamlessly.
- If rewriting a file from scratch (due to heavy damage), utilize `git mv` or a clean clone (`cp`) to maintain history, and immediately re-apply the `# --- DEVIATION:` subsets manually with precision. 

## 4. Preservation of Integrity Checks
Do not strip or relax native integrity checks simply because it is a "test environment". 
- Validations such as `gzip -t` verifications, associative array definitions (`declare -A`), and memory constraints (`JAVA_MEM`) must remain fully intact within the Test sandbox to accurately simulate production resilience.

## 5. Bioinformatician Perspective: Subset/Mock Parameter Appropriateness
When evaluating the pipeline through the lens of a **Subset Mock** (e.g., small read counts mapped to single chromosomes like Human Chr21 or Dog Chr38), production parameters may be mathematically or statistically incompatible. The expert role must assess and apply logical overrides:
- **Genome Indexing**: STAR indexing for small genomes MUST use `--genomeSAindexNbases` calculated as `min(14, log2(GenomeLength)/2 - 1)`.
- **Statistical Power**: Small datasets lack the power of production runs. Significance thresholds (p-values) and minimum gene set sizes (minGSSize) MUST be adjusted in test scripts to ensure downstream code blocks (like Enrichment or Heatmaps) actually execute and produce visual output for verification.
- **Dynamic Routing**: Avoid hardcoded metadata keys or dosage strings in test visualizations to ensure the test suite is resilient to changes in the sample table.
- **Strict Documentation**: All such adjustments **MUST** be explicitly categorized under a `# --- DEVIATION:` block as per Principle 2. Avoid "silent" relaxing of parameters.

## 6. HPC Operations (HPC-01)
Operating on the high-performance computing environment (HPC-01) requires extreme caution and adherence to Slurm-based scheduling.

### 6.1 Security and Scrutiny
- **MANDATORY REVIEW**: All HPC scripts MUST be heavily scrutinized before execution. This is to prevent personal info leaks, accidental resource wastage, and unforeseen system conflicts.
- **RESTRICTED ACCESS**: Strictly adhere to the following directory permissions:
    - `/mnt/22T/phongdinh`: Personal storage (Micromamba packages, etc.).
    - `/mnt/22T/vuphong`: Shared directory for project collaborations.
    - `/home/phongdinh`: Home directory—use ONLY for tools, scripts, and environment management.
    - **PROHIBITED**: Do not access or store data in any other locations unless explicitly authorized.

### 6.2 Job Submission (Slurm)
- **Batch vs. Shell**: Scripts that require task initialization (substantial computational load) MUST be formatted as an `sbatch` script with appropriate headers. Scripts that do not involve task initialization or significant resource use should remain as regular `.sh` files.
- **Interactive Exception**: Only use `srun` or direct bash execution if the resource consumption is negligible compared to the system's total capacity (88 CPUs, 500GB RAM).
- **Package Management**: Micromamba is pre-initialized on this account. **DO NOT** attempt to automate package installations or environment creations. If a script requires a missing dependency, notify the user to install it manually.

### 6.3 Data Management
- **Storage Limits**: Never store large datasets in the home directory (`~`). Utilize the partitioned `/mnt/` volumes (10T2, 12T, 18T, 22T) for active data processing.
