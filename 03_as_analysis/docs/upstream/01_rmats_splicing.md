# 06 — rMATS-turbo Alternative Splicing Pipeline

This stage performs differential alternative splicing analysis using **rMATS-turbo**. It is designed for high-performance execution on the cluster using an internal task manager.

## Overview

The rMATS pipeline compares two groups of BAM files (e.g., Treatment vs. Control) to identify five types of splicing events:
- **SE**: Skipped Exon
- **RI**: Retained Intron
- **MXE**: Mutually Exclusive Exons
- **A3SS**: Alternative 3' Splice Site
- **A5SS**: Alternative 5' Splice Site

## Workflow

### 1. BAM List Generation
Before running rMATS, the script executes `utils/generate_bam_lists.py`. This utility parses the project metadata and generates the necessary input files. For a detailed breakdown of the logic, see [Utility Dissection](../../../utils/utils_dissection.md#3-generate_bam_listspy).

### 2. Parallel Execution (Picard-Style)
The script `06_rmats_splicing.sh` manages multiple comparisons concurrently.
- It detects available CPUs (via `SLURM_CPUS_PER_TASK`).
- It runs a set number of comparisons in parallel (`PARALLEL_TASKS`), assigning a subset of threads to each.
- This ensures maximum throughput without exceeding resource allocations.

### 3. Splicing Analysis (prep & post)
For each comparison, rMATS is run in two stages:
- **Prep**: Processes BAM files and annotations.
- **Post**: Performs the statistical comparison.
- **Parameters**:
  - `libType`: `fr-firststrand` (default for this project).
  - `readLength`: 150bp.
  - `novelSS`: Enabled to detect unannotated splice sites.

### 4. Significance Filtering
Following the rMATS run, an `awk` filter is applied to the `.JCEC.txt` files to extract high-confidence events:
- **FDR** <= 0.05
- **|DeltaPSI|** (IncLevelDifference) >= 0.1
- **Average Counts** (IJC + SJC) >= 10 across replicates.

## Output Structure

Results are organized in `_data/rMATS-turbo/results/{CellLine}/{Comparison}/`:
- `filtered/`: Contains the `.filtered.txt` files for each AS type.
- `raw/`: Original rMATS output files.
- `intermediate/`: Internal rMATS files used for `prep` and `post` tasks.

## Usage

```bash
# Ensure the mamba environment is active
mamba activate cancer_rnaseq

# Run the pipeline (usually via sbatch)
bash 01_core_alignment/scripts/06_rmats_splicing.sh
```
