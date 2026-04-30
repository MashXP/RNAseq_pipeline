# Test Script Deviations

This document tracks all fundamental divergences between the **Test Suite** scripts and the **Production** scripts. These deviations are necessary to accommodate small mock datasets (e.g., Chromosome 21 subsets) and local execution environments.

## 01_core_alignment (Test)

### 01_test_genome_prep.sh
- **Pathing**: Outputs to `_data/genome_test`, `_data/index_test`, and `_data/fastq_trimmed_test`.
- **Targeting**: Downloads Chromosome 21 (Human) and Chromosome 38 (Dog) FASTA instead of primary assemblies.
- **Filtering**: Filters Ensembl GTF to keep only Chr 21/38 annotations to match FASTA.
- **Resources**: Adjusted `runThreadN` to 10 for local machine specs.
- **STAR Parameters**: Appended `--genomeSAindexNbases 11` for small subset genomes.

### 02_test_star_align.sh
- **Pathing**: Redirects output to `_data/bam_test` and reads from `_data/fastq_trimmed_test`.
- **Resources**: Fallback changed from 32 to 10 threads to match local machine specs.

### 03_test_alignment_qc.sh
- **Pathing**: Redirects results to `_data/qc_test`.
- **Memory**: Lower memory allocation for local test environment (`JAVA_MEM_PER_SAMPLE=4`).
- **RefFlat**: Selects test subset RefFlat (`chr21.refFlat`, `chr38.refFlat`) instead of primary assembly.

### 04_test_multiqc.sh
- **Pathing**: Points to test output directories for BAM, QC, and Counts.
- **Naming**: Report filename includes `_test` suffix.

---

## 02_dge_analysis (Test)

### 01_test_quantification.sh (Upstream)
- **Pathing**: Outputs to `_data/counts_test`.
- **Resources**: Fallback changed from 16 to 10 threads.
- **GTF**: Selects test subset GTF instead of primary assembly.

### 01_test_data_prep.R (Downstream)
- **Pathing**: Reads from `../_data/counts_test/` and saves to `./.RData_test/`.
- **Filtering**: Relaxed raw count filtering (`rowSums(counts_matrix >= 1) >= 2`) for small mock datasets.

### 02_test_deseq2_dge.R (Downstream)
- **Pathing**: Redirects results to `../results_test/` and loads from `./.RData_test/`.

### 03_test_enrichment.R (Downstream)
- **Pathing**: Redirects results to `../results_test/` and loads from `./.RData_test/`.
- **Thresholds**: Relaxed raw p-value ORA threshold for small mock dataset.
- **Robustness**: Wraps GSEA in `tryCatch` to handle sparse chromosome subsets without crashing.

### Visualizations (04-16_test_XX)
- **Pathing**: All scripts redirect results to `../results_test/` and load from `./.RData_test/`.
- **Transformation**: Use `varianceStabilizingTransformation` directly for small mock datasets (blind = FALSE).
- **Venn/Volcano**: Relaxed thresholds (e.g., `padj < 0.2`) to ensure figures are generated for verification.
- **Heatmaps**: Always take top 6 by raw p-value to ensure figures are generated.

---

> [!NOTE]
> These deviations are flagged in the code with `# --- DEVIATION: [rationale]`. 
> When re-implementing test scripts, ensure these overrides are preserved to maintain compatibility with the Chr21/38 mock data.
