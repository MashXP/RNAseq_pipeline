# Bioinformatics Guide & Benchmarks

This document serves as a reference for bioinformatics best practices, quality thresholds, and statistical integrity standards established for the RNA-seq pipeline.

## Quality Benchmarks

### FASTQ Initial Check
*   **Format**: Headers should start with `@`, followed by sequence and quality.
*   **Quality Encoding**: Standard modern Illumina is **Phred+33**.
*   **Observation**: `N` bases at the start of reads are a common Illumina artifact and are handled by trimming.

### Trimming (Trimmomatic)
*   **Both Surviving (PE)**: **70% - 90%** is considered excellent.
*   **Dropped Reads**: Should be **< 5%**.
*   **Indicator**: If "Both Surviving" is < 50%, you likely have high adapter contamination or severe library quality issues.

### Alignment (STAR)
*   **Production (Full Genome)**: 
    *   **70% - 95%** Uniquely Mapped reads is considered good.
    *   **< 50%** suggests potential sample contamination, incorrect species reference, or severe degradation.
*   **Test Suite (Chr21 Only)**:
    *   **1.5% - 3.0%** Unique mapping is **EXPECTED**.
    *   Since Chr21 is ~1.5% of the genome, only ~1.5% of genome-wide reads should map to it. 85%+ should be categorized as "unmapped: too short".

### Quantification (featureCounts)
*   **Assignment Rate**: **60% - 80%** of mapped reads should assign to features (genes).
*   **Zero counts**: If most genes have 0 counts, check if your GTF `gene_id` matches the BAM file attributes.

### QC (Picard)
*   **Read Distribution**: For RNA-seq, the majority of reads (>80%) should map to **Exons**. High Intronic/Intergenic mapping indicates DNA contamination.
*   **Gene Body Coverage**: Coverage should be relatively flat. A sharp peak at the 3' end indicates RNA degradation.

## Statistical Integrity
*   **Ref Level**: Always ensure `DMSO_Romi` (or primary control) is set as the reference level using `relevel(condition, ref = "DMSO_Romi")` in R.
*   **GSEA Ranking**: Use `sign(log2FoldChange) * -log10(pvalue)` to preserve both direction and significance.
*   **P-values**: Always use **unshrunken** p-values for ranking, as shrinkage models (like apeglm) often focus on LFC and do not provide a direct p-value output.
