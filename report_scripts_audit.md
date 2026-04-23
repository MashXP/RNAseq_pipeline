# Pipeline Documentation vs Script Audit Report

This report verifies that the logic, thresholds, and flow documented in `docs/upstream/` and `docs/downstream/` accurately reflect the actual implementation in `scripts_upstream/` and `scripts_downstream/`.

## 1. Upstream Audit

### upstream_README
- **Documentation Verified**: `docs/upstream/README.md`
- **Script Verified**: Directory `scripts_upstream/utils/`
- **Audit Findings**:
  - **Stale Reference**: The README documented `gtf2bed.py` as a shared utility logic script residing in `scripts_upstream/utils/`. However, `gtf2bed.py` only exists in `test_upstream/utils/`.
- **Status**: **FAIL (Documentation Drift)**. Fixed by removing the stale reference from the README to perfectly mirror the production directory.

### 01_genome_prep
- **Documentation Verified**: `docs/upstream/01_genome_prep.md`
- **Script Verified**: `scripts_upstream/01_genome_prep.sh`
- **Audit Findings**:
  - **Safety Flags**: `set -e` and `set -o pipefail` are accurately documented and implemented.
  - **Genomic Sources**: Links perfectly match Ensembl Release 113 for Human and Canine.
  - **Picard refFlat Transformation**: The `awk` parsing logic to convert the 12-column GenePred output into the 11-column refFlat format matches character-for-character.
  - **Trimming Logic (Trimmomatic)**: The JVM memory heuristic (`MEM_GB - 12`) and strict Illumina clip settings (`SLIDINGWINDOW:4:15 MINLEN:36`) are faithfully implemented in the script.
- **Status**: **PASS**. Perfect alignment.

### 02_star_align
- **Documentation Verified**: `docs/upstream/02_star_align.md`
- **Script Verified**: `scripts_upstream/02_star_align.sh`
- **Audit Findings**:
  - **Resource Auto-Detection**: Dynamic allocation of `THREADS` via SLURM and `SORT_RAM_BYTES` calculation (`MEM_GB * 700000000`) perfectly matches the script.
  - **Integrity Validation**: Defensive logic checking compressed FASTQ integrity (`gzip -t`) before alignment is correctly documented.
  - **STAR Parameters**: The complex mapping flags (`--readFilesCommand gunzip -c`, `--outSAMtype BAM SortedByCoordinate`, `--quantMode GeneCounts`) are strictly implemented.
  - **Success Metric Extraction**: `grep` extraction of uniquely mapped reads is accurate.
- **Status**: **PASS**. Perfect alignment.

### 03_alignment_qc
- **Documentation Verified**: `docs/upstream/03_alignment_qc.md`
- **Script Verified**: `scripts_upstream/03_alignment_qc.sh`
- **Audit Findings**:
  - **Resource & Concurrency Management**: Calculation for `PARALLEL_SAMPLES=4` and JVM allocation per thread (`128 / 4 - 4`) mathematically tracks with the script logic.
  - **Dynamic refFlat Selection**: Python `parse_samples.py` combined with grep successfully pulls species info and dynamically points Picard to either the Human or Canine `.refFlat` map.
  - **Picard Parameters**: The execution uses the documented `--STRAND=SECOND_READ_TRANSCRIPTION_STRAND` (matching TruSeq kits) and correctly nulls ribosomal intervals.
  - **Background Job Control**: The `while` loop monitoring `jobs -r` combined with `wait` at the end functions perfectly as a traffic guard.
- **Status**: **PASS**. Fully matches documentation.

### 04_quantification
- **Documentation Verified**: `docs/upstream/04_quantification.md`
- **Script Verified**: `scripts_upstream/04_quantification.sh`
- **Audit Findings**:
  - **Group Iteration**: Utilizing python `parse_samples.py` to identify unique groups and map them to the proper Species GTF is successfully implemented.
  - **FeatureCounts Call**: The strict strand-specific counting logic (`-p -s 2 -t exon -g gene_id`) is correctly applied to group-concatenated BAM strings.
  - **Biotype Classification**: Secondary featureCounts call utilizing `-g gene_biotype` and the downstream `biotype_to_multiqc.py` parsing function identically matches the documented logic.
- **Status**: **PASS**. Logic is robust and identical to the documentation.

### 05_multiqc
- **Documentation Verified**: `docs/upstream/05_multiqc.md`
- **Script Verified**: `scripts_upstream/05_multiqc.sh`
- **Audit Findings**:
  - **Aggregation Call**: The target paths (`BAM_DIR`, `QC_DIR`, `QC_RAW_DIR`, `COUNTS_DIR`) and execution flags (`-o`, `-f`, `-n`, `-c`) match the documentation perfectly.
  - **Discrepancy (Minor/Cosmetic)**: The bash script contains a comment `# 3. RSeQC logs`, but the markdown documentation accurately notes that `RSeQC ... is no longer used`. The code execution itself is correct, so this is just a stale comment in the bash script.
- **Status**: **PASS**. Execution maps exactly to the documented intent.


### utils_dissection
- **Documentation Verified**: `docs/upstream/utils_dissection.md`
- **Script Verified**: `scripts_upstream/utils/parse_samples.py` and `scripts_upstream/utils/biotype_to_multiqc.py`
- **Audit Findings**:
  - **parse_samples.py**: The documentation presents a simplified CSV parsing logic (`next(reader)`). The actual python script contains more robust logic (hunting for the `File1` header and skipping empty rows). The core intent matches perfectly.
  - **biotype_to_multiqc.py**: The documentation shows explicit column dropping (`df.drop(columns=...)`), while the script uses dynamic slicing (`df.columns[5:]`). Again, the script is simply more robust, but the output exactly matches the documented `# plot_type: 'bargraph'` MultiQC requirement.
- **Status**: **PASS (with note)**. The documentation provides simplified pseudocode that correctly explains the logic, while the actual scripts are fully robust.

## 2. Downstream Audit

### downstream_README
- **Documentation Verified**: `docs/downstream/README.md`
- **Script Verified**: Directory `scripts_downstream/`
- **Audit Findings**:
  - **Filename Link**: The README contained a stale link to `00_bridge.md`.
- **Status**: **FAIL (Documentation Drift)**. Fixed by updating the link to `01_bridge_data_prep.md`.

### downstream_libraries
- **Documentation Verified**: `docs/downstream/libraries.md`
- **Script Verified**: Entire `scripts_downstream/` directory
- **Audit Findings**:
  - **Patchwork Rationale**: The documentation listed "Dotplots" as an example of `patchwork` stitching. As discovered in the `04_07` and `04_08` audit, dotplots actually use `facet_wrap`.
- **Status**: **FAIL (Documentation Drift)**. Fixed by replacing "Dotplots" with "NES plots" in the rationale example.

### 01_bridge_data_prep
- **Documentation Verified**: `docs/downstream/01_bridge_data_prep.md`
- **Script Verified**: `scripts_downstream/01_bridge_data_prep.R`
- **Audit Findings**:
  - **Filename Update**: The documentation file has been renamed from `00_bridge.md` to `01_bridge_data_prep.md` to precisely match the script.
  - **Metadata Re-alignment**: The factor leveling logic (`c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM")`) matches identically.
  - **Count Extraction & Filtering**: Slicing the matrix by `rownames(metadata)` and pre-filtering `rowSums(counts_matrix >= 10) >= 3` are perfectly implemented.
- **Status**: **PASS**. Fully matches documentation.

### 01_data_prep
- **Documentation Verified**: `docs/downstream/01_data_prep.md`
- **Script Verified**: `scripts_downstream/01_data_prep.R`
- **Audit Findings**:
  - **Metadata Parsing**: The CSV loading, NA filtering, and specific Group filtering match perfectly.
  - **Desktop-Quality Naming**: `str_replace` and `row_number()` string construction to build the `display_name` string match the code exactly.
  - **Count Matrix Cleaning & Alignment**: `rename_with` combined with `str_remove` to clean the STAR BAM paths from the column headers matches precisely. The column reordering (`counts_matrix[, metadata$OriginalSample]`) is accurate.
  - **Serialization**: Saving logic matches.
- **Status**: **PASS**. Fully matches documentation.

### 02_deseq2_dge
- **Documentation Verified**: `docs/downstream/02_deseq2_dge.md`
- **Script Verified**: `scripts_downstream/02_deseq2_dge.R`
- **Audit Findings**:
  - **High-Performance Gene Mapping**: The `awk` command embedded in `system()` effectively pulls gene_id and gene_name directly from the `.gtf` file just as documented.
  - **Reference Levels**: The control references (`DMSO_Romi` and `H9`/`UL1`) are explicitly set using `relevel()`.
  - **Multifactorial Model**: The `~ cell_line + condition` design is implemented accurately inside an `if` block checking for cell_line existence.
  - **Shrinkage & Subsetting**: `apeglm` shrinkage and loop-based subsetting for consistency comparisons are perfectly aligned.
- **Status**: **PASS**. Fully matches documentation.

### 03_enrichment
- **Documentation Verified**: `docs/downstream/03_enrichment.md`
- **Script Verified**: `scripts_downstream/03_enrichment.R`
- **Audit Findings**:
  - **Species Org.DB**: Accurate switching between `org.Hs.eg.db`/`hsa` and `org.Cf.eg.db`/`cfa`.
  - **ORA Thresholds**: Extracting `sig_genes` using `padj < 0.05, abs(log2FoldChange) > 2` matches the strict documented logic.
  - **GSEA Ranking**: Extracting the Wald-statistic (`stat`) and sorting decreasingly matches the "Gold Standard" claim.
  - **Robust Error Handling**: The `tryCatch` wrappers successfully catch failures in `GSEA` or `gseGO` and output `[CAUTION]` without crashing the loop.
- **Status**: **PASS**. Fully matches documentation.

### 04_01_pca
- **Documentation Verified**: `docs/downstream/04_01_pca.md`
- **Script Verified**: `scripts_downstream/04_01_pca.R`
- **Audit Findings**:
  - **VST Transformation**: `vst(dds, blind = TRUE)` aligns exactly.
  - **PCA Calculations**: Base DESeq2 `plotPCA` with `intgroup` mapping accurately tracks the documentation.
  - **Subset Safety Check**: `if (ncol(vsd_cell) >= 2)` perfectly implements the robustness check for cell-line specific plots.
  - **Discrepancy (Shared Dose Palette)**: The documentation shows a programmatic palette definition using `RColorBrewer`. However, the actual script defines the palette via explicit hex codes (`#E41A1C`, `#377EB8`, `grey85`). The explicit hex mapping in the script is actually better for enforcing absolute consistency, so the script is superior to the doc snippet.
- **Status**: **PASS (with note)**. Script logic is robust; color palette definition in the docs is slightly stale compared to the explicit hex mapping in the codebase.

### 04_02_volcano
- **Documentation Verified**: `docs/downstream/04_02_volcano.md`
- **Script Verified**: `scripts_downstream/04_02_volcano.R`
- **Audit Findings**:
  - **NA Filtering & Infinity Capping**: Removing `NA` values and capping `-log10(0)` to `500` via `case_when` precisely match the script logic.
  - **Significance Coding**: Log2FC thresholds (>= 2.0, <= -2.0) and padj (< 0.05) classifications accurately categorize genes into "up", "down", and "not_significant".
  - **Labeling Heuristic**: Slicing the top 5 genes sorted by `padj` and `desc(abs(log2FoldChange))` strictly mimics the documentation.
  - **Annotations & Ggrepel**: On-plot counts (`annotate`) and non-overlapping label properties (`geom_text_repel` with `max.overlaps = 20`) align exactly.
- **Status**: **PASS**. Fully matches documentation.

### 04_03_venn
- **Documentation Verified**: `docs/downstream/04_03_venn.md`
- **Script Verified**: `scripts_downstream/04_03_venn.R`
- **Audit Findings**:
  - **Directional Filtering**: The function list for `All`, `Up`, and `Down` padj/LFC thresholds (abs(log2FoldChange) > 2) is exactly as documented.
  - **Label Formatting & Clipping**: `str_replace_all`, `expansion(mult = .4)`, and `coord_cartesian(clip = "off")` perfectly implement the label visibility logic.
  - **Statistical Rigor**: `phyper` (p_overlap), Jaccard Index, and Representation Factor math are fully and accurately implemented and saved to `04_03_venn_rigor_stats.csv`.
  - **Shared Genes Extraction**: Extracting `shared_ids` via `Reduce(intersect, sig_list)` matches exactly.
  - **Stacking**: `wrap_plots(venn_plots, ncol = 1)` correctly enforces the vertical stack strategy.
- **Status**: **PASS**. Fully matches documentation.

### 04_04_heatmap_pathway
- **Documentation Verified**: `docs/downstream/04_04_heatmap_pathway.md`
- **Script Verified**: `scripts_downstream/04_04_heatmap_pathway.R`
- **Audit Findings**:
  - **Robustness**: Empty set checking `is.null(hallmark_res)` and capping at `head(..., 6)` pathways match the docs exactly.
  - **Leading Edge Extraction**: Slicing the string (`str_split`) by `/` and grabbing the intersection matches perfectly.
  - **Heuristics**: The gene limitation (`head(12)`) and duplicate handling (`paste0(..., "__", i)`) perfectly match the logic presented.
  - **Z-Score**: The scaling (`t(scale(t(...)))`) and capping (`pmin(pmax(..., -2), 2)`) align perfectly.
  - **Row Splitting**: `row_split = pathway_factor` acts as documented.
  - **Dynamic Sizing**: The `png` height algorithm `250 + nrow(mat_pathway)*25` implements the dynamic scaling discussed in the docs.
- **Status**: **PASS**. Fully matches documentation.

### 04_05_heatmap_variable
- **Documentation Verified**: `docs/downstream/04_05_heatmap_variable.md`
- **Script Verified**: `scripts_downstream/04_05_heatmap_variable.R`
- **Audit Findings**:
  - **Variance Extraction**: Identifying the top 50 via `rowVars` is correctly implemented. Z-score standardisation and capping (+/- 2) align perfectly.
  - **Column Splitting**: Splitting by `conditions_var` with `column_gap = unit(0, "mm")` enforces the documented visual grouping.
  - **Canvas Math**: The algorithm to manually calculate `canvas_h_var` in pixels (`px_per_mm_v`, `mm_rows_v`, etc.) perfectly tracks the script logic.
  - **Discrepancy (Color Mapping)**: The documentation cites `"navy"` and `"firebrick3"` for the `colorRamp2`, whereas the script uses explicit hex codes `"#3B4CC0"` and `"#B40426"`. The hex codes are preferred for precise publication figures.
- **Status**: **PASS (with note)**. Script logic is perfectly aligned; color definition in the docs is slightly simplified.

### 04_06_enrichment_nes
- **Documentation Verified**: `docs/downstream/04_06_enrichment_nes.md`
- **Script Verified**: `scripts_downstream/04_06_enrichment_nes.R`
- **Audit Findings**:
  - **Null Skipping & Top 10 Heuristic**: Checking `!is.null` and `nrow > 0`, followed by `slice_max(..., n = 10)` matches exactly.
  - **Biological Response Logic**: Generating "Activated" and "Suppressed" based on `NES > 0` matches exactly.
  - **Coloring**: The `#D62728` and `#1F77B4` assignments match exactly.
  - **Discrepancy 1 (Faceting)**: The documentation explicitly states that `facet_grid(Response ~ ., scales = "free_y", space = "free_y")` is used to create a clean "Mirror" effect splitting the plot. **The script completely lacks this faceting logic**, relying solely on `reorder(Description, NES)`.
  - **Discrepancy 2 (Scaling)**: The documentation claims the plot dynamically scales height using `max(6, n_rows * 0.38)`. The script actually uses hardcoded heights (`height = 7` or `10`).
- **Status**: **FAIL (Documentation Drift)**. The script is missing the praised `facet_grid` functionality and dynamic sizing described in the documentation.

### 04_07_enrichment_dotplot
- **Documentation Verified**: `docs/downstream/04_07_enrichment_dotplot.md`
- **Script Verified**: `scripts_downstream/04_07_enrichment_dotplot.R`
- **Audit Findings**:
  - **Severe Discrepancy (Aesthetic Mapping)**: The documentation states dots are colored by `pvalue` using a continuous gradient, and sized by `setSize`. The actual script builds a "Mirror Dotplot" where the X-axis dictates Activated/Suppressed, color is mapped statically to `Response` (Red/Blue), and `alpha` is mapped to `-log10(padj)`.
  - **Severe Discrepancy (Stitching)**: The documentation claims `patchwork::wrap_plots` is used for stitching multi-dose comparisons. The script actually relies entirely on `facet_wrap(~ Contrast, nrow = 1)` to generate the combined visual.
  - **Discrepancy (Target Database)**: The doc references extracting `res$gsea_go`. The script extracts `res$gsea_hallmark`.
  - **Discrepancy (String Wrapping)**: The doc praises the use of `str_wrap` to format long pathway names. The script instead uses a custom `clean_pathway_name` function leveraging `tools::toTitleCase` without any line-break wrapping.
  - **Discrepancy (Overcrowding Limit)**: The doc claims there is a fallback to cap the Y-axis at 40 pathways if it gets too crowded. This logic is completely absent from the script, which merely relies on `slice_min(n=12)` per contrast.
- **Status**: **FAIL (Severe Documentation Drift)**. The script's implementation has diverged completely from the documentation's description of a standard ORA/GO dotplot, evolving into a faceted Hallmark Mirror Dotplot without the corresponding documentation updates.

### 04_08_ora_dotplot
- **Documentation Verified**: `docs/downstream/04_08_ora_dotplot.md`
- **Script Verified**: `scripts_downstream/04_08_ora_dotplot.R`
- **Audit Findings**:
  - **Ratio Calculation**: Parsing the `GeneRatio` fraction string (`strsplit`) into a float metric aligns exactly.
  - **Aesthetic Mapping**: `aes(x = Ratio, y = Description, color = p.adjust, size = Count)` is accurately implemented (unlike 04_07).
  - **String Wrapping**: `str_wrap(Description, width = 35)` is accurately implemented.
  - **Discrepancy (Stitching)**: The documentation claims `patchwork::wrap_plots` is used. The script actually uses `facet_wrap(~ Contrast, nrow = 1)`.
  - **Discrepancy (Overcrowding Limit)**: The documentation details a fallback to cap unique pathways at 40 globally. This logic is completely absent from the script.
- **Status**: **FAIL (Documentation Drift)**. While closer than 04_07, it still suffers from the same `facet_wrap` vs `patchwork` divergence and lacks the documented Y-axis limiter.

### 04_09_upset_consistency
- **Documentation Verified**: `docs/downstream/04_09_upset_consistency.md`
- **Script Verified**: `scripts_downstream/04_09_upset_consistency.R`
- **Audit Findings**:
  - **Dynamic Detection**: The `str_extract` Regex identifying cell-line subsets matches precisely, along with the `< 2` error handling.
  - **List Compilation**: The `romi_sets` logic appending `get_sig_ids` lists matches.
  - **UpSet Parameters**: `upset(fromList(all_sets), ...)` matches exactly, including the high-resolution `text.scale`, `point.size`, and `line.size` parameters. Note: Script contains advanced grid viewport rendering to fix clipping, which goes beyond the docs but is functionally superior.
- **Status**: **PASS**. Fully matches documentation.

### 04_10_correlation_plots
- **Documentation Verified**: `docs/downstream/04_10_correlation_plots.md`
- **Script Verified**: `scripts_downstream/04_10_correlation_plots.R`
- **Audit Findings**:
  - **Dynamic Detection**: `< 2` error handling for subset checks matches exactly.
  - **Comparing Cell-Line Results**: The `inner_join` on `Geneid` filtering `log2FoldChange` to `lfc1`/`lfc2` matches exactly.
  - **Consensus Definition**: The `case_when` assigning "Consensus (Sig in both)" and "Specific (Sig in one)" using `padj < 0.05` matches perfectly.
  - **Correlation Calculation**: Computing `cor(..., method = "pearson")` and displaying it dynamically in the subtitle is accurately implemented.
  - **Plotting**: The `ggplot` structure with `geom_smooth(method = "lm")` matches perfectly.
- **Status**: **PASS**. Fully matches documentation.

