# Downstream Dissection: 04_01_pca.R

PCA is the first visualization we look at. It tells us if our treatment groups are actually different and if our replicates "cluster" together like they should.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData` (Contains the `dds` object).
- **Processing**: Variance Stabilizing Transformation (VST), Principal Component calculation.
- **Output**: 
    - **QC Figures**: `../results/[Group]/figures/04_01_pca_combined.png` and `04_01_pca_[CellLine].png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for visualization. See [**libraries.md**](libraries.md) for full technical justifications.
- `DESeq2`: For the `vst()` transformation.
- `ggplot2`: For high-quality PCA plotting.
- `tidyverse`: For data piping and manipulation.
- `RColorBrewer`: For consistent treatment color palettes.

---

## 1. Variance Stabilizing Transformation (VST)
```r
vsd <- vst(dds, blind = TRUE)
```
- **The Job**: Compresses the dynamic range of your gene counts.
- **The Reasoning**: In raw RNA-seq data, a few extremely high-count genes can "drown out" the rest of the data on a PCA plot. **VST** makes the variance constant across all genes, ensuring that the PCA reflects the behavior of the *whole* transcriptome, not just the top 5 genes.

---

## 2. Species-Level Combined PCA
```r
pca_data <- plotPCA(vsd, intgroup = c("cell_line", "condition"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition, shape = cell_line)) +
  geom_point(size = 4, alpha = 0.8)
```
- **The Job**: Extracts the percentage of variance for PC1 and PC2 directly from the DESeq2 calculation, then plots every sample from the species on a single graph.
- **The Reasoning**: 
    - **Shape (Cell Line)**: You should see clear separation between your healthy and aggressive lines (e.g., H9 vs SUPM2). 
    - **Color (Condition)**: You should see "shifts" in position caused by Kromastat and Romidepsin. If the Romi samples are right on top of the DMSO samples, the drug didn't work.

---

## 3. Shared Dose Palette
```r
dose_palette <- c(
  "DMSO_Romi"      = "grey85",
  "Romi_6nM"       = "#E41A1C", # Brighter Red
  "DMSO_Kromastat" = "grey70",
  "Kromastat_6nM"  = "#377EB8"  # Brighter Blue
)
```
- **The Job**: Defines a "Master Color Set" for your treatments.
- **Why it matters**: This ensures that `Kromastat_6nM` is the same color in every single plot (PCA, Volcano, Heatmap). This consistency makes your final research presentation much easier for your mentor to follow.

---

## 4. Sub-PCA & Robustness (Cell-Line Specific)
```r
vsd_cell <- vsd[, colData(vsd)$cell_line == cl]
if (ncol(vsd_cell) >= 2) {
  pca_cell_data <- plotPCA(vsd_cell, intgroup = "condition", returnData = TRUE)
}
```
- **The Job**: Generates a "Zoomed In" PCA for each individual cell line, but safely skips lines that don't have enough replicates to form a valid matrix.
- **The Reasoning**: Sometimes the difference between Two Cell Lines is so huge (PC1 > 90%) that you can't see the drug effect (which might only be 2% of the variance). By running sub-PCAs, we "zoom in" to see exactly how Kromastat is shifting the aggressive cancer line specifically. The `ncol >= 2` check prevents R from crashing if a cell line subset is unexpectedly empty or contains only one sample.

---

## 5. Metadata Integration
- **Input**: Loads the `dds` object from Step 02.
- **Output**: 
    - `04_01_pca_combined.png`: The big-picture overview.
    - `04_01_pca_[CellLine].png`: The high-resolution focus on drug efficacy within specific models.

---

## 6. How to Interpret Your PCA
Principal Component Analysis (PCA) is a "sanity check" for your entire experiment.

1. **PC1 (X-axis) and PC2 (Y-axis)**: The percentages in the labels (e.g., PC1: 65%) represent how much of the "story" of your transcriptome is captured by that axis. PC1 is the most significant source of variation.
2. **Replicate Clustering**: Points from the same condition (e.g., three dots of the same color) should be close together. If one dot is far away from its partners, it may be a "technical outlier" that needs investigation.
3. **Biological Separation**: If the Treatment dots are far away from the DMSO dots along PC1, your drug had a massive, global impact on the cell.
4. **Species vs Treatment**: In combined plots, you will often see PC1 splitting the cell lines (H9 vs SUPM2) and PC2 splitting the treatment (DMSO vs Romi). This is normal; it means the biological difference between cell types is greater than the impact of the drug.
