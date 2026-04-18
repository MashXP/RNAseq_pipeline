# Downstream Dissection: 04_05_heatmap_variable.R

While the "Pathway Heatmap" looks at specific biological groups, this script looks at the **50 Most Variable Genes** across your entire experiment. It is the best way to see the "Raw Power" of the drugs without any biological bias.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `scripts_downstream/.RData/[Group]/02_deseq_results.RData`.
    - **Species Org.DB**: `org.Hs.eg.db` or `org.Cf.eg.db`.
- **Processing**: VST transformation, identification of top 50 genes by `rowVars`, Z-score standardization, column-splitting by `condition`.
- **Output**: 
    - **Global Heatmap**: `../results/[Group]/figures/04_heatmap_top_variable.png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for visualization. See [**libraries.md**](libraries.md) for full technical justifications.
- `ComplexHeatmap`: Essential for treatment-based column splitting.
- `circlize`: For high-contrast Z-score color mapping.
- `DESeq2`: For VST transformation.
- `tidyverse`: Data cleaning and mapping.

---

## 1. Variance Extraction & Standardization
```r
top_genes_idx <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat_var       <- assay(vsd)[top_genes_idx, ]

# Z-score standardization and capping
mat_var       <- t(scale(t(mat_var))) 
mat_var       <- pmin(pmax(mat_var, -2), 2)
```
- **The Job**: Identifies the 50 genes that show the largest changes in expression across all your samples, then standardizes them to a -2 to +2 scale.
- **The Reasoning**: This is a "blind" search. It doesn't care about your metadata; it just cares about where the action is. The Z-score standardization is critical because raw count variances can span thousands of units; standardizing ensures that both high-abundance and low-abundance variable genes are equally visible on the final color map.

---

## 2. Column Splitting (By Treatment)
```r
ht_var <- Heatmap(
  mat_var,
  column_split = factor(conditions_var, levels = unique(conditions_var)),
  column_gap = unit(0, "mm")
)
```
- **The Job**: Groups the columns of the heatmap by their Treatment (DMSO, Romi, Krom).
- **Why it matters**: It makes the heatmap highly organized. You can clearly see a "Block" of Red (High expression) in the drug groups and a "Block" of Navy (Low expression) in the DMSO groups. This proves the drug effect is consistent across replicates.

---

## 3. High-Contrast Color Mapping
```r
col_fun_var <- colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))
```
- **The Job**: Maps the Z-scores to a "Fire and Ice" color scale.
- **The Reasoning**: Standard R colors can be muddy. By using explicit hex codes like `#3B4CC0` (deep blue) and `#B40426` (deep red), we create a high-contrast image that looks professional in a journal article or a mentor's slide deck.

---

## 4. Precise Canvas Math (Formatting)
```r
px_per_mm_v <- 150 / 25.4
mm_rows_v    <- nrow(mat_var) * 6
canvas_h_var <- round((mm_title + ... + mm_rows_v) * px_per_mm_v)
```
- **The Job**: Manually calculates the exact pixel height and width of the final PNG image.
- **The Reasoning**: Standard `ggsave` often squashes heatmaps or makes the text unreadable. This code ensures that every row of the heatmap is exactly 6mm tall, preventing any blurring or distortion of your results.

---

## 5. How to Interpret Your Variable Gene Heatmap
This heatmap provides an **Unbiased Global View** of your experiment. Unlike pathway heatmaps, it shows the "top 50 action items" regardless of their known biological function.

1. **Horizontal Strips (Gene Signatures)**: Look for genes that are **Red** (High) in one treatment group and **Blue** (Low) in another. This defines the "Unique Signature" of that drug.
2. **Column Consistency**: All replicates (e.g., three DMSO columns) should look nearly identical. If one replicate looks like a treatment sample, you may have a "sample swap" or technical error.
3. **Conserved Responses**: If you see a gene that is Red across both Kromastat and Romidepsin, but Blue in DMSO, you have found a **Conserved Mechanism of Action**. This is high-value data for proving the drugs behave similarly.
4. **Primary Drivers**: Because these are the 50 most variable genes in the *entire* study, they represent the "Primary Drivers" of the cellular response. Check the symbols (e.g., *CDKN1A*) to see if they match your expected therapeutic targets.
