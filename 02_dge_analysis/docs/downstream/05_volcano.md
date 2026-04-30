# Downstream Dissection: 04_02_volcano.R

The Volcano plot is the "Money Shot" of RNA-seq. It identifies the most heavily changed and most statistically confident genes (like your drug targets) in a single visual.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData` (Contains results_list with pre-mapped symbols).
- **Processing**: P-value capping at 500, fold-change significance classification (LFC > 2.0), and top-gene labeling.
- **Output**: 
    - **DGE Figures**: `../results/[Group]/figures/04_02_volcano_[Contrast].png` and `04_02_volcano_combined.png`.
    - **Summary Tables**: `../results/[Group]/tables/04_02_dge_summary_stats.csv` and `04_02_top_dge_genes.csv`.
    - **Conserved Targets**: `../results/[Group]/tables/04_02_conserved_targets.csv` (H9 ∩ SUPM2 overlap).

---

## 0.1 Library Rationales
This script utilizes the following libraries for visualization. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggplot2`: Base engine for plotting.
- `ggrepel`: Critical for non-overlapping gene labels.
- `patchwork`: For stitching multiple volcano plots into a single report.
- `tidyverse`: Data manipulation and text formatting.

---

## 1. Contrast Iteration & NA Filtering
```r
for (contrast in names(results_list)) {
  plot_df <- results_list[[contrast]]$df
  
  plot_df <- plot_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
```
- **The Job**: Loops through all treatment comparisons and strips out any genes that lack statistical data.
- **The Reasoning**: Sometimes DESeq2 filters out low-count genes and sets their `padj` to `NA`. If these are left in, `ggplot` will throw warnings or errors. This step ensures a clean, error-free dataset for plotting.

---

## 2. Handling Mathematical Infinities
```r
plot_df <- plot_df %>%
  mutate(
    neg_log10_padj = case_when(
      padj == 0 ~ 500,
      -log10(padj) > 500 ~ 500,
      TRUE ~ -log10(padj)
    )
  )
```
- **The Job**: Prevents "infinite" p-values from breaking the plot.
- **The Reasoning**: Sometimes DESeq2 is so confident in a gene (e.g., a massive drug response) that it reports a p-value of 0. Mathematically, the `-log10(0)` is infinite. This block caps the visual significance at 500, ensuring the plot remains beautiful and displayable while still showing that the gene is "off the charts" significant.

---

## 3. Color Coding (Significance Status)
```r
sig_status = case_when(
  padj < 0.05 & log2FoldChange >= 2.0  ~ "up",
  padj < 0.05 & log2FoldChange <= -2.0 ~ "down",
  TRUE ~ "not_significant"
)
```
- **The Job**: Categorizes every gene into your project's standardized significance bins.
- **The Reasoning**:
    - **Up (Red)**: Genes activated by the drug (Krom/Romi).
    - **Down (Blue)**: Genes suppressed by the drug.
    - **NS (Grey)**: The "Background" transcriptome that didn't change.

---

## 4. Top Gene Labeling (Heuristic)
```r
label_df <- plot_df %>%
  filter(sig_status != "not_significant") %>%
  group_by(sig_status) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 5)
```
- **The Job**: Automatically identifies the 5 "Most Important" genes on each side (Up/Down) to label with text.
- **The Reasoning**: A plot with 20,000 labels is unreadable. Your mentor wants to see the "Big Hits." This logic picks the genes with the lowest p-values AND the highest fold changes, ensuring the most biologically relevant targets (like cell cycle regulators or apoptosis initiators) are named on the graph.

---

## 5. On-Plot Summary Annotations
```r
annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", down_count), ...)
annotate("text", x = Inf, y = Inf, label = paste0("Up: ", up_count), ...)
```
- **The Job**: Prints the total number of significant genes directly in the top corners of the plot.
- **Why it matters**: It provides an "At-A-Glance" summary. Instead of looking at a separate Excel table to find how many genes Kromastat changed, you can read the total count directly from the image.

---

## 6. Visual Rigor (`ggrepel`)
```r
p <- p + geom_text_repel(
    data = label_df, aes(label = gene_label),
    max.overlaps = 20, force = 2, clip = "off"
)
```
- **The Job**: Uses smart placement to ensure gene names don't overlap with each other or the data points.
- **The Reasoning**: By setting `clip = "off"`, we prevent labels from being cut off if they land near the edge of the image. This ensures every gene name is 100% legible for your final presentation.

---

## 7. Authoritative DGE Exports & Conserved Targets
```r
write.csv(final_summary_df, file.path(res_dir, "tables/04_02_dge_summary_stats.csv"))
write.csv(final_conserved_df, file.path(res_dir, "tables/04_02_conserved_targets.csv"))
```
- **The Job**: Exports three critical tables for your report:
    1. **Summary Stats**: Total Up/Down counts per contrast.
    2. **Top Genes**: The top 5 Up/Down genes for each comparison.
    3. **Conserved Targets**: The genes that are significant in BOTH the healthy line (H9) and the cancer line (SUPM2) for a given treatment.
- **The Reasoning**: The "Conserved Targets" table is particularly high-value as it highlights genes that respond to the drug regardless of the cell line—the "Core Mechanism."

## 8. How to Interpret Your Volcano Plot
The Volcano plot helps you balance **magnitude** (how much it changed) with **significance** (how sure we are).

1. **The X-axis (log2FoldChange)**: Represents the effect size. 
    - **0**: No change.
    - **2**: 4x expression (Upregulated).
    - **-2**: 4x suppression (Downregulated).
    - **>5**: Massive activation (often seen in stress response or apoptosis).
2. **The Y-axis (-log10 Adjusted p-value)**: Represents statistical confidence. 
    - The higher the point, the more "sure" we are that the change is real and not just a random fluctuation.
    - Points above the horizontal line (usually $padj < 0.05$) are statistically significant.
3. **The Upper-Right Quadrant**: These are the **Upregulated** genes you likely care about (potential drug targets or activated pathways).
4. **The Upper-Left Quadrant**: These are the **Downregulated** genes (genes being suppressed by the treatment).
5. **Labels**: The genes labeled in the top corners are your "Top Hits"—they have the best combination of huge change and high confidence.
