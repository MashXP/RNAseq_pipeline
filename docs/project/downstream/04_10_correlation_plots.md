# Downstream Dissection: 04_10_correlation_plots.R

This is the "Proof of Mechanism" script. While other plots show that genes changed, this plot proves that they changed in the **exact same direction** and with the **same strength** across different biological models.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData` (Contains cell-line specific DGE results).
- **Processing**: Inner-joining LFCs between two cell lines, Consensus classification (LFC > 2.0), Pearson correlation calculation.
- **Output**: 
    - **Consistency Figure**: `../results/[Group]/figures/04_10_lfc_correlation.png` (Side-by-side Romi vs Krom with top gene labels).

---

## 0.1 Library Rationales
This script utilizes the following libraries for directional analysis. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggplot2`: For LFC-LFC scatter plotting.
- `ggrepel`: For non-overlapping gene labels.
- `patchwork`: For horizontal stitching of drug-specific correlation plots.
- `tidyverse`: For `inner_join` and data filtering.

## 1. Dynamic Cell-Line Detection & Robustness
```r
keys <- names(results_list)
known_cl <- c("H9", "SUPM2", "CNK89", "UL1")
cell_lines <- intersect(known_cl, unique(str_extract(keys, "^[A-Za-z0-9]+(?=_)")))

if (length(cell_lines) < 2) {
  message("\n[CRITICAL] Error: Could not find cell-line subset results...")
  stop("Missing subset results.")
}
```
- **The Job**: Uses regex to automatically detect which cell lines were processed in your dataset by looking at the prefixes of the contrast names.
- **The Reasoning**: A correlation requires exactly two things to compare. If the pipeline was run on a single cell line, LFC correlation is mathematically impossible. This block protects the script from crashing silently by verifying that at least two cell line subsets exist before attempting the math.

---

## 2. Comparing Cell-Line Results
```r
merged <- inner_join(
  df1 %>% select(Geneid, gene_label, lfc1 = log2FoldChange, padj1 = padj),
  df2 %>% select(Geneid, lfc2 = log2FoldChange, padj2 = padj),
  by = "Geneid"
)
```
- **The Job**: Merges the results from two different cell lines (e.g., H9 vs SUPM2) into a single table, aligning them by Gene ID.
- **The Reasoning**: This creates a "pairwise comparison." We can now look at `Gene X` and see exactly what Kromastat did to it in the healthy line (`lfc1`) versus the cancer line (`lfc2`).

---

## 3. Defining "Consensus" & Automatic Labeling
```r
consensus_top <- merged %>% 
  filter(status == "Consensus (Sig in both)") %>%
  mutate(abs_lfc = (abs(lfc1) + abs(lfc2))/2) %>%
  slice_max(order_by = abs_lfc, n = 5)

specific_top <- merged %>%
  filter(status == "Specific (Sig in one)") %>%
  mutate(abs_lfc = pmax(abs(lfc1), abs(lfc2))) %>%
  slice_max(order_by = abs_lfc, n = 5)
```
- **The Job**: Identifies which genes are "Universal Hits" (Consensus) vs "Model-Specific Hits" and automatically selects the top 5 from each category for labeling.
- **Why it matters**: A "Red" dot on this plot (Consensus) is a high-value target. Labeling the top 5 consensus genes allows your mentor to instantly see the most robust targets (e.g., *HMGCS1* or *SQLE*) that respond to the drug regardless of the cell type.

---

## 4. The Pearson Correlation (R-value)
```r
corr_val <- cor(merged$lfc1, merged$lfc2, method = "pearson", use = "complete.obs")
...
subtitle = paste(cl1, "vs", cl2, "| R =", round(corr_val, 3))
```
- **The Job**: Calculates a mathematical score (from -1 to 1) for the similarity of the two datasets.
- **The Reasoning**: 
    - **R > 0.8**: Extremely high consistency. The drug works identically in both models.
    - **R < 0.3**: Low consistency. The drug behaves differently in the cancer line compared to the healthy line. 
    This number is the most rigorous proof of "Mechanism Conservation" for your mentor's final report.

---

## 5. LFC-LFC Plotting with Ggrepel
```r
ggplot(merged, aes(x = lfc1, y = lfc2, color = status)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_text_repel(data = label_df, aes(label = label), size = 3, max.overlaps = 20)
```
- **The Job**: 
    1. Plots the Fold Change of one cell line on the X-axis and the other on the Y-axis.
    2. Adds a "Trend Line" to visualize the relationship.
    3. Overlays the names of the most impactful genes using `geom_text_repel`.
- **The Reasoning**: This creates a "Diagonal Cloud." The inclusion of `ggrepel` labels transforms a simple scatter plot into a biological map, allowing researchers to identify specific gene targets at a glance.

---

## 6. Summary Output
- **Title**: Describes the "Mechanism Consistency."
- **Subtitle**: Includes a legend for the cell line mapping: `(I): Indolent | (A): Aggressive`.
- **Result**: `04_10_lfc_correlation.png` — a final, high-impact proof that your drugs have a reliable and predictable biological effect across different biological systems.

---

## 7. How to Interpret Your Correlation Plot
This plot is the ultimate test of **Biological Reproducibility**.

1. **The R-Value (Pearson Correlation)**:
    - **R > 0.7**: Strong evidence that the drug's mechanism is **Conserved**. The drug does the same thing to both cell lines.
    - **R < 0.4**: Evidence of **Cell-Line Specificity**. The drug affects the cancer line differently than the healthy line. 
2. **Consensus Genes (Red Dots)**: These are the "Universal Hits" that passed the **|LFC| > 2.0** and **padj < 0.05** thresholds in BOTH cell lines. The labeled genes are the **Highest-Confidence Targets**.
3. **The Diagonal Cloud**: 
    - Most points should cluster along the diagonal dashed line. This means that if a gene was upregulated in Model A, it was also upregulated in Model B.
4. **Outliers (Dots far from the line)**: These represent "Model-Specific" biology. For example, a gene that is inhibited in cancer but not in healthy cells would appear far away from the diagonal.
