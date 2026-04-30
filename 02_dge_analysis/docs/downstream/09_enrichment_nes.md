# Downstream Dissection: 04_06_enrichment_nes.R

This script creates the "Executive Summary" of your functional analysis. It takes the complex GSEA results and turns them into a simple, high-impact Barplot showing which biological themes were **Activated** or **Suppressed** by the drug.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Enrichment RData**: `./.RData/[Group]/03_enrichment_results.RData`.
- **Processing**: Selection of top 10 activated and top 10 suppressed pathways based on absolute NES.
- **Output**: 
    - **NES Figures**: `../results/[Group]/figures/04_06_hallmark_nes_[Contrast].png` and `04_06_hallmark_nes_combined.png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for visualization. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggplot2`: For faceted barplotting.
- `patchwork`: For horizontal stitching of comparison plots.
- `tidyverse`: Data filtering and factor leveling.

## 1. Contrast Iteration & Null Skipping
```r
for (contrast in names(enrichment_results_all)) {
  hallmark_res <- enrichment_results_all[[contrast]]$gsea_hallmark
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    p <- make_nes_plot(hallmark_res, contrast)
  }
}
```
- **The Job**: Loops through every drug contrast (e.g., Kromastat vs DMSO) and extracts the Hallmark GSEA results.
- **The Reasoning**: This includes a critical "Null Check". If a specific treatment fails to yield any significantly enriched pathways, the script will silently skip it rather than crashing. This robustness ensures the pipeline completes successfully even with sparse data.

---

## 2. Defining Biological Response
```r
hall_df <- as.data.frame(hallmark_res) %>%
  mutate(Response = ifelse(NES > 0, "Activated", "Suppressed"))
```
- **The Job**: Categorizes every pathway based on its Normalized Enrichment Score (NES).
- **The Reasoning**: Positive NES means the genes in that pathway were generally "turned up" (Activated) by Kromastat or Romidepsin. Negative NES means they were "turned down" (Suppressed). This is the key metric your mentor will look for in the results.

---

## 3. The "Top 10" Heuristic
```r
hall_df <- hall_df %>%
  group_by(Response) %>%
  slice_max(order_by = abs(NES), n = 10) %>%
  ungroup()
```
- **The Job**: Selects only the 10 most extreme (highest/lowest) pathways in each direction.
- **The Reasoning**: If GSEA finds 50 significant pathways, the plot becomes unreadable. By focusing on the "Top 10," we ensure that only the most impactful biological changes (like Cell Cycle arrest or Inflammation) are presented in your final figure.

---

## 4. Logical Ordering
```r
ggplot(hall_df, aes(x = NES, y = reorder(Description, NES), fill = Response)) +
  geom_col(width = 0.8)
```
- **The Job**: Plots the pathways and automatically sorts them by their NES score (`reorder(Description, NES)`).
- **Why it matters**: This creates a clean gradient effect. The Activated pathways naturally group at the top, and the Suppressed pathways group at the bottom based on their mathematical scores. This visual organization makes the "story" of the drug immediately obvious.

---

## 5. Industry-Standard Coloring
```r
scale_fill_manual(
  values = c("Activated" = "#D62728", "Suppressed" = "#1F77B4"),
  name = "Direction"
)
```
- **The Job**: Applies your project's signature Red/Blue color scheme to the bars.
- **The Reasoning**: By reusing the exact same colors from the Volcano plots, you build a cohesive visual identity for your results. Red always means "Increase," Blue always means "Decrease."

---

## 6. Robust Conditional Stitching
```r
if (length(nes_plots_curated) > 1) {
  combined_nes <- wrap_plots(nes_plots_curated, nrow = 1)
} else {
  message("Skipping Combined NES plot: Fewer than 2 curated contrasts with NES data.")
}
```
- **The Job**: Safely attempts to combine curated drug-specific barplots into one figure.
- **The Reasoning**: The script uses a `COMBINED_EXCLUDE` list to remove non-informative QC contrasts (e.g., DMSO vs DMSO) from the final faceted comparison. This ensures the "Combined" plot focuses only on the high-value treatment comparisons that add biological insight.

---

## 7. Fixed Canvas Sizing
```r
ggsave(..., width = 8, height = 7, limitsize = FALSE)
```
- **The Job**: Applies a standardized aspect ratio to the individual NES plots.
- **The Reasoning**: Because we specifically slice the "Top 10" activated and suppressed pathways, the maximum number of bars on the plot is predictable (up to 20). Using a fixed height of `7` ensures uniform figure sizes across all comparisons for your final report without the need for dynamic stretching.

---

## 8. How to Interpret Your NES Barplot
The NES barplot is your high-level functional dashboard. It shows you the direction and strength of biological changes.

1. **Normalized Enrichment Score (NES)**: 
    - **Positive NES (Red Bars)**: These pathways are **Activated**. The genes in this set are predominantly upregulated by the drug.
    - **Negative NES (Blue Bars)**: These pathways are **Suppressed**. The genes in this set are predominantly downregulated.
    - **Magnitude**: A higher absolute NES (e.g., 3.0 vs 1.5) indicates a much stronger and more coordinated shift in that biological process.
2. **Biological Significance**: 
    - If you see `E2F_TARGETS` and `MYC_TARGETS` in the **Suppressed** (Blue) section, it strongly indicates that the drug is stopping cell proliferation.
    - If you see `P53_PATHWAY` or `APOPTOSIS` in the **Activated** (Red) section, it suggests the drug is triggering cancer cell death.
3. **FDR (False Discovery Rate)**: All pathways shown have a p-adjust < 0.05. This means you can be 95% confident that these functional shifts are not due to random chance.
