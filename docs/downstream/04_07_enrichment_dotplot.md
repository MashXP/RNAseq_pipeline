# Downstream Dissection: 04_07_enrichment_dotplot.R

The Dotplot is the best way to compare your two drugs (Kromastat vs Romidepsin) side-by-side. It shows you not just *if* a pathway is significant, but how many genes it involves and how strong the statistical evidence is.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Enrichment RData**: `scripts_downstream/.RData/[Group]/03_enrichment_results.RData` (GSEA Hallmark results).
- **Processing**: Selection of Top 12 pathways per contrast, title-casing for labels, mirror dot-encoding (X-axis=direction, color=Response, size=Leading Edge, alpha=padj), and faceted stitching.
- **Output**: 
    - **GSEA Figures**: `../results/[Group]/figures/04_gsea_dotplot_[Contrast].png` and `04_gsea_dotplot_combined.png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for comparison plotting. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggplot2`: For dot-encoding and faceting (`facet_wrap`).
- `tidyverse`: For pathway selection and factor level sorting.
- `tools`: For programmatic `toTitleCase` string cleaning.

---

## 1. Contrast Mapping (`map2_dfr`)
```r
gsea_combined <- map2_dfr(
  enrichment_results_all, names(enrichment_results_all),
  function(res, contrast_name) {
    gsea_res <- res$gsea_hallmark
    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) { ... }
  }
)
```
- **The Job**: Iterates through all drug contrasts and securely compiles their GSEA results into one giant table.
- **The Reasoning**: `map2_dfr` is a robust `tidyverse` function that automatically stitches lists of dataframes together while gracefully ignoring `NULL` (empty) results, preventing catastrophic pipeline failures if one drug dose didn't work.

---

## 2. Targeted Pathway Selection
```r
top_paths_per_contrast <- gsea_combined_curated %>%
  group_by(Contrast) %>%
  slice_min(order_by = p.adjust, n = 12) %>%
  ungroup() %>%
  pull(Description) %>%
  unique()
```
- **The Job**: Identifies the Top 12 pathways for *every* comparison based on the most rigorous statistical metric (`padj`).
- **The Reasoning**: If Romidepsin has 100 significant pathways but Kromastat only has 5, a standard global "Top 20" list would only show Romidepsin results. This logic ensures both drugs contribute fairly to the final comparison plot while keeping the Y-axis concise.

---

## 3. Clean Naming (`clean_pathway_name`)
```r
clean_pathway_name <- function(pathway_name) {
  x <- sub("^HALLMARK_", "", pathway_name)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}
```
- **The Job**: Removes the ugly `HALLMARK_` prefix and converts ALL_CAPS technical names into Title Case.
- **Why it matters**: The default database strings are aesthetically poor for publication. Converting them to readable Title Case ensures the plot looks professional in a journal article or a mentor's slide deck.

---

## 4. The "Mirror Dotplot" (Aesthetic Mapping)
```r
ggplot(dotplot_df, aes(x = PlotX, y = Description)) +
  geom_point(aes(size = Count, color = Response, alpha = neg_log10_padj))
```
- **The Job**: Encodes three different variables into every dot:
    - **X-Axis (`PlotX`)**: Places Activated pathways on the left (-1) and Suppressed on the right (1).
    - **Color (`Response`)**: Solid Red for Activated, Solid Blue for Suppressed.
    - **Alpha (`neg_log10_padj`)**: Transparency shows statistical confidence. Faded dots are less significant.
    - **Size (`Count`)**: Shows how many genes are in the Leading Edge.
- **The Reasoning**: This highly advanced "Mirror Dotplot" packs an incredible amount of information into one view, instantly separating the biological direction from the statistical strength.

---

## 5. Vertical Separators
```r
geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5)
```
- **The Job**: Draws a clear line down the center between the "Activated" and "Suppressed" columns.
- **The Reasoning**: This simple visual "fence" makes it much easier for a reader to navigate the grid and instantly see which pathways were turned on vs. turned off by the drug.

---

## 6. Faceted Multi-Dose Stitching
```r
combined_p <- ggplot(...) +
  facet_wrap(~ Contrast, nrow = 1)
```
- **The Job**: Uses ggplot's native faceting to place the drug comparisons side-by-side using a single global Y-axis.
- **The Reasoning**: This is the "A/B Comparison." You can look across the rows to see if `Kromastat` is hitting the same biological targets (at the same strength) as `Romidepsin`. Faceting ensures perfect Y-axis alignment across all plots.

---

## 7. How to Interpret Your Enrichment Dotplot
Dotplots pack four different dimensions of data into one view.

1. **X-Axis (PlotX)**: Splits the dots into **Activated** (left) and **Suppressed** (right). 
2. **Dot Color (Response)**: Matches the X-Axis direction (Red = Activated, Blue = Suppressed).
3. **Dot Transparency (-log10 padj)**: Represents the statistical confidence. Faded dots are less significant; solid, opaque dots are extremely significant.
4. **Dot Size (Leading Edge Count)**: Represents how many genes from your dataset belong to that pathway's leading edge. Larger dots indicate a more global impact.
5. **Row-wise Comparison**: By looking across the "staircase" of dots, you can see if different drug doses affect the same pathways. If the dots stay in the same row across facets, the drug mechanism is consistent.
