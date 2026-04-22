# Downstream Dissection: 04_08_ora_dotplot.R

While GSEA looks at the "Global" trend of all genes, the ORA Dotplot focuses strictly on your **Significant Hits** (the genes that passed your p-value and Fold Change cutoffs).

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Enrichment RData**: `./.RData/[Group]/03_enrichment_results.RData` (ORA GO results).
- **Processing**: GeneRatio parsing, top 15 pathway selection per contrast, dot-encoding (color=padj, size=Count).
- **Output**: 
    - **ORA Figures**: `../results/[Group]/figures/04_08_ora_dotplot_[Contrast].png` and `04_08_ora_dotplot_combined.png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for comparison plotting. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggplot2`: For dot-encoding and faceting (`facet_wrap`).
- `tidyverse`: For GeneRatio parsing and `str_wrap`.

---

## 1. Contrast Mapping & Ratio Calculation
```r
ora_combined <- map2_dfr(
  enrichment_results_all, names(enrichment_results_all),
  function(res, contrast_name) {
    ora_res <- res$ora_go
    if (!is.null(ora_res) && nrow(as.data.frame(ora_res)) > 0) {
      as.data.frame(ora_res) %>%
        mutate(
          Contrast = contrast_name,
          Ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        )
    }
  }
)
```
- **The Job**: Iterates through all drug contrasts and securely compiles their ORA results, while converting the biological "fraction" (e.g., 10 genes out of 100 in the pathway) into a decimal number (0.1).
- **The Reasoning**: This `Ratio` is the primary metric for ORA. It tells you what percentage of a biological pathway was "hit" by your drug. The script also employs a `COMBINED_EXCLUDE` list to filter out QC contrasts (e.g., DMSO vs DMSO) from the final faceted comparison, keeping the visualization clinically relevant.

---

## 2. Balanced Pathway Selection
```r
top_paths <- ora_combined %>%
  group_by(Contrast) %>%
  slice_min(order_by = p.adjust, n = 15) %>%
  ungroup() %>% pull(Description) %>% unique()
```
- **The Job**: Identifies the Top 15 pathways for *every* comparison.
- **The Reasoning**: Ensures both strong drugs (Romidepsin) and milder drugs (Kromastat) are represented fairly on the final comparison plot.

---

---

## 3. Faceted Multi-Dose Stitching
```r
combined_p <- ggplot(...) +
  facet_wrap(~ Contrast, nrow = 1)
```
- **The Job**: Uses ggplot's native faceting to place the drug comparisons side-by-side using a single global Y-axis.
- **The Reasoning**: This provides a side-by-side view to confirm if Kromastat achieves the same biological pathway saturation as Romidepsin. Faceting ensures perfect Y-axis alignment across all plots.

---

## 4. How to Interpret Your ORA Dotplot
While GSEA shows global trends, ORA focuses strictly on your **Significant Hits** (the genes that passed your p-value and Fold Change cutoffs).

1. **The X-Axis (Ratio)**: Represents the **Gene Ratio**.
    - If a pathway has 100 genes, and your Ratio is **0.2**, it means you "hit" 20 of them (20%).
    - A high Ratio is often more biologically meaningful than a low p-value, as it indicates a large chunk of that specific pathway was affected.
2. **Dot Size (Count)**: The total number of significant genes in that pathway.
3. **Dot Color (padj)**: The statistical significance. **Red** is better (more confident).
4. **Interpreting the Gaps**: If you see a pathway in Romidepsin but not in Kromastat, it means that specific biological process is **Unique** to Romidepsin's mechanism of action at that dose.
