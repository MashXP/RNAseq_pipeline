# Downstream Dissection: 04_04_heatmap_pathway.R

This script generates the most detailed visual in the pipeline: the **Leading Edge Heatmap**. It doesn't just show that a pathway is "on"—it shows you exactly *which* genes are driving that signal.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData` (Counts and LFC results).
    - **Enrichment RData**: `./.RData/[Group]/03_enrichment_results.RData` (GSEA pathways).
    - **Species Org.DB**: `org.Hs.eg.db` or `org.Cf.eg.db`.
- **Processing**: VST transformation, leading-edge gene extraction from GSEA, Z-score standardization, value capping at +/- 1.
- **Output**: 
    - **Heatmap Figures**: `../results/[Group]/figures/04_04_heatmap_pathway_[Contrast].png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for advanced visualization. See [**libraries.md**](libraries.md) for full technical justifications.
- `ComplexHeatmap`: Essential for pathway-based row splitting.
- `circlize`: For precise Z-score color mapping.
- `DESeq2`: For VST transformation.
- `tidyverse`: Data cleaning and string manipulation.
- `clusterProfiler`: For `bitr` gene mapping and GSEA processing.

---

## 1. Pathway Prioritization & Robustness
```r
for (contrast in contrasts_to_plot) {
  if (is.null(hallmark_res) || nrow(as.data.frame(hallmark_res)) == 0) {
    message("  Skipping: No significant Hallmark pathways found.")
    next
  }
  top_pathways <- head(df_halo, 6)
}
```
- **The Job**: Iterates through each drug contrast and safely skips any that lack significant pathway data, then selects a maximum of 6 top pathways.
- **The Reasoning**: This prevents script crashes when a drug fails to enrich any pathways. Limiting to the Top 6 pathways ensures the final heatmap is not overwhelmingly tall and focuses only on the most dominant biological themes.

---

## 2. Extracting the "Leading Edge"
```r
genes_path <- str_split(top_pathways$core_enrichment[i], "/")[[1]]
genes_path <- intersect(genes_path, rownames(mat_full))
```
- **The Job**: GSEA identifies a specific subset of genes that account for the enrichment signal, known as the "Leading Edge." 
- **The Reasoning**: A pathway might contain 200 genes, but maybe only 15 are actually moving. This code extracts those 15 "drivers" so we don't clutter the heatmap with dead weight.

---

## 3. Gene Limitation Heuristics
```r
if (length(genes_path) > 12) {
  genes_path <- res_lfc[genes_path, ] %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(12) %>% pull(gene)
}
# Allow duplicate genes across pathways
rownames(mat_slice) <- paste0(genes_path, "__", i)
```
- **The Job**: If a pathway has more than 12 leading-edge genes, it keeps only the Top 12 (ranked by absolute fold change). It also ensures genes can appear in multiple pathways by appending a unique suffix.
- **The Reasoning**: If `MYC_TARGETS` has 150 leading-edge genes, it would dominate the entire heatmap and push all other pathways off the page. Capping at 12 preserves equal visual weight for all pathways. The suffix trick prevents R from crashing when the same driver gene (e.g., `CDKN1A`) appears in both `P53_PATHWAY` and `APOPTOSIS`.

---

## 4. Z-Score Standardization
```r
mat_pathway <- t(scale(t(mat_pathway)))
mat_pathway <- pmin(pmax(mat_pathway, -1), 1)
```
- **The Job**: Converts raw expression values into "standard deviations from the mean."
- **The Reasoning**: 
    - **Scale**: If Gene A has 10,000 counts and Gene B has 10 counts, you can't see them on the same map. Scaling (Z-score) puts them on the same -2 to +2 scale.
    - **Capping**: We cap the values at +/- 1 so that one "outlier" sample doesn't wash out the subtle differences in the other samples. This also intensifies the color contrast for publication.

---

## 5. Row Splitting (By Pathway)
```r
ht <- Heatmap(
  mat_pathway,
  row_split = pathway_factor,
  row_title_rot = 0,
  cluster_rows = FALSE
)
```
- **The Job**: Groups the rows of the heatmap into distinct "boxes" based on their biological pathway (e.g., E2F Targets, DNA Repair).
- **Why it matters**: This is a premium visualization technique. Instead of one giant wall of colors, it creates an organized report where you can see exactly which genes belong to which functional group.

---

## 6. Metadata Annotations (The Headers)
```r
col_meta$condition <- factor(col_meta$condition, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM"))
col_order <- order(col_meta$cell_line, col_meta$condition)
...
column_split = data.frame(cell_lines, drug_groups)
```
- **The Job**: 
    1. Enforces a strict order: Romidepsin (Control + Treated) first, then Kromastat.
    2. Splits the columns visually by both Cell Line and Drug Group.
- **The Reasoning**: In a large 24-sample study, you need to be able to tell at a glance which columns are `Kromastat_6nM` and which are `DMSO`. This header provides the necessary context and ensures the "Conserved" pattern is easy to spot across cell lines.

---

## 7. Dynamic Sizing
```r
png(..., width = 1200, height = 200 + nrow(mat_pathway)*25, res = 150)
```
- **The Job**: Automatically calculates how "tall" the image should be based on how many genes were found, and places a packed legend list on the LEFT.
- **The Reasoning**: This prevents the "Squashed Heatmap" problem. By moving the legend to the left and calculating height dynamically, we ensure every row label is readable and the overall figure follows professional multi-panel standards.

---

## 8. How to Interpret Your Pathway Heatmap
Heatmaps provide a granular look at the expression signatures within specific biological processes.

1. **The Color Scale (Z-score)**: 
    - **Red (Positive)**: Gene expression is **above average** for that row.
    - **Blue (Negative)**: Gene expression is **below average** for that row.
    - **White (Zero)**: Gene expression is at the mean level.
    - *Note: We use Z-scores to compare genes with different raw magnitudes on the same scale.*
2. **Column Grouping**: Look for clear blocks of color. In a successful drug study, you should see a solid block of Blue in DMSO samples turning into a solid block of Red in Treatment samples (or vice versa).
3. **Row Splitting**: Genes are grouped by their **Biological Pathway**. If one "box" (e.g., DNA Repair) is entirely Red while another box (e.g., Cell Cycle) is entirely Blue, you can conclude that the drug is differentially regulating these distinct systems.
4. **Leading Edge Drivers**: The genes listed on the right are the "drivers" identified by GSEA. If you see specific symbols like *CDKN1A* or *CCNE2* in these blocks, it confirms that the drug is hitting those specific regulatory checkpoints.
