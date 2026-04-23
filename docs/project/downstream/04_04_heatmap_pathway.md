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
for (contrast in names(enrichment_results_all)) {
  res <- enrichment_results_all[[contrast]]$gsea_hallmark
  df_res <- as.data.frame(res) %>% filter(p.adjust < 0.05)
  
  # Select Top Pathways (Activated and Suppressed)
  top_act <- head(df_res %>% filter(NES > 0) %>% arrange(desc(NES)), 3)
  top_sup <- head(df_res %>% filter(NES < 0) %>% arrange(NES), 2)
  sel_paths <- rbind(top_act, top_sup)
}
```
- **The Job**: Iterates through each drug contrast and selects the top 3 activated and top 2 suppressed Hallmark pathways.
- **The Reasoning**: This ensures a balanced view of the biological response. Instead of just showing the "Top 5" (which might all be upregulated), we explicitly pull both directions to see what the drug activates AND what it shuts down.

---

## 2. Extracting the "Leading Edge"
```r
genes <- unique(unlist(strsplit(sel_paths$core_enrichment[i], "/")))
genes <- intersect(genes, rownames(dds))
genes <- head(genes, 12) # Mentor logic: Top 12 genes per pathway
```
- **The Job**: GSEA identifies a specific subset of genes that account for the enrichment signal, known as the "Leading Edge." We cap this at the top 12 genes per pathway.
- **The Reasoning**: If a pathway has 200 genes, it would dominate the entire heatmap. Capping at 12 (ranked by their contribution to the enrichment score) preserves equal visual weight for all pathways and highlights the "Core Drivers" as per mentor guidelines.

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
colnames(mat) <- colnames(mat) %>%
  str_remove_all("human_|canine_") %>% 
  str_replace_all("DMSO", "D") %>%
  str_replace_all("Romidepsin|Romi", "R") %>%
  str_replace_all("Kromastat|Kroma", "K")
```
- **The Job**: Automatically shortens sample names for the plot.
- **The Reasoning**: In a large 24-sample study, long names like `SUPM2_DMSO_Rep1` cause labels to overlap and become unreadable. This "Thesis-Ready" cleanup ensures every sample is identifiable at a glance (e.g., `SUPM2_D_1`).

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
