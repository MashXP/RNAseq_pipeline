# Downstream Dissection: 04_12_ortholog_gene_overlap.R

This script identifies **Direct Gene Conservation**. While GSEA (03) and Alluvials (11) show conservation of biological *themes*, this script proves that the *exact same genes* are being hit by the drug in both Human and Canine models.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Human DGE Tables**: `../results/human/tables/02_dge_*.csv`
    - **Canine DGE Tables**: `../results/canine/tables/02_dge_*.csv`
    - **Genome Annotations**: GTF files for both species (to verify gene symbols).
- **Processing**: Inner join on Gene Symbols, significance filtering (padj < 0.05, |LFC| > 2), and correlation mapping.
- **Output**: 
    - **Scatter Plots**: `../results/comparative/figures/04_12_ortholog_scatter_[Contrast].png`
    - **Panoramic Heatmap**: `../results/comparative/figures/04_12_conserved_ortholog_heatmap.png`
    - **Summary Table**: `../results/comparative/tables/04_12_ortholog_overlap_summary.csv`

---

## 1. Symbol-Based Ortholog Mapping
```r
ortholog_tbl <- inner_join(human_map, canine_map, by = "gene_name", suffix = c("_human", "_canine"))
```
- **The Job**: Creates a lookup table of genes that exist in both the Human and Canine genomes using the Gene Symbol (e.g., *TP53*, *MYC*) as the unique identifier.
- **The Reasoning**: While complex ortholog mapping tools exist, gene symbols remain the most biologically interpretable bridge. By matching symbols, we ensure that the "Conserved Targets" we find are actionable and well-understood in the literature.

---

## 2. Cross-Species LFC Correlation
```r
merged <- inner_join(human_res, canine_res, by = "gene_label") %>%
  mutate(
    status = case_when(
      h_padj < 0.05 & c_padj < 0.05 & h_lfc >= 2.0 & c_lfc >= 2.0 ~ "shared_up",
      h_padj < 0.05 & c_padj < 0.05 & h_lfc <= -2.0 & c_lfc <= -2.0 ~ "shared_down",
      TRUE ~ "other"
    )
  )
```
- **The Job**: Combines the Log2FoldChange data for Human and Canine for the same drug treatment and identifies genes that go **UP in both** or **DOWN in both**.
- **The Reasoning**: This is the "Validation Pillar." If a drug suppresses *CCNE2* in Human and *CCNE2* in Canine, we have high confidence that the drug's mechanism of action is evolutionarily conserved.

---

## 3. High-Resolution Scatter Plots
- **The Job**: Visualizes every matched gene on a grid (Human X-axis vs. Canine Y-axis).
- **Why it matters**: A high correlation (points clustering along a diagonal line) proves that the drug is behaving identically in both species. The points highlighted in Red (Shared Up) and Blue (Shared Down) represent the "Universal Mechanism" of your therapeutic agent.

---

## 4. Top-40 Conserved Heatmap
```r
heatmap_mat <- as.matrix(heatmap_data)
pheatmap(heatmap_mat, cluster_rows = TRUE, cluster_cols = FALSE, ...)
```
- **The Job**: Pulls the 40 genes with the highest combined effect size across both species and plots them in a single panoramic heatmap.
- **The Reasoning**: This provides a "High-Value Target List" for the publication. It shows the drug's strongest, most consistent effects across the species barrier in a single, high-impact image.

---

## 5. Summary Statistics Table
- **The Job**: Exports a final count of how many genes were Shared Up, Shared Down, or divergent.
- **The Reasoning**: This table provides the hard numbers for the "Results" section of the paper (e.g., "We identified 142 genes consistently downregulated by Romidepsin across both Human and Canine models").
