# Ortholog NES Overlap (04_13)

## Overview
The `04_13_ortholog_nes_overlap.R` script provides a high-level functional comparison between species. While `04_12` looks at individual genes, this script looks at the **Pathway Activity (NES)** to determine if the same biological processes (e.g., MTORC1 Signaling, MYC Targets) are consistently regulated across Human and Canine models.

---

## Technical Specifications
- **Script**: [04_13_ortholog_nes_overlap.R](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/scripts_downstream/04_13_ortholog_nes_overlap.R)
- **Input**: 
    - **Human Hallmark GSEA**: `../results/human/tables/03_gsea_hallmark_*.csv`
    - **Canine Hallmark GSEA**: `../results/canine/tables/03_gsea_hallmark_*.csv`
- **Output**: 
    - **Comparison Bar Plots**: `../results/comparative/figures/04_13_nes_bar_[ID].png`
    - **Overlap Data Tables**: `../results/comparative/tables/04_13_nes_overlap_[ID].csv`

---

## Analysis Logic

### 1. Species-Agnostic GSEA Alignment
```r
merged <- inner_join(h_df, c_df, by = "ID") %>%
  mutate(ID = str_replace(ID, "HALLMARK_", "")) %>%
  mutate(sig_status = case_when(
    h_padj < 0.05 & c_padj < 0.05 ~ "Both Significant",
    h_padj < 0.05 ~ "Human Only",
    c_padj < 0.05 ~ "Canine Only",
    TRUE ~ "NS"
  ))
```
- **The Job**: Merges the Hallmark GSEA results from both species and classifies pathways based on their cross-species statistical significance.
- **The Reasoning**: This identifies "Universal Pathways." If MTORC1 Signaling has a high absolute NES and is significant in both species, it is a conserved drug-responsive module.

### 2. Side-by-Side NES Bar Plots
Visualizes the Normalized Enrichment Score (NES) for each pathway side-by-side for Human and Canine.
- **Ordering**: Pathways are ordered by the Human NES for a clean, descending visual "story."
- **Coloring**: Static species colors (Human = Pink, Canine = Green) are used to clearly separate the data points.
- **Why it matters**: This plot is the "Final Proof" of functional conservation. It demonstrates that the drug's impact on biological themes is reproducible across species, even if the underlying genes might vary slightly.

---

## Execution
Run from the project root:
```bash
./run down
```
