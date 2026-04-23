# Unified Alluvial Analysis (04_11)

## Overview
The `04_11_alluvial_plot.R` script is the primary visualization engine for cross-species functional conservation. It uses `ggalluvial` to track the "flow" of Hallmark pathway significance between Human and Canine transcriptomes, enabling both fine-grained comparison of specific cell-line pairs and high-level global views of the entire study.

---

## Technical Specifications
- **Script**: [04_11_alluvial_plot.R](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/scripts_downstream/04_11_alluvial_plot.R)
- **Input**: 
    - **Human GSEA Results**: `../results/human/tables/03_gsea_hallmark_*.csv`
    - **Canine GSEA Results**: `../results/canine/tables/03_gsea_hallmark_*.csv`
- **Output**: 
    - **Figures**: `../results/comparative/figures/04_11_alluvial_[Name].png`
    - **Tables**: `../results/comparative/tables/04_11_alluvial_gene_flow_[summary/detailed].csv`

---

## Analysis Modes

### 1. Individual Comparative Plots
Generates 2-way alluvials between corresponding cell lines (e.g., Human H9 vs. Canine UL1).
- **Naming**: `04_11_alluvial_[Drug]_[Human]_vs_[Canine].png`

### 2. Global Panoramic Flows
Generates multi-axis flows (4 axes) to connect the entire experimental design.
- **Phenotypic Chains**: Visualizes how the drug response evolves across the "Indolent-to-Aggressive" spectrum.
- **Drug-Species Bridge**: `Human:Romi -> Human:Kroma -> Canine:Romi -> Canine:Kroma`. Proves mechanistic conservation.
- **Sub-grouping**: Includes specific flows for **Aggressive** and **Indolent** cell line pairs.
- **Naming**: `04_11_alluvial_hallmark_[Drug].png` or `04_11_alluvial_[aggressive/indolent].png`

### 3. Gene-Level Flow Analysis
Tracks individual gene response classes (shared_up, specific_to_human, etc.) across matched pairs.
- **Naming**: `04_11_alluvial_gene_flow.png`

---

## Visual Design Principles
- **Color Palette**: Synced with NES barplots (`Activated` = Red, `Suppressed` = Blue, `NS` = Grey).
- **Readability**: 
    - Thin white outlines around flows to prevent "muddy" overlaps.
    - Semi-transparent strata boxes (`alpha = 0.3`) for a premium look.
    - Centered, bolded labels for categorical clarity.

---

## Execution
Run from the project root:
```bash
./run down
```
