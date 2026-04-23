# Unified Alluvial Analysis (04_11)

## Overview
The `04_11_alluvial_plot.R` script is the primary visualization engine for cross-species functional conservation. It uses `ggalluvial` to track the "flow" of Hallmark pathway significance between Human and Dog transcriptomes, enabling both fine-grained comparison of specific cell-line pairs and high-level panoramic views of the entire study.

---

## Technical Specifications
- **Script**: [04_11_alluvial_plot.R](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/scripts_downstream/04_11_alluvial_plot.R)
- **Input**: 
    - **Human GSEA Results**: `../results/human/tables/03_gsea_hallmark_*.csv`
    - **Dog GSEA Results**: `../results/dog/tables/03_gsea_hallmark_*.csv`
- **Output**: 
    - **Figures**: `../results/comparative/figures/04_11_alluvial_[comp/master]_[Name].png`
    - **Tables**: `../results/comparative/tables/04_11_alluvial_[comp/master]_[Name].csv`

---

## Analysis Modes

### 1. Individual Comparative Plots
Generates 2-way alluvials between corresponding cell lines (e.g., Human H9 vs. Dog UL1).
- **Metric**: Displays **Global Conservation** (all pathways) and **Signal Conservation** (significant in human).
- **Labels**: Strata boxes show the status name and the raw count `(n)`.
- **Naming**: `04_11_alluvial_comp_[Drug]_[Human]_vs_[Dog].png`

### 2. Master Panoramic Flows (Huge Alluvials)
Generates multi-axis flows (4+ axes) to connect the entire experimental design.
- **Phenotypic Chains**: `H9 -> SUPM2 -> UL1 -> CNK89`. Visualizes how the drug response evolves across the "Indolent-to-Aggressive" spectrum in both species.
- **Drug-Species Bridge**: `Human:Romi -> Human:Kroma -> Dog:Romi -> Dog:Kroma`. Proves mechanisitic conservation between the two therapeutic agents across the species barrier.
- **Naming**: `04_11_alluvial_master_[romidepsin/kromastat/bridge].png`

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
mamba run -n cancer_rnaseq Rscript scripts_downstream/04_11_alluvial_plot.R
```
