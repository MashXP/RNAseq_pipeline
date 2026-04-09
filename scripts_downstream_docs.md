# Downstream Analysis Documentation (R)

This document provides a technical breakdown of the R scripts used for differential gene expression (DGE) and functional enrichment.

> [!NOTE]
> The R environment and all necessary Bioconductor packages are pre-configured in the `cancer_rnaseq` Mamba environment.
---

## 01_data_prep.R (in scripts_downstream/)

**Purpose:** Cleans the count matrix and integrates sample metadata.

**Key Actions:**
- **Metadata Cleaning:** Extracts sample details from `MiaPAca-2_Sample_Data_Table.csv`.
- **Gene Filtering:** Removes genes with low counts (threshold: 10 reads in >= 3 samples).
- **Output:** Saves a processed RData object for the next step.

---

## 02_deseq2_dge.R (in scripts_downstream/)

**Purpose:** Identifies differentially expressed genes (DEGs) across multiple doses.

**Key Actions:**
- **Normalization:** Estimates size factors and dispersion.
- **Model Fitting:** Fits a Negative Binomial GLM (~ condition).
- **LFC Shrinkage:** Applies `apeglm` to all treatment doses.
- **Contrasts:** Automatically loops through all concentrations (`0.5nM`, `1nM`, `5nM`) vs. `DMSO (NG)`.

---

## 03_enrichment.R (in scripts_downstream/)

**Purpose:** Performs functional annotation and Gene Set Enrichment Analysis (GSEA).

**Key Actions:**
- **ORA (Over-Representation Analysis):** Runs GO and KEGG for each dose using significantly up/down genes.
- **GSEA (Gene Set Enrichment Analysis):** Uses ranked gene lists (LFC) to identify activated/suppressed biological pathways across the full expression spectrum.
- **Output:** Generates CSV tables for ORA, GSEA (Full), and subsets for Activated vs. Suppressed gene sets.

---

## 04_visualization.R (in scripts_downstream/)

**Purpose:** Generates comprehensive visual reports for the multi-dose study.

**Key Actions:**
- **PCA Plot:** Visualizes sample clustering across all treatment levels.
- **Volcano Plots:** Generated for every dose vs. DMSO comparison.
- **Venn Diagram:** Shows overalpping DEGs between different drug concentrations.
- **Pathway Heatmap:** Extracts genes from the top-enriched pathways (identified in `03_enrichment.R`) and displays their expression patterns across all study groups.
- **Interactive Reports:** Outputs are optimized for inclusion in terminal summaries and report documents.
