# Downstream Analysis Documentation (R)

This document provides a technical breakdown of the R scripts used for differential gene expression (DGE) and functional enrichment.

> [!NOTE]
> The R environment and all necessary Bioconductor packages are pre-configured in the `cancer_rnaseq` environment.
---

## 01_data_prep.R (in scripts_downstream/)

**Purpose:** Cleans the count matrix and integrates sample metadata.

**Key Actions:**
- **Metadata Cleaning:** Extracts sample details from the provided Sample Data Table (CSV).
- **Gene Filtering:** Removes genes with low counts (threshold: 10 reads in >= 3 samples).
- **Output:** Saves a processed RData object (`01_data_prep.RData`) for the next step.

---

## 02_deseq2_dge.R (in scripts_downstream/)

**Purpose:** Identifies differentially expressed genes (DEGs) across multiple groups.

**Key Actions:**
- **Normalization:** Estimates size factors and dispersion.
- **Model Fitting:** Fits a Negative Binomial GLM (~ condition).
- **LFC Shrinkage:** Applies `apeglm` to all treatment groups.
- **Contrasts:** Automatically loops through all treatment vs. control comparisons.

---

## 03_enrichment.R (in scripts_downstream/)

**Purpose:** Performs functional annotation and Gene Set Enrichment Analysis (GSEA).

**Key Actions:**
- **ORA (Over-Representation Analysis):** Runs GO and KEGG for each group using significantly up/down genes.
- **GSEA (Gene Set Enrichment Analysis):** Uses ranked gene lists (LFC) to identify biological pathways across the full expression spectrum.
- **Output:** Generates CSV tables and RData objects (`03_enrichment_results.RData`).

---

## 04_Visualization Suite (Modular scripts)

The visualization step is modularized into specialized scripts to allow for targeted execution and easier maintenance.

### 04_01_pca.R
- **Goal:** Principal Component Analysis.
- **Description:** Visualizes sample clustering and batch effects based on top variable genes.

### 04_02_volcano.R
- **Goal:** Volcano Plots.
- **Description:** Stitches per-dose Volcano plots showing log2 Fold Change vs. Adjusted P-value. Includes custom coloring for up/downregulated genes.

### 04_03_venn.R
- **Goal:** DEG Overlap Venn Diagram.
- **Description:** Calculates the intersection of significantly regulated genes across different drug concentrations.

### 04_04_heatmap_pathway.R
- **Goal:** Pathway-specific Gene Heatmap.
- **Description:** Extracts genes from the top-enriched GO pathways and displays Z-score normalized expression patterns across all samples using `ComplexHeatmap`.

### 04_05_heatmap_variable.R
- **Goal:** Top Variable Genes Heatmap.
- **Description:** Displays the top 50 most variable genes across the entire dataset to highlight overall expression variance.

### 04_06_enrichment_nes.R
- **Goal:** Hallmark NES Barplots.
- **Description:** Visualizes the Normalized Enrichment Score (NES) for MSigDB Hallmark gene sets, highlighting activated vs. suppressed pathways.

### 04_07_enrichment_dotplot.R
- **Goal:** GSEA Dotplots.
- **Description:** Generates dotplots summarizing GSEA GO enrichment results across all treatment doses for comparative analysis.

### 04_08_ora_dotplot.R
- **Goal:** ORA Dotplots.
- **Description:** Visualizes Over-Representation Analysis (ORA) results, showing the Ratio of significant genes found in specific GO pathways.
