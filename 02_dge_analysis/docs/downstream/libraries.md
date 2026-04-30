# Downstream R Libraries: Rationales

This document explains the technical and scientific reasoning behind the specific R packages used in the downstream analysis pipeline.

---

## 1. Core Data Science
### `tidyverse`
- **Rationale**: The `tidyverse` is a collection of packages (`dplyr`, `ggplot2`, `readr`, `stringr`) that share a common philosophy and syntax.
- **Why it matters**: It allows for "piped" workflows (`%>%`), making the code much easier for a human to read. It is the industry standard for cleaning and transforming messy biological metadata.

---

## 2. Differential Expression
### `DESeq2`
- **Rationale**: The gold standard for RNA-seq analysis. It uses a **Negative Binomial distribution** to model the discrete nature of gene counts.
- **Why it matters**: 
    - **Wald Test**: Provides high statistical confidence for small sample sizes.
    - **Shrinkage (APEGLM)**: Prevents "false positive" large fold changes in low-count genes. Without this, your Volcano plots would be dominated by noise.

---

## 3. Functional Biology (Enrichment)
### `clusterProfiler`
- **Rationale**: A universal interface for GO, KEGG, and GSEA analysis.
- **Why it matters**: It is highly reliable and integrates perfectly with `DESeq2` outputs. It handles the complex math of Gene Set Enrichment Analysis (GSEA) with high performance.

### `msigdbr`
- **Rationale**: Direct interface to the MSigDB (Molecular Signatures Database).
- **Why it matters**: Instead of manually downloading GMT files, this package pulls the latest **Hallmark** gene sets directly into R, ensuring your "Biological Story" is based on the most current curated data.

---

## 4. Visual Excellence (Heatmaps)
### `ComplexHeatmap`
- **Rationale**: Far more powerful than the standard `pheatmap` or `heatmap.2` packages.
- **Why it matters**: 
    - **Row Splitting**: Essential for our "Pathway Heatmaps" where we need to separate genes by biological process.
    - **Annotations**: Allows us to add multiple color-coded bars at the top (Condition, Cell Line) to provide context for every column.

### `circlize`
- **Rationale**: Provides sophisticated color-mapping functions.
- **Why it matters**: It allows us to define "Exact" Z-score colors (e.g., Navy to White to Firebrick), ensuring that our heatmaps are both beautiful and scientifically accurate.

---

## 5. Layout & Presentation
### `patchwork`
- **Rationale**: The simplest and most robust way to combine multiple `ggplot2` objects.
- **Why it matters**: It allows us to "Stitch" Romidepsin and Kromastat results side-by-side (e.g., in Volcano or NES plots) with perfect alignment.

### `ggrepel`
- **Rationale**: Automatically calculates non-overlapping positions for text labels.
- **Why it matters**: In a Volcano plot with 20,000 points, simple labeling is a mess. `ggrepel` ensures your top gene names (like your drug targets) are perfectly legible without blocking the data points.

---

## 6. Set Analysis
### `ggVennDiagram`
- **Rationale**: Provides aesthetically pleasing Venn diagrams using the `ggplot2` engine.
- **Why it matters**: It supports gradient fills based on gene counts, making the "Intensity" of overlap visually obvious at a glance.

---

## 7. Comparative Flow
### `ggalluvial`
- **Rationale**: Specifically designed to track categorical changes between groups.
- **Why it matters**: This is the engine for our **Alluvial Plots (04_11)**. It allows us to visualize how the statistical significance of a pathway "flows" from Human to Canine, proving evolutionary conservation of the drug's effect.

### `UpSetR`
- **Rationale**: A matrix-based visualization for set intersections.
- **Why it matters**: Venn diagrams fail when you have more than 3 groups. For your multi-cell line study, `UpSetR` is the only way to clearly show which genes are shared across 4 or more experimental sets.
