# Contrast Logic and Subset Selection

In complex transcriptomic studies (e.g., 2 Drugs × 2 Cell Lines), a single comparison is insufficient. We use a multi-contrast approach to "triangulate" biological truth and identify high-confidence gene subsets.

## 1. Why Multiple Contrasts?
The primary goal of running multiple contrasts is **Comprehensive Subset Selection**. By looking at the data through different "lenses," we can filter out technical noise and focus on robust biological signals.

### A. The "Efficacy" Lens (Treatment vs. DMSO)
*   **Goal**: Identify genes that respond to the drug.
*   **Subset**: The "Drug Responsive" transcriptome.
*   **Example**: `Romi_6nM_vs_DMSO_Romi`.

### B. The "Specificity" Lens (Cell Line Subsets)
*   **Goal**: Determine if a drug's effect is universal or dependent on the biological background (Healthy vs. Cancer).
*   **Subset**: "Context-Dependent" targets.
*   **Example**: `H9_Romi vs DMSO` vs. `SUPM2_Romi vs DMSO`.

### C. The "Shared Core" Lens (Global/Pooled Contrasts)
*   **Goal**: Identify the fundamental mechanism of action that transcends cell line identity.
*   **Subset**: The "Conserved Core Mechanism."
*   **Example**: The intersection of Romidepsin and Kromastat targets.

---

## 2. Anatomy of a Heatmap Contrast
When we generate a heatmap for a specific contrast (e.g., `DMSO_Kromastat_vs_DMSO_Romi`), we are combining two distinct selection processes:

1.  **Gene Selection (The Filter)**: The genes displayed are the **Leading Edge** genes derived from that specific contrast's GSEA/DGE results.
2.  **Expression Visualization (The Context)**: The colors represent the **Z-scores across ALL samples** in the study.

### Why do we show all samples?
Showing all samples in a heatmap—even those not involved in the original contrast—allows for **Cross-Validation**. 
*   If a gene is selected because it varies between two DMSO controls, but it shows no response to the Drugs in the other columns, we can flag it as a **Technical Artifact**.
*   If a gene varies between controls AND responds to drugs, it suggests a **Baseline Interaction** that warrants deeper investigation.

---

## 3. Summary: The Pipeline as a Sieve
Think of the multiple contrasts as a series of sieves:
*   **Sieve 1 (DGE)**: Removes low-count/non-significant noise.
*   **Sieve 2 (Intersections/Venn)**: Removes genes that aren't consistently regulated across replicates or cell lines.
*   **Sieve 3 (Functional Enrichment)**: Removes genes that don't belong to a recognized biological pathway.

The final genes that survive all these "sieves" (the top shared targets in your report) are the most biologically certain and therapeutically relevant candidates.
