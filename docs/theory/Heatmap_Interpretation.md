# Interpreting Heatmap Quality: The "Visual Sanity Check"

When viewing a Z-score heatmap, the spatial distribution of color provides immediate feedback on the technical quality of the transcriptomic experiment. Before interpreting the biology, you must verify the technical integrity.

## 1. Vertical Consistency (The Replicate Rule)
In a high-quality experiment, all replicates of the same condition (e.g., Rep 1, 2, and 3 of H9+Romi) should look nearly identical.
*   **Success**: The three columns form a solid, uniform vertical block of color.
*   **Failure (Patchiness)**: If Rep 1 is Red and Rep 2 is Blue for the same gene, it indicates high technical noise or a sample swap.

## 2. Homogeneity in QC Contrasts (The "Boring" Heatmap)
For control-vs-control heatmaps where no change is expected (e.g., `DMSO_Kromastat_vs_DMSO_Romi`):
*   **Success**: The heatmap is **homogenous and "boring."** No clear blocks of red or blue exist between the two groups. This proves that your baseline is stable and the two vehicle controls are equivalent.
*   **Failure (Divergence)**: A clear divide between the two controls means the baseline is shifting, which would invalidate the drug-effect analysis.

## 3. The "Staircase" or "Gradient" Warning
If the color intensity seems to fade or strengthen sequentially from Rep 1 to Rep 3 across many genes:
*   **Indication**: This often indicates a **technical gradient** (e.g., a processing delay or a position effect on the sequencing flow cell) rather than true biology.

## 4. Horizontal Banding (Pathway Robustness)
In pathway-grouped heatmaps, look for horizontal consistency across different genes in the same pathway.
*   **Success**: Most genes in the "MYC Targets" block shift to Blue together under treatment.
*   **Failure**: Only one or two genes shift, while the rest stay neutral. This suggests the "enrichment" may be driven by outliers rather than a true pathway-wide response.
