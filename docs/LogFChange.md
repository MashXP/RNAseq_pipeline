# Biological Analysis Memo: Effect-Size Thresholding

**Date**: April 15, 2026  
**Subject**: Elevation of log2FoldChange (LFC) Stringency to 2.0  
**Status**: Investigatory / Recommended

## Observation
Current analysis of Romidepsin and Kromastat treatment in SUPM2 and H9 lymphoma models using a standard $|log2FC| \ge 1.0$ threshold ($FDR < 0.05$) yields an extremely high number of differentially expressed genes (DEGs), reaching nearly 9,000 in some contrasts. This represent ~40-50% of the active transcriptome, which significantly dilutes biological signal during GSEA and functional enrichment.

## Rationale for Switching to 2.0
Following expert mentorship advice, we investigated shifting to a $|log2FC| \ge 2.0$ (4-fold change) threshold. 

1. **Functional Precision**: Prioritizes genes with massive regulatory shifts, likely to be the primary drivers of therapeutic efficacy (P53/MYC antagonism).
2. **Noise Suppression**: HDAC inhibitors (Romidepsin) induce widespread transcriptomic "ripples." A 2.0 threshold filters out secondary/tertiary effects to focus on primary targets.
3. **Statistical Integrity**: All counts are based on **apeglm-moderated** LFC values, which already suppress noise in low-count genes. A 2.0 threshold on moderated values is extremely robust.

## Impact Audit (|LFC| $\ge$ 2.0)

| Cell Line | Comparison | DEGs (1.0) | DEGs (2.0) | Reduction |
| :--- | :--- | :---: | :---: | :---: |
| **SUPM2** | Romi vs DMSO | 8,863 | 4,497 | ~49% |
| **SUPM2** | Krom vs DMSO | 5,439 | 2,057 | ~62% |
| **H9** | Romi vs DMSO | 7,835 | 3,757 | ~52% |
| **H9** | Krom vs DMSO | 4,326 | 1,435 | ~67% |
| **UL1** | Romi vs DMSO | 4,803 | 2,168 | ~55% |
| **CNK89** | Krom vs DMSO | 593 | 82 | ~86% |

## Recommendation
Transition the [Main Analysis Report](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/docs/MAIN_ANALYSIS.md) and all visualization scripts (`04_02_volcano.R`, `04_04_heatmap_pathway.R`) to the **2.0 threshold** for the final peer-reviewed version.
