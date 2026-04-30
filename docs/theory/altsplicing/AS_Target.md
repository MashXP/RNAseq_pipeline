# Alternative Splicing: Selection Criteria for Sashimi Plots

This document outlines the priority criteria for selecting alternative splicing (AS) events from rMATS-turbo results for downstream visualization and biological validation.

## 1. Effect Size: The "Big Change" Rule
The primary metric for biological impact is the **`IncLevelDifference`** (also known as **DeltaPSI** or **ΔΨ**).
- **Threshold**: Prioritize events with `|IncLevelDifference| > 0.2`.
- **Top Tier**: Events with `|IncLevelDifference| > 0.5` represent massive isoform switches (e.g., from 80% inclusion to 30% inclusion) and are high-priority candidates for publication.

## 2. Statistical Robustness: FDR
Ensure the event is statistically significant across replicates.
- **Threshold**: `FDR <= 0.05` is mandatory.
- **Note**: A very low FDR (e.g., < 0.001) combined with a high DeltaPSI indicates a highly consistent and strong biological signal.

## 3. Visual Clarity: The "Clean Loop" Rule
Sashimi plots rely on read junctions. For a plot to look "convincing" in a paper:
- **Read Depth**: Check the **IJC** (Inclusion Junction Counts) and **SJC** (Skipping Junction Counts).
- **Recommendation**: Target events where the total junction reads across replicates are **> 30-50**. Low-coverage events (e.g., < 10 reads) result in thin, jittery arcs that may be dismissed as noise.

## 4. Ease of Interpretation: AS Type
- **Skipped Exon (SE)**: High Priority. These are the most common and the easiest to explain (e.g., "Drug X causes exon 4 to be skipped, likely truncating the protein").
- **Mutually Exclusive Exons (MXE)**: Medium-High Priority. Very interesting for protein domain switching.
- **Alternative 3'/5' Splice Sites (A3SS/A5SS)**: Medium Priority. Often lead to subtle frame shifts or small amino acid changes.
- **Retained Intron (RI)**: Lower Priority. Often lead to nonsense-mediated decay (NMD) rather than a functional protein isoform.

## 5. Functional Context: "Double-Hit" Candidates
Cross-reference the **`geneSymbol`** with your Differential Gene Expression (DGE) results:
- **"Double-Hit"**: If a gene is both **Differentially Expressed (DEG)** and **Differentially Spliced (DAS)**, it suggests a multi-layered regulatory response to the treatment.
- **Gene Function**: Prioritize genes known to be involved in the treatment pathway (e.g., Chromatin remodeling, Cancer cell cycle).

---

## Top Verified Candidates (Project Snapshot)

Based on the criteria above, the following events have been identified as high-priority targets for initial Sashimi plot generation.

| Gene | Type | Cell Line | Comparison | ΔΨ | ID | Selection Rationale |
| :--- | :--- | :--- | :--- | :---: | :--- | :--- |
| **CTSS** | SE | UL1 | Romi vs DMSO | -0.85 | 94880 | **Top Canine Target**: Massive effect size in canine study. |
| **IL31RA** | SE | SUPM2 | Romi vs DMSO | -0.78 | 114255 | High magnitude switch in human lymphoma cells. |
| **ARL4A** | SE | SUPM2 | Romi vs DMSO | +0.77 | 88115 | Strong inclusion signal; solid read depth. |
| **ALDOA** | SE | SUPM2 | Romi vs DMSO | +0.74 | 159236 | **Best Visual Candidate**: Huge shift + >2000 reads! |
| **SEMA4D** | SE | H9 | Romi vs DMSO | +0.74 | 218658 | Highly consistent target across multiple comparisons. |
| **CTSB** | MXE | SUPM2 | Romi vs DMSO | +0.76 | 29982 | Best Mutually Exclusive Exon candidate. |

---
*Last Updated: 2026-04-29*
*Reference: rMATS-turbo Pipeline Summary*
