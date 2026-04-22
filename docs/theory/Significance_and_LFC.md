# Theory: Statistical Significance and Effect Size in RNA-Seq

This document outlines the relationship between Adjusted P-values (padj) and Log2 Fold Changes (LFC) in the context of differential gene expression analysis.

## 1. Adjusted P-value (padj)
The `padj` (or False Discovery Rate, FDR) represents the probability that a gene is identified as differentially expressed purely by chance (a False Positive).

*   **Standard Threshold**: `padj < 0.05`. This means there is a 95% confidence that the observed change is not due to random noise.
*   **Stringency**: In large-scale drug studies, we often observe extremely low p-values (e.g., `1e-100` or `0.0`). These values indicate that the change is highly reproducible across all biological replicates.
*   **Max padj (Intersection Analysis)**: When looking for genes shared by two different drugs (e.g., Romidepsin AND Kromastat), we use the **Maximum padj** of the two. This represents the "Weakest Link" significance, ensuring the gene is reliable in both experiments.

## 2. Log2 Fold Change (LFC)
The `LFC` represents the magnitude of the change in gene expression.

*   **Definition**: A Log2 conversion of the ratio of expression between treated and control samples.
    *   `LFC = 1.0`: 2-fold increase (Up-regulated).
    *   `LFC = -1.0`: 2-fold decrease (Down-regulated).
    *   `LFC = 2.0`: 4-fold increase.
*   **Our Threshold**: We use `|LFC| > 2.0` (4-fold change) to filter for biologically impactful shifts, moving beyond simple statistical significance.

## 3. The Relationship: Significance vs. Magnitude
While often correlated, `padj` and `LFC` measure different things:

| Scenario | LFC (Magnitude) | padj (Significance) | Interpretation |
| :--- | :--- | :--- | :--- |
| **High LFC, High padj** | Large change | Not significant | Likely a outlier; inconsistent between replicates. |
| **Low LFC, Low padj** | Small change | Highly significant | Very consistent change, but perhaps not biologically impactful. |
| **High LFC, Low padj** | Large change | Highly significant | **Golden Target**: Both impactful and robustly consistent. |

### The "Volcano Effect"
In a Volcano Plot, we visualize this relationship. Genes "erupting" toward the top corners of the plot have both high magnitude (x-axis) and high significance (y-axis).
