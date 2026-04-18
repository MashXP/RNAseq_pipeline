# Intuition: ORA vs. GSEA

This document explains the conceptual difference between Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) as implemented in this pipeline.

---

## 1. ORA (Over-Representation Analysis)
**The "VIP List" Approach**

*   **Mechanism:** You define a "Significant List" using a hard cutoff (e.g., $p_{adj} < 0.05$ and $|log2FC| > 1$). You then compare this list against a background (all genes on the chip/genome) to see if specific pathways appear more often than random chance.
*   **Intuition:** Imagine a bag of 10,000 marbles (genes). 100 are "red" (belong to Pathway A). You grab 50 marbles blindly (your significant genes). If you find 20 red marbles in your hand, that pathway is "Over-Represented."
*   **Pros:** Fast, easy to understand, focuses on the "big winners."
*   **Cons:** Highly dependent on the arbitrary cutoff. Ignores genes that barely missed the cutoff (e.g., $p = 0.051$).

**What to look for in the graph:**
*   **Gene Ratio:** How much of the pathway you "captured."
*   **Count:** Absolute number of genes from your list in that pathway.

---

## 2. GSEA (Gene Set Enrichment Analysis)
**The "Group Trend" Approach**

*   **Mechanism:** No cutoffs used. You rank **every single gene** detected in the experiment from "Most Upregulated" to "Most Downregulated." You then walk down this ranked list and see where the members of a specific pathway cluster.
*   **Intuition:** Imagine a marathon race where every gene is a runner. If all runners from "Team Mitochondria" finish in the Top 10%, the team is "Enriched" at the top (Activated). If they all lag at the back, they are "Suppressed."
*   **Pros:** Detects subtle but consistent shifts. Even if no single gene is "significant," a pathway can be enriched if all its members move 10% in the same direction.
*   **Cons:** Computationally intensive; harder to interpret if the ranking metric is noisy.

**What to look for in the graph:**
*   **NES (Normalized Enrichment Score):** Direction and strength. Positive = Activated; Negative = Suppressed.
*   **Leading Edge:** The subset of genes that actually contribute to the enrichment score.

---

## 3. Comparison Summary

| Feature | ORA | GSEA |
| :--- | :--- | :--- |
| **Data Used** | Only "Top" genes (Subset) | All genes (Ranked list) |
| **Analogy** | A VIP Guest List | A Team Ranking in a Race |
| **Main Metric** | Gene Ratio / p-value | NES (Enrichment Score) |
| **Best For** | Finding obvious biological hits | Finding subtle, coordinated shifts |
| **Threshold** | Sensitive to $p$ cutoff | Threshold-free |

---

## 4. Application in this Project
*   **Use ORA** to identify the most dramatic processes being turned ON or OFF by Romidepsin.
*   **Use GSEA** to see if broader biological "hallmarks" (like the entire G2M Cell Cycle checkpoint) are shifting, even if some individual genes have high variance.
