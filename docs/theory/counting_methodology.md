# Bioinformatics Fundamentals: Methodology Memo

This document serves as a summary of core bioinformatics concepts encountered during the project, specifically regarding the statistical basis of RNA-seq analysis.

---

## 1. What is Dispersion?
In RNA-seq, "Dispersion" is the technical term for **Biological Noise**. 

When we sequence three replicates of the same treatment, they never have identical gene counts (e.g., 98, 105, 100). This "spread" is the dispersion. 
*   **High Dispersion**: The gene is naturally "flickering" and unstable.
*   **Low Dispersion**: The gene is precise and stable across replicates.

## 2. Why do we need a "Dispersion Model"?
We typically have many genes (20,000+) but very few replicates (n=3). We cannot trust the noise calculation of just 3 samples alone. DESeq2 solves this using **Information Sharing (Shrinkage)**:
1.  It looks at the noise levels of *all* genes together.
2.  It builds a "model" of how much noise is expected for a gene at a certain expression level.
3.  It "shrinks" the individual gene's noise estimate toward the global average.

## 3. The Power of the "Full Model" (Global Context)
This project compared two ways of building that noise model:
*   **Subset Model**: Only lets the math "see" a few samples (e.g., 12 samples). The noise estimates are less stable because the "average" is based on less data.
*   **Full Model (Our Approach)**: Lets the math see the whole experiment (e.g., 24-48 samples). 

> [!TIP]
> **Why it matters**: A gene that looks like a "Top Gene" in a subset might actually be revealed as "just naturally noisy" when you look at the full dataset. Conversely, a stable signal like **RCC2** is only prioritized when the model has enough data (the Full Model) to realize how uniquely stable it is across all conditions.

## 4. Why use LFC Shrinkage (Apeglm vs. Normal)?
In RNA-seq, genes with low counts are very "jittery." If a gene has 1 count in DMSO and 2 counts in Romi, that is technically a **2x fold change**, but it’s completely unreliable noise. **Shrinkage** identifies these noisy genes and pulls their fold-change toward zero.

*   **`normal` (Legacy)**: The original method used in your mentor's scripts. It treats all genes the same way. It’s reliable but can sometimes "dampen" the signal of truly important genes.
*   **`apeglm` (Modern)**: The method used in your pipeline. It uses a more sophisticated mathematical distribution (Cauchy) that is much better at:
    1.  Preserving large, real fold changes.
    2.  Aggressively silencing low-count noise.

## 5. Which factor matters most?
Our verification testing (**verify_mentor_discrepancy.R**) proved that while both factors contribute, the **Dispersion Pool (Full vs. Subset)** is the dominant driver of your results:

1.  **Shrinkage (Minimal Impact on Ranking)**: When we tested your mentor's "Subset" model using your modern `apeglm` shrinkage, we still got a **5/5 match** with their Top 5 list. This means changing the shrinkage method alone didn't change the top genes.
2.  **Model Pool (Major Impact on Ranking)**: When we switched from a "Subset" model (12 samples) to your "Full" model (24 samples), the Top 5 ranking **immediately shifted** (only a 2/5 match).

**Conclusion**: The decision to analyze all samples as one "Full Model" is the single most important factor that distinguishes your high-confidence results (including the rise of **RCC2**) from the legacy reports.

## 6. Venn Diagram Principles: Directionality & Precision
Comparing your drug treatments requires consistent logic. We have implemented two key principles:

1.  **Why we use Ensembl IDs for overlaps**: Different genes can sometimes share the same Symbol (e.g., pseucanineenes). By using Ensembl IDs for the Venn circle calculations, we ensure 100% mathematical precision. We only map back to Symbols for the final report to ensure legibility.
2.  **The Importance of Directionality (Up vs. Down)**: 
    *   An **"All DEGs"** Venn can be misleading. A gene might show up in the "Shared" overlap even if it's going **UP** in drug A and **DOWN** in drug B.
    *   By separating **"Up-regulated"** and **"Down-regulated"** Venns, we isolate the truly **Conserved Mechanisms**—the specific biological pathways that both drugs are activating or silencing in the same way.

## 7. Advanced Consistency Metrics: UpSet and LFC Correlation
As the complexity of the experiment grows (2 drugs across 2 cell lines), standard Venn diagrams become unreadable. We have implemented two premium tools to resolve this:

1.  **UpSet Analysis (Intersectional Matrix)**: 
    - Unlike Venn circles, UpSet plots use a matrix to show every combination of overlaps.
    - **How to Read**: The horizontal bars show total DEGs per group. The vertical bars show the size of specific intersections. The "dots" below show which groups participate in that intersection.
    - **The Story**: We use this to identify the **"Conserved Mechanism"**: the genes that respond to both drugs in both cell lines.
2.  **LFC Correlation (Biological Proof)**: 
    - Rather than a simple "In/Out" list, we plot the **magnitude** of change (Log2 Fold Change) in H9 vs. SUPM2.
    - **The Pearson Correlation ($R$)**: A high $R$ value (e.g., > 0.8) is definitive proof that the drug mechanism is stable regardless of the cell line's unique genetic background.

## 8. GSEA Statistical Rigor: The Wald Stat
To match the highest standards of bioinformatics storytelling, we utilize the **DESeq2 Wald Statistic (`stat`)** as the primary ranking metric for Gene Set Enrichment Analysis (GSEA).

### 8.1 Why not rank by Fold Change alone?
Ranking by Log2FC alone can be noisy because low-count genes often have inflated, unstable fold changes. 

### 8.2 The power of the `stat` column
The `stat` column equals the **Log2FC divided by its Standard Error**.
- **High Error (Noisy gene)** → Low `stat` → Lower priority in GSEA.
- **Low Error (Stable gene)** → High `stat` → Higher priority in GSEA.

This ensures that the enriched pathways are driven by your most **reliable** and **robust** genes, not just the loudest noise.

### 8.3 "Leading Edge" Interpretation
A pathway enrichment score is driven by a subset of genes called the **Leading Edge**.
- **Biological Drivers**: These are the "VIP" genes in a pathway (e.g., the specific Cyclins driving a "Cell Cycle" enrichment).
- **Visualization**: We use **Pathway-Grouped Heatmaps** to show the actual expression of these leading-edge genes, transforming a simple dotplot into a biological narrative.
