# Z-Score vs. Log Fold Change in Heatmap Visualization

When visualizing RNA-seq results, two common metrics are used for heatmap coloring: **Z-scores** and **Log Fold Change (LogFC)**. Choosing the right one depends on the biological question being asked.

## 1. Z-Score (Normalized Counts)
**Definition**: A Z-score measures how many standard deviations a sample’s expression value is from the mean of all samples in the experiment.

$$z = \frac{x - \mu}{\sigma}$$

### When to use:
*   **Replicate Verification**: When you need to show that biological replicates (e.g., Rep 1, 2, 3) are consistent.
*   **Baseline Visibility**: When you want to see the expression level in the control group (DMSO) relative to the treatment group.
*   **GSEA Leading Edge**: Standard for showing the expression profile of genes driving a pathway enrichment across the whole study.

### Pros:
*   Shows the "raw" behavior of the data across all samples.
*   Highlights outliers or technical artifacts in specific replicates.
*   Clear visual representation of "The Big Flip" (e.g., Control is all Blue, Treatment is all Red).

---

## 2. Log Fold Change (LogFC)
**Definition**: The ratio of expression between two conditions, usually expressed in $log_2$.

### When to use:
*   **Drug-to-Drug Comparisons**: When you want to compare the *magnitude* of the effect of 10 different drugs in a single summary view.
*   **Multi-Contrast Overviews**: When the goal is to show which genes are up-regulated across different comparisons.

### Pros:
*   Very concise; one column represents one entire contrast.
*   Directly visualizes the "Efficacy" or "Potency" of a treatment.

### Cons:
*   Hides replicate variance (it averages them together).
*   Loses the baseline context (you don't know if the gene was highly expressed or lowly expressed in the control, only that it changed).

---

## Summary Table

| Feature | Z-Score (VST) | Log Fold Change (LogFC) |
| :--- | :--- | :--- |
| **Input Data** | Per-sample counts | Comparison results |
| **Replicate View** | Visible (e.g., 3 columns per group) | Hidden (1 column per contrast) |
| **Goal** | Biological Robustness | Magnitude of Effect |
| **Best For** | Proving a pathway shift is real | Comparing drug potency |
