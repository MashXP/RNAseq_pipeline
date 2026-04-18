# Pipeline Design Logic: Full Model vs. Subset Model

This document captures the strategic rationale for using different statistical models within the RNA-seq pipeline.

---

## 1. The Full Model (The "Big Picture")
*   **What it does**: It lets the computer look at all 24 samples (H9 + SUPM2) at once to calculate the "Noise" (Dispersion).
*   **Why use it?**: This is your **Discovery Tool**. Because the math has more data points, it is much more confident.
*   **Example**: **RCC2** was "found" by the Full Model. In a subset, its signal was too quiet, but in the Full Model, the math realized, "Wait, this gene is moving perfectly across EVERY sample."
*   **When to use it**: Use this for your **Main Conclusion**. It proves that the drug works on a "Human Cancer" level, not just one specific petri dish.

## 2. The Subset Model (The "Validation Tool")
*   **What it does**: It isolates just one cell line (e.g., H9 only, 12 samples).
*   **Why use it?**: This is your **Rigor Check**. Even if a drug works on "Human Cancer" in general, you need to prove it doesn't fail in a specific cell line.
*   **Example**: If you only had the Full Model, you might miss the fact that "Drug A" is 10x stronger in H9 than in SUPM2. Subsets reveal these differences.
*   **When to use it**: Use this for **Consistency Analysis** (like your UpSet and Correlation plots). It proves your results aren't just an "average" that is being driven by only one cell line.
