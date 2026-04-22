# Downstream Dissection: 04_09_upset_consistency.R

This script answers the ultimate "Rigor" question: **"Are these results consistent across different models?"** It compares the drug effects in your healthy and aggressive cancer lines to find "Universal Targets."

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData` (Contains results_list with cell-line subsets).
- **Processing**: Extraction of significant gene IDs per cell-line/drug combination, UpSet intersection matrix calculation.
- **Output**: 
    - **Consistency Figure**: `../results/[Group]/figures/04_09_upset_consistency.png`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for intersection analysis. See [**libraries.md**](libraries.md) for full technical justifications.
- `UpSetR`: For matrix-based set intersections (better than Venn for 4+ sets).
- `tidyverse`: For extracting significant ID lists and string cleaning.
- `grid`: For adding custom title text to the UpSet plot.

## 1. Dynamic Cell-Line Detection & Robustness
```r
keys <- names(results_list)
known_cl <- c("H9", "SUPM2", "CNK89", "UL1")
cell_lines <- intersect(known_cl, unique(str_extract(keys, "^[A-Za-z0-9]+(?=_)")))

if (length(cell_lines) < 2) {
  message("\n[CRITICAL] Error: Could not find cell-line subset results...")
  stop("Missing subset results.")
}
```
- **The Job**: Uses regex (String Extraction) to automatically detect which cell lines were processed in your dataset by looking at the prefixes of the contrast names (e.g., extracting `H9` from `H9_Romi_vs_DMSO`).
- **The Reasoning**: This makes the script perfectly modular. Whether you are running the pipeline on just 2 cell lines or a massive 10-cell-line cohort, the script automatically adjusts. The `< 2` check prevents the script from crashing silently, alerting you immediately if you forgot to generate sub-comparisons in Step 02.

---

## 2. Multi-Set "UpSet" Logic
```r
library(UpSetR)
...
print(upset(fromList(all_sets), nsets = length(all_sets), order.by = "freq"))
```
- **The Job**: Generates an UpSet plot, which is a modern, high-power alternative to a Venn diagram.
- **The Reasoning**: Venn diagrams become a messy "flower" if you try to compare more than 3 things. Since you are comparing 2 drugs (Krom/Romi) across 2+ cell lines, we have 4+ sets. UpSet handles this by showing a "Matrix" of dots and bars, making it instantly clear how many genes are shared across all 4 groups.

---

## 3. Global Consistency Detection
```r
romi_case <- "Romi_6nM_vs_DMSO_Romi"
romi_sets <- list()
for (cl in cell_lines) {
  romi_sets[[paste0(cl, " Romi")]] <- get_sig_ids(...)
}
```
- **The Job**: Groups the significant gene hits from every cell line for a specific drug.
- **Why it matters**: If you find 50 genes that are changed by Kromastat in *all* cell lines (the "intersection"), you have found the drug's **Core Mechanism**. These 50 genes are much more important than genes that only changed in one specific line.

---

## 4. High-Resolution Legibility
```r
text.scale = c(1.5, 1.2, 1.2, 1.2, 1.2, 1.2),
point.size = 3.5, 
line.size = 1.0
```
- **The Job**: Increases the size of all labels and data points on the final plot.
- **The Reasoning**: Default UpSet plots often have tiny text that is impossible to read in a PowerPoint slide. This ensures the labels (like `H9 Romi`) and the numbers of shared genes are bold and clear for your mentor's final audit.

---

## 5. Directional Consistency
- **Input**: Loads the `results_list` from the DGE step (Step 02).
- **Threshold**: Uses a strict `padj < 0.05` and `|log2FC| > 2` to ensure we are only comparing high-confidence biological shifts.
- **Output**: `04_09_upset_consistency.png` — a comprehensive map of how well your experimental groups overlap.

---

## 6. How to Interpret Your UpSet Plot
The UpSet plot is used to identify **Consistently Regulated Genes** across multiple experimental conditions.

1. **Vertical Bars (Intersection Size)**: The height of the top bars represents the number of genes in a specific intersection.
2. **The Dot Matrix (Bottom)**: 
    - **A Single Dot**: Represents genes that are **unique** to that specific group (e.g., genes changed in SUPM2 but not H9).
    - **Connected Dots**: Represents genes that are **shared** between those specific groups.
3. **The "All-Points" Connection**: Look for the bar where all dots in a column are connected. This bar represents your **Universal Drug Signature**—genes that respond to the treatment regardless of the cell line or model.
4. **Set Size (Left Horizontal Bars)**: Shows the total number of significant genes for each individual group before looking at overlaps. This helps you see at a glance which cell line was more responsive to the drug.
