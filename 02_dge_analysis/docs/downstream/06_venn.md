# Downstream Dissection: 04_03_venn.R

The Venn diagram answers one of your most important project questions: **"Do Kromastat and Romidepsin change the same genes?"**

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData`.
- **Processing**: Extraction of significant IDs (padj < 0.05, abs(LFC) > 2) for All/Up/Down directions, intersection of sets.
- **Output**: 
    - **Overlap Figures**: `../results/[Group]/figures/04_03_venn_[Direction].png` and `04_03_venn_combined.png`.
    - **Shared Gene Lists**: `../results/[Group]/tables/venn_lists/shared_genes_[Direction].csv`.
    - **Master Shared Targets**: `../results/[Group]/tables/04_03_shared_targets_stats.csv`.
    - **Rigor Stats Table**: `../results/[Group]/tables/04_03_venn_rigor_stats.csv`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for set analysis. See [**libraries.md**](libraries.md) for full technical justifications.
- `ggVennDiagram`: Interactive and aesthetic Venn diagrams.
- `tidyverse`: Data filtering and list manipulation.
- `patchwork`: For vertical stitching of directional Venns.

---

## 1. Directional Filtering Logic
```r
directions <- list(
  All  = function(lfc, padj) padj < 0.05 & abs(lfc) > 2,
  Up   = function(lfc, padj) padj < 0.05 & lfc > 2,
  Down = function(lfc, padj) padj < 0.05 & lfc < -2
)
```
- **The Job**: This creates three separate Venn diagrams: 
    - **All**: Every gene shifted by the drugs.
    - **Up**: Genes activated by both.
    - **Down**: Genes suppressed by both.
- **The Reasoning**: Sometimes two drugs shift the same 100 genes, but one drug pushes them "Up" while the other pushes them "Down." By looking at directional Venns, we prove that the drugs have the **same biological mechanism**.

---

## 2. Preventing "Label Clipping"
```r
clean_names <- names(sig_list) %>% 
  str_replace_all("_vs_", "\nvs\n") %>%
  str_replace_all("_", " ")

scale_x_continuous(expand = expansion(mult = .4)) + 
coord_cartesian(clip = "off")
```
- **The Job**: Reformats the long sample names (e.g., `Kromastat_6nM_vs_DMSO_Kromastat`) to fit on the screen without being cut off.
- **The Reasoning**: This represents high-quality figure preparation for a thesis or publication. By expanding the axes (`mult = .4`) and turning off "clipping," we ensure that every word of the comparison title is visible.

---

## 3. Robust Multi-Directional Loop
```r
for (dir_name in names(directions)) {
  # Extract gene sets per direction
  sig_list <- lapply(results_list[priority_contrasts], function(x) {
    if (is.null(x$df)) return(character(0))
    res_df <- x$df
    pass_filter <- directions[[dir_name]](res_df$log2FoldChange, res_df$padj)
    res_df$Geneid[pass_filter]
  })
  
  # Filter out empty sets
  sig_list <- sig_list[sapply(sig_list, length) > 0]
  
  if (length(sig_list) < 2) {
    message("  -- Skipping ", dir_name, " Venn: Fewer than 2 priority contrasts have significant genes.")
    next
  }
}
```
- **The Job**: This loop automatically processes "All," "Up," and "Down" gene sets, while gracefully skipping any direction that doesn't have enough data to form a comparison.
- **The Reasoning**: RNA-seq can be unpredictable. If Romidepsin suppresses 500 genes but Kromastat suppresses 0, a Venn diagram is mathematically impossible. This logic prevents the script from crashing and provides a clear "Skip" message, allowing the rest of your pipeline to finish successfully.

---

## 4. Statistical Rigor (Hypergeometric Testing)
```r
p_overlap <- phyper(n_shared - 1, n_romi, total_genes - n_romi, n_kroma, lower.tail = FALSE)
jaccard <- n_shared / length(union_ids)
rep_factor <- n_shared / expected_shared
```
- **The Job**: Calculates three high-impact metrics to prove the significance of the overlap:
    - **P-value**: The mathematical probability that the observed overlap happened by chance.
    - **Jaccard Index**: Quantifies the "Similarity" between sets (0 = completely different, 1 = identical).
    - **Representation Factor**: Shows if the overlap is "richer" than expected (e.g., 7.2x).
- **Why it matters**: A Venn diagram shows *that* an overlap exists; this math proves that the overlap is **biologically meaningful**.

---

## 5. Automated "Rigor Summary" Export
```r
write.csv(rigor_summary, 
          file = paste0(res_dir, "/tables/04_03_venn_rigor_stats.csv"),
          row.names = FALSE)
```
- **The Job**: Consolidates the statistical metrics for All, Up, and Down directions into a single reference table.
- **The Reasoning**: This provides a ready-made table for the "Results" section of a paper or thesis, moving the data from a "vague visualization" to "rigorous statistical evidence."

---

## 6. Shared Gene Extraction (Mentor Style)
```r
shared_ids <- Reduce(intersect, sig_list)

romi_df <- results_list[[priority_contrasts[1]]]$df %>% 
  filter(Geneid %in% shared_ids) %>%
  select(Geneid, gene_label, log2FoldChange, padj)

kroma_df <- results_list[[priority_contrasts[2]]]$df %>% 
  filter(Geneid %in% shared_ids) %>%
  select(Geneid, log2FoldChange, padj)

shared_table <- inner_join(romi_df, kroma_df, by = "Geneid", suffix = c("_Romi", "_Kroma")) %>%
  mutate(
    direction = dir_name, 
    max_padj = pmax(padj_Romi, padj_Kroma, na.rm = TRUE),
    avg_abs_lfc = (abs(log2FoldChange_Romi) + abs(log2FoldChange_Kroma)) / 2
  ) %>%
  arrange(max_padj, desc(avg_abs_lfc))
```
- **The Job**: Mathematically identifies the specific genes that are "in the middle" of the Venn circles and exports them to a CSV.
- **Why it matters**: A picture is nice, but your mentor will want the *list* of names. This exports the specific genes that both Kromastat and Romidepsin hit. These are your "Candidate Mechanisms."

---

## 7. Vertical Stack Strategy
```r
p_combined <- wrap_plots(venn_plots, ncol = 1) + 
  plot_annotation(title = paste0("Directional Venn Comparisons: ", group_name))
```
- **The Job**: Combines the All/Up/Down Venns into one tall, vertical image.
- **The Reasoning**: Venn diagrams take up a lot of horizontal space. By stacking them vertically, we preserve the legibility of the gene counts and labels, making the final `04_03_venn_combined.png` a perfect one-page figure for your report.

---

## 8. How to Interpret Your Venn Diagram
The Venn diagram is the ultimate tool for comparing the "Biological Fingerprint" of two or more drugs.

1. **The Overlap (Center)**: This is your **Conserved Mechanism**. The genes in this section respond to *both* drugs. A large overlap in the "Up" Venn confirms that the drugs are likely activating the same therapeutic pathways (e.g., P53 induction).
2. **The Unique Areas (Outer Circles)**: These represent **Drug-Specific Effects**. If Romidepsin has 1,000 unique genes but Kromastat only has 10, it suggests Romidepsin is far more potent (or potentially more toxic).
3. **Directional Consistency Check**: 
    - Compare the "All" Venn to the "Up" and "Down" Venns. 
    - **Success**: Most shared genes in "All" are also shared in "Up" or "Down." This proves the drugs move genes in the *same* direction.
    - **Discovery**: If genes are shared in "All" but *not* in "Up" or "Down," it means the two drugs are moving the same genes in **opposite directions**—a very rare and exciting biological finding!
4. **Gene Lists**: Don't just look at the numbers! Check the [**shared_genes_up.csv**](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/tables/venn_lists/) files. If you see the same gene symbols appearing in the overlap across multiple cell lines, you have found a **Universal Biomarker**.
