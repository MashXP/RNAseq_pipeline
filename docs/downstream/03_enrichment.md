# Downstream Dissection: 03_enrichment.R

This script answers the question: "What are these genes actually doing?" It translates a list of significant genes into biological processes (GO), metabolic pathways (KEGG), and clinical hallmarks (MSigDB).

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Results RData**: `scripts_downstream/.RData/[Group]/02_deseq_results.RData` (Contains results_list).
    - **Species Org.DB**: `org.Hs.eg.db` or `org.Cf.eg.db`.
- **Processing**: ORA (GO/KEGG) for significant hits, GSEA (Hallmark/GO) based on Wald-statistic ranking.
- **Output**: 
    - **Enrichment CSVs**: `../results/[Group]/tables/03_[go/kegg/gsea]_...csv`.
    - **Final Enrichment RData**: `scripts_downstream/.RData/[Group]/03_enrichment_results.RData`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for functional annotation. See [**libraries.md**](libraries.md) for full technical justifications.
- `clusterProfiler`: Principal package for GO/KEGG/GSEA math.
- `msigdbr`: Pulls latest hallmark sets from MSigDB.
- `tidyverse`: Data manipulation.
- `org.Hs.eg.db` / `org.Cf.eg.db`: Essential species-specific annotation.

---

## 1. Species-Specific Org.DB Selection
```r
if (species_name == "Human") {
  library(org.Hs.eg.db)
  kegg_org <- 'hsa'
} else if (species_name == "Dog") {
  library(org.Cf.eg.db)
  kegg_org <- 'cfa'
}
```
- **The Job**: Loads the correct biological database for the species being analyzed.
- **The Reasoning**: Genes mean different things in different organisms. For your dual-species project, we must use `hsa` (Homo sapiens) and `cfa` (Canis familiaris) to ensure our pathway maps are species-accurate.

---

## 2. ORA (Over-Representation Analysis)
```r
sig_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  rownames()
...
ego <- enrichGO(gene = gene_map$ENSEMBL, OrgDb = org_db, ont = "BP", ...)
```
- **The Job**: Asks, "Are there more genes in 'Pathway X' than we would expect by random chance among my significant hits?"
- **The Reasoning**: This is the "Classic" approach. It uses your most confident genes (Fold Change > 1) to find clear, obvious biological shifts.

---

## 3. High-Performance GSEA Ranking
```r
res_clean <- res_df[!is.na(res_df$stat), ]
ranked_genes <- res_clean$stat
names(ranked_genes) <- rownames(res_clean)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
```
- **The Job**: Creates a "Ranked List" of every single gene in the dataset, from most activated to most suppressed, using the **Wald-Statistic (stat)**.
- **Why it matters**: This is your mentor's "Gold Standard." Unlike ORA, GSEA doesn't ignore "almost significant" genes. By ranking by the Wald-statistic, we incorporate both the size of the change AND the statistical certainty into one number.

---

## 4. GSEA Hallmark Logic
```r
egsea_hallmark <- GSEA(ranked_genes, TERM2GENE = h_t2g, minGSSize = 10, maxGSSize = 500, ...)
```
- **The Job**: Directly checks the **MSigDB Hallmark** gene sets.
- **The Reasoning**: Hallmark sets (produced by the Broad Institute) are curated specifically to represent the 50 most clearly defined biological "themes." For a cancer study like yours, this often provides the most "interpretable" results for a publication.

---

## 5. Robust Error Handling (`tryCatch`)
```r
ego_gsea <- tryCatch({
  gseGO(...)
}, error = function(e) {
  message("  [CAUTION] GSEA GO failed: ", e$message)
  return(NULL)
})
```
- **The Job**: Wraps complex mathematical functions in a "Safety Net."
- **The Reasoning**: GSEA sometimes fails if a dataset is too small or has zero overlap with certain pathways. Instead of the whole pipeline crashing, `tryCatch` skips the failed test and continues with the next one, ensuring you get as many results as possible.

---

## 6. Biological Result Serialization
```r
save(enrichment_results_all, file = paste0("./.RData/", group_name, "/03_enrichment_results.RData"))
```
- **The Job**: Packages all Functional Biology results into an `.RData` file.
- **The Reasoning**: Enrichment analysis can take several minutes per sample. Saving these results once allows the next visualization scripts (**04.x**) to generate infinite plots without needing to re-calculate the math every time.
