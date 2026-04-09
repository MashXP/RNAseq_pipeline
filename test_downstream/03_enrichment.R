# [[test_downstream/03_enrichment.R]]
# TEST VERSION: Simplified functional enrichment.
# Logic: Maps genes to Entrez and runs basic GO ORA.
# Divergence: RELAXED P-VALUE (pvalue < 1.0) to verify technical mapping logic on small data.

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# 1. Load test DGE results
load("./.RData/02_test_deseq_results.RData")

dir.create("tables_test", showWarnings = FALSE)

# 2. Simplified logic for testing (just mapping and trying one enrichment)
message("--- [TEST] Running Enrichment Check ---")

res <- as.data.frame(results(dds))
# LOWERED THRESHOLD FOR TEST
sig_genes <- res %>%
  filter(pvalue < 0.5) %>% 
  rownames()

if (length(sig_genes) > 2) {
  message("Found ", length(sig_genes), " potential genes for test mapping.")
  
  # Try mapping
  gene_map <- bitr(sig_genes, fromType = "ENSEMBL",
                   toType = c("ENTREZID", "SYMBOL"),
                   OrgDb = org.Hs.eg.db)
  
  # Try basic GO
  ego <- enrichGO(gene          = gene_map$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1, # Accept everything for test
                  qvalueCutoff  = 1)
  
  write.csv(as.data.frame(ego), file = "tables_test/03_test_go_report.csv")
} else {
  message("[TEST] Not enough genes even for test thresholds.")
}

message("[TEST] Enrichment script finished.")
