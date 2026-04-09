# 03_enrichment.R
# Goal: Functional enrichment (GO/KEGG) for MiaPAca-2 DGE results.

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# 1. Load DGE results
load("./.RData/02_deseq_results.RData")

enrichment_results_all <- list()

for (dose in names(results_list)) {
  message("--- Running Enrichment for: ", dose, " ---")
  
  res_shrunk <- results_list[[dose]]$shrunk
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  
  # 2. ORA (Over-Representation Analysis)
  sig_genes <- as.data.frame(res_shrunk) %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    rownames()
  
  if (length(sig_genes) > 5) {
    # Map Ensembl to Entrez for KEGG
    gene_map <- bitr(sig_genes, fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org.Hs.eg.db)
    
    # GO ORA
    ego <- enrichGO(gene          = gene_map$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    readable      = TRUE)
    
    # KEGG ORA
    ekegg <- enrichKEGG(gene         = gene_map$ENTREZID,
                        organism     = 'hsa')
    
    write.csv(as.data.frame(ego), file = paste0("../results/tables/03_go_ora_", safe_name, ".csv"))
    write.csv(as.data.frame(ekegg), file = paste0("../results/tables/03_kegg_ora_", safe_name, ".csv"))
  } else {
    message("Too few significant genes for ORA in ", dose)
    ego <- NULL
    ekegg <- NULL
  }

  # 3. GSEA (Gene Set Enrichment Analysis)
  # Prepare ranked gene list (LFC based)
  ranked_genes <- res_shrunk$log2FoldChange
  names(ranked_genes) <- rownames(res_shrunk)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  message("Running GSEA...")
  ego_gsea <- gseGO(geneList      = ranked_genes,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = "BP",
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    pvalueCutoff  = 0.05,
                    verbose       = FALSE,
                    OrgDb         = org.Hs.eg.db)
  
  # Save GSEA results (including activated/suppressed)
  if (!is.null(ego_gsea)) {
    gsea_df <- as.data.frame(ego_gsea)
    write.csv(gsea_df, file = paste0("../results/tables/03_gsea_functional_", safe_name, ".csv"))
    
    # Identify Activated (NES > 0) and Suppressed (NES < 0)
    activated <- gsea_df %>% filter(NES > 0)
    suppressed <- gsea_df %>% filter(NES < 0)
    
    write.csv(activated, file = paste0("../results/tables/03_gsea_activated_", safe_name, ".csv"))
    write.csv(suppressed, file = paste0("../results/tables/03_gsea_suppressed_", safe_name, ".csv"))
  }

  enrichment_results_all[[dose]] <- list(ora_go = ego, ora_kegg = ekegg, gsea_go = ego_gsea)
}

# 4. Save results
save(enrichment_results_all, file = "./.RData/03_enrichment_results.RData")

message("Enrichment analysis complete. All results saved to ../results/tables/")
