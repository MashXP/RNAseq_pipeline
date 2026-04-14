# [[scripts_downstream/03_enrichment.R]]
# Goal: Functional enrichment (GO/KEGG) for specific group and species.

library(clusterProfiler)
library(tidyverse)
library(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 03_enrichment.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

# Load specific org.db based on species
if (species_name == "Human") {
  library(org.Hs.eg.db)
  org_db <- org.Hs.eg.db
  msig_species <- "Homo sapiens"
  kegg_org <- 'hsa'
} else if (species_name == "Dog") {
  # BiocManager::install("org.Cf.eg.db") needs to be ensured by users
  library(org.Cf.eg.db)
  org_db <- org.Cf.eg.db
  msig_species <- "Canis lupus familiaris"
  kegg_org <- 'cfa'
} else {
  stop("Unsupported species: ", species_name)
}

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# Prepare MSigDB Hallmark gene sets (ENSG keys)
# Note: Dog MSigDB might map to symbols or entrez, but we will keep as Ensembl if available
h_t2g <- msigdbr(species = msig_species, category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_results_all <- list()

for (contrast in names(results_list)) {
  message("--- Running Enrichment for: ", contrast, " ---")
  
  res_shrunk <- results_list[[contrast]]$shrunk
  safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
  
  # 2. ORA (Over-Representation Analysis)
  sig_genes <- as.data.frame(res_shrunk) %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    rownames()
  
  if (length(sig_genes) > 5) {
    # Map Ensembl to Entrez for KEGG
    gene_map <- bitr(sig_genes, fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org_db)
    
    # GO ORA
    ego <- enrichGO(gene          = gene_map$ENSEMBL,
                    OrgDb         = org_db,
                    keyType       = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    readable      = FALSE)
    
    # KEGG ORA
    ekegg <- enrichKEGG(gene         = gene_map$ENTREZID,
                        organism     = kegg_org)
    
    write.csv(as.data.frame(ego), file = paste0(res_dir, "/tables/03_go_ora_", safe_name, ".csv"))
    write.csv(as.data.frame(ekegg), file = paste0(res_dir, "/tables/03_kegg_ora_", safe_name, ".csv"))
  } else {
    message("Too few significant genes for ORA in ", contrast)
    ego <- NULL
    ekegg <- NULL
  }

  # 3. GSEA (Gene Set Enrichment Analysis)
  # rank = sign(log2FC) * -log10(pvalue)
  message("Preparing ranked list for GSEA...")

  res_df <- as.data.frame(results_list[[contrast]]$res)
  res_df$ranking_metric <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)

  # Cap Infinite values
  max_rank <- max(res_df$ranking_metric[is.finite(res_df$ranking_metric)], na.rm = TRUE)
  if (is.infinite(max_rank)) max_rank <- 100
  res_df$ranking_metric[is.infinite(res_df$ranking_metric)] <- max_rank * sign(res_df$log2FoldChange[is.infinite(res_df$ranking_metric)])
  res_df <- res_df[!is.na(res_df$ranking_metric), ]

  ranked_genes <- res_df$ranking_metric
  names(ranked_genes) <- rownames(res_df)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  message("Running GSEA (fgsea engine)...")
  # Standard GSEA wrapper for production
  ego_gsea <- tryCatch({
    gseGO(geneList     = ranked_genes,
          OrgDb        = org_db,
          keyType      = "ENSEMBL",
          ont          = "BP",
          minGSSize    = 10,
          maxGSSize    = 500,
          pvalueCutoff = 0.05,
          verbose      = FALSE)
  }, error = function(e) {
    message("  [CAUTION] GSEA GO failed for ", contrast, ": ", e$message)
    return(NULL)
  })
  
  # Save GSEA results (including activated/suppressed)
  if (!is.null(ego_gsea) && nrow(as.data.frame(ego_gsea)) > 0) {
    gsea_df <- as.data.frame(ego_gsea)
    write.csv(gsea_df, file = paste0(res_dir, "/tables/03_gsea_functional_", safe_name, ".csv"))
    
    # Identify Activated (NES > 0) and Suppressed (NES < 0)
    activated <- gsea_df %>% filter(NES > 0)
    suppressed <- gsea_df %>% filter(NES < 0)
    
    write.csv(activated, file = paste0(res_dir, "/tables/03_gsea_activated_", safe_name, ".csv"))
    write.csv(suppressed, file = paste0(res_dir, "/tables/03_gsea_suppressed_", safe_name, ".csv"))
  }

  # 4. GSEA Hallmark
  message("Running GSEA Hallmark for: ", contrast)
  # Standard GSEA wrapper for production
  egsea_hallmark <- tryCatch({
    GSEA(ranked_genes, TERM2GENE = h_t2g, minGSSize = 10, pvalueCutoff = 0.05)
  }, error = function(e) {
    message("  [CAUTION] GSEA Hallmark failed for ", contrast, ": ", e$message)
    return(NULL)
  })
  if (!is.null(egsea_hallmark) && nrow(as.data.frame(egsea_hallmark)) > 0) {
    write.csv(as.data.frame(egsea_hallmark),
              file = paste0(res_dir, "/tables/03_gsea_hallmark_", safe_name, ".csv"))
  }

  enrichment_results_all[[contrast]] <- list(
    ora_go = ego,
    ora_kegg = ekegg,
    gsea_go = ego_gsea,
    gsea_hallmark = egsea_hallmark
  )
}

# 4. Save results
save(enrichment_results_all, file = paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

message("Enrichment analysis complete for ", group_name, ". All results saved to ", res_dir, "/tables/")
