# [[scripts_downstream/03_enrichment.R]]
# Goal: Functional enrichment (GO/KEGG/Hallmark).
# Upgrade: Uses Wald-statistic ranking (stat) for GSEA and extracts leading edge genes.

library(clusterProfiler)
library(tidyverse)
library(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 03_enrichment.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

# Load specific org.db based on species
if (species_name == "Human") {
  library(org.Hs.eg.db)
  org_db <- org.Hs.eg.db
  msig_species <- "Homo sapiens"
  kegg_org <- 'hsa'
} else if (species_name == "Canine") {
  library(org.Cf.eg.db)
  org_db <- org.Cf.eg.db
  msig_species <- "Canis lupus familiaris"
  kegg_org <- 'cfa'
} else {
  stop("Unsupported species: ", species_name)
}

res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# Prepare MSigDB Hallmark gene sets (ENSG keys)
h_t2g <- msigdbr(species = msig_species, category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)

enrichment_results_all <- list()

for (contrast in names(results_list)) {
  message("\n--- Running Enrichment for: ", contrast, " ---")
  
  res_obj <- results_list[[contrast]]$res
  res_df <- as.data.frame(res_obj)
  safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
  
  # 2. ORA (Over-Representation Analysis)
  sig_genes <- res_df %>%
    filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
    rownames()
  
  ego <- NULL
  ekegg <- NULL
  if (length(sig_genes) > 5) {
    gene_map <- tryCatch({
      bitr(sig_genes, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org_db)
    }, error = function(e) return(NULL))
    
    if (!is.null(gene_map)) {
      ego <- tryCatch({
        enrichGO(gene = gene_map$ENSEMBL, OrgDb = org_db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
      }, error = function(e) return(NULL))
      
      ekegg <- tryCatch({
        enrichKEGG(gene = gene_map$ENTREZID, organism = kegg_org)
      }, error = function(e) return(NULL))
    }
  }

  # 3. GSEA (Gene Set Enrichment Analysis)
  # Raking: Using the Wald-statistic (stat column) as per Mentor's gold standard
  message("Preparing ranked list (stat-based)...")
  
  # Clean result: remove NA stats
  res_clean <- res_df[!is.na(res_df$stat), ]
  ranked_genes <- res_clean$stat
  names(ranked_genes) <- rownames(res_clean)
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  # A. GSEA Hallmark (Highest statistical priority)
  message("Running GSEA Hallmark...")
  egsea_hallmark <- tryCatch({
    GSEA(ranked_genes, TERM2GENE = h_t2g, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  }, error = function(e) {
    message("  [CAUTION] GSEA Hallmark failed: ", e$message)
    return(NULL)
  })

  # B. GSEA GO
  message("Running GSEA GO...")
  ego_gsea <- tryCatch({
    gseGO(geneList = ranked_genes, OrgDb = org_db, keyType = "ENSEMBL", ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  }, error = function(e) {
    message("  [CAUTION] GSEA GO failed: ", e$message)
    return(NULL)
  })
  
  # Save ORA results if they exist
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file = paste0(res_dir, "/tables/03_go_ora_", safe_name, ".csv"), row.names = FALSE)
  }
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write.csv(as.data.frame(ekegg), file = paste0(res_dir, "/tables/03_kegg_ora_", safe_name, ".csv"), row.names = FALSE)
  }

  # Save GSEA results (including activated/suppressed)
  if (!is.null(ego_gsea) && nrow(as.data.frame(ego_gsea)) > 0) {
    gsea_df <- as.data.frame(ego_gsea)
    write.csv(gsea_df, file = paste0(res_dir, "/tables/03_gsea_functional_", safe_name, ".csv"), row.names = FALSE)
  }

  if (!is.null(egsea_hallmark) && nrow(as.data.frame(egsea_hallmark)) > 0) {
    write.csv(as.data.frame(egsea_hallmark),
              file = paste0(res_dir, "/tables/03_gsea_hallmark_", safe_name, ".csv"), row.names = FALSE)
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

message("\n[OK] Enrichment complete. Stat-based GSEA results saved.")
