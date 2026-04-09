# 04_visualization.R
# Goal: Generate PCA, Volcano, and Heatmap plots for MiaPAca-2 results.

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)

# 1. Load data
load("./.RData/02_deseq_results.RData")
load("./.RData/03_enrichment_results.RData")

library(ggVennDiagram)

# 2. PCA Plot
# Transform counts (VST)
vsd <- vst(dds, blind=FALSE)
p_pca <- plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  ggtitle("PCA Plot: MiaPAca-2 Drug Response")

ggsave("../results/figures/04_pca_plot.png", p_pca, width = 8, height = 6)

# 3. Multi-dose Volcano Plots
for (dose in names(results_list)) {
  res_shrunk <- results_list[[dose]]$shrunk
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  
  p_volcano <- EnhancedVolcano(res_shrunk,
                               lab = rownames(res_shrunk),
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = paste(dose, 'vs DMSO (NG)'),
                               pCutoff = 0.05,
                               FCcutoff = 1.0)
  
  ggsave(paste0("../results/figures/04_volcano_", safe_name, ".png"), 
         p_volcano, width = 10, height = 8)
}

# 4. Venn Diagram (3.1a)
# Extract significant gene sets for each dose
sig_list <- lapply(results_list, function(x) {
  rownames(subset(as.data.frame(x$shrunk), padj < 0.05 & abs(log2FoldChange) > 1))
})

# Filter out empty sets if any
sig_list <- sig_list[sapply(sig_list, length) > 0]

if (length(sig_list) >= 2) {
  p_venn <- ggVennDiagram(sig_list) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle("Overlap of DEGs across Kromastat Doses")
  
  ggsave("../results/figures/04_venn_diagram_doses.png", p_venn, width = 8, height = 8)
}

# 5. Pathway-specific Heatmap (3.1b)
# Extract genes from the top enriched GO terms (Biological Process)
# We'll use the results from the 5nM group as primary pathway source
top_go <- enrichment_results_all[["5nM"]]$ora_go
if (!is.null(top_go)) {
  # Get top 3 pathways
  pathway_genes <- top_go@result %>%
    filter(p.adjust < 0.05) %>%
    head(3) %>%
    pull(geneID) %>%
    str_split("/") %>%
    unlist() %>%
    unique()
  
  if (length(pathway_genes) > 2) {
    # Subset VSD to these genes
    pathway_genes_ensembl <- pathway_genes # enrichGO with readable=TRUE might produce symbols, let's check
    # If readable=TRUE was used, they are symbols. We need to match rownames of vsd (ENSEMBL)
    # Actually enrichGO with keyType="ENSEMBL" and readable=TRUE returns symbols in geneID or ENSEMBL?
    # Usually it replaces the ID. I'll try to find them in rownames.
    
    mat_pathway <- assay(vsd)[rownames(vsd) %in% pathway_genes_ensembl, ]
    
    if (nrow(mat_pathway) > 5) {
      mat_pathway <- mat_pathway - rowMeans(mat_pathway)
      df_annot <- as.data.frame(colData(dds)[,c("condition")])
      colnames(df_annot) <- "Concentration"
      
      png("../results/figures/04_heatmap_pathway_genes.png", width = 1000, height = 800)
      pheatmap(mat_pathway, 
               annotation_col = df_annot,
               main = "Expression of Top Enriched Pathway DEGs",
               show_rownames = (nrow(mat_pathway) < 50), # Hide if too many
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
      dev.off()
    }
  }
}

# 6. Default Heatmap (Top 50 Variable)
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat_var  <- assay(vsd)[top_genes, ]
mat_var  <- mat_var - rowMeans(mat_var)

png("../results/figures/04_heatmap_top_variable.png", width = 800, height = 600)
pheatmap(mat_var, annotation_col = as.data.frame(colData(dds)[,"condition", drop=F]), 
         main = "Top 50 Most Variable Genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

message("Visualization complete. All plots saved to ../results/figures/")
