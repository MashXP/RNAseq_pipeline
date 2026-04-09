# [[test_downstream/04_visualization.R]]
# TEST VERSION: Visual verification (PCA, Volcano, Venn, Heatmap).
# Logic: Generates production-style plots with test data.
# Divergence: High thresholds (p < 0.5) to ensure some plot features are visible.

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)
library(ggVennDiagram)

# 1. Load data
load("./.RData/02_test_deseq_results.RData")
# Load enrichment if file exists
if (file.exists("./.RData/03_test_enrichment_results.RData")) {
  load("./.RData/03_test_enrichment_results.RData")
}

dir.create("plots_test", showWarnings = FALSE)

# 2. PCA Plot
vsd <- vst(dds, blind=FALSE)
p_pca <- plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  ggtitle("[TEST] PCA Plot: MiaPAca-2 (Chr21 Only)")

ggsave("plots_test/04_pca_plot.png", p_pca, width = 8, height = 6)

# 3. Volcano Plot Loop
for (dose in names(results_list)) {
  res_shrunk <- results_list[[dose]]$shrunk
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  
  message("[TEST] Generating Volcano for ", dose)
  p_volcano <- EnhancedVolcano(res_shrunk,
                               lab = rownames(res_shrunk),
                               x = 'log2FoldChange',
                               y = 'pvalue', # Use raw pvalue for test visibility
                               title = paste('[TEST]', dose),
                               pCutoff = 0.5, # HIGHLY RELAXED FOR TEST
                               FCcutoff = 0.5)
  
  ggsave(paste0("plots_test/04_volcano_", safe_name, ".png"), 
         p_volcano, width = 10, height = 8)
}

# 4. Venn Diagram (Simplified overlap check)
sig_list <- lapply(results_list, function(x) {
  rownames(subset(as.data.frame(x$shrunk), pvalue < 0.5))
})

if (length(sig_list) >= 2) {
  p_venn <- ggVennDiagram(sig_list) +
    ggtitle("[TEST] Degree of Overlap in DEGs")
  ggsave("plots_test/04_venn_diagram.png", p_venn, width = 8, height = 8)
}

# 5. Default Heatmap (Top 30 Variable genes)
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat_var  <- assay(vsd)[top_genes, ]
mat_var  <- mat_var - rowMeans(mat_var)

png("plots_test/04_heatmap_top_variable.png", width = 800, height = 600)
pheatmap(mat_var, annotation_col = as.data.frame(colData(dds)[,"condition", drop=F]), 
         main = "[TEST] Top 30 Most Variable Genes (Chr21)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

message("[TEST] Visualization complete. Results in plots_test/")
