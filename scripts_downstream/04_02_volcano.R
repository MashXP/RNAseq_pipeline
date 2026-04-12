# [[scripts_downstream/04_02_volcano.R]]
# Goal: Generate Volcano plots for all doses in specific group.

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(tidyverse)
library(clusterProfiler)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_02_volcano.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

if (species_name == "Human") {
  library(org.Hs.eg.db)
  org_db <- org.Hs.eg.db
} else if (species_name == "Dog") {
  library(org.Cf.eg.db)
  org_db <- org.Cf.eg.db
} else {
  stop("Unsupported species: ", species_name)
}

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# 1.5 ID Mapping (ENSG -> Symbol)
message("Mapping gene IDs...")
gene_ids <- rownames(dds)
gene_map <- bitr(gene_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org_db) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
row.names(gene_map) <- gene_map$ENSEMBL

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 3. Multi-dose Volcano Plots -- individual per dose + stitched combined
# Loop through contrasts for all doses
volcano_plots <- list()
for (dose in names(results_list)) {
  message("Generating Volcano for ", dose)
  res_shrunk <- results_list[[dose]]$shrunk
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  
  # Map ENSG IDs to Gene Symbols for labels
  ensg_ids <- rownames(res_shrunk)
  sym_labels <- gene_map[ensg_ids, "SYMBOL"]
  sym_labels <- ifelse(is.na(sym_labels), ensg_ids, sym_labels) # fallback to ENSG
  
  # Direction-based custom coloring: Red=Up, Blue=Down, Grey=Not Significant
  res_df <- as.data.frame(res_shrunk)
  pCut <- 0.05; fcCut <- 1.0
  keyvals <- case_when(
    res_df$log2FoldChange >  fcCut & !is.na(res_df$padj) & res_df$padj < pCut ~ "#D62728",
    res_df$log2FoldChange < -fcCut & !is.na(res_df$padj) & res_df$padj < pCut ~ "#1F77B4",
    TRUE ~ "grey60"
  )
  names(keyvals) <- case_when(
    keyvals == "#D62728" ~ "Upregulated",
    keyvals == "#1F77B4" ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
  
  # Gene count summary for subtitle
  n_up   <- sum(keyvals == "#D62728", na.rm = TRUE)
  n_down <- sum(keyvals == "#1F77B4", na.rm = TRUE)
  n_ns   <- sum(keyvals == "grey60",  na.rm = TRUE)
  count_subtitle <- paste0("Up: ", n_up, "  |  Down: ", n_down, "  |  Not Significant: ", n_ns)
  
  p_volcano <- EnhancedVolcano(res_shrunk,
                               lab = sym_labels,
                               x = 'log2FoldChange',
                               y = 'padj',
                               title = paste(dose, 'vs DMSO'),
                               subtitle = count_subtitle,
                               ylab = "-log10(Adjusted p-value)",
                               pCutoff = pCut,
                               FCcutoff = fcCut,
                               colCustom = keyvals,
                               colAlpha = 0.8)
  
  volcano_plots[[dose]] <- p_volcano
  out_vol <- paste0(res_dir, "/figures/04_volcano_", safe_name, ".png")
  ggsave(out_vol, p_volcano, width = 10, height = 8, bg = "white")
  message("[OK] Volcano saved  -> 04_volcano_", safe_name, ".png",
          "  |  Up: ", n_up, "  Down: ", n_down, "  NS: ", n_ns)
}

# Stitch all volcano plots side-by-side
if (length(volcano_plots) > 1) {
  combined_volcano <- wrap_plots(volcano_plots, nrow = 1) +
    plot_annotation(
      title = paste0("Volcano Plots: ", group_name, " vs DMSO"),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  ggsave(paste0(res_dir, "/figures/04_volcano_combined.png"),
         combined_volcano, width = 10 * length(volcano_plots), height = 9, bg = "white")
  message("[OK] Combined volcano -> 04_volcano_combined.png",
          "  (", length(volcano_plots), " panels, ", 10 * length(volcano_plots), "x9 in)")
}
