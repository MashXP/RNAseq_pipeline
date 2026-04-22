# [[scripts_downstream/04_04_heatmap_pathway.R]]
# Goal: Generate Pathway-Grouped Heatmaps (Leading Edge Analysis).
# Strategy: Visualize the specific genes driving the top Hallmark GSEA results.

library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(clusterProfiler)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_04_heatmap_pathway.R <Group> [Species]")
}
group_name <- args[1]
species_name <- if (length(args) >= 2) args[2] else str_to_title(group_name)

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
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# Transform counts (VST)
vsd <- vst(dds, blind=FALSE)
mat_full <- assay(vsd)

# Mapping for labels
gene_ids <- rownames(vsd)
gene_map <- bitr(gene_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org_db) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
rownames(gene_map) <- gene_map$ENSEMBL

# 2. Iterate through primary contrasts to generate heatmaps
# We prioritize Hallmark GSEA results for the heatmaps
contrasts_to_plot <- names(enrichment_results_all)

for (contrast in contrasts_to_plot) {
  message("\n--- Generating Leading-Edge Heatmap for: ", contrast, " ---")
  
  res_obj <- enrichment_results_all[[contrast]]
  hallmark_res <- res_obj$gsea_hallmark
  
  if (is.null(hallmark_res) || nrow(as.data.frame(hallmark_res)) == 0) {
    message("  Skipping: No significant Hallmark pathways found.")
    next
  }
  
  # Select Top Pathways (e.g., 3 Activated, 3 Suppressed)
  df_halo <- as.data.frame(hallmark_res) %>%
    filter(p.adjust < 0.05) %>%
    mutate(abs_NES = abs(NES)) %>%
    arrange(desc(abs_NES))
  
  if (nrow(df_halo) == 0) {
    message("  Skipping: No pathways with padj < 0.05.")
    next
  }
  
  # Pick top pathways to avoid an infinite heatmap
  top_pathways <- head(df_halo, 6)
  
  expanded_matrix_list <- list()
  annotation_row_list  <- list()
  
  for (i in seq_len(nrow(top_pathways))) {
    path_raw  <- top_pathways$ID[i]
    path_clean <- str_to_title(gsub("_", " ", gsub("HALLMARK_", "", path_raw)))
    
    # Extract Leading Edge (from core_enrichment column)
    # core_enrichment is a list of / separated Ensembl IDs
    genes_path <- str_split(top_pathways$core_enrichment[i], "/")[[1]]
    genes_path <- intersect(genes_path, rownames(mat_full))
    
    # Limit to top 15 genes per pathway by absolute LFC to keep it clean
    # extracted from our results_list
    res_lfc <- as.data.frame(results_list[[contrast]]$shrunk)
    genes_path <- intersect(genes_path, rownames(res_lfc))
    
    if (length(genes_path) > 12) {
      genes_path <- res_lfc[genes_path, ] %>%
        mutate(gene = rownames(.)) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        head(12) %>%
        pull(gene)
    }
    
    if (length(genes_path) > 0) {
      mat_slice <- mat_full[genes_path, , drop = FALSE]
      # Unique rownames to allow same gene in multiple pathways
      rownames(mat_slice) <- paste0(genes_path, "__", i) 
      
      expanded_matrix_list[[path_clean]] <- mat_slice
      
      annot_slice <- data.frame(Pathway = rep(path_clean, length(genes_path)))
      rownames(annot_slice) <- rownames(mat_slice)
      annotation_row_list[[path_clean]] <- annot_slice
    }
  }
  
  if (length(expanded_matrix_list) > 0) {
    mat_pathway  <- do.call(rbind, expanded_matrix_list)
    df_annot_row <- do.call(rbind, annotation_row_list)
    
    # 4. Shorten Sample Names for the plot (Thesis-Ready Aesthetics)
    colnames(mat_pathway) <- colnames(mat_pathway) %>%
      str_remove_all("human_|dog_") %>% 
      str_replace_all("DMSO", "D") %>%
      str_replace_all("Romidepsin|Romi", "R") %>%
      str_replace_all("Kromastat|Kroma", "K") %>%
      str_replace_all("6nM", "") %>%
      str_replace_all("(?i)replicate|rep", "") %>%
      str_replace_all("_+", "_") %>%
      str_remove("^_|_$")
    
    # Standardize (Z-score)
    mat_pathway <- t(scale(t(mat_pathway)))
    mat_pathway <- pmin(pmax(mat_pathway, -1), 1) # Cap at +/- 1 for contrast
    
    # Prep Labels
    orig_ensg <- sub("__.*$", "", rownames(mat_pathway))
    symbols_display <- gene_map[orig_ensg, "SYMBOL"]
    labels_display <- ifelse(is.na(symbols_display), orig_ensg, symbols_display)
    
    # Sorting and Factors
    pathway_factor <- factor(df_annot_row$Pathway, levels = unique(df_annot_row$Pathway))
    
    # Colors
    col_fun <- colorRamp2(c(-1, 0, 1), c("#3B4CC0", "white", "#B40426"))
    n_pw <- nlevels(pathway_factor)
    pw_palette <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(min(n_pw, 8), "Dark2"))(n_pw),
      levels(pathway_factor)
    )
    
    # Metadata for annotation (REORDERED for cell line grouping)
    col_meta <- as.data.frame(colData(vsd))
    if (!"cell_line" %in% colnames(col_meta)) col_meta$cell_line <- species_name
    
    # Define Column Order: Cell Line first, then Condition (Explicitly Romi first)
    col_meta$condition <- factor(col_meta$condition, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM"))
    col_order <- order(col_meta$cell_line, col_meta$condition)
    mat_pathway <- mat_pathway[, col_order]
    col_meta_ordered <- col_meta[col_order, ]
    
    conditions <- col_meta_ordered$condition
    cell_lines <- col_meta_ordered$cell_line
    drug_groups <- factor(ifelse(grepl("Romi", conditions), "Romidepsin", "Kromastat"), 
                          levels = c("Romidepsin", "Kromastat"))
    
    dose_palette <- c(
      "DMSO_Romi"      = "grey85",
      "Romi_6nM"       = "#E41A1C", # Brighter Red
      "DMSO_Kromastat" = "grey70",
      "Kromastat_6nM"  = "#377EB8"  # Brighter Blue
    )

    ha_list <- list(condition = dose_palette)
    ha_list$cell_line <- setNames(RColorBrewer::brewer.pal(max(3, length(unique(cell_lines))), "Accent")[seq_along(unique(cell_lines))], unique(cell_lines))
    
    top_ha <- HeatmapAnnotation(df = col_meta_ordered[, names(ha_list), drop=FALSE], 
                                col = ha_list,
                                show_legend = FALSE)
    
    # The Heatmap
    # Wrap title if too long
    clean_contrast <- str_replace_all(contrast, "_", " ")
    wrapped_title <- str_wrap(paste0(species_name, ": ", clean_contrast), width = 50)
    
    ht <- Heatmap(
      mat_pathway,
      name = "z-score",
      col = col_fun,
      row_split = pathway_factor,
      row_title_side = "right",
      row_title_rot = 0,
      row_title_gp = gpar(fontsize = 9, fontface = "bold"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_split = data.frame(cell_lines, drug_groups),
      column_gap = unit(c(2, 10, 2), "mm"),
      show_row_names = FALSE,
      right_annotation = rowAnnotation(
        id = anno_text(labels_display, gp = gpar(fontsize = 8)),
        Pathway = pathway_factor,
        col = list(Pathway = pw_palette),
        show_legend = FALSE
      ),
      top_annotation = top_ha,
      column_title = wrapped_title,
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 10),
      show_heatmap_legend = FALSE # We will create it manually for custom ordering
    )
    
    # 5. Create Custom Legends for Ordered Stacking
    lgd_shorthand = Legend(labels = c("D = DMSO", "R = Romidepsin", "K = Kromastat"), 
                           title = "Shorthand", 
                           graphics = list(function(x, y, w, h) {}, 
                                           function(x, y, w, h) {}, 
                                           function(x, y, w, h) {}))
    
    lgd_condition = Legend(title = "Condition", 
                           at = names(ha_list$condition), 
                           legend_gp = gpar(fill = ha_list$condition))
    
    lgd_list = list(lgd_shorthand, lgd_condition)
    
    if (!is.null(cell_lines)) {
      lgd_cell = Legend(title = "Cell Line", 
                        at = names(ha_list$cell_line), 
                        legend_gp = gpar(fill = ha_list$cell_line))
      lgd_list = c(lgd_list, list(lgd_cell))
    }
    
    lgd_z = Legend(title = "z-score", 
                   col_fun = col_fun, 
                   at = c(-1, 0, 1))
    lgd_list = c(lgd_list, list(lgd_z))
    
    # Pack them vertically
    pd = packLegend(list = lgd_list, direction = "vertical")
    
    safe_target <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
    png(paste0(res_dir, "/figures/04_04_heatmap_pathway_", safe_target, ".png"), 
        width = 2000, height = 250 + nrow(mat_pathway)*25, res = 150)
    
    # Draw with the packed legend list on the left
    draw(ht, 
         annotation_legend_list = pd,
         heatmap_legend_side = "left", 
         annotation_legend_side = "left",
         padding = unit(c(5, 25, 5, 10), "mm")) # bottom, left, top, right
    dev.off()
    
    message("  [OK] Saved heatmap -> 04_04_heatmap_pathway_", safe_target, ".png")
  }
}

message("\n[OK] Leading-edge heatmap generation complete.")
