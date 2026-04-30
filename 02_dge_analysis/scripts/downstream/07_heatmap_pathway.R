# [[02_dge_analysis/scripts/downstream/07_heatmap_pathway.R]]
# Goal: Exact replication of mentor's pheatmap styling for grouped pathways.

library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(clusterProfiler)

# --- Dynamic Pathing Setup ---
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) { curr <- dirname(curr) }
  project_root <- curr
}
data_dir <- file.path(project_root, "_data")
results_root <- file.path(project_root, "results")
rdata_dir <- file.path(data_dir, "RData", "02_dge_analysis")
# -----------------------------


args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 07_heatmap_pathway.R <Group>")
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

if (species_name == "Human") {
  library(org.Hs.eg.db)
  org_db <- org.Hs.eg.db
} else {
  library(org.Cf.eg.db)
  org_db <- org.Cf.eg.db
}

res_dir <- file.path(results_root, "02_dge_analysis", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# --- Visualization Tweaks (Edit these to adjust figure sizes) ---
SINGLE_WIDTH <- 1200
SHARED_WIDTH <- 1750
rdata_dir <- file.path(data_dir, "RData", "02_dge_analysis")
# -----------------------------------------------------------------

# 1. Load data
load(file.path(rdata_dir, group_name, "02_deseq_results.RData"))
load(file.path(rdata_dir, group_name, "03_enrichment_results.RData"))

# Mapping for labels
gene_map <- bitr(rownames(dds), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org_db) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
rownames(gene_map) <- gene_map$ENSEMBL

# Colors
col_fun <- colorRamp2(c(-1, 0, 1), c("#3B4CC0", "white", "#B40426"))

message("\n--- Generating Mentor-Style Pathway Heatmaps ---")

for (contrast in names(enrichment_results_all)) {
  res <- enrichment_results_all[[contrast]]$gsea_hallmark
  if (is.null(res)) next
  
  df_res <- as.data.frame(res) %>% filter(p.adjust < 0.05)
  if (nrow(df_res) == 0) next
  
  # Select Top Pathways (Activated and Suppressed)
  top_act <- head(df_res %>% filter(NES > 0) %>% arrange(desc(NES)), 3)
  top_sup <- head(df_res %>% filter(NES < 0) %>% arrange(NES), 2)
  sel_paths <- rbind(top_act, top_sup)
  
  # Build heatmap_df
  gene_rows <- list()
  for (i in seq_len(nrow(sel_paths))) {
    path_name <- sel_paths$Description[i]
    path_clean <- str_to_title(gsub("_", " ", gsub("HALLMARK_", "", path_name)))
    
    # Extract leading edge and map to valid ENSEMBL
    genes <- unique(unlist(strsplit(sel_paths$core_enrichment[i], "/")))
    genes <- intersect(genes, rownames(dds))
    genes <- head(genes, 12) # Mentor logic: Top 12 genes per pathway
    
    if (length(genes) > 0) {
      gene_rows[[length(gene_rows) + 1]] <- data.frame(
        Geneid = genes,
        gene_label = gene_map[genes, "SYMBOL"],
        pathway_display = path_clean,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(gene_rows) == 0) next
  heatmap_df <- do.call(rbind, gene_rows)
  heatmap_df <- heatmap_df[!duplicated(heatmap_df$Geneid), ]
  
  # Fetch DESeq2 stats to replicate mentor's exact sorting (Fixes the internal mosaic)
  res_df <- as.data.frame(results_list[[contrast]]$shrunk)
  res_df$Geneid <- rownames(res_df)
  heatmap_df <- merge(heatmap_df, res_df[, c("Geneid", "padj", "log2FoldChange")], by="Geneid", all.x=TRUE)
  
  # Sort: Pathway first, then significance, then effect size
  heatmap_df <- heatmap_df[order(heatmap_df$pathway_display, heatmap_df$padj, -abs(heatmap_df$log2FoldChange)), ]
  
  # Parse contrast to strictly isolate cell lines and conditions
  parts <- str_split(contrast, "_vs_")[[1]]
  treat_full <- parts[1]
  control_full <- parts[2]
  
  cell_line_filter <- ""
  cl_list <- unique(as.character(colData(dds)$cell_line))
  for (cl_name in cl_list) {
    if (grepl(paste0("^", cl_name, "_"), treat_full)) {
      cell_line_filter <- cl_name
      break
    }
  }
  
  conds <- c(gsub(paste0("^", cell_line_filter, "_"), "", treat_full),
             gsub(paste0("^", cell_line_filter, "_"), "", control_full))
  
  all_samples <- rownames(colData(dds))
  sample_ids <- all_samples[colData(dds)$condition %in% conds]
  
  if (cell_line_filter != "") {
    sample_ids <- sample_ids[grepl(cell_line_filter, sample_ids)]
  }
  
  dds_local <- dds[, sample_ids]
  design(dds_local) <- ~ 1
  vsd_local <- vst(dds_local, blind = FALSE)
  
  # Extract matrix in the EXACT sorted order
  mat <- assay(vsd_local)[heatmap_df$Geneid, ]
  rownames(mat) <- heatmap_df$gene_label
  
  mat <- t(scale(t(mat))) # Scale rows
  mat <- pmin(pmax(mat, -1), 1) # Cap at +/- 1 for contrast
  
  # 4. Shorten Sample Names for the plot (Thesis-Ready Aesthetics)
  colnames(mat) <- colnames(mat) %>%
    str_remove_all("human_|canine_") %>% 
    str_replace_all("DMSO", "D") %>%
    str_replace_all("Romidepsin", "R") %>%
    str_replace_all("Kromastat", "K") %>%
    str_replace_all("6nM", "") %>%
    str_replace_all("(?i)replicate|rep", "") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_|_$")
  
  # Annotations and Meta
  col_meta <- as.data.frame(colData(vsd_local)[, c("condition", "cell_line"), drop = FALSE])
  colnames(col_meta) <- c("Treatment", "Cell_Line")
  
  # Sorting and Factors
  pathway_factor <- factor(heatmap_df$pathway_display, levels = unique(heatmap_df$pathway_display))
  
  # Define Column Order
  col_order <- order(col_meta$Treatment)
  mat <- mat[, col_order]
  col_meta_ordered <- col_meta[col_order, , drop = FALSE]
  
  conditions <- col_meta_ordered$Treatment
  cell_lines <- col_meta_ordered$Cell_Line
  
  # Colors
  n_pw <- nlevels(pathway_factor)
  pw_palette <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(min(n_pw, 8), "Dark2"))(n_pw),
    levels(pathway_factor)
  )
  
  dose_palette <- c(
    "DMSO_Romidepsin"      = "grey85",
    "Romidepsin_6nM"       = "#E41A1C", # Brighter Red
    "DMSO_Kromastat" = "grey70",
    "Kromastat_6nM"  = "#377EB8"  # Brighter Blue
  )
  
  ha_list <- list(Treatment = dose_palette[levels(factor(conditions))])
  # Dynamic Cell Line coloring
  cl_levels <- levels(factor(cell_lines))
  cl_colors <- setNames(colorRampPalette(c("#7FC97F", "#BEAED4"))(length(cl_levels)), cl_levels)
  ha_list$Cell_Line <- cl_colors
  
  top_ha <- HeatmapAnnotation(df = col_meta_ordered[, c("Treatment", "Cell_Line"), drop=FALSE], 
                              col = ha_list,
                              show_legend = FALSE,
                              annotation_name_gp = gpar(fontsize = 15, fontface = "bold"))
  
  # Prepare Row Labels
  labels_display <- heatmap_df$gene_label
  
  # Wrap title if too long
  clean_contrast <- str_replace_all(contrast, "_", " ")
  wrapped_title <- str_wrap(paste0(species_name, ": ", clean_contrast), width = 50)
  
  ht <- Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    row_split = pathway_factor,
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 13.5, fontface = "bold"),
    cluster_rows = FALSE,
    cluster_columns = FALSE, # We manually ordered them
    column_split = factor(conditions, levels = unique(conditions)),
    column_gap = unit(2, "mm"),
    show_row_names = FALSE,
    right_annotation = rowAnnotation(
      id = anno_text(labels_display, gp = gpar(fontsize = 12)),
      Pathway = pathway_factor,
      col = list(Pathway = pw_palette),
      show_legend = FALSE,
      annotation_name_rot = 45,
      annotation_name_side = "bottom",
      annotation_name_gp = gpar(fontsize = 15, fontface = "bold"),
      show_annotation_name = c(id = FALSE, Pathway = TRUE)
    ),
    top_annotation = top_ha,
    column_title = wrapped_title,
    column_title_gp = gpar(fontsize = 18, fontface = "bold"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 15),
    show_heatmap_legend = FALSE # Custom legend handling
  )
  
  # 5. Create Custom Legends for Ordered Stacking
  lgd_shorthand = Legend(labels = c("D = DMSO (Control)", "R = Romidepsin", "K = Kromastat"), 
                         title = "Shorthand", 
                         graphics = list(function(x, y, w, h) {}, 
                                         function(x, y, w, h) {}, 
                                         function(x, y, w, h) {}))
  
  lgd_condition = Legend(title = "Treatment", 
                         at = names(ha_list$Treatment), 
                         legend_gp = gpar(fill = ha_list$Treatment))
  
  lgd_cell = Legend(title = "Cell Line", 
                    at = names(ha_list$Cell_Line), 
                    legend_gp = gpar(fill = ha_list$Cell_Line))
  
  lgd_z = Legend(title = "z-score", 
                 col_fun = col_fun, 
                 at = c(-1, 0, 1))
                 
  lgd_list = list(lgd_shorthand, lgd_condition, lgd_cell, lgd_z)
  
  pd = packLegend(list = lgd_list, direction = "vertical")
  
  safe_target <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
  png(paste0(res_dir, "/figures/07_heatmap_pathway_", safe_target, ".png"), 
      width = SINGLE_WIDTH, height = 250 + nrow(mat)*25, res = 150)
  
  draw(ht, 
       annotation_legend_list = pd,
       heatmap_legend_side = "left", 
       annotation_legend_side = "left",
       padding = unit(c(5, 25, 5, 10), "mm")) # bottom, left, top, right
  dev.off()
  
  message("  [OK] Saved heatmap -> 07_heatmap_pathway_", safe_target, ".png")
}

# --- Shared Legend Components (Reused Below) ---
lgd_shorthand = Legend(labels = c("D = DMSO (Control)", "R = Romidepsin", "K = Kromastat"), 
                       title = "Shorthand", 
                       graphics = list(function(x, y, w, h) {}, 
                                       function(x, y, w, h) {}, 
                                       function(x, y, w, h) {}))

lgd_condition = Legend(title = "Treatment", 
                       at = c("DMSO_Romidepsin", "Romidepsin_6nM", "DMSO_Kromastat", "Kromastat_6nM"), 
                       legend_gp = gpar(fill = c("grey85", "#E41A1C", "grey70", "#377EB8")))

lgd_cell = Legend(title = "Cell Line", 
                  at = unique(as.character(colData(dds)$cell_line)), 
                  legend_gp = gpar(fill = colorRampPalette(c("#7FC97F", "#BEAED4"))(length(unique(colData(dds)$cell_line)))))

lgd_z = Legend(title = "z-score", 
               col_fun = col_fun, 
               at = c(-1, 0, 1))

lgd_list = list(lgd_shorthand, lgd_condition, lgd_cell, lgd_z)
pd_common = packLegend(list = lgd_list, direction = "vertical")

# --- Generating Shared Leading-Edge Heatmaps (Mentor's Logic Parity) ---
if ("cell_line" %in% colnames(colData(dds))) {
  message("\n--- Generating Shared Leading-Edge Heatmaps (Mentor's Logic) ---")
  cell_lines <- unique(as.character(colData(dds)$cell_line))
  
  for (cl in cell_lines) {
    c1 <- paste0(cl, "_Romidepsin_6nM_vs_DMSO_Romidepsin")
    c2 <- paste0(cl, "_Kromastat_6nM_vs_DMSO_Kromastat")
    
    if (!all(c(c1, c2) %in% names(results_list))) next
    
    # 1. Identify Shared Significant Pathways
    p1 <- as.data.frame(enrichment_results_all[[c1]]$gsea_hallmark)
    p2 <- as.data.frame(enrichment_results_all[[c2]]$gsea_hallmark)
    
    if (nrow(p1) == 0 || nrow(p2) == 0) next
    
    shared_pws <- intersect(p1$Description[p1$p.adjust < 0.05], 
                            p2$Description[p2$p.adjust < 0.05])
    
    if (length(shared_pws) == 0) {
      message("  [SKIP] No shared significant pathways for ", cl)
      next
    }
    
    # 2. Extract Shared Leading-Edge Genes using Mentor's Scoring
    df1 <- results_list[[c1]]$df
    df2 <- results_list[[c2]]$df
    rownames(df1) <- df1$Geneid
    rownames(df2) <- df2$Geneid
    
    selected_genes_list <- list()
    for (pw in shared_pws) {
      genes1 <- unlist(strsplit(p1$core_enrichment[p1$Description == pw], "/"))
      genes2 <- unlist(strsplit(p2$core_enrichment[p2$Description == pw], "/"))
      
      # Use intersection of leading edge (with fallback to union like mentor)
      shared_le <- intersect(genes1, genes2)
      if (length(shared_le) == 0) shared_le <- union(genes1, genes2)
      
      # Filter to valid IDs and calculate mentor's selection score
      valid_genes <- intersect(shared_le, rownames(df1))
      if (length(valid_genes) == 0) next
      
      tmp_df <- df1[valid_genes, ] %>%
        left_join(df2[valid_genes, ] %>% dplyr::select(Geneid, log2FoldChange, padj), 
                  by = "Geneid", suffix = c("_R1", "_K1")) %>%
        mutate(
          same_dir = sign(log2FoldChange_R1) == sign(log2FoldChange_K1),
          shared_sig = padj_R1 < 0.05 & padj_K1 < 0.05 & abs(log2FoldChange_R1) > 2 & abs(log2FoldChange_K1) > 2,
          # Mentor's exact scoring logic:
          score = abs(log2FoldChange_R1) + abs(log2FoldChange_K1) + 
                  ifelse(shared_sig & same_dir, 4, 0) + 
                  ifelse(same_dir, 1.5, 0),
          pathway_display = str_to_title(gsub("_", " ", gsub("HALLMARK_", "", pw)))
        ) %>%
        arrange(desc(score)) %>%
        head(12) # Mentor's limit per pathway
      
      selected_genes_list[[pw]] <- tmp_df
    }
    
    if (length(selected_genes_list) == 0) next
    shared_df <- do.call(rbind, selected_genes_list)
    shared_df <- shared_df[!duplicated(shared_df$Geneid), ]
    
    # Final pool capped at 50 and sorted for visualization
    shared_df <- head(shared_df[order(shared_df$pathway_display, -shared_df$score), ], 50)
    pathway_factor_shared <- factor(shared_df$pathway_display, levels = unique(shared_df$pathway_display))
    
    # 4. Extract matrix for all relevant samples
    relevant_conds <- c("DMSO_Romidepsin", "Romidepsin_6nM", "DMSO_Kromastat", "Kromastat_6nM")
    sample_ids <- colnames(dds)[colData(dds)$cell_line == cl & colData(dds)$condition %in% relevant_conds]
    
    dds_cl <- dds[, sample_ids]
    design(dds_cl) <- ~ 1
    vsd_cl <- vst(dds_cl, blind = FALSE)
    
    mat <- assay(vsd_cl)[shared_df$Geneid, ]
    rownames(mat) <- shared_df$gene_label
    mat <- t(scale(t(mat)))
    mat <- pmin(pmax(mat, -1), 1)
    
    # Shorten Names
    colnames(mat) <- colnames(mat) %>%
      str_remove_all("human_|canine_") %>% 
      str_replace_all("DMSO", "D") %>%
      str_replace_all("Romidepsin", "R") %>%
      str_replace_all("Kromastat", "K") %>%
      str_replace_all("6nM", "") %>%
      str_replace_all("(?i)replicate|rep", "") %>%
      str_replace_all("_+", "_") %>%
      str_remove("^_|_$")
      
    # Annotations
    col_meta <- as.data.frame(colData(vsd_cl)[, "condition", drop = FALSE])
    colnames(col_meta) <- "Treatment"
    col_meta$Cell_Line <- cl
    
    # Mentor's order: Treatment first, then DMSO
    # We'll use: Romidepsin, DMSO_Romidepsin, Kromastat, DMSO_Kromastat (Adjusted to match image logic)
    image_order <- c("Romidepsin_6nM", "DMSO_Romidepsin", "Kromastat_6nM", "DMSO_Kromastat")
    col_order <- order(factor(col_meta$Treatment, levels = image_order))
    mat <- mat[, col_order]
    col_meta_ordered <- col_meta[col_order, , drop = FALSE]
    
    # Dynamic palette for the shared heatmap section
    cl_all <- unique(as.character(colData(dds)$cell_line))
    ha_list <- list(
      Treatment = c("DMSO_Romidepsin" = "grey85", "Romidepsin_6nM" = "#E41A1C", "DMSO_Kromastat" = "grey70", "Kromastat_6nM" = "#377EB8"),
      Cell_Line = setNames(colorRampPalette(c("#7FC97F", "#BEAED4"))(length(cl_all)), cl_all)
    )
    
    top_ha <- HeatmapAnnotation(df = col_meta_ordered[, c("Treatment", "Cell_Line"), drop=FALSE], 
                                col = ha_list,
                                show_legend = FALSE,
                                annotation_name_gp = gpar(fontsize = 15, fontface = "bold"))
    
    n_pw_s <- nlevels(pathway_factor_shared)
    pw_palette_s <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(min(n_pw_s, 8), "Set3"))(n_pw_s),
      levels(pathway_factor_shared)
    )
    
    wrapped_title <- str_wrap(paste0(species_name, " ", cl, ": Shared Leading-Edge Overlap"), width = 40)
    
    ht <- Heatmap(
      mat,
      name = "z-score",
      col = col_fun,
      row_split = pathway_factor_shared,
      row_title_side = "right",
      row_title_rot = 0,
      row_title_gp = gpar(fontsize = 13.5, fontface = "bold"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_split = factor(col_meta_ordered$Treatment, levels = image_order),
      column_gap = unit(2, "mm"),
      show_row_names = FALSE,
      right_annotation = rowAnnotation(
        id = anno_text(shared_df$gene_label, gp = gpar(fontsize = 12)),
        Pathway = pathway_factor_shared,
        col = list(Pathway = pw_palette_s),
        show_legend = FALSE,
        annotation_name_rot = 45,
        annotation_name_side = "bottom",
        annotation_name_gp = gpar(fontsize = 15, fontface = "bold"),
        show_annotation_name = c(id = FALSE, Pathway = TRUE)
      ),
      top_annotation = top_ha,
      column_title = wrapped_title,
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_rot = 45,
      column_names_gp = gpar(fontsize = 15),
      show_heatmap_legend = FALSE
    )
    
    png(paste0(res_dir, "/figures/07_heatmap_shared_deg_", cl, ".png"), 
        width = SHARED_WIDTH, height = 250 + nrow(mat)*25, res = 150)
    
    draw(ht, 
         annotation_legend_list = pd_common,
         heatmap_legend_side = "left", 
         annotation_legend_side = "left",
         padding = unit(c(5, 25, 5, 10), "mm"))
    dev.off()
    
    message("  [OK] Saved Shared Leading-Edge heatmap -> 07_heatmap_shared_deg_", cl, ".png")
  }
}

message("\n[OK] Heatmap generation complete.")
