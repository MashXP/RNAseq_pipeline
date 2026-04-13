# [[scripts_downstream/04_04_heatmap_pathway.R]]
# Goal: Generate Pathway-specific Heatmap for specific group.

library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(clusterProfiler)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_04_heatmap_pathway.R <Group> <Species>")
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
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# 1.5 ID Mapping (ENSG -> Symbol)
message("Mapping gene IDs...")
gene_ids <- rownames(dds)
gene_map <- bitr(gene_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org_db) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
row.names(gene_map) <- gene_map$ENSEMBL

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# Transform counts (VST)
vsd <- vst(dds, blind=FALSE)

# 1.6 Create clean Sample Names for display
colData(vsd)$display_name <- rownames(colData(vsd))
colnames(vsd) <- colData(vsd)$display_name

# Order VSD columns by condition
vsd <- vsd[, order(vsd$condition)]

# Shared dose colour palette
cond_levels   <- unique(as.character(colData(vsd)$condition))
dose_palette_shared <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(min(length(cond_levels), 8), "Set2"))(length(cond_levels)),
  cond_levels
)

# 5. Pathway-specific Heatmap
enrichment_doses <- names(enrichment_results_all)
target_dose <- if("5nM" %in% enrichment_doses) "5nM" else tail(enrichment_doses, 1)

# Robust selection fallback (ORA GO -> GSEA GO -> GSEA Hallmark)
res_obj <- enrichment_results_all[[target_dose]]
top_go_df <- if(!is.null(res_obj$ora_go) && nrow(as.data.frame(res_obj$ora_go)) > 0) {
  message("Using ORA GO for heatmap.")
  as.data.frame(res_obj$ora_go)
} else if (!is.null(res_obj$gsea_go) && nrow(as.data.frame(res_obj$gsea_go)) > 0) {
  message("ORA GO empty. Using GSEA GO for heatmap.")
  as.data.frame(res_obj$gsea_go)
} else if (!is.null(res_obj$gsea_hallmark) && nrow(as.data.frame(res_obj$gsea_hallmark)) > 0) {
  message("ORA/GSEA GO empty. Using GSEA Hallmark for heatmap.")
  as.data.frame(res_obj$gsea_hallmark)
} else {
  NULL
}

if (!is.null(top_go_df)) {
  top_pathways <- top_go_df %>%
    # Standard significance filtering for production
    filter(p.adjust < 0.05) %>%
    head(6)

  expanded_matrix_list <- list()
  annotation_row_list  <- list()

  for (i in seq_len(nrow(top_pathways))) {
    path_name  <- top_pathways$Description[i]
    safe_path  <- str_trunc(str_replace_all(path_name, "[^a-zA-Z0-9]", "_"), 30)

    # Dynamic column detection (ORA uses geneID, GSEA uses core_enrichment)
    gene_col <- if("geneID" %in% colnames(top_pathways)) "geneID" else "core_enrichment"

    if (!gene_col %in% colnames(top_pathways) || is.na(top_pathways[[gene_col]][i])) {
      message("Skipping pathway '", path_name, "': No gene column found.")
      next
    }

    genes_path <- str_split(top_pathways[[gene_col]][i], "/")[[1]]
    genes_path <- intersect(genes_path, rownames(vsd))
    message("Pathway '", path_name, "': ", length(genes_path), " genes found in VSD.")

    if (length(genes_path) > 0) {
      mat_slice <- assay(vsd)[genes_path, , drop = FALSE]
      rownames(mat_slice) <- paste0(genes_path, "__", safe_path)
      expanded_matrix_list[[path_name]] <- mat_slice
      annot_slice <- data.frame(Pathway = rep(path_name, length(genes_path)))
      rownames(annot_slice) <- rownames(mat_slice)
      annotation_row_list[[path_name]] <- annot_slice
    }
  }

  if (length(expanded_matrix_list) > 0) {
    mat_pathway  <- do.call(rbind, expanded_matrix_list)
    df_annot_row <- do.call(rbind, annotation_row_list)
    mat_pathway  <- mat_pathway - rowMeans(mat_pathway)
    mat_pathway  <- pmin(pmax(mat_pathway, -2), 2)

    orig_ensg       <- str_split_i(rownames(mat_pathway), "__", 1)
    symbols_display <- gene_map[orig_ensg, "SYMBOL"]
    labels_display  <- ifelse(is.na(symbols_display), orig_ensg, symbols_display)

    wrapped_pw     <- str_wrap(df_annot_row$Pathway, width = 20)
    pathway_factor <- factor(wrapped_pw, levels = unique(wrapped_pw))

    col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

    pathway_names <- levels(pathway_factor)
    n_pw          <- length(pathway_names)
    pw_palette    <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(min(n_pw, 8), "Set1"))(n_pw),
      pathway_names
    )

    conditions   <- colData(vsd)$condition
    top_ha <- HeatmapAnnotation(
      Dose = as.character(conditions),
      col  = list(Dose = dose_palette_shared),
      show_annotation_name = FALSE,
      show_legend          = FALSE
    )

    right_ha <- rowAnnotation(
      id      = anno_text(labels_display, gp = gpar(fontsize = 7)),
      Pathway = pathway_factor,
      gap     = unit(4, "mm"),
      col     = list(Pathway = pw_palette),
      show_annotation_name = FALSE,
      show_legend          = c(id = FALSE, Pathway = FALSE)
    )

    ht_pathway <- Heatmap(
      mat_pathway,
      name              = "log2FC",
      col               = col_fun,
      row_split         = pathway_factor,
      row_title_side    = "right",
      row_title_rot     = 0,
      row_title_gp      = gpar(fontsize = 8, fontface = "bold"),
      row_gap           = unit(0, "mm"),
      cluster_rows      = FALSE,
      cluster_row_slices = FALSE,
      cluster_columns   = FALSE,
      column_split      = factor(as.character(conditions), levels = unique(as.character(conditions))),
      column_title_gp   = gpar(fontsize = 9, fontface = "bold"),
      column_title_rot  = 0,
      column_gap        = unit(0, "mm"),
      top_annotation    = top_ha,
      right_annotation  = right_ha,
      show_row_names    = FALSE,
      show_column_names = TRUE,
      column_names_rot  = 45,
      column_names_gp   = gpar(fontsize = 10),
      use_raster        = FALSE,
      height            = unit(nrow(mat_pathway) * 6, "mm"),
      width             = unit(ncol(mat_pathway) * 10, "mm"),
      heatmap_legend_param = list(
        title          = "log2(FC)",
        title_position = "leftcenter-rot",
        at             = c(-2, -1, 0, 1, 2),
        labels         = c("-2", "-1", "0", "1", "2"),
        direction      = "vertical",
        legend_height  = unit(4, "cm")
      )
    )

    px_per_mm    <- 150 / 25.4
    mm_legend    <- 0
    mm_title     <- 22
    mm_col_title <- 12
    mm_dose_bar  <- 8
    mm_col_names <- 40
    mm_padding   <- 5
    mm_rows      <- nrow(mat_pathway) * 6
    canvas_h_px  <- round((mm_legend + mm_title + mm_col_title + mm_dose_bar + mm_col_names + mm_padding + mm_rows) * px_per_mm)

    png(paste0(res_dir, "/figures/04_heatmap_pathway_genes.png"),
        width = 1280, height = canvas_h_px, res = 150)
    ht_opt(TITLE_PADDING = unit(c(2, 8), "mm"))
    draw(ht_pathway,
         heatmap_legend_side    = "left",
         annotation_legend_side = "top",
         merge_legend           = FALSE,
         align_heatmap_legend   = "heatmap_top",
         column_title           = paste0(group_name, ": Genes Grouped by Pathway"),
         column_title_gp        = gpar(fontsize = 12, fontface = "bold"))
    decorate_annotation("Pathway", {
      grid.text("Pathway",
                x    = unit(0, "npc"),
                y    = unit(1, "npc") + unit(3, "mm"),
                just = c("left", "bottom"),
                gp   = gpar(fontsize = 9, fontface = "bold"))
    })
    dev.off()
    message("[OK] Pathway heatmap -> 04_heatmap_pathway_genes.png",
            "  (", nrow(mat_pathway), " rows x ", ncol(mat_pathway), " cols",
            "  |  ", nlevels(pathway_factor), " pathways",
            "  |  canvas: 1280x", canvas_h_px, "px)")
  } else {
    message("Skipping Pathway heatmap: No qualifying genes for top pathways.")
  }
} else {
  message("Skipping Pathway heatmap: No ORA GO enrichment results found.")
}
