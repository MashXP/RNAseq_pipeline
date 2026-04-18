# [[scripts_downstream/04_05_heatmap_variable.R]]
# Goal: Generate Top 50 Most Variable Genes Heatmap for specific group.

library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(clusterProfiler)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_05_heatmap_variable.R <Group> [Species]")
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

# Transform counts (VST)
vsd <- vst(dds, blind=FALSE)

# 1.6 Create clean Sample Names for display
# --- DEVIATION: Standardize naming with pathway heatmap logic (04_04)
display_names <- rownames(colData(vsd)) %>%
  str_remove_all("human_|dog_") %>% 
  str_replace_all("DMSO", "D") %>%
  str_replace_all("Romidepsin|Romi", "R") %>%
  str_replace_all("Kromastat|Kroma", "K") %>%
  str_replace_all("6nM", "") %>%
  str_replace_all("(?i)replicate|rep", "") %>%
  str_replace_all("_+", "_") %>%
  str_remove("^_|_$")

colData(vsd)$display_name <- display_names
colnames(vsd) <- display_names

# Order VSD columns by cell_line then condition
vsd <- vsd[, order(vsd$cell_line, vsd$condition)]

# Color palettes (Intuitive Distinct Palette)
conditions <- colData(vsd)$condition
cell_lines <- colData(vsd)$cell_line

dose_palette <- c(
  "DMSO_Romi"      = "grey85",
  "Romi_6nM"       = "#E41A1C", # Brighter Red
  "DMSO_Kromastat" = "grey70",
  "Kromastat_6nM"  = "#377EB8"  # Brighter Blue
)

cl_palette <- setNames(
  RColorBrewer::brewer.pal(max(3, length(unique(cell_lines))), "Accent")[seq_along(unique(cell_lines))], 
  unique(cell_lines)
)

# 6. Default Heatmap (Top 50 Variable) with ComplexHeatmap
top_genes_idx <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat_var       <- assay(vsd)[top_genes_idx, ]
# Intensify colors: Full Z-score scaling (center and scale by SD)
mat_var       <- t(scale(t(mat_var))) 
mat_var       <- pmin(pmax(mat_var, -2), 2)

symbols_var       <- gene_map[rownames(mat_var), "SYMBOL"]
rownames(mat_var) <- ifelse(is.na(symbols_var), rownames(mat_var), symbols_var)

col_fun_var <- colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))

conditions_var <- as.character(colData(vsd)$condition)
cell_lines_var <- as.character(colData(vsd)$cell_line)

top_ha_var <- HeatmapAnnotation(
  Dose     = conditions_var,
  CellLine = cell_lines_var,
  col      = list(
    Dose     = dose_palette,
    CellLine = cl_palette
  ),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
  show_legend          = TRUE
)

ht_var <- Heatmap(
  mat_var,
  name              = "Z-score",
  col               = col_fun_var,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  column_split      = factor(conditions_var, levels = unique(conditions_var)),
  column_title_gp   = gpar(fontsize = 9, fontface = "bold"),
  column_title_rot  = 0,
  column_gap        = unit(0, "mm"),
  top_annotation    = top_ha_var,
  show_row_names    = TRUE,
  row_names_gp      = gpar(fontsize = 8),
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 10),
  use_raster        = FALSE,
  height            = unit(nrow(mat_var) * 6, "mm"),
  width             = unit(ncol(mat_var) * 10, "mm"),
  heatmap_legend_param = list(
    title          = "Z-score",
    title_position = "leftcenter-rot",
    at             = c(-2, -1, 0, 1, 2),
    labels         = c("-2", "-1", "0", "1", "2"),
    direction      = "vertical",
    legend_height  = unit(4, "cm")
  )
)

px_per_mm_v    <- 150 / 25.4
mm_title_v     <- 22
mm_col_title_v <- 12
mm_dose_bar_v  <- 12  # Increased for dual annotation
mm_col_names_v <- 40
mm_padding_v   <- 10
mm_rows_v      <- nrow(mat_var) * 6
canvas_h_var   <- round((mm_title_v + mm_col_title_v + mm_dose_bar_v + mm_col_names_v + mm_padding_v + mm_rows_v) * px_per_mm_v)

png(paste0(res_dir, "/figures/04_05_heatmap_top_variable.png"),
    width = 1800, height = canvas_h_var, res = 150)
ht_opt(TITLE_PADDING = unit(c(2, 8), "mm"))
draw(ht_var,
     heatmap_legend_side    = "left",
     annotation_legend_side = "top",
     merge_legend           = FALSE,
     align_heatmap_legend   = "heatmap_top",
     column_title           = "Top 50 Most Variable Genes",
     column_title_gp        = gpar(fontsize = 12, fontface = "bold"))
dev.off()
message("[OK] Top-var heatmap  -> 04_05_heatmap_top_variable.png",
        "  (", nrow(mat_var), " rows x ", ncol(mat_var), " cols",
        "  |  canvas: 1280x", canvas_h_var, "px)")
