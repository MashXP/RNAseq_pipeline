# [[scripts_downstream/04_01_pca.R]]
# Goal: Generate PCA plots at the species level (combined) and cell-line level.
# Aligned with mentor's multifactorial design principles.

library(DESeq2)
library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_01_pca.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
# We load the results from step 02 as it contains the dds object
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("--- Generating Multifactorial PCA for ", group_name, " ---")

# 2. Transform counts (VST) 
# blind = TRUE is standard for QC PCA plots per mentor's script
vsd <- vst(dds, blind = TRUE)

# 3. Combined PCA (Species Level)
pca_data <- plotPCA(vsd, intgroup = c("cell_line", "condition"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Shared dose colour palette (Intuitive Distinct Palette)
dose_palette <- c(
  "DMSO_Romidepsin"      = "grey85",
  "Romidepsin_6nM"       = "#E41A1C", # Brighter Red
  "DMSO_Kromastat" = "grey70",
  "Kromastat_6nM"  = "#377EB8"  # Brighter Blue
)

# Combined Plot logic
p_combined <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = cell_line)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = dose_palette) +
  theme_bw(base_size = 21) +
  labs(
    title = paste0("Species-Level PCA: ", species_name),
    subtitle = "Color = Treatment | Shape = Cell Line",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(res_dir, "figures/04_01_pca_combined.png"), p_combined, width = 11, height = 7, bg = "white")
message("[OK] Saved Combined PCA: 04_01_pca_combined.png")

# 4. Individual Cell-Line PCAs
cell_lines <- unique(as.character(colData(vsd)$cell_line))

for (cl in cell_lines) {
  message("Processing sub-PCA for cell line: ", cl)
  
  vsd_cell <- vsd[, colData(vsd)$cell_line == cl]
  
  # plotPCA requires at least 2 samples, usually 3+ for meaningful plot
  if (ncol(vsd_cell) >= 2) {
    pca_cell_data <- plotPCA(vsd_cell, intgroup = "condition", returnData = TRUE)
    percent_var_cell <- round(100 * attr(pca_cell_data, "percentVar"))
    
    p_cell <- ggplot(pca_cell_data, aes(PC1, PC2, color = condition)) +
      geom_point(size = 5) +
      scale_color_manual(values = dose_palette) +
      theme_bw(base_size = 21) +
      labs(
        title = paste0("Cell-Line PCA: ", cl),
        subtitle = paste0("Species: ", species_name),
        x = paste0("PC1: ", percent_var_cell[1], "% variance"),
        y = paste0("PC2: ", percent_var_cell[2], "% variance")
      ) +
      theme(panel.grid.minor = element_blank())
    
    filename_cl <- paste0("04_01_pca_", cl, ".png")
    ggsave(file.path(res_dir, "figures", filename_cl), p_cell, width = 9, height = 6, bg = "white")
    message("  [OK] Saved Sub-PCA: ", filename_cl)
  }
}

message("PCA Step Complete for ", group_name)
