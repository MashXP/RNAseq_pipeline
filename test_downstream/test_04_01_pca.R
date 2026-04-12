# [[test_downstream/test_04_01_pca.R]]
# Goal: Generate PCA plot for specific group. (TEST)

library(DESeq2)
library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript test_04_01_pca.R <Group> <Species>")
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

# --- DEVIATION: routing results to test directory
res_dir <- paste0("../results_test/", group_name)

# 1. Load data
# --- DEVIATION: routing input from test directory
load(paste0("./.RData_test/", group_name, "/02_deseq_results.RData"))

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

p_pca <- plotPCA(vsd, intgroup = "condition") +
  scale_color_manual(values = dose_palette_shared) +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey93", linewidth = 0.2)
  ) +
  ggtitle(paste0("PCA Plot: ", group_name, " Drug Response"))

ggsave(paste0(res_dir, "/figures/04_pca_plot.png"), p_pca, width = 8, height = 6, bg = "white")
message("[OK] PCA plot saved -> 04_pca_plot.png  (8x6 in, ", length(cond_levels), " conditions)")
