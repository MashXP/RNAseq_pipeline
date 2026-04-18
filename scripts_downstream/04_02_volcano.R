# [[scripts_downstream/04_02_volcano.R]]
# Goal: Generate Volcano plots using prior DGE results.
# FIXED: NA handling for padj and robust labeling for sparse results.

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_02_volcano.R <Group> [Species]")
}
group_name <- args[1]
species_name <- if (length(args) >= 2) args[2] else str_to_title(group_name)

res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load prior DGE data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("--- Generating Informative Volcano Plots for ", group_name, " ---")

# User's Color Palette
color_up   <- "#D62728"
color_down <- "#1F77B4"
color_ns   <- "grey70"

volcano_plots <- list()

for (contrast in names(results_list)) {
  message("Processing Volcano for: ", contrast)
  
  # results_list[[contrast]]$df contains pre-mapped gene symbols from 02_deseq2_dge.R
  plot_df <- results_list[[contrast]]$df
  
  # 2. Robust Significance Classification
  # Filter out rows with NA in log2FoldChange or padj for plotting stability
  plot_df <- plot_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  
  padj_cutoff <- 0.05
  lfc_cutoff <- 1.0
  
  plot_df <- plot_df %>%
    mutate(
      # Handle padj=0 and apply a visual cap at 500 for clarity
      neg_log10_padj = case_when(
        padj == 0 ~ 500,
        -log10(padj) > 500 ~ 500,
        TRUE ~ -log10(padj)
      ),
      sig_status = case_when(
        padj < padj_cutoff & log2FoldChange >= lfc_cutoff  ~ "up",
        padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down",
        TRUE ~ "not_significant"
      )
    )
  
  up_count <- sum(plot_df$sig_status == "up", na.rm = TRUE)
  down_count <- sum(plot_df$sig_status == "down", na.rm = TRUE)
  
  # 3. Mentor's Labeling Heuristic (Fixing tie-breaks for padj = 0)
  # Label top 5 genes per side, prioritized by padj then abs(LFC)
  label_per_side <- 5
  label_df <- plot_df %>%
    filter(sig_status != "not_significant") %>%
    group_by(sig_status) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = label_per_side) %>%
    ungroup()
  
  # 4. Construct Plot
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = sig_status)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("up" = color_up, "down" = color_down, "not_significant" = color_ns)) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40", alpha = 0.7) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40", alpha = 0.7) +
    coord_cartesian(clip = "off") + # Prevent labels from being cut off
    theme_bw(base_size = 14) +
    labs(
      title = str_replace_all(contrast, "_", " "),
      subtitle = paste0("Species: ", species_name),
      x = "log2foldChange",
      y = "-log10 adjusted p-value"
    ) +
    annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", down_count), 
             hjust = -0.2, vjust = 1.5, color = color_down, fontface = "bold", size = 5) +
    annotate("text", x = Inf, y = Inf, label = paste0("Up: ", up_count), 
             hjust = 1.2, vjust = 1.5, color = color_up, fontface = "bold", size = 5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm") # Extra margin for annotations
    )
  
  # Only add labels if there are significant genes
  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = gene_label),
      size = 3.5,
      max.overlaps = 20,
      box.padding = 0.6,
      point.padding = 0.4,
      segment.alpha = 0.5,
      show.legend = FALSE,
      force = 2
    )
  }
  
  volcano_plots[[contrast]] <- p
  
  safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
  ggsave(paste0(res_dir, "/figures/04_02_volcano_", safe_name, ".png"), p, width = 8, height = 7, dpi = 300, bg = "white")
}

# 5. Combined Plot
if (length(volcano_plots) > 1) {
  p_combined <- wrap_plots(volcano_plots, ncol = 2) + 
    plot_annotation(title = paste0("Volcano Plot Summary: ", group_name))
  
  ggsave(paste0(res_dir, "/figures/04_02_volcano_combined.png"), p_combined, 
         width = 16, height = 8 * ceiling(length(volcano_plots)/2), dpi = 300, bg = "white")
}

message("[OK] Volcano plots successfully completed for ", group_name)
