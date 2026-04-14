# [[scripts_downstream/04_06_enrichment_nes.R]]
# Goal: Generate Hallmark NES Barplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_06_enrichment_nes.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 7. Hallmark NES Barplots -- individual per dose + stitched combined
message("Building Hallmark NES plots...")

# Helper: build a single NES barplot for one dose
make_nes_plot <- function(hallmark_res, dose_label) {
  hall_df <- as.data.frame(hallmark_res) %>%
    mutate(Response = ifelse(NES > 0, "Activated", "Suppressed")) %>%
    group_by(Response) %>%
    slice_max(order_by = abs(NES), n = 10) %>%
    ungroup() %>%
    arrange(Response, ifelse(Response == "Activated", NES, -NES)) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  ggplot(hall_df, aes(x = NES, y = Description, fill = Response)) +
    geom_col(width = 0.75) +
    scale_fill_manual(
      values = c("Activated" = "#D62728", "Suppressed" = "#1F77B4"),
      name = "Direction"
    ) +
    facet_grid(Response ~ ., scales = "free_y", space = "free_y") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
    theme_bw(base_size = 11) +
    theme(
      strip.text.y = element_text(angle = -90, face = "bold"),
      strip.background = element_rect(fill = "grey88"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0("Contrast: ", dose_label),
      x = "Normalized Enrichment Score (NES)", y = NULL
    )
}

nes_plots <- list()
for (contrast in names(enrichment_results_all)) {
  hallmark_res <- enrichment_results_all[[contrast]]$gsea_hallmark
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
    p <- make_nes_plot(hallmark_res, contrast)
    nes_plots[[contrast]] <- p
    
    # Save individual plot
    n_rows <- nrow(as.data.frame(hallmark_res))
    ggsave(
      paste0(res_dir, "/figures/04_hallmark_nes_", safe_name, ".png"),
      p, width = 10, height = max(6, n_rows * 0.38)
    )
  }
}

# Save stitched combined NES plot
if (length(nes_plots) > 1) {
  combined_nes <- wrap_plots(nes_plots, nrow = 1) +
    plot_annotation(title = paste0("Hallmark NES: ", group_name),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  ggsave(paste0(res_dir, "/figures/04_hallmark_nes_combined.png"),
         combined_nes, width = 10 * length(nes_plots), height = 12)
  message("[OK] Combined NES plot -> 04_hallmark_nes_combined.png")
} else {
  message("Skipping Combined NES plot: Fewer than 2 doses with NES data.")
}
