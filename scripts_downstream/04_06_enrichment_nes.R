# [[scripts_downstream/04_06_enrichment_nes.R]]
# Goal: Generate Hallmark NES Barplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_06_enrichment_nes.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 7. Hallmark NES Barplots -- individual per dose + stitched combined
message("Building Hallmark NES plots...")

# Helper: clean Hallmark names for professional display (Mentor style)
clean_pathway_name <- function(pathway_name) {
  x <- sub("^HALLMARK_", "", pathway_name)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# Helper: build a single NES barplot for one dose
make_nes_plot <- function(hallmark_res, dose_label) {
  is_direct <- grepl("vs_Kromastat", dose_label)
  
  hall_df <- as.data.frame(hallmark_res) %>%
    mutate(
      Response = if(is_direct) {
                   ifelse(NES > 0, "Romi higher", "Krom higher")
                 } else {
                   ifelse(NES > 0, "Activated", "Suppressed")
                 },
      Description = vapply(Description, clean_pathway_name, character(1))
    ) %>%
    group_by(Response) %>%
    slice_max(order_by = abs(NES), n = 10) %>%
    ungroup()
  
  ggplot(hall_df, aes(x = NES, y = reorder(Description, NES), fill = Response)) +
    geom_col(width = 0.8) +
    scale_fill_manual(
      values = c("Activated" = "#D62728", "Suppressed" = "#1F77B4",
                 "Romi higher" = "#D62728", "Krom higher" = "#1F77B4"),
      name = if(is_direct) "Preference" else "Direction"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 9, face = "bold"),
      axis.title.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10))
    ) +
    labs(
      title = if(grepl("^(H9|SUPM2)_", dose_label)) {
                cl_prefix <- str_extract(dose_label, "^(H9|SUPM2)")
                rest <- sub("^(H9|SUPM2)_", "", dose_label)
                paste0(cl_prefix, ": ", str_replace_all(rest, "_", " "))
              } else {
                paste0("Global: ", str_replace_all(dose_label, "_", " "))
              },
      subtitle = if(is_direct) str_wrap("Positive NES favors Romidepsin; negative NES favors Kromastat", width = 50) else NULL,
      x = "Normalized Enrichment Score (NES)", y = "Hallmark Pathway"
    )
}

nes_plots <- list()
for (contrast in names(enrichment_results_all)) {
  hallmark_res <- enrichment_results_all[[contrast]]$gsea_hallmark
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
    p <- make_nes_plot(hallmark_res, contrast)
    nes_plots[[contrast]] <- p

    # Save individual plot for ALL contrasts (including QC and global pooled)
    ggsave(
      paste0(res_dir, "/figures/04_06_hallmark_nes_", safe_name, ".png"),
      p, width = 8, height = 7, limitsize = FALSE
    )
  }
}

# Curated combined NES plot -- exclude statistically weaker and QC-only contrasts:
#   - DMSO_Kromastat_vs_DMSO_Romi: vehicle QC baseline, near-empty by design
#   - Romi_6nM_vs_DMSO_Romi: global pooled across both cell lines, weaker signal
#   - Kromastat_6nM_vs_DMSO_Kromastat: global pooled, same reason
# These are preserved as individual saves above for reference.
COMBINED_EXCLUDE <- c(
  "DMSO_Kromastat_vs_DMSO_Romi",
  "Romi_6nM_vs_DMSO_Romi",
  "Kromastat_6nM_vs_DMSO_Kromastat",
  "Romi_6nM_vs_Kromastat_6nM"
)
nes_plots_curated <- nes_plots[!names(nes_plots) %in% COMBINED_EXCLUDE]

if (length(nes_plots_curated) > 1) {
  combined_nes <- wrap_plots(nes_plots_curated, nrow = 1) +
    plot_annotation(title = paste0("Hallmark NES: ", group_name),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  ggsave(paste0(res_dir, "/figures/04_06_hallmark_nes_combined.png"),
         combined_nes, width = 6 * length(nes_plots_curated), height = 10, limitsize = FALSE)
  message("[OK] Combined NES plot -> 04_06_hallmark_nes_combined.png")
} else {
  message("Skipping Combined NES plot: Fewer than 2 curated contrasts with NES data.")
}
