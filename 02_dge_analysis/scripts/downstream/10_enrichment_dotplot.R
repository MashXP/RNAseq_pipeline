# [[02_dge_analysis/scripts/downstream/10_enrichment_dotplot.R]]
# Goal: Generate GSEA Hallmark Dotplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

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
if (length(args) == 0) {
  stop("Usage: Rscript 10_enrichment_dotplot.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- file.path(results_root, "02_dge_analysis", group_name)

# 1. Load data
load(file.path(rdata_dir, group_name, "02_deseq_results.RData"))
load(file.path(rdata_dir, group_name, "03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 8. GSEA Dotplots -- individual per dose + single faceted combined
message("Building GSEA dotplots...")

# Helper: clean Hallmark names for professional display (Mentor style)
clean_pathway_name <- function(pathway_name) {
  x <- sub("^HALLMARK_", "", pathway_name)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# Collect GSEA results (Hallmark) across all contrasts
gsea_combined <- map2_dfr(
  enrichment_results_all,
  names(enrichment_results_all),
  function(res, contrast_name) {
    gsea_res <- res$gsea_hallmark
    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      as.data.frame(gsea_res) %>%
        dplyr::select(Description, NES, pvalue, p.adjust, core_enrichment) %>%
        mutate(
          Contrast = contrast_name,
          Response = ifelse(NES > 0, "activated", "suppressed"),
          # Mirror logic: -1 for activated (Left), 1 for suppressed (Right)
          PlotX = ifelse(NES > 0, -1, 1),
          # Leading edge count (biological relevance)
          Count = vapply(str_split(core_enrichment, "/"), length, integer(1)),
          neg_log10_padj = -log10(p.adjust),
          Description = vapply(Description, clean_pathway_name, character(1))
        )
    }
  }
)

# Helper: build a single per-contrast dotplot (independent Y-axis, for individual saves)
make_gsea_plot_single <- function(df_dose, dose_label, top_n = 10) {
  hall_df <- df_dose %>%
    group_by(Response) %>%
    slice_min(order_by = p.adjust, n = top_n) %>%
    ungroup() %>%
    mutate(neg_log10_padj = -log10(p.adjust))

  ggplot(hall_df, aes(x = PlotX, y = reorder(Description, NES))) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_point(aes(size = Count, color = Response, alpha = neg_log10_padj), shape = 16) +
    scale_x_continuous(
      breaks = c(-1, 1),
      labels = c("Activated", "Suppressed"),
      limits = c(-1.5, 1.5)
    ) +
    scale_color_manual(
      values = c("activated" = "#D62728", "suppressed" = "#1F77B4"),
      name = "Direction"
    ) +
    scale_alpha_continuous(name = "-log10(padj)", range = c(0.4, 1)) +
    scale_size_continuous(name = "Leading Edge", range = c(3, 10)) +
    theme_bw(base_size = 16.5) +
    theme(
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 13.5, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16.5)
    ) +
    labs(title = paste0("GSEA Hallmark: ", dose_label), x = NULL, y = NULL) +
    guides(
      color = guide_legend(order = 1, title = "Direction", override.aes = list(size = 5)),
      alpha = guide_legend(order = 2, title = "-log10(padj)"),
      size  = guide_legend(order = 3, title = "Leading Edge")
    )
}

if (!is.null(gsea_combined) && nrow(gsea_combined) > 0) {

  # --- INDIVIDUAL PLOTS (per contrast, independent Y-axis) ---
  for (contrast in names(enrichment_results_all)) {
    df_dose <- gsea_combined %>% filter(Contrast == contrast)
    if (nrow(df_dose) > 0) {
      safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
      p <- make_gsea_plot_single(df_dose, contrast)
      ggsave(
        paste0(res_dir, "/figures/10_gsea_dotplot_", safe_name, ".png"),
        p, width = 8, height = 9, limitsize = FALSE
      )
    }
  }

  # --- COMBINED FACET PLOT (single ggplot, global Y-axis, mentor style) ---
  # Curated combined dotplot -- exclude statistically weaker and QC-only contrasts:
  #   - DMSO_Kromastat_vs_DMSO_Romidepsin: vehicle QC baseline, near-empty by design
  #   - Romidepsin_6nM_vs_DMSO_Romidepsin: global pooled across both cell lines, weaker signal
  #   - Kromastat_6nM_vs_DMSO_Kromastat: global pooled, same reason
  # These are preserved as individual saves above for reference.
  COMBINED_EXCLUDE <- c(
    "DMSO_Kromastat_vs_DMSO_Romidepsin",
    "Romidepsin_6nM_vs_DMSO_Romidepsin",
    "Kromastat_6nM_vs_DMSO_Kromastat"
  )
  gsea_combined_curated <- gsea_combined %>%
    filter(!Contrast %in% COMBINED_EXCLUDE)

  # Select top_n pathways per contrast by padj (not balanced, to avoid empty rows)
  TOP_N <- 12
  top_paths_per_contrast <- gsea_combined_curated %>%
    group_by(Contrast) %>%
    slice_min(order_by = p.adjust, n = TOP_N) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()

  # Build combined data: only pathways in the global selection
  dotplot_df <- gsea_combined_curated %>%
    filter(Description %in% top_paths_per_contrast)

  # Global Y-axis: order by NES across all contrasts (descending = activated on top)
  global_order <- dotplot_df %>%
    group_by(Description) %>%
    summarise(mean_NES = mean(NES), .groups = "drop") %>%
    arrange(mean_NES) %>%
    pull(Description)

  dotplot_df <- dotplot_df %>%
    mutate(Description = factor(Description, levels = global_order))

  combined_p <- ggplot(dotplot_df, aes(x = PlotX, y = Description)) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_point(aes(size = Count, color = Response, alpha = neg_log10_padj), shape = 16) +
    scale_x_continuous(
      breaks = c(-1, 1),
      labels = c("Activated", "Suppressed"),
      limits = c(-1.5, 1.5)
    ) +
    scale_color_manual(
      values = c("activated" = "#D62728", "suppressed" = "#1F77B4"),
      name = "Direction"
    ) +
    scale_alpha_continuous(name = "-log10(padj)", range = c(0.4, 1)) +
    scale_size_continuous(name = "Leading Edge", range = c(3, 10)) +
    facet_wrap(~ Contrast, nrow = 1) +
    theme_bw(base_size = 16.5) +
    theme(
      strip.text = element_text(face = "bold", size = 13.5),
      strip.background = element_rect(fill = "grey92"),
      axis.text.x = element_text(size = 13.5, face = "bold"),
      axis.text.y = element_text(size = 13.5, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 19.5)
    ) +
    labs(
      title = paste0("GSEA Hallmark Comparison: ", species_name),
      x = NULL, y = "Hallmark Pathway"
    ) +
    guides(
      color = guide_legend(order = 1, title = "Direction", override.aes = list(size = 5)),
      alpha = guide_legend(order = 2, title = "-log10(padj)"),
      size  = guide_legend(order = 3, title = "Leading Edge")
    )

  n_contrasts <- length(unique(dotplot_df$Contrast))
  n_paths     <- length(global_order)

  ggsave(
    paste0(res_dir, "/figures/10_gsea_dotplot_combined.png"),
    combined_p,
    width  = 5 * n_contrasts + 3,
    height = max(8, n_paths * 0.45),
    limitsize = FALSE
  )
  message("[OK] Combined Hallmark GSEA dotplot -> 10_gsea_dotplot_combined.png")

} else {
  message("Skipping GSEA dotplots: No Hallmark GSEA results found.")
}
