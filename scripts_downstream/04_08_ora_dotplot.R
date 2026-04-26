# [[scripts_downstream/04_08_ora_dotplot.R]]
# Goal: Generate ORA Dotplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_08_ora_dotplot.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 9. ORA Dotplots -- individual per dose + faceted combined
message("Building ORA dotplots...")

# Collect ORA results across all contrasts
ora_combined <- map2_dfr(
  enrichment_results_all,
  names(enrichment_results_all),
  function(res, contrast_name) {
    ora_res <- res$ora_go
    if (!is.null(ora_res) && nrow(as.data.frame(ora_res)) > 0) {
      as.data.frame(ora_res) %>%
        dplyr::select(Description, pvalue, p.adjust, Count, GeneRatio) %>%
        mutate(
          Contrast = contrast_name,
          # Parse GeneRatio (e.g., "10/100" -> 0.1)
          Ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
          # Wrap long descriptions once, consistently
          Description = str_wrap(Description, width = 35)
        )
    }
  }
)

# Helper: build a single ORA dotplot (local top_n, independent Y-axis)
make_ora_plot_single <- function(df_dose, dose_label, top_n = 15) {
  hall_df <- df_dose %>%
    slice_min(order_by = p.adjust, n = top_n) %>%
    arrange(p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))

  ggplot(hall_df, aes(x = Ratio, y = Description, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#D62728", high = "#1F77B4",
      name = "padj"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x = element_text(size = 13.5),
      axis.text.y = element_text(size = 13.5, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16.5)
    ) +
    labs(title = paste0("ORA GO: ", dose_label),
         x = "Gene Ratio", y = NULL) +
    guides(
      color = guide_colorbar(order = 1, title = "padj"),
      size  = guide_legend(order = 2, title = "Count")
    )
}

if (!is.null(ora_combined) && nrow(ora_combined) > 0) {

  # --- INDIVIDUAL PLOTS (per contrast, independent Y-axis) ---
  for (contrast in names(enrichment_results_all)) {
    df_dose <- ora_combined %>% filter(Contrast == contrast)
    if (nrow(df_dose) > 0) {
      safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
      p <- make_ora_plot_single(df_dose, contrast)
      ggsave(
        paste0(res_dir, "/figures/04_08_ora_dotplot_", safe_name, ".png"),
        p, width = 8, height = 9, limitsize = FALSE
      )
    }
  }

  # --- COMBINED FACET PLOT (single ggplot, global Y-axis) ---
  # Curated combined ORA plot -- exclude statistically weaker and QC-only contrasts:
  #   - DMSO_Kromastat_vs_DMSO_Romidepsin: vehicle QC baseline, near-empty by design
  #   - Romidepsin_6nM_vs_DMSO_Romidepsin: global pooled across both cell lines, weaker signal
  #   - Kromastat_6nM_vs_DMSO_Kromastat: global pooled, same reason
  # These are preserved as individual saves above for reference.
  COMBINED_EXCLUDE <- c(
    "DMSO_Kromastat_vs_DMSO_Romidepsin",
    "Romidepsin_6nM_vs_DMSO_Romidepsin",
    "Kromastat_6nM_vs_DMSO_Kromastat"
  )
  ora_combined_curated <- ora_combined %>%
    filter(!Contrast %in% COMBINED_EXCLUDE)

  # Select top_n pathways per curated contrast by padj
  TOP_N <- 12
  top_paths <- ora_combined_curated %>%
    group_by(Contrast) %>%
    slice_min(order_by = p.adjust, n = TOP_N) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()

  # Build combined data: only selected pathways
  dotplot_df <- ora_combined_curated %>%
    filter(Description %in% top_paths)

  # Global Y-axis: order by mean padj across contrasts (most significant on top)
  global_order <- dotplot_df %>%
    group_by(Description) %>%
    summarise(mean_padj = mean(p.adjust), .groups = "drop") %>%
    arrange(desc(mean_padj)) %>%
    pull(Description)

  dotplot_df <- dotplot_df %>%
    mutate(Description = factor(Description, levels = global_order))

  n_contrasts <- length(unique(dotplot_df$Contrast))
  n_paths     <- length(global_order)

  combined_p <- ggplot(dotplot_df, aes(x = Ratio, y = Description,
                                        color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#D62728", high = "#1F77B4",
      name = "padj"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    facet_wrap(~ Contrast, nrow = 1) +
    theme_bw(base_size = 16.5) +
    theme(
      strip.text = element_text(face = "bold", size = 13.5),
      strip.background = element_rect(fill = "grey92"),
      axis.text.x = element_text(size = 13.5),
      axis.text.y = element_text(size = 13.5, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 19.5)
    ) +
    labs(
      title = paste0("GO ORA Comparison: ", species_name),
      x = "Gene Ratio", y = "GO Term"
    ) +
    guides(
      color = guide_colorbar(order = 1, title = "padj"),
      size  = guide_legend(order = 2, title = "Count")
    )

  ggsave(
    paste0(res_dir, "/figures/04_08_ora_dotplot_combined.png"),
    combined_p,
    width  = 5 * n_contrasts + 3,
    height = max(8, n_paths * 0.45),
    limitsize = FALSE
  )
  message("[OK] Combined ORA dotplot saved -> 04_08_ora_dotplot_combined.png")

} else {
  message("Skipping ORA dotplots: No significant results found.")
}
