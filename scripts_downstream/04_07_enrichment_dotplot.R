# [[scripts_downstream/04_07_enrichment_dotplot.R]]
# Goal: Generate GSEA Dotplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_07_enrichment_dotplot.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData/", group_name, "/03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 8. GSEA Dotplots -- individual per dose + stitched combined
message("Building GSEA dotplots...")

# Collect GSEA results across all contrasts
gsea_combined <- map2_dfr(
  enrichment_results_all,
  names(enrichment_results_all),
  function(res, contrast_name) {
    gsea_res <- res$gsea_go
    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      as.data.frame(gsea_res) %>%
        dplyr::select(Description, NES, pvalue, p.adjust, setSize) %>%
        mutate(
          Contrast = contrast_name,
          Response = ifelse(NES > 0, "activated", "suppressed")
        )
    }
  }
)

# Helper: build a single GSEA dotplot for one dose
make_gsea_plot <- function(df_dose, dose_label, top_paths_wrapped) {
  df_dose <- df_dose %>%
    filter(Description %in% top_paths_wrapped) %>%
    mutate(
      Description = factor(Description, levels = rev(top_paths_wrapped)),
      Response = factor(Response, levels = c("activated", "suppressed"))
    )
  
  ggplot(df_dose, aes(x = Response, y = Description,
                      color = pvalue, size = setSize)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#D62728", high = "#1F77B4",
      name = "pvalue"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 11)) +
    geom_vline(xintercept = 1.5, color = "grey40", linewidth = 0.6) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9, hjust = 1),
      strip.background = element_rect(fill = "#2F6F8E"),
      strip.text = element_text(color = "white", face = "bold", size = 11),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    ) +
    labs(title = paste0("Enrichment GSEA: ", dose_label),
         x = NULL, y = NULL) +
    guides(
      color = guide_colorbar(order = 1, title = "pvalue"),
      size  = guide_legend(order = 2, title = "Count")
    )
}

if (!is.null(gsea_combined) && nrow(gsea_combined) > 0) {
  # Balance Selection: Take Top 10 Activated and Top 10 Suppressed per Contrast
  # To ensure weaker signals (like SUPM2 activation) are not drowned out by global significance
  top_paths_balanced <- gsea_combined %>%
    group_by(Contrast, Response) %>%
    slice_min(order_by = pvalue, n = 10) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # Limit to a reasonable number for the Y-axis if the set is too large
  if (length(top_paths_balanced) > 40) {
     # If too many, prioritize by global significance but keep at least some from each contrast
     top_paths <- gsea_combined %>%
        filter(Description %in% top_paths_balanced) %>%
        group_by(Description) %>%
        summarise(min_p = min(pvalue)) %>%
        arrange(min_p) %>%
        head(40) %>%
        pull(Description)
  } else {
     top_paths <- top_paths_balanced
  }
  top_paths_wrapped <- str_wrap(top_paths, width = 35)
  
  gsea_combined <- gsea_combined %>%
    mutate(Description = str_wrap(Description, width = 35))
  
  gsea_plots <- list()
  for (contrast in names(enrichment_results_all)) {
    df_dose <- gsea_combined %>% filter(Contrast == contrast)
    if (nrow(df_dose) > 0) {
      safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
      p <- make_gsea_plot(df_dose, contrast, top_paths_wrapped)
      gsea_plots[[contrast]] <- p
      
      # Save individual plot
      ggsave(
        paste0(res_dir, "/figures/04_gsea_dotplot_", safe_name, ".png"),
        p, width = 7, height = max(6, length(top_paths) * 0.45)
      )
    }
  }
  
  # Save stitched combined GSEA plot
  if (length(gsea_plots) > 1) {
    combined_gsea <- wrap_plots(gsea_plots, nrow = 1) +
      plot_annotation(
        title = "GSEA: Multi-Dose Comparison",
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
      )
    ggsave(
      paste0(res_dir, "/figures/04_gsea_dotplot_combined.png"),
      combined_gsea,
      width = 8 * length(gsea_plots) + 1,
      height = max(8, length(top_paths) * 0.45)
    )
    message("[OK] Combined GSEA dotplot -> 04_gsea_dotplot_combined.png")
  }
} else {
  message("Skipping GSEA dotplots: No GSEA GO results found.")
}
