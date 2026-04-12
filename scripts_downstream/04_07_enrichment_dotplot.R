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

# Collect GSEA results across all doses
gsea_combined <- map2_dfr(
  enrichment_results_all,
  names(enrichment_results_all),
  function(res, dose_name) {
    gsea_res <- res$gsea_go
    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      as.data.frame(gsea_res) %>%
        dplyr::select(Description, NES, pvalue, p.adjust, setSize) %>%
        mutate(
          Dose = dose_name,
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
  # Get globally consistent top pathways (union across all doses)
  top_paths <- gsea_combined %>%
    group_by(Description) %>%
    summarise(min_pval = min(pvalue)) %>%
    arrange(min_pval) %>%
    head(25) %>%
    pull(Description)
  top_paths_wrapped <- str_wrap(top_paths, width = 35)
  
  gsea_combined <- gsea_combined %>%
    mutate(Description = str_wrap(Description, width = 35))
  
  gsea_plots <- list()
  for (dose in names(enrichment_results_all)) {
    df_dose <- gsea_combined %>% filter(Dose == dose)
    if (nrow(df_dose) > 0) {
      safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
      p <- make_gsea_plot(df_dose, dose, top_paths_wrapped)
      gsea_plots[[dose]] <- p
      
      # Save individual plot
      ggsave(
        paste0(res_dir, "/figures/04_gsea_dotplot_", safe_name, ".png"),
        p, width = 7, height = max(6, length(top_paths) * 0.45)
      )
    }
  }
  
  # Save stitched combined GSEA plot
  if (length(gsea_plots) > 1) {
    combined_gsea <- wrap_plots(gsea_plots, nrow = 1, guides = "collect") +
      plot_annotation(
        title = "GSEA: Multi-Dose Comparison",
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
      )
    ggsave(
      paste0(res_dir, "/figures/04_gsea_dotplot_combined.png"),
      combined_gsea,
      width = 7 * length(gsea_plots) + 2,
      height = max(8, length(top_paths) * 0.45)
    )
    message("[OK] Combined GSEA dotplot -> 04_gsea_dotplot_combined.png")
  }
} else {
  message("Skipping GSEA dotplots: No GSEA GO results found.")
}
