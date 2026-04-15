# [[test_downstream/test_04_08_ora_dotplot.R]]
# Goal: Generate ORA Dotplots for specific group.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript test_04_08_ora_dotplot.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

# --- DEVIATION: routing results to test directory
res_dir <- paste0("../results_test/", group_name)

# 1. Load data
# --- DEVIATION: routing input from test directory
load(paste0("./.RData_test/", group_name, "/02_deseq_results.RData"))
load(paste0("./.RData_test/", group_name, "/03_enrichment_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 9. ORA Dotplots -- individual per dose + stitched combined
message("Building ORA dotplots...")

# Collect ORA results across all contrasts
ora_combined <- map2_dfr(
  enrichment_results_all,
  names(enrichment_results_all),
  function(res, contrast_name) {
    ora_res <- res$go_ora
    if (!is.null(ora_res) && nrow(as.data.frame(ora_res)) > 0) {
      as.data.frame(ora_res) %>%
        dplyr::select(Description, pvalue, p.adjust, Count, GeneRatio) %>%
        mutate(
          Contrast = contrast_name,
          # Parse GeneRatio (e.g., "10/100" -> 0.1)
          Ratio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        )
    }
  }
)

# Helper: build a single ORA dotplot
make_ora_plot <- function(df_dose, dose_label, top_paths_wrapped) {
  df_dose <- df_dose %>%
    filter(Description %in% top_paths_wrapped) %>%
    mutate(
      Description = factor(Description, levels = rev(top_paths_wrapped))
    )
  
  ggplot(df_dose, aes(x = Ratio, y = Description,
                      color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#E41A1C", high = "#377EB8",
      name = "padj"
    ) +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    ) +
    labs(title = paste0("ORA Enrichment: ", dose_label),
         x = "Gene Ratio", y = NULL)
}

if (!is.null(ora_combined) && nrow(ora_combined) > 0) {
  # Select Top 15 pathways per contrast to ensure diversity
  top_paths <- ora_combined %>%
    group_by(Contrast) %>%
    slice_min(order_by = p.adjust, n = 15) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # Limit to 40 for readability
  if (length(top_paths) > 40) {
     top_paths <- ora_combined %>%
        filter(Description %in% top_paths) %>%
        group_by(Description) %>%
        summarise(min_p = min(p.adjust)) %>%
        arrange(min_p) %>%
        head(40) %>%
        pull(Description)
  }
  
  top_paths_wrapped <- str_wrap(top_paths, width = 40)
  
  ora_combined <- ora_combined %>%
    mutate(Description = str_wrap(Description, width = 40))
  
  ora_plots <- list()
  for (contrast in names(enrichment_results_all)) {
    df_dose <- ora_combined %>% filter(Contrast == contrast)
    if (nrow(df_dose) > 0) {
      safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
      p <- make_ora_plot(df_dose, contrast, top_paths_wrapped)
      ora_plots[[contrast]] <- p
      
      ggsave(
        paste0(res_dir, "/figures/04_ora_dotplot_", safe_name, ".png"),
        p, width = 8, height = max(6, length(top_paths) * 0.4)
      )
    }
  }
  
  # Stitched plot
  if (length(ora_plots) > 1) {
    combined_ora <- wrap_plots(ora_plots, nrow = 1) +
      plot_annotation(title = "GO ORA: Pathway Comparison")
    
    ggsave(
      paste0(res_dir, "/figures/04_ora_dotplot_combined.png"),
      combined_ora,
      width = 8 * length(ora_plots),
      height = max(8, length(top_paths) * 0.4)
    )
    message("[OK] Combined ORA dotplot saved.")
  }
} else {
  message("Skipping ORA dotplots: No significant results found.")
}
