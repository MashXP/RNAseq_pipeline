# [[02_dge_analysis/scripts/downstream/09_enrichment_nes.R]]
# Goal: Generate Hallmark NES Barplots for specific group.

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
  stop("Usage: Rscript 09_enrichment_nes.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

# Helper: format drug names for professional display
format_drug_label <- function(name) {
  name %>%
    str_replace_all("_", " ")
}

res_dir <- file.path(results_root, "02_dge_analysis", group_name)

# 1. Load data
load(file.path(rdata_dir, group_name, "02_deseq_results.RData"))
load(file.path(rdata_dir, group_name, "03_enrichment_results.RData"))

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
                   ifelse(NES > 0, "Romidepsin higher", "Kromastat higher")
                 } else {
                   ifelse(NES > 0, "Activated", "Suppressed")
                 },
      Description = vapply(Description, clean_pathway_name, character(1))
    ) %>%
    group_by(Response) %>%
    slice_max(order_by = abs(NES), n = 10) %>%
    ungroup() %>%
    mutate(Response = factor(Response, levels = c("Activated", "Romidepsin higher", "Suppressed", "Kromastat higher")))

  # Dynamic height calculation based on number of pathways
  n_rows <- nrow(hall_df)
  plot_height <- max(6, n_rows * 0.38)

  p <- ggplot(hall_df, aes(x = NES, y = reorder(Description, NES), fill = Response)) +
    geom_col(width = 0.8) +
    facet_grid(Response ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c("Activated" = "#D62728", "Suppressed" = "#1F77B4",
                 "Romidepsin higher" = "#D62728", "Kromastat higher" = "#1F77B4"),
      name = if(is_direct) "Preference" else "Direction"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 13.5, face = "bold"),
      axis.title.x = element_text(size = 15),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 15, hjust = 0.5, margin = margin(b = 10)),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    labs(
      title = {
                cl_list <- unique(as.character(metadata$cell_line))
                cl_pattern <- paste0("^(", paste(cl_list, collapse="|"), ")_")
                if(grepl(cl_pattern, dose_label)) {
                  cl_prefix <- str_extract(dose_label, paste0("^(", paste(cl_list, collapse="|"), ")"))
                  rest <- sub(paste0("^", cl_prefix, "_"), "", dose_label)
                  paste0(cl_prefix, ": ", format_drug_label(rest))
                } else {
                  paste0("Global: ", format_drug_label(dose_label))
                }
              },
      subtitle = if(is_direct) "Positive NES favors Romidepsin;\nnegative NES favors Kromastat" else NULL,
      x = "Normalized Enrichment Score (NES)", y = "Hallmark Pathway"
    )
  
  return(list(plot = p, height = plot_height))
}

nes_plots <- list()
for (contrast in names(enrichment_results_all)) {
  hallmark_res <- enrichment_results_all[[contrast]]$gsea_hallmark
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    safe_name <- str_replace_all(contrast, "[^a-zA-Z0-9]", "_")
    res <- make_nes_plot(hallmark_res, contrast)
    p <- res$plot
    h <- res$height
    nes_plots[[contrast]] <- p

    # Save individual plot for ALL contrasts (use expanded name for consistency)
    safe_title <- str_replace_all(format_drug_label(contrast), "[^a-zA-Z0-9]", "_")
    ggsave(
      paste0(res_dir, "/figures/09_hallmark_nes_", safe_title, ".png"),
      p, width = 8, height = h, limitsize = FALSE
    )
  }
}

# Curated combined NES plot -- exclude statistically weaker and QC-only contrasts:
#   - DMSO_Kromastat_vs_DMSO_Romidepsin: vehicle QC baseline, near-empty by design
#   - Romidepsin_6nM_vs_DMSO_Romidepsin: global pooled across both cell lines, weaker signal
#   - Kromastat_6nM_vs_DMSO_Kromastat: global pooled, same reason
# These are preserved as individual saves above for reference.
COMBINED_EXCLUDE <- c(
  "DMSO_Kromastat_vs_DMSO_Romidepsin",
  "Romidepsin_6nM_vs_DMSO_Romidepsin",
  "Kromastat_6nM_vs_DMSO_Kromastat",
  "Romidepsin_6nM_vs_Kromastat_6nM"
)
nes_plots_curated <- nes_plots[!names(nes_plots) %in% COMBINED_EXCLUDE]

if (length(nes_plots_curated) > 1) {
  combined_nes <- wrap_plots(nes_plots_curated, nrow = 1) +
    plot_annotation(title = paste0("Hallmark NES: ", group_name),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  ggsave(paste0(res_dir, "/figures/09_hallmark_nes_combined.png"),
         combined_nes, width = 6 * length(nes_plots_curated), height = 10, limitsize = FALSE)
  message("[OK] Combined NES plot -> 09_hallmark_nes_combined.png")
} else {
  message("Skipping Combined NES plot: Fewer than 2 curated contrasts with NES data.")
}

# 8. Cell-Line Overlap NES Barplots (Dynamic detection)
message("Building Cell-Line Overlap NES plots...")

# Function to build side-by-side barplot for detected cell lines
make_nes_overlap_plot <- function(drug_label, contrasts, cl_list) {
  # Extract data frames
  df_list <- lapply(contrasts, function(ct) {
    # Extract cell line from the beginning of the contrast name
    cl_match <- str_extract(ct, paste0("^(", paste(cl_list, collapse="|"), ")"))
    if (is.na(cl_match)) return(NULL)
    
    hallmark_res <- enrichment_results_all[[ct]]$gsea_hallmark
    if (is.null(hallmark_res) || nrow(as.data.frame(hallmark_res)) == 0) return(NULL)
    
    as.data.frame(hallmark_res) %>%
      mutate(CellLine = cl_match, Description = vapply(Description, clean_pathway_name, character(1))) %>%
      select(Description, NES, p.adjust, CellLine)
  })

  # Remove NULLs and bind
  merged_df <- bind_rows(df_list)
  if (nrow(merged_df) == 0) return(NULL)
  
  # Identify top pathways for the union of all cell lines (Top 10 abs NES each)
  top_pathways <- merged_df %>%
    group_by(CellLine) %>%
    slice_max(order_by = abs(NES), n = 10) %>%
    pull(Description) %>%
    unique()

  plot_df <- merged_df %>%
    filter(Description %in% top_pathways) %>%
    complete(Description, CellLine, fill = list(NES = 0, p.adjust = 1)) %>%
    mutate(Description = factor(Description, levels = rev(sort(unique(Description)))))

  # Calculate dynamic height
  n_paths <- length(unique(plot_df$Description))
  plot_height <- max(7, n_paths * 0.45)

  p <- ggplot(plot_df, aes(x = NES, y = reorder(Description, NES), fill = CellLine)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    # Use a discrete color palette for cell lines
    scale_fill_brewer(palette = "Set1") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 13.5, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Pathway NES Comparison: ", format_drug_label(drug_label)),
      subtitle = paste0("Cell-Line Overlap (", species_name, ")"),
      x = "Normalized Enrichment Score (NES)",
      y = "Hallmark Pathway"
    )

  return(list(plot = p, height = plot_height))
}

# 8.1 Identify Drug-specific contrast groups automatically
all_contrasts <- names(enrichment_results_all)
cl_list <- unique(as.character(metadata$cell_line))

# Find treatments (e.g., Romidepsin_6nM_vs_DMSO_Romidepsin) by stripping cell line prefix
treatment_map <- list()
for (ct in all_contrasts) {
  for (cl in cl_list) {
    if (grepl(paste0("^", cl, "_"), ct)) {
      treatment <- sub(paste0("^", cl, "_"), "", ct)
      treatment_map[[treatment]] <- c(treatment_map[[treatment]], ct)
    }
  }
}


# 8.2 Generate Overlap Plots for each treatment found in multiple cell lines
for (tr in names(treatment_map)) {
  cts <- unique(treatment_map[[tr]])
  if (length(cts) >= 2) {
    res <- make_nes_overlap_plot(tr, cts, cl_list)
    if (!is.null(res)) {
      # Use expanded drug name for filename but keep it safe
      safe_tr <- str_replace_all(format_drug_label(tr), "[^a-zA-Z0-9]", "_")
      ggsave(
        paste0(res_dir, "/figures/09_hallmark_nes_overlap_", safe_tr, ".png"),
        res$plot, width = 13, height = res$height, bg = "white"
      )
      message("[OK] Overlap NES plot -> 09_hallmark_nes_overlap_", safe_tr, ".png")
    }
  }
}
