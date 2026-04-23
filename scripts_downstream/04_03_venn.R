# [[scripts_downstream/04_03_venn.R]]
# Goal: Generate Directional DEG Overlap Venn Diagram and export shared gene lists.

library(DESeq2)
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_03_venn.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)
# Ensure output directories exist
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(res_dir, "tables/venn_lists"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("--- Generating Rigorous Venn Analysis for ", group_name, " ---")

# 2. Logic Setup
directions <- list(
  All  = function(lfc, padj) padj < 0.05 & abs(lfc) > 2,
  Up   = function(lfc, padj) padj < 0.05 & lfc > 2,
  Down = function(lfc, padj) padj < 0.05 & lfc < -2
)

priority_contrasts <- c(
  "Romi_6nM_vs_DMSO_Romi",
  "Kromastat_6nM_vs_DMSO_Kromastat"
)

# Premium colors: Up=Red (#D62728), Down=Blue (#1F77B4), All=Orange (#E69F00)
dir_colors <- c(All = "#E69F00", Up = "#D62728", Down = "#1F77B4")

# Storage for combined plot and tables
venn_plots <- list()
shared_targets_list <- list()

# 3. Process Each Direction
rigor_summary <- data.frame()

for (dir_name in names(directions)) {
  message("Processing Direction: ", dir_name)
  
  # Extract gene sets
  sig_list <- lapply(results_list[priority_contrasts], function(x) {
    if (is.null(x$df)) return(character(0))
    res_df <- x$df
    pass_filter <- directions[[dir_name]](res_df$log2FoldChange, res_df$padj)
    res_df$Geneid[pass_filter]
  })
  
  # Filter out empty sets
  sig_list <- sig_list[sapply(sig_list, length) > 0]
  
  if (length(sig_list) < 2) {
    message("  -- Skipping ", dir_name, " Venn: Fewer than 2 priority contrasts have significant genes.")
    next
  }
  
  # A. Generate Visualization
  clean_names <- names(sig_list) %>% 
    str_replace_all("_vs_", "\nvs\n") %>%
    str_replace_all("_", " ")
  names(sig_list) <- clean_names
  
  p_venn <- ggVennDiagram(sig_list) +
    scale_fill_gradient(low = "white", high = dir_colors[dir_name]) + 
    scale_x_continuous(expand = expansion(mult = .4)) + 
    scale_y_continuous(expand = expansion(mult = .1)) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(dir_name, " DEGs"), fill = "Count") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.margin = margin(10, 10, 10, 10))

  venn_plots[[dir_name]] <- p_venn
  
  # B. Export Shared Lists and Rigor Stats
  shared_ids <- Reduce(intersect, sig_list)
  union_ids <- Reduce(union, sig_list)
  total_genes <- nrow(results_list[[priority_contrasts[1]]]$df)
  
  n_shared <- length(shared_ids)
  n_romi <- length(sig_list[[1]])
  n_kroma <- length(sig_list[[2]])
  
  p_overlap <- phyper(n_shared - 1, n_romi, total_genes - n_romi, n_kroma, lower.tail = FALSE)
  jaccard <- n_shared / length(union_ids)
  expected_shared <- (n_romi * n_kroma) / total_genes
  rep_factor <- n_shared / expected_shared
  
  rigor_summary <- rbind(rigor_summary, data.frame(
    Direction = dir_name,
    Set1_Romi = n_romi,
    Set2_Kroma = n_kroma,
    Shared = n_shared,
    P_Overlap = p_overlap,
    Jaccard_Index = jaccard,
    Representation_Factor = rep_factor
  ))

  if (length(shared_ids) > 0) {
    romi_df <- results_list[[priority_contrasts[1]]]$df %>% 
      filter(Geneid %in% shared_ids) %>%
      select(Geneid, gene_label, log2FoldChange, padj)
    
    kroma_df <- results_list[[priority_contrasts[2]]]$df %>% 
      filter(Geneid %in% shared_ids) %>%
      select(Geneid, log2FoldChange, padj)
    
    shared_table <- inner_join(romi_df, kroma_df, by = "Geneid", suffix = c("_Romi", "_Kroma")) %>%
      mutate(
        direction = dir_name, 
        max_padj = pmax(padj_Romi, padj_Kroma, na.rm = TRUE),
        avg_abs_lfc = (abs(log2FoldChange_Romi) + abs(log2FoldChange_Kroma)) / 2
      ) %>%
      arrange(max_padj, desc(avg_abs_lfc))
    
    shared_targets_list[[dir_name]] <- shared_table
  }
}

# Save the master shared targets table
if (length(shared_targets_list) > 0) {
  final_shared_df <- bind_rows(shared_targets_list)
  write.csv(final_shared_df, 
            file = file.path(res_dir, "tables/04_03_shared_targets_stats.csv"),
            row.names = FALSE)
  message("[OK] Master shared targets table saved -> 04_03_shared_targets_stats.csv")
}

# Save the rigor summary table
write.csv(rigor_summary, 
          file = paste0(res_dir, "/tables/04_03_venn_rigor_stats.csv"),
          row.names = FALSE)
message("[OK] Rigor summary table saved -> 04_03_venn_rigor_stats.csv")

# 4. Save Individual and Combined Plots
if (length(venn_plots) > 0) {
  # Individual save
  for (name in names(venn_plots)) {
    ggsave(paste0(res_dir, "/figures/04_03_venn_", tolower(name), ".png"), venn_plots[[name]],
           width = 8, height = 7, bg = "white")
  }
  
  # Combined save (Vertical stack is better for long labels)
  p_combined <- wrap_plots(venn_plots, ncol = 1) + 
    plot_annotation(
      title = paste0("Directional Venn Comparisons: ", group_name),
      subtitle = "Criteria: padj < 0.05 & |log2FoldChange| > 2",
      theme = theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
  
  # Save with a portrait aspect ratio for the vertical stack
  ggsave(paste0(res_dir, "/figures/04_03_venn_combined.png"), p_combined,
         width = 10, height = 8 * length(venn_plots), dpi = 300, bg = "white")
  
  message("[OK] Combined Venn diagram saved -> 04_03_venn_combined.png")
}

message("[OK] Venn analysis complete.")
