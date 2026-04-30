# [[02_dge_analysis/scripts/downstream/05_volcano.R]]
# Goal: Generate Volcano plots using prior DGE results.
# FIXED: NA handling for padj and robust labeling for sparse results.

library(DESeq2)
library(ggplot2)
library(ggrepel)
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
  stop("Usage: Rscript 05_volcano.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- file.path(results_root, "02_dge_analysis", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load prior DGE data
load(file.path(rdata_dir, group_name, "02_deseq_results.RData"))

message("--- Generating Informative Volcano Plots for ", group_name, " ---")

# User's Color Palette
color_up   <- "#D62728"
color_down <- "#1F77B4"
color_ns   <- "grey70"

volcano_plots <- list()
summary_stats <- list() # Store counts for formal export
all_top_genes <- list() # Store top 10 genes per contrast

for (contrast in names(results_list)) {
  message("Processing Volcano for: ", contrast)
  
  # results_list[[contrast]]$df contains pre-mapped gene symbols from 02_deseq2_dge.R
  plot_df <- results_list[[contrast]]$df
  
  # 2. Robust Significance Classification
  plot_df <- plot_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  
  padj_cutoff <- 0.05
  lfc_cutoff <- 2.0
  
  plot_df <- plot_df %>%
    mutate(
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
  
  # Store stats for the summary table
  summary_stats[[contrast]] <- data.frame(
    Contrast = contrast,
    Up = up_count,
    Down = down_count,
    Total = up_count + down_count
  )

  # 3. Top Genes Extraction (5 Up / 5 Down for formal export)
  top_up <- plot_df %>%
    filter(sig_status == "up", !is.na(gene_name)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 5) %>%
    mutate(Contrast = contrast, Direction = "Up")
    
  top_down <- plot_df %>%
    filter(sig_status == "down", !is.na(gene_name)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 5) %>%
    mutate(Contrast = contrast, Direction = "Down")
    
  all_top_genes[[contrast]] <- bind_rows(top_up, top_down) %>%
    select(Contrast, Direction, gene_name, log2FoldChange, padj)

  # 4. Mentor's Labeling Heuristic
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
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 18) +
    labs(
      title = {
                cl_list <- unique(as.character(metadata$cell_line))
                cl_pattern <- paste0("^(", paste(cl_list, collapse="|"), ")_")
                if(grepl(cl_pattern, contrast)) {
                  cl_prefix <- str_extract(contrast, paste0("^(", paste(cl_list, collapse="|"), ")"))
                  rest <- sub(paste0("^", cl_prefix, "_"), "", contrast)
                  paste0(cl_prefix, ": ", str_replace_all(rest, "_", " "))
                } else {
                  paste0("Global: ", str_replace_all(contrast, "_", " "))
                }
              },
      subtitle = if(grepl("vs_Kromastat", contrast)) "Positive log2FC means higher in Romidepsin;\nnegative log2FC means higher in Kromastat" else NULL,
      x = "log2foldChange",
      y = "-log10 adjusted p-value"
    )
    
    is_drug_drug <- grepl("vs_Kromastat", contrast) && !grepl("DMSO", contrast)
    
    p <- p +
    annotate("text", x = -Inf, y = Inf, 
             label = if(grepl("vs_Kromastat", contrast)) paste0("Higher in Kromastat: ", down_count) else paste0("Downregulated: ", down_count), 
             hjust = -0.05, vjust = -1.2, color = color_down, fontface = "bold", size = if(is_drug_drug) 5.5 else 7.5) +
    annotate("text", x = Inf, y = Inf, 
             label = if(grepl("vs_Kromastat", contrast)) paste0("Higher in Romidepsin: ", up_count) else paste0("Upregulated: ", up_count), 
             hjust = 1.05, vjust = -1.2, color = color_up, fontface = "bold", size = if(is_drug_drug) 5.5 else 7.5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, 
                                margin = margin(b = if(grepl("vs_Kromastat", contrast)) 10 else 50)),
      plot.subtitle = element_text(hjust = 0.5, size = if(is_drug_drug) 13.5 else 18, margin = margin(b = 45)),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2.5, 1, 1, 1, "cm")
    )
  
  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = gene_label),
      size = 5.25,
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
  ggsave(paste0(res_dir, "/figures/05_volcano_", safe_name, ".png"), p, width = 9, height = 7, dpi = 300, bg = "white")
}

# 5. Export Master Summary Tables
message("--- Exporting Authoritative DGE Summary Tables ---")
dir.create(file.path(res_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

# 5.1 Global Counts
final_summary_df <- bind_rows(summary_stats)
write.csv(final_summary_df, file.path(res_dir, "tables/05_dge_summary_stats.csv"), row.names = FALSE)

# 5.2 Top 10 Genes per Contrast
final_top_genes_df <- bind_rows(all_top_genes)
write.csv(final_top_genes_df, file.path(res_dir, "tables/05_top_dge_genes.csv"), row.names = FALSE)

# 5.3 Conserved Targets (H9 ∩ SUPM2)
message("--- Identifying Conserved Targets across Cell Lines ---")
conserved_list <- list()
treatments <- c("Romidepsin_6nM_vs_DMSO_Romidepsin", "Kromastat_6nM_vs_DMSO_Kromastat", "Romidepsin_6nM_vs_Kromastat_6nM")

for (tr in treatments) {
  cl_list <- unique(as.character(metadata$cell_line))
  if (length(cl_list) >= 2) {
    cl1_key <- paste0(cl_list[1], "_", tr)
    cl2_key <- paste0(cl_list[2], "_", tr)
    
    if (cl1_key %in% names(results_list) && cl2_key %in% names(results_list)) {
      cl1_df <- results_list[[cl1_key]]$df %>% 
        filter(padj < 0.05, abs(log2FoldChange) > 2.0) %>%
        select(gene_name, log2FoldChange, padj)
        
      cl2_df <- results_list[[cl2_key]]$df %>% 
        filter(padj < 0.05, abs(log2FoldChange) > 2.0) %>%
        select(gene_name, log2FoldChange, padj)
        
      conserved <- inner_join(cl1_df, cl2_df, by = "gene_name", suffix = paste0("_", cl_list[1:2])) %>%
        filter(!is.na(gene_name)) %>%
        mutate(Treatment = tr)
        
      conserved_list[[tr]] <- conserved
    }
  }
}

if (length(conserved_list) > 0) {
  final_conserved_df <- bind_rows(conserved_list)
  write.csv(final_conserved_df, file.path(res_dir, "tables/05_conserved_targets.csv"), row.names = FALSE)
  message("[OK] Conserved targets table saved.")
}

# 6. Combined Plot
if (length(volcano_plots) > 1) {
  p_combined <- wrap_plots(volcano_plots, ncol = 2) + 
    plot_annotation(title = paste0("Volcano Plot Summary: ", group_name))
  
  ggsave(paste0(res_dir, "/figures/05_volcano_combined.png"), p_combined, 
         width = 16, height = 8 * ceiling(length(volcano_plots)/2), dpi = 300, bg = "white")
}

message("[OK] Volcano plots and summary table completed for ", group_name)
