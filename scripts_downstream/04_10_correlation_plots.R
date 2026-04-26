# [[scripts_downstream/04_10_correlation_plots.R]]
# Goal: Generate LFC-LFC correlation plots for consistency proofing.

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_10_correlation_plots.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# Identify Cell Line Subsets dynamically from the results keys
keys <- names(results_list)
# We assume cell line subsets are the ones that are NOT the global pooled contrasts
global_contrasts <- c("Romidepsin_6nM_vs_DMSO_Romidepsin", "Kromastat_6nM_vs_DMSO_Kromastat", "Romidepsin_6nM_vs_Kromastat_6nM", "DMSO_Kromastat_vs_DMSO_Romidepsin")
cell_lines <- unique(str_extract(keys[!keys %in% global_contrasts], "^[A-Za-z0-9]+(?=_)"))
cell_lines <- cell_lines[!is.na(cell_lines)]

if (length(cell_lines) < 2) {
  message("\n[CRITICAL] Error: Could not find cell-line subset results (e.g., H9_..., SUPM2_...).")
  message(">> Did you remember to re-run '02_deseq2_dge.R' after the upgrade?")
  stop("Missing subset results.")
}

cl1 <- cell_lines[1]
cl2 <- cell_lines[2]

message("--- Generating LFC Correlation Report (", cl1, " vs ", cl2, ") ---")

# Threshold definition for coloring
lfc_thresh <- 2.0
padj_thresh <- 0.05

create_corr_plot <- function(drug_name, stub) {
  # Join results from both cell lines
  df1 <- results_list[[paste0(cl1, "_", stub)]]$df
  df2 <- results_list[[paste0(cl2, "_", stub)]]$df
  
  if (is.null(df1) || is.null(df2)) return(NULL)
  
  # All genes present in both cell lines are retained to show the unbiased global LFC trend.
  # Significance status is encoded by color. 
  # Criteria: padj < 0.05 AND |LFC| > 2.0
  merged <- inner_join(
    df1 %>% select(Geneid, gene_label, lfc1 = log2FoldChange, padj1 = padj),
    df2 %>% select(Geneid, lfc2 = log2FoldChange, padj2 = padj),
    by = "Geneid"
  ) %>%
    mutate(
      is_sig1 = (padj1 < padj_thresh & abs(lfc1) > lfc_thresh),
      is_sig2 = (padj2 < padj_thresh & abs(lfc2) > lfc_thresh),
      status = case_when(
        (is_sig1 & is_sig2) ~ "Consensus (Sig in both)",
        (is_sig1 | is_sig2) ~ "Specific (Sig in one)",
        TRUE ~ "Non-significant"
      )
    )
  
  # Calculate Pearson Correlation
  corr_val <- cor(merged$lfc1, merged$lfc2, method = "pearson", use = "complete.obs")
  
  # Annotations: Top 5 Consensus, Top 5 Specific
  consensus_top <- merged %>% 
    filter(status == "Consensus (Sig in both)") %>%
    mutate(abs_lfc = (abs(lfc1) + abs(lfc2))/2) %>%
    slice_max(order_by = abs_lfc, n = 5) %>%
    mutate(label = gene_label)
    
  specific_top <- merged %>%
    filter(status == "Specific (Sig in one)") %>%
    mutate(abs_lfc = pmax(abs(lfc1), abs(lfc2))) %>%
    slice_max(order_by = abs_lfc, n = 5) %>%
    mutate(label = ifelse(is_sig1, paste0(gene_label, " (I)"), paste0(gene_label, " (A)")))
    
  label_df <- bind_rows(consensus_top, specific_top)

  ggplot(merged, aes(x = lfc1, y = lfc2, color = status)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_text_repel(data = label_df, aes(label = label), size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(
      "Consensus (Sig in both)" = "#D62728", # Premium Red
      "Specific (Sig in one)" = "#1F77B4",   # Premium Blue
      "Non-significant" = "grey80"
    )) +
    theme_bw(base_size = 18) +
    labs(
      title = paste0(drug_name, " Mechanism Consistency"),
      subtitle = paste0(cl1, " vs ", cl2, " | R = ", round(corr_val, 3), " | |LFC| > ", lfc_thresh),
      x = paste0("log2FoldChange (", cl1, ")"),
      y = paste0("log2FoldChange (", cl2, ")"),
      color = "Significance"
    )
}

p_romi <- create_corr_plot("Romidepsin", "Romidepsin_6nM_vs_DMSO_Romidepsin")
p_krom <- create_corr_plot("Kromastat", "Kromastat_6nM_vs_DMSO_Kromastat")

if (!is.null(p_romi) && !is.null(p_krom)) {
  p_combined <- (p_romi + p_krom) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Transcriptome Response Consistency: ", group_name),
      subtitle = paste0("(I): Indolent (", cl1, ") | (A): Aggressive (", cl2, ")"),
      theme = theme(
        plot.title = element_text(size = 27, face = "bold", hjust = 0.5, margin = margin(b = 10, t = 10)),
        plot.subtitle = element_text(size = 18, face = "italic", hjust = 0.5, margin = margin(b = 10))
      )
    ) & theme(legend.position = "bottom")
  
  ggsave(file.path(res_dir, "figures/04_10_lfc_correlation.png"), p_combined, 
         width = 14, height = 9, dpi = 300, bg = "white")
  
  message("[OK] Correlation plot saved -> 04_10_lfc_correlation.png")
}
