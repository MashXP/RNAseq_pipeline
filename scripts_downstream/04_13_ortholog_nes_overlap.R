# [[scripts_downstream/04_13_ortholog_nes_overlap.R]]
# Goal: Human-Canine functional enrichment (GSEA NES) comparison.
# Implements cross-species pathway activity validation using side-by-side bar plots.

library(tidyverse)
library(ggplot2)

# Paths
project_dir <- ".."
human_results_dir  <- file.path(project_dir, "results", "human", "tables")
canine_results_dir <- file.path(project_dir, "results", "canine", "tables")
output_dir         <- file.path(project_dir, "results", "comparative")
plots_dir          <- file.path(output_dir, "figures")
tables_dir         <- file.path(output_dir, "tables")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Define Comparison Pairs
# Format: list(label, human_file, canine_file)
comparison_list <- list(
  # a/ Pooled Data
  list(id = "Romidepsin", label = "Globally Pooled: Romidepsin", 
       h_file = "03_gsea_hallmark_Romidepsin_6nM_vs_DMSO_Romidepsin.csv", 
       c_file = "03_gsea_hallmark_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  list(id = "Kromastat", label = "Globally Pooled: Kromastat", 
       h_file = "03_gsea_hallmark_Kromastat_6nM_vs_DMSO_Kromastat.csv", 
       c_file = "03_gsea_hallmark_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
       
  # b/ Specific Cell Line Matches
  list(id = "H9_CNK89_Kromastat", label = "H9 (Human) vs CNK89 (Canine) - Kromastat", 
       h_file = "03_gsea_hallmark_H9_Kromastat_6nM_vs_DMSO_Kromastat.csv", 
       c_file = "03_gsea_hallmark_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  list(id = "H9_CNK89_Romidepsin", label = "H9 (Human) vs CNK89 (Canine) - Romidepsin", 
       h_file = "03_gsea_hallmark_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.csv", 
       c_file = "03_gsea_hallmark_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  list(id = "SUPM2_UL1_Kromastat", label = "SUPM2 (Human) vs UL1 (Canine) - Kromastat", 
       h_file = "03_gsea_hallmark_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.csv", 
       c_file = "03_gsea_hallmark_UL1_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  list(id = "SUPM2_UL1_Romidepsin", label = "SUPM2 (Human) vs UL1 (Canine) - Romidepsin", 
       h_file = "03_gsea_hallmark_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.csv", 
       c_file = "03_gsea_hallmark_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.csv")
)

# 3. Plotting Function
plot_nes_overlap <- function(pair) {
  h_path <- file.path(human_results_dir, pair$h_file)
  c_path <- file.path(canine_results_dir, pair$c_file)
  
  if (!file.exists(h_path) || !file.exists(c_path)) {
    message("  [SKIP] Files not found for: ", pair$label)
    return(NULL)
  }
  
  h_df <- read_csv(h_path, show_col_types = FALSE) %>%
    select(ID, h_NES = NES, h_padj = p.adjust)
  c_df <- read_csv(c_path, show_col_types = FALSE) %>%
    select(ID, c_NES = NES, c_padj = p.adjust)
    
  merged <- inner_join(h_df, c_df, by = "ID") %>%
    mutate(ID = str_replace(ID, "HALLMARK_", "")) %>%
    mutate(sig_status = case_when(
      h_padj < 0.05 & c_padj < 0.05 ~ "Both Significant",
      h_padj < 0.05 ~ "Human Only",
      c_padj < 0.05 ~ "Canine Only",
      TRUE ~ "NS"
    ))

  # Reshape for Bar Plot
  plot_df <- merged %>%
    select(ID, Human = h_NES, Canine = c_NES, sig_status) %>%
    pivot_longer(cols = c(Human, Canine), names_to = "Species", values_to = "NES") %>%
    # Sort by Human NES for better visualization
    mutate(ID = factor(ID, levels = merged$ID[order(merged$h_NES)]))

  p <- ggplot(plot_df, aes(x = ID, y = NES, fill = Species)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Human" = "#E377C2", "Canine" = "#2CA02C")) +
    theme_minimal(base_size = 18) +
    labs(title = paste0("Pathway NES Comparison: ", pair$label),
         x = "Hallmark Pathway",
         y = "Normalized Enrichment Score (NES)") +
    theme(
      axis.text.y = element_text(size = 10.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(10, 160, 10, 10)
    )

  ggsave(file.path(plots_dir, paste0("04_13_nes_bar_", pair$id, ".png")), p, width = 14, height = 12, bg = "white", dpi = 300)
  write_csv(merged, file.path(tables_dir, paste0("04_13_nes_overlap_", pair$id, ".csv")))
}

# 4. Execute Analysis
message("\n--- Generating Species-Agnostic NES Overlap Plots ---")
for (pair in comparison_list) {
  message("Processing: ", pair$label)
  plot_nes_overlap(pair)
}

message("\n[OK] NES Overlap Analysis Completed.")
