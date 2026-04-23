# [[scripts_downstream/04_11_alluvial_plot.R]]
# Goal: Unified Alluvial Analysis engine for Human-Canine comparative transcriptomics.
# Generates both individual 2-way comparisons and master multi-axis panoramic flows.

library(ggalluvial)
library(tidyverse)

# Path Setup (Relative to scripts_downstream/)
human_tables_dir <- "../results/human/tables"
canine_tables_dir   <- "../results/canine/tables"
output_dir       <- "../results/comparative/figures"
table_dir        <- "../results/comparative/tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Colors for statuses (Synced with NES plots)
status_colors <- c(
  "Activated" = "#D62728",       # Red
  "Suppressed" = "#1F77B4",      # Blue
  "Non-significant" = "#D1D1D1"  # Grey
)

# Helper: Clean Hallmark names
clean_pathway_name <- function(pathway_name) {
  x <- sub("^HALLMARK_", "", pathway_name)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# Helper: Read and harmonize GSEA status
read_gsea_status <- function(file_path, padj_cutoff = 0.05) {
  if (!file.exists(file_path)) {
    message("Warning: File not found: ", file_path)
    return(NULL)
  }
  
  df <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(ID = clean_pathway_name(ID)) %>%
    mutate(status = case_when(
      p.adjust < padj_cutoff & NES > 0 ~ "Activated",
      p.adjust < padj_cutoff & NES < 0 ~ "Suppressed",
      TRUE ~ "Non-significant"
    )) %>%
    select(ID, status)
  
  return(df)
}

# Master Dataset Loading
nodes <- list(
  h_h9_r    = file.path(human_tables_dir, "03_gsea_hallmark_H9_Romi_6nM_vs_DMSO_Romi.csv"),
  h_supm2_r = file.path(human_tables_dir, "03_gsea_hallmark_SUPM2_Romi_6nM_vs_DMSO_Romi.csv"),
  h_h9_k    = file.path(human_tables_dir, "03_gsea_hallmark_H9_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  h_supm2_k = file.path(human_tables_dir, "03_gsea_hallmark_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  d_ul1_r   = file.path(canine_tables_dir, "03_gsea_hallmark_UL1_Romi_6nM_vs_DMSO_Romi.csv"),
  d_cnk89_r = file.path(canine_tables_dir, "03_gsea_hallmark_CNK89_Romi_6nM_vs_DMSO_Romi.csv"),
  d_ul1_k   = file.path(canine_tables_dir, "03_gsea_hallmark_UL1_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  d_cnk89_k = file.path(canine_tables_dir, "03_gsea_hallmark_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.csv")
)

message("Merging datasets for unified analysis...")
master_df <- NULL
for (node_name in names(nodes)) {
  temp_df <- read_gsea_status(nodes[[node_name]])
  if (is.null(temp_df)) next
  colnames(temp_df)[2] <- node_name
  if (is.null(master_df)) { master_df <- temp_df } else { master_df <- full_join(master_df, temp_df, by = "ID") }
}
master_df[is.na(master_df)] <- "Non-significant"
total_n <- nrow(master_df)

# Unified Plotting Function
make_alluvial_plot <- function(df, axes_cols, labels, title, filename_base, is_master = FALSE) {
  plot_data <- df %>%
    select(ID, all_of(axes_cols)) %>%
    group_by_at(axes_cols) %>%
    summarise(Freq = n(), .groups = 'drop')
  
  # Calculate conservation if exactly 2 axes
  subtitle_str <- "Conserved Hallmark Pathways across conditions"
  if (length(axes_cols) == 2) {
    cons_n <- sum(df[[axes_cols[1]]] == df[[axes_cols[2]]])
    cons_pct <- round(100 * cons_n / nrow(df), 1)
    
    h_sig <- df %>% filter(!!sym(axes_cols[1]) != "Non-significant")
    if (nrow(h_sig) > 0) {
      sig_cons_n <- sum(h_sig[[axes_cols[1]]] == h_sig[[axes_cols[2]]])
      sig_cons_pct <- round(100 * sig_cons_n / nrow(h_sig), 1)
      subtitle_str <- paste0("Global Conservation: ", cons_pct, "% | Signal Conservation: ", sig_cons_pct, "%")
    }
  }

  # Programmatically build the aesthetic mapping
  mapping <- aes(y = Freq)
  for (i in seq_along(axes_cols)) {
    mapping[[paste0("axis", i)]] <- sym(axes_cols[i])
  }
  
  p <- ggplot(plot_data, mapping) +
    geom_alluvium(aes(fill = !!sym(axes_cols[1])), width = 1/4, alpha = 0.75, color = "white", linewidth = 0.2) +
    geom_stratum(width = 1/4, fill = "white", color = "black", linewidth = 0.7, alpha = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(paste0(stratum, "\n(", count, ")"))), 
              size = 3.5, fontface = "bold", color = "black") +
    scale_x_discrete(limits = labels, expand = c(.12, .12)) +
    scale_fill_manual(values = status_colors) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "none"
    ) +
    labs(title = title, subtitle = subtitle_str, caption = paste0("N = ", nrow(df), " Hallmark Pathways"))

  ggsave(file.path(output_dir, paste0(filename_base, ".png")), p, width = 12, height = 8, bg = "white")
  
  # Export table for this specific plot
  write_csv(df %>% select(ID, all_of(axes_cols)), 
            file.path(table_dir, paste0(filename_base, ".csv")))
}

# --- Part 1: Individual Comparative Plots ---
message("Generating individual comparisons...")
comparisons <- list(
  list(cols = c("h_h9_r", "d_ul1_r"), labels = c("Human:H9", "Canine:UL1"), title = "Romidepsin Response (H9 vs UL1)", name = "04_11_alluvial_comp_Romi_H9_vs_UL1"),
  list(cols = c("h_supm2_r", "d_cnk89_r"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Romidepsin Response (SUPM2 vs CNK89)", name = "04_11_alluvial_comp_Romi_SUPM2_vs_CNK89"),
  list(cols = c("h_h9_k", "d_ul1_k"), labels = c("Human:H9", "Canine:UL1"), title = "Kromastat Response (H9 vs UL1)", name = "04_11_alluvial_comp_Kroma_H9_vs_UL1"),
  list(cols = c("h_supm2_k", "d_cnk89_k"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Kromastat Response (SUPM2 vs CNK89)", name = "04_11_alluvial_comp_Kroma_SUPM2_vs_CNK89")
)

for (comp in comparisons) {
  make_alluvial_plot(master_df, comp$cols, comp$labels, comp$title, comp$name)
}

# --- Part 2: Master Panoramic Flows ---
message("Generating master panoramic flows...")
# Master Romidepsin Chain
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_r", "h_supm2_r", "d_ul1_r", "d_cnk89_r"),
  labels    = c("H9 (Romi)", "SUPM2 (Romi)", "UL1 (Romi)", "CNK89 (Romi)"),
  title     = "Master Romidepsin Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_master_romidepsin"
)

# Master Kromastat Chain
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_k", "h_supm2_k", "d_ul1_k", "d_cnk89_k"),
  labels    = c("H9 (Kroma)", "SUPM2 (Kroma)", "UL1 (Kroma)", "CNK89 (Kroma)"),
  title     = "Master Kromastat Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_master_kromastat"
)

# Drug Bridge (Aggressive Cancer)
make_alluvial_plot(
  master_df,
  axes_cols = c("h_supm2_r", "h_supm2_k", "d_cnk89_r", "d_cnk89_k"),
  labels    = c("Human:Romi", "Human:Kroma", "Canine:Romi", "Canine:Kroma"),
  title     = "Master Drug-Species Bridge (Aggressive/Cancer Lines)",
  filename_base = "04_11_alluvial_master_bridge_aggressive"
)

message("[OK] Unified Alluvial Engine Completed.")
