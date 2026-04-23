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
  list(cols = c("h_h9_r", "d_ul1_r"), labels = c("Human:H9", "Canine:UL1"), title = "Romidepsin Response (H9 vs UL1)", name = "04_11_alluvial_Romi_H9_vs_UL1"),
  list(cols = c("h_supm2_r", "d_cnk89_r"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Romidepsin Response (SUPM2 vs CNK89)", name = "04_11_alluvial_Romi_SUPM2_vs_CNK89"),
  list(cols = c("h_h9_k", "d_ul1_k"), labels = c("Human:H9", "Canine:UL1"), title = "Kromastat Response (H9 vs UL1)", name = "04_11_alluvial_Kroma_H9_vs_UL1"),
  list(cols = c("h_supm2_k", "d_cnk89_k"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Kromastat Response (SUPM2 vs CNK89)", name = "04_11_alluvial_Kroma_SUPM2_vs_CNK89")
)

for (comp in comparisons) {
  make_alluvial_plot(master_df, comp$cols, comp$labels, comp$title, comp$name)
}

# --- Part 2: Master Panoramic Flows ---
message("Generating master panoramic flows...")
# Global Romidepsin Chain
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_r", "h_supm2_r", "d_ul1_r", "d_cnk89_r"),
  labels    = c("H9 (Romi)", "SUPM2 (Romi)", "UL1 (Romi)", "CNK89 (Romi)"),
  title     = "Global Romidepsin Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_hallmark_Romi"
)

# Global Kromastat Chain
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_k", "h_supm2_k", "d_ul1_k", "d_cnk89_k"),
  labels    = c("H9 (Kroma)", "SUPM2 (Kroma)", "UL1 (Kroma)", "CNK89 (Kroma)"),
  title     = "Global Kromastat Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_hallmark_Kroma"
)

# Global Drug Bridge (Aggressive Cancer)
make_alluvial_plot(
  master_df,
  axes_cols = c("h_supm2_r", "h_supm2_k", "d_cnk89_r", "d_cnk89_k"),
  labels    = c("Human:Romi", "Human:Kroma", "Canine:Romi", "Canine:Kroma"),
  title     = "Global Drug-Species Bridge (Aggressive Cell Lines)",
  filename_base = "04_11_alluvial_aggressive"
)

# Global Drug Bridge (Indolent/Healthy-ish Lines)
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_r", "h_h9_k", "d_ul1_r", "d_ul1_k"),
  labels    = c("Human:Romi", "Human:Kroma", "Canine:Romi", "Canine:Kroma"),
  title     = "Global Drug-Species Bridge (Indolent CellLines)",
  filename_base = "04_11_alluvial_indolent"
)

# --- Gene-Level Alluvial Flow Analysis (Ortholog Specific) ---
# We compare specific matched cell line pairs for each drug
# Romidepsin: H9 vs UL1 (Aggressive) | SUPM2 vs CNK89 (Very Aggressive)
# Kromastat:  H9 vs UL1 | SUPM2 vs CNK89

message("\n--- Generating Gene-Level Alluvial Flow Plots ---")

alluvial_pairs <- list(
  list(drug = "Romi", human_cl = "H9", canine_cl = "UL1", h_file = "02_dge_H9_Romi_6nM_vs_DMSO_Romi.csv", c_file = "02_dge_UL1_Romi_6nM_vs_DMSO_Romi.csv"),
  list(drug = "Romi", human_cl = "SUPM2", canine_cl = "CNK89", h_file = "02_dge_SUPM2_Romi_6nM_vs_DMSO_Romi.csv", c_file = "02_dge_CNK89_Romi_6nM_vs_DMSO_Romi.csv"),
  list(drug = "Krom", human_cl = "H9", canine_cl = "UL1", h_file = "02_dge_H9_Kromastat_6nM_vs_DMSO_Kromastat.csv", c_file = "02_dge_UL1_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  list(drug = "Krom", human_cl = "SUPM2", canine_cl = "CNK89", h_file = "02_dge_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.csv", c_file = "02_dge_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.csv")
)

all_alluvial_data <- list()

for (pair in alluvial_pairs) {
  h_path <- file.path(human_tables_dir, pair$h_file)
  c_path <- file.path(canine_tables_dir, pair$c_file)
  
  if (!file.exists(h_path) || !file.exists(c_path)) next
  
  h_df <- read_csv(h_path, show_col_types = FALSE) %>%
    mutate(h_state = case_when(
      padj < 0.05 & log2FoldChange >= 2.0 ~ "Up",
      padj < 0.05 & log2FoldChange <= -2.0 ~ "Down",
      TRUE ~ "NotSig"
    )) %>%
    select(gene_label, h_state, h_log2FC = log2FoldChange) %>% # Keep fold change
    distinct(gene_label, .keep_all = TRUE)
    
  c_df <- read_csv(c_path, show_col_types = FALSE) %>%
    mutate(c_state = case_when(
      padj < 0.05 & log2FoldChange >= 2.0 ~ "Up",
      padj < 0.05 & log2FoldChange <= -2.0 ~ "Down",
      TRUE ~ "NotSig"
    )) %>%
    select(gene_label, c_state, c_log2FC = log2FoldChange) %>% # Keep fold change
    distinct(gene_label, .keep_all = TRUE)
    
  pair_df <- inner_join(h_df, c_df, by = "gene_label") %>%
    filter(h_state != "NotSig" | c_state != "NotSig") %>% # Only significant in at least one
    mutate(
      flow_class = case_when(
        h_state == "Up" & c_state == "Up" ~ "concordant_up",
        h_state == "Down" & c_state == "Down" ~ "concordant_down",
        h_state != "NotSig" & c_state != "NotSig" & h_state != c_state ~ "discordant",
        h_state != "NotSig" & c_state == "NotSig" ~ "human_only",
        h_state == "NotSig" & c_state != "NotSig" ~ "dog_only"
      ),
      drug = pair$drug,
      human_cl = pair$human_cl,
      canine_cl = pair$canine_cl
    )
    
  all_alluvial_data[[length(all_alluvial_data) + 1]] <- pair_df
}

if (length(all_alluvial_data) > 0) {
  full_alluvial_df <- bind_rows(all_alluvial_data)

  p_alluvial <- ggplot(full_alluvial_df,
         aes(axis1 = human_cl, axis2 = h_state, axis3 = c_state, axis4 = canine_cl)) +
    geom_alluvium(aes(fill = flow_class), width = 1/8, alpha = 0.95, color = NA) + # Increased opacity, removed borders
    geom_stratum(width = 1/8, fill = "white", color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
    scale_x_discrete(limits = c("Human cell line", "Human state", "Dog state", "Dog cell line"), expand = c(.1, .1)) +
    scale_fill_manual(values = c(
      "concordant_up" = "#D62728", 
      "concordant_down" = "#1F77B4", 
      "discordant" = "#9467BD", 
      "human_only" = "#E377C2", 
      "dog_only" = "#2CA02C"
    )) +
    facet_wrap(~drug) + # Fixed scale for species-to-species magnitude comparison
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = "Human-Dog ortholog response flows for matched cell-line pairs",
         fill = "Flow class")

  ggsave(file.path(output_dir, "04_11_alluvial_gene_flow.png"), p_alluvial, width = 14, height = 8, bg = "white", dpi = 300)

  # --- Part 4: Export Gene Flow Data Tables ---
  message("Exporting Gene Flow data tables...")
  
  # 4.1 Summary Counts per Drug/Classification
  flow_summary <- full_alluvial_df %>%
    group_by(drug, human_cl, canine_cl, flow_class) %>%
    summarise(gene_count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = flow_class, values_from = gene_count, values_fill = 0) %>%
    mutate(Comparison = paste0(human_cl, " vs ", canine_cl, " ", drug)) %>%
    select(Comparison, everything())
    
  write_csv(flow_summary, file.path(table_dir, "04_11_alluvial_gene_flow_summary.csv"))
  
  # 4.2 Detailed Sorted Gene List (Mentor's Style)
  # Priority: Concordant first, then species-specific
  flow_priority <- c("concordant_up" = 1, "concordant_down" = 2, "discordant" = 3, "human_only" = 4, "dog_only" = 5)
  
  gene_report <- full_alluvial_df %>%
    mutate(
      priority = flow_priority[flow_class],
      Comparison = paste0(human_cl, " vs ", canine_cl, " ", drug)
    ) %>%
    arrange(drug, priority, desc(abs(h_log2FC))) %>% # Sort by impact
    select(
      Drug = drug,
      Comparison,
      `Human gene` = gene_label,
      `Dog gene` = gene_label,
      `Human log2FC` = h_log2FC,
      `Dog log2FC` = c_log2FC,
      FlowClass = flow_class
    )
    
  write_csv(gene_report, file.path(table_dir, "04_11_alluvial_gene_flow_detailed.csv"))
  message("  [OK] Saved Tables: 04_11_alluvial_gene_flow_summary.csv, 04_11_alluvial_gene_flow_detailed.csv")
}

message("[OK] Unified Alluvial Engine Completed.")
