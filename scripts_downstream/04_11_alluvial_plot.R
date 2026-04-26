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

# Colors for statuses and flows
status_colors <- c(
  "Activated" = "#D62728",       # Red
  "Suppressed" = "#1F77B4",      # Blue
  "Non-significant" = "#D1D1D1"  # Grey
)

flow_colors <- c(
  "concordant_activated" = "#D62728",    # Red
  "concordant_suppressed" = "#1F77B4",  # Blue
  "discordant" = "#9467BD",             # Purple
  "human_only" = "#E377C2",             # Pink
  "canine_only" = "#2CA02C",            # Green
  "Non-significant" = "#D1D1D1"         # Grey
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
  h_h9_r    = file.path(human_tables_dir, "03_gsea_hallmark_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  h_supm2_r = file.path(human_tables_dir, "03_gsea_hallmark_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  h_h9_k    = file.path(human_tables_dir, "03_gsea_hallmark_H9_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  h_supm2_k = file.path(human_tables_dir, "03_gsea_hallmark_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  d_ul1_r   = file.path(canine_tables_dir, "03_gsea_hallmark_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  d_cnk89_r = file.path(canine_tables_dir, "03_gsea_hallmark_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
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
make_alluvial_plot <- function(df, axes_cols, labels, title, filename_base, ref_cols = NULL) {
  is_two_way <- length(axes_cols) == 2
  
  if (is_two_way) {
    # Add flow classification for 2-way pathway comparison
    plot_df <- df %>%
      mutate(
        flow_class = case_when(
          !!sym(axes_cols[1]) == "Activated" & !!sym(axes_cols[2]) == "Activated" ~ "concordant_activated",
          !!sym(axes_cols[1]) == "Suppressed" & !!sym(axes_cols[2]) == "Suppressed" ~ "concordant_suppressed",
          !!sym(axes_cols[1]) != "Non-significant" & !!sym(axes_cols[2]) != "Non-significant" & !!sym(axes_cols[1]) != !!sym(axes_cols[2]) ~ "discordant",
          !!sym(axes_cols[1]) != "Non-significant" & !!sym(axes_cols[2]) == "Non-significant" ~ "human_only",
          !!sym(axes_cols[1]) == "Non-significant" & !!sym(axes_cols[2]) != "Non-significant" ~ "canine_only",
          TRUE ~ "Non-significant"
        )
      )
    
    plot_data <- plot_df %>%
      filter(flow_class != "Non-significant") %>%
      group_by_at(c(axes_cols, "flow_class")) %>%
      summarise(Freq = n(), .groups = 'drop')
      
    fill_col <- "flow_class"
    fill_scale <- scale_fill_manual(values = flow_colors, name = "Pathway Flow")
  } else if (!is.null(ref_cols)) {
    # Multi-axis plot with cross-species flow classification based on ref_cols
    plot_df <- df %>%
      mutate(
        flow_class = case_when(
          !!sym(ref_cols[1]) == "Activated" & !!sym(ref_cols[2]) == "Activated" ~ "concordant_activated",
          !!sym(ref_cols[1]) == "Suppressed" & !!sym(ref_cols[2]) == "Suppressed" ~ "concordant_suppressed",
          !!sym(ref_cols[1]) != "Non-significant" & !!sym(ref_cols[2]) != "Non-significant" & !!sym(ref_cols[1]) != !!sym(ref_cols[2]) ~ "discordant",
          !!sym(ref_cols[1]) != "Non-significant" & !!sym(ref_cols[2]) == "Non-significant" ~ "human_only",
          !!sym(ref_cols[1]) == "Non-significant" & !!sym(ref_cols[2]) != "Non-significant" ~ "canine_only",
          TRUE ~ "Non-significant"
        )
      )
    
    plot_data <- plot_df %>%
      group_by_at(c(axes_cols, "flow_class")) %>%
      summarise(Freq = n(), .groups = 'drop')
    
    fill_col <- "flow_class"
    fill_scale <- scale_fill_manual(values = flow_colors, name = "Cross-Species Flow (Romidepsin ref)")
  } else {
    plot_df <- df
    plot_data <- plot_df %>%
      group_by_at(axes_cols) %>%
      summarise(Freq = n(), .groups = 'drop')
    
    fill_col <- axes_cols[1]
    fill_scale <- scale_fill_manual(values = status_colors, name = "Initial State")
  }

  # Calculate conservation if exactly 2 axes
  subtitle_str <- "Conserved Hallmark Pathways across conditions"
  if (is_two_way) {
    cons_n <- sum(plot_df[[axes_cols[1]]] == plot_df[[axes_cols[2]]])
    cons_pct <- round(100 * cons_n / nrow(plot_df), 1)
    
    h_sig <- plot_df %>% filter(!!sym(axes_cols[1]) != "Non-significant")
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
  
  # Dynamic expand: 4-axis needs more horizontal padding than 2-axis
  x_expand <- if (length(axes_cols) <= 2) c(.18, .18) else c(.25, .25)
  # Dynamic width: 4-axis plots need more canvas
  plot_width <- if (length(axes_cols) <= 2) 16 else 22

  p <- ggplot(plot_data, mapping) +
    geom_alluvium(aes(fill = !!sym(fill_col)), width = 1/4, alpha = 0.75, color = "white", linewidth = 0.2) +
    geom_stratum(width = 1/4, fill = "white", color = "black", linewidth = 0.7, alpha = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(paste0(stratum, "\n(", count, ")"))), 
              size = 5.25, fontface = "bold", color = "black") +
    scale_x_discrete(limits = labels, expand = x_expand) +
    fill_scale +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 21) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 15, face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(10, 40, 10, 40)
    ) +
    labs(title = title, subtitle = subtitle_str, caption = paste0("N = ", nrow(plot_df), " Hallmark Pathways"))

  ggsave(file.path(output_dir, paste0(filename_base, ".png")), p, width = plot_width, height = 9, bg = "white", dpi = 300)
  
  # Export table for this specific plot
  write_csv(plot_df %>% select(ID, all_of(axes_cols)), 
            file.path(table_dir, paste0(filename_base, ".csv")))
}

# --- Part 1: Individual Comparative Plots ---
message("Generating individual comparisons...")
comparisons <- list(
  list(cols = c("h_h9_r", "d_ul1_r"), labels = c("Human:H9", "Canine:UL1"), title = "Romidepsin Response (H9 vs UL1)", name = "04_11_alluvial_Romidepsin_H9_vs_UL1"),
  list(cols = c("h_supm2_r", "d_cnk89_r"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Romidepsin Response (SUPM2 vs CNK89)", name = "04_11_alluvial_Romidepsin_SUPM2_vs_CNK89"),
  list(cols = c("h_h9_k", "d_ul1_k"), labels = c("Human:H9", "Canine:UL1"), title = "Kromastat Response (H9 vs UL1)", name = "04_11_alluvial_Kromastat_H9_vs_UL1"),
  list(cols = c("h_supm2_k", "d_cnk89_k"), labels = c("Human:SUPM2", "Canine:CNK89"), title = "Kromastat Response (SUPM2 vs CNK89)", name = "04_11_alluvial_Kromastat_SUPM2_vs_CNK89")
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
  labels    = c("H9 (Romidepsin)", "SUPM2 (Romidepsin)", "UL1 (Romidepsin)", "CNK89 (Romidepsin)"),
  title     = "Global Romidepsin Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_hallmark_Romidepsin"
)

# Global Kromastat Chain
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_k", "h_supm2_k", "d_ul1_k", "d_cnk89_k"),
  labels    = c("H9 (Kromastat)", "SUPM2 (Kromastat)", "UL1 (Kromastat)", "CNK89 (Kromastat)"),
  title     = "Global Kromastat Phenotypic Flow (Healthy -> Cancer | Human -> Canine)",
  filename_base = "04_11_alluvial_hallmark_Kromastat"
)

# Global Drug Bridge (Aggressive Cancer)
make_alluvial_plot(
  master_df,
  axes_cols = c("h_supm2_r", "h_supm2_k", "d_cnk89_r", "d_cnk89_k"),
  ref_cols  = c("h_supm2_r", "d_cnk89_r"),
  labels    = c("Human:Romidepsin", "Human:Kromastat", "Canine:Romidepsin", "Canine:Kromastat"),
  title     = "Global Drug-Species Bridge (Aggressive Cell Lines)",
  filename_base = "04_11_alluvial_aggressive"
)

# Global Drug Bridge (Indolent/Healthy-ish Lines)
make_alluvial_plot(
  master_df,
  axes_cols = c("h_h9_r", "h_h9_k", "d_ul1_r", "d_ul1_k"),
  ref_cols  = c("h_h9_r", "d_ul1_r"),
  labels    = c("Human:Romidepsin", "Human:Kromastat", "Canine:Romidepsin", "Canine:Kromastat"),
  title     = "Global Drug-Species Bridge (Indolent Cell Lines)",
  filename_base = "04_11_alluvial_indolent"
)

# --- Gene-Level Alluvial Flow Analysis (Ortholog Specific) ---
# We compare specific matched cell line pairs for each drug
# Romidepsin: H9 vs UL1 (Aggressive) | SUPM2 vs CNK89 (Very Aggressive)
# Kromastat:  H9 vs UL1 | SUPM2 vs CNK89

message("\n--- Generating Gene-Level Alluvial Flow Plots ---")

alluvial_pairs <- list(
  list(drug = "Romidepsin", human_cl = "H9", canine_cl = "UL1", h_file = "02_dge_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.csv", c_file = "02_dge_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
  list(drug = "Romidepsin", human_cl = "SUPM2", canine_cl = "CNK89", h_file = "02_dge_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.csv", c_file = "02_dge_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.csv"),
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
      padj < 0.05 & log2FoldChange >= 2.0 ~ "Activated",
      padj < 0.05 & log2FoldChange <= -2.0 ~ "Suppressed",
      TRUE ~ "Non-significant"
    )) %>%
    select(gene_label, h_state, h_log2FC = log2FoldChange) %>% # Keep fold change
    distinct(gene_label, .keep_all = TRUE)
    
  c_df <- read_csv(c_path, show_col_types = FALSE) %>%
    mutate(c_state = case_when(
      padj < 0.05 & log2FoldChange >= 2.0 ~ "Activated",
      padj < 0.05 & log2FoldChange <= -2.0 ~ "Suppressed",
      TRUE ~ "Non-significant"
    )) %>%
    select(gene_label, c_state, c_log2FC = log2FoldChange) %>% # Keep fold change
    distinct(gene_label, .keep_all = TRUE)
    
  pair_df <- inner_join(h_df, c_df, by = "gene_label") %>%
    filter(h_state != "Non-significant" | c_state != "Non-significant") %>% # Only significant in at least one
    mutate(
      flow_class = case_when(
        h_state == "Activated" & c_state == "Activated" ~ "concordant_activated",
        h_state == "Suppressed" & c_state == "Suppressed" ~ "concordant_suppressed",
        h_state != "Non-significant" & c_state != "Non-significant" & h_state != c_state ~ "discordant",
        h_state != "Non-significant" & c_state == "Non-significant" ~ "human_only",
        h_state == "Non-significant" & c_state != "Non-significant" ~ "canine_only"
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
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4.5, fontface = "bold") +
    scale_x_discrete(limits = c("Human cell line", "Human state", "Canine state", "Canine cell line"), expand = c(.18, .18)) +
    scale_fill_manual(values = c(
      "concordant_activated" = "#D62728", 
      "concordant_suppressed" = "#1F77B4", 
      "discordant" = "#9467BD", 
      "human_only" = "#E377C2", 
      "canine_only" = "#2CA02C"
    )) +
    facet_wrap(~drug) + # Fixed scale for species-to-species magnitude comparison
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = "Human-Canine ortholog response flows for matched cell-line pairs",
         fill = "Flow class")

  ggsave(file.path(output_dir, "04_11_alluvial_gene_flow.png"), p_alluvial, width = 18, height = 10, bg = "white", dpi = 300)

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
  flow_priority <- c("concordant_activated" = 1, "concordant_suppressed" = 2, "discordant" = 3, "human_only" = 4, "canine_only" = 5)
  
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
      `Canine gene` = gene_label,
      `Human log2FC` = h_log2FC,
      `Canine log2FC` = c_log2FC,
      FlowClass = flow_class
    )
    
  write_csv(gene_report, file.path(table_dir, "04_11_alluvial_gene_flow_detailed.csv"))
  message("  [OK] Saved Tables: 04_11_alluvial_gene_flow_summary.csv, 04_11_alluvial_gene_flow_detailed.csv")
}

message("[OK] Unified Alluvial Engine Completed.")
