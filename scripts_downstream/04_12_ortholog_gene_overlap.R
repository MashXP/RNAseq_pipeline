# [[scripts_downstream/04_12_ortholog_gene_overlap.R]]
# Goal: Human-Canine ortholog comparison of differential gene expression.
# Implements mentor's gene-specific overlap logic for cross-species validation.

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Paths
project_dir <- ".."
human_results_dir  <- file.path(project_dir, "results", "human", "tables")
canine_results_dir <- file.path(project_dir, "results", "canine", "tables")
output_dir         <- file.path(project_dir, "results", "comparative")
plots_dir          <- file.path(output_dir, "figures")
tables_dir         <- file.path(output_dir, "tables")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Build Ortholog Map (Symbol-based Fallback)
message("Building ortholog map from local GTF files...")

extract_gtf_attribute <- function(attr_text, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  value <- sub(paste0(".*", pattern, ".*"), "\\1", attr_text)
  value[!grepl(pattern, attr_text)] <- NA_character_
  value
}

read_gtf_gene_map <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF file not found: ", gtf_path)
  }
  gtf_df <- read.delim(gtf_path, header = FALSE, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE)
  colnames(gtf_df)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gene_df <- gtf_df[gtf_df$feature == "gene", "attribute", drop = FALSE]
  
  out <- data.frame(
    gene_id = extract_gtf_attribute(gene_df$attribute, "gene_id"),
    gene_name = extract_gtf_attribute(gene_df$attribute, "gene_name"),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(gene_id), !is.na(gene_name), gene_name != "") %>%
    distinct(gene_name, .keep_all = TRUE)
  return(out)
}

human_gtf  <- file.path(project_dir, "_data", "genome", "Human", "Homo_sapiens.GRCh38.113.gtf")
canine_gtf <- file.path(project_dir, "_data", "genome", "Canine", "Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf")

human_map  <- read_gtf_gene_map(human_gtf)
canine_map <- read_gtf_gene_map(canine_gtf)

# Join on gene_name (the common bridge)
ortholog_tbl <- inner_join(human_map, canine_map, by = "gene_name", suffix = c("_human", "_canine"))
message("Matched ", nrow(ortholog_tbl), " orthologs by gene symbol.")

# 2. Comparison Logic
# We compare matched contrasts between species
comparisons <- data.frame(
  label = c("Romi 6nM", "Kromastat 6nM"),
  human_stub  = c("02_dge_Romi_6nM_vs_DMSO_Romi.csv", "02_dge_Kromastat_6nM_vs_DMSO_Kromastat.csv"),
  canine_stub = c("02_dge_Romi_6nM_vs_DMSO_Romi.csv", "02_dge_Kromastat_6nM_vs_DMSO_Kromastat.csv")
)

all_merged_results <- list()

for (i in 1:nrow(comparisons)) {
  label <- comparisons$label[i]
  message("\n--- Comparing cross-species response for: ", label, " ---")
  
  human_file  <- file.path(human_results_dir, comparisons$human_stub[i])
  canine_file <- file.path(canine_results_dir, comparisons$canine_stub[i])
  
  if (!file.exists(human_file) || !file.exists(canine_file)) {
    message("  Warning: One or both species files missing for ", label)
    next
  }
  
  human_res <- read_csv(human_file, show_col_types = FALSE) %>%
    select(gene_label, log2FoldChange, padj) %>%
    rename(h_lfc = log2FoldChange, h_padj = padj)
    
  canine_res <- read_csv(canine_file, show_col_types = FALSE) %>%
    select(gene_label, log2FoldChange, padj) %>%
    rename(c_lfc = log2FoldChange, c_padj = padj)
    
  merged <- inner_join(human_res, canine_res, by = "gene_label") %>%
    mutate(
      status = case_when(
        h_padj < 0.05 & c_padj < 0.05 & h_lfc >= 2.0 & c_lfc >= 2.0 ~ "shared_up",
        h_padj < 0.05 & c_padj < 0.05 & h_lfc <= -2.0 & c_lfc <= -2.0 ~ "shared_down",
        TRUE ~ "other"
      ),
      comparison = label
    )
  
  all_merged_results[[label]] <- merged
  
  # A. Scatter Plot
  p <- ggplot(merged, aes(x = h_lfc, y = c_lfc, color = status)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("shared_up" = "#D62728", "shared_down" = "#1F77B4", "other" = "grey80")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    theme_bw() +
    labs(
      title = paste0("Ortholog Overlap: ", label),
      subtitle = "Criteria: padj < 0.05 & |log2FC| > 2",
      x = "Human log2FC",
      y = "Canine log2FC"
    )
    
  safe_name <- str_replace_all(label, " ", "_")
  ggsave(file.path(plots_dir, paste0("04_12_ortholog_scatter_", safe_name, ".png")), p, width = 8, height = 7, bg = "white")
  
  # B. Export Tables
  write_csv(merged, file.path(tables_dir, paste0("04_12_ortholog_merged_", safe_name, ".csv")))
}

# 3. Conserved Ortholog Heatmap (Panoramic)
message("\n--- Generating Conserved Ortholog Heatmap ---")
master_df <- bind_rows(all_merged_results) %>%
  filter(status != "other") %>%
  mutate(mean_abs_lfc = (abs(h_lfc) + abs(c_lfc)) / 2) %>%
  arrange(desc(mean_abs_lfc))

if (nrow(master_df) >= 2) {
  # Select Top 40 genes for the heatmap
  top_genes <- master_df %>%
    distinct(gene_label, .keep_all = TRUE) %>%
    head(40) %>%
    pull(gene_label)
    
  heatmap_data <- master_df %>%
    filter(gene_label %in% top_genes) %>%
    select(gene_label, h_lfc, c_lfc, comparison) %>%
    pivot_longer(cols = c(h_lfc, c_lfc), names_to = "species", values_to = "lfc") %>%
    mutate(Species = ifelse(species == "h_lfc", "Human", "Canine")) %>%
    select(-species) %>%
    pivot_wider(names_from = c(comparison, Species), values_from = lfc) %>%
    column_to_rownames("gene_label")
    
  heatmap_mat <- as.matrix(heatmap_data)
  heatmap_mat[is.na(heatmap_mat)] <- 0
  
  png(file.path(plots_dir, "04_12_conserved_ortholog_heatmap.png"), width = 1200, height = 1600, res = 150)
  pheatmap(
    heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(255),
    border_color = NA,
    main = "Top Conserved Orthologs (log2FC)"
  )
  dev.off()
}

# 4. Global Summary Stats
summary_stats <- bind_rows(all_merged_results) %>%
  group_by(comparison, status) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = status, values_from = count, values_fill = 0)

write_csv(summary_stats, file.path(tables_dir, "04_12_ortholog_overlap_summary.csv"))

message("[OK] Gene-specific ortholog overlap completed.")
