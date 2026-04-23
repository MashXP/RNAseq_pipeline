# [[scripts_downstream/04_09_upset_consistency.R]]
# Goal: Generate UpSet plots for cross-cellline consistency analysis.

library(DESeq2)
library(ggplot2)
library(UpSetR)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 04_09_upset_consistency.R <Group> [Species]")
}
group_name <- tolower(args[1])
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# 2. Extract Significant Gene Lists
# We want to compare Romidepsin consistency and Kromastat consistency
# Threshold: padj < 0.05 & |log2FoldChange| > 2
get_sig_ids <- function(res_entry) {
  if (is.null(res_entry$df)) return(character(0))
  df <- res_entry$df
  df$Geneid[df$padj < 0.05 & abs(df$log2FoldChange) > 2 & !is.na(df$padj)]
}

# 2. Identify Cell Line Subsets
# Identify Cell Line Subsets dynamically from the results keys
keys <- names(results_list)
# We assume cell line subsets are the ones that are NOT the global pooled contrasts
global_contrasts <- c("Romi_6nM_vs_DMSO_Romi", "Kromastat_6nM_vs_DMSO_Kromastat", "Romi_6nM_vs_Kromastat_6nM", "DMSO_Kromastat_vs_DMSO_Romi")
cell_lines <- unique(str_extract(keys[!keys %in% global_contrasts], "^[A-Za-z0-9]+(?=_)"))
cell_lines <- cell_lines[!is.na(cell_lines)]

if (length(cell_lines) < 2) {
  message("\n[CRITICAL] Error: Could not find cell-line subset results (e.g., H9_..., SUPM2_...).")
  message(">> Did you remember to re-run '02_deseq2_dge.R' after the upgrade?")
  message(">> Current keys in results_list: ", paste(head(keys), collapse=", "))
  stop("Missing subset results.")
}

message("--- Generating UpSet Consistency Report for ", group_name, " (", paste(cell_lines, collapse=", "), ") ---")

# A. Romidepsin Consistency
romi_case <- "Romi_6nM_vs_DMSO_Romi"
romi_sets <- list()
for (cl in cell_lines) {
  romi_sets[[paste0(cl, " Romi")]] <- get_sig_ids(results_list[[paste0(cl, "_", romi_case)]])
}

# B. Kromastat Consistency
krom_case <- "Kromastat_6nM_vs_DMSO_Kromastat"
krom_sets <- list()
for (cl in cell_lines) {
  krom_sets[[paste0(cl, " Krom")]] <- get_sig_ids(results_list[[paste0(cl, "_", krom_case)]])
}

# C. Combine for a 4-way comparison if possible
all_sets <- c(romi_sets, krom_sets)

# Remove empty sets to avoid UpSet errors
all_sets <- all_sets[sapply(all_sets, length) > 0]

if (length(all_sets) >= 2) {
  message("Saving 4-way UpSet Plot...")
  
  png(file.path(res_dir, "figures/04_09_upset_consistency.png"), 
      width = 1800, height = 900, res = 150)
  
  # UpSetR uses grid graphics. par(oma) is often ignored.
  # We use a grid viewport to shift the entire plot to the right, creating a large left margin.
  library(grid)
  grid.newpage()
  pushViewport(viewport(x = 0.50, width = 0.9)) # Shift center to 0.54 (from 0.5) to add moderate left padding
  
  # Standard UpSet Plot with Red/Blue premium theme
  # We use draw.plot = FALSE then grid.draw to avoid newpage issues if necessary, 
  # but print() usually works within a viewport if we don't call grid.newpage() inside it.
  upset_plot <- upset(fromList(all_sets), 
         nsets = length(all_sets), 
         order.by = "freq", 
         main.bar.color = "#D62728", # Red for intersections
         sets.bar.color = "#1F77B4", # Blue for sets
         matrix.color = "#1F77B4",
         text.scale = c(1.5, 1.2, 1.2, 1.0, 1.2, 1.2), # 4th is Set Size axis labels
         point.size = 3.5, 
         line.size = 1.0)
  
  print(upset_plot, newpage = FALSE)
  popViewport()
  
  
  grid.text(paste0("Consistency Matrix: ", group_name), 
            x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
  
  dev.off()
  message("[OK] UpSet plot saved -> 04_09_upset_consistency.png")
} else {
  message("Skipping UpSet plot: Not enough significant genes in subsets.")
}
