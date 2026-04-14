# [[scripts_downstream/04_03_venn.R]]
# Goal: Generate DEG Overlap Venn Diagram for specific group.

library(DESeq2)
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_03_venn.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

res_dir <- paste0("../results/", group_name)

# 1. Load data
load(paste0("./.RData/", group_name, "/02_deseq_results.RData"))

# Ensure output directories exist
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 4. Venn Diagram (Refined)
# Extract significant gene sets for each contrast
sig_list_all <- lapply(results_list, function(x) {
  rownames(subset(as.data.frame(x$shrunk), padj < 0.05 & abs(log2FoldChange) > 1))
})

# Prioritize the 3 primary biological contrasts for a clear 3-way Venn
# but also handle cases where they might be missing
priority_contrasts <- c(
  "Romi_6nM_vs_DMSO_Romi",
  "Kromastat_6nM_vs_DMSO_Kromastat",
  "Romi_6nM_vs_Kromastat_6nM"
)

sig_list_focused <- sig_list_all[names(sig_list_all) %in% priority_contrasts]

# Filter out empty sets to avoid ambiguous "empty" circles
sig_list_final <- sig_list_focused[sapply(sig_list_focused, length) > 0]

if (length(sig_list_final) >= 2) {
  # Shorten names for the plot labels if needed
  names(sig_list_final) <- names(sig_list_final) %>%
    str_replace_all("_vs_", " vs ")
  
  p_venn <- ggVennDiagram(sig_list_final) +
    scale_fill_gradient(low = "white", high = "red") +
    scale_x_continuous(expand = expansion(mult = .2)) +
    scale_y_continuous(expand = expansion(mult = .1)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("DEG Overlap: ", group_name, " Primary Comparisons"),
      subtitle = "Significant Genes: padj < 0.05 & |log2FC| > 1",
      fill = "Gene Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 20, 10, 20)
    )

  ggsave(paste0(res_dir, "/figures/04_venn_diagram_primary.png"), p_venn,
         width = 10, height = 8, bg = "white")
  message("[OK] Refined Venn diagram saved -> 04_venn_diagram_primary.png",
          "  (", length(sig_list_final), " sets included)")
} else {
  message("Skipping Venn diagram: Fewer than 2 primary contrasts had significant genes.")
}
