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

# 4. Venn Diagram (3.1a)
# Extract significant gene sets for each dose
sig_list <- lapply(results_list, function(x) {
  rownames(subset(as.data.frame(x$shrunk), padj < 0.05 & abs(log2FoldChange) > 1))
})

# Filter out empty sets if any
sig_list <- sig_list[sapply(sig_list, length) > 0]

if (length(sig_list) >= 2) {
  p_venn <- ggVennDiagram(sig_list) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(paste0("DEG Overlap: ", group_name, " Treatments vs DMSO"))
  
  ggsave(paste0(res_dir, "/figures/04_venn_diagram_doses.png"), p_venn,
         width = 8, height = 8, bg = "white")
  message("[OK] Venn diagram  -> 04_venn_diagram_doses.png",
          "  (", length(sig_list), " dose sets)")
} else {
  message("Skipping Venn diagram: Fewer than 2 doses with significant genes.")
}
