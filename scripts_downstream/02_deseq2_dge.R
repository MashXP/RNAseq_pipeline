# [[scripts_downstream/02_deseq2_dge.R]]
# Goal: Run Differential Gene Expression analysis with DESeq2 for specific group.

library(DESeq2)
library(tidyverse)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_deseq2_dge.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

# Ensure output directories exist
res_dir <- paste0("../results/", group_name)
dir.create(paste0(res_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(res_dir, "/tables"), showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load(paste0("./.RData/", group_name, "/01_processed_counts.RData"))

# 2. Construct DESeqDataSet
metadata$condition <- relevel(metadata$condition, ref = "DMSO")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata,
                              design = ~ condition)

# 3. Process
dds <- DESeq(dds)

# 4. Extract results for all doses vs DMSO
doses <- c("Romi (6 nM)", "GNL", "NanoRomi (6 nM)")
results_list <- list()

for (dose in doses) {
  message("Processing contrast: ", dose, " vs DMSO")
  
  res <- results(dds, contrast = c("condition", dose, "DMSO"))
  
  all_coefs <- resultsNames(dds)
  clean_dose <- make.names(dose)
  clean_ref  <- make.names("DMSO")
  
  target_pattern <- paste0("condition_", clean_dose, "_vs_", clean_ref)
  coef_name <- all_coefs[grepl(target_pattern, all_coefs, fixed = TRUE)]
  
  if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    message("Warning: apeglm coefficient not found exactly. Falling back to 'normal' shrinkage.")
    res_shrunk <- lfcShrink(dds, contrast = c("condition", dose, "DMSO"), type = "normal")
  }
  
  results_list[[dose]] <- list(res = res, shrunk = res_shrunk)
  
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  write.csv(as.data.frame(res_shrunk), 
            file = paste0(res_dir, "/tables/02_dge_", safe_name, "_vs_DMSO.csv"))
}

# 5. Save comprehensive results
save(dds, results_list, file = paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("DGE analysis complete for ", group_name, ". Results saved to ", res_dir, "/tables/")
