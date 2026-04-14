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
metadata$condition <- relevel(metadata$condition, ref = "DMSO_Romi")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata,
                              design = ~ condition)

# 3. Process
dds <- DESeq(dds)

# 4. Extract results for specific comparisons
# Mapping: list(c(Treatment, Control))
comparisons <- list(
  c("Romi_6nM",       "DMSO_Romi"),
  c("Kromastat_6nM",  "DMSO_Kromastat"),
  c("Romi_6nM",       "Kromastat_6nM"),
  c("DMSO_Kromastat", "DMSO_Romi")
)

results_list <- list()

for (comp in comparisons) {
  dose <- comp[1]
  ref  <- comp[2]
  contrast_name <- paste0(dose, "_vs_", ref)
  message("Processing contrast: ", contrast_name)
  
  res <- results(dds, contrast = c("condition", dose, ref))
  
  all_coefs <- resultsNames(dds)
  clean_dose <- make.names(dose)
  clean_ref  <- make.names(ref)
  
  # Try to find exactly matching coefficient for apeglm shrinkage
  target_pattern <- paste0("condition_", clean_dose, "_vs_", clean_ref)
  coef_name <- all_coefs[grepl(target_pattern, all_coefs, fixed = TRUE)]
  
  if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    message("Warning: apeglm coefficient not found exactly for ", contrast_name, ". Falling back to 'normal' shrinkage.")
    res_shrunk <- lfcShrink(dds, contrast = c("condition", dose, ref), type = "normal")
  }
  
  results_list[[contrast_name]] <- list(res = res, shrunk = res_shrunk)
  
  safe_contrast <- str_replace_all(contrast_name, "[^a-zA-Z0-9]", "_")
  write.csv(as.data.frame(res_shrunk), 
            file = paste0(res_dir, "/tables/02_dge_", safe_contrast, ".csv"))
}

# 5. Save comprehensive results
save(dds, results_list, file = paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("DGE analysis complete for ", group_name, ". Results saved to ", res_dir, "/tables/")
