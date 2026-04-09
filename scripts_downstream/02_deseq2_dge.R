# 02_deseq2_dge.R
# Goal: Run Differential Gene Expression analysis with DESeq2 for MiaPAca-2 cancer cells.

library(DESeq2)
library(tidyverse)
library(ggplot2)

# Ensure output directories exist
dir.create("../results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("../results/tables", showWarnings = FALSE, recursive = TRUE)

# 1. Load data
load("./.RData/01_processed_counts.RData")

# 2. Construct DESeqDataSet
# Ensure reference level is DMSO (NG)
metadata$condition <- relevel(metadata$condition, ref = "DMSO (NG)")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata,
                              design = ~ condition)

# 3. Process
dds <- DESeq(dds)

# 4. Extract results for all doses vs DMSO
doses <- c("0.5nM (500pM)", "1nM", "5nM")
results_list <- list()

for (dose in doses) {
  message("Processing contrast: ", dose, " vs DMSO (NG)")
  
  # Extract results
  res <- results(dds, contrast = c("condition", dose, "DMSO (NG)"))
  
  # LFC Shrinkage (using apeglm)
  # Dynamic coefficient name generation
  clean_dose <- make.names(dose)
  clean_ref  <- make.names("DMSO (NG)")
  coef_name  <- paste0("condition_", clean_dose, "_vs_", clean_ref)
  
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # Store and Save
  results_list[[dose]] <- list(res = res, shrunk = res_shrunk)
  
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  write.csv(as.data.frame(res_shrunk), 
            file = paste0("../results/tables/02_dge_", safe_name, "_vs_DMSO.csv"))
}

# 5. Save comprehensive results
save(dds, results_list, file = "./.RData/02_deseq_results.RData")

message("DGE analysis complete. Results for all doses saved to ../results/tables/")
