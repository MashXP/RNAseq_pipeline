# [[test_downstream/02_deseq2_dge.R]]
# TEST VERSION: Differential Gene Expression.
# Logic: Runs DESeq2 and generates results_list (matches production structure).
# Divergence: RELAXED P-VALUE (pvalue < 0.5) to ensure test visualization visibility.

library(DESeq2)
library(tidyverse)

# 1. Load data
load("./.RData/01_test_processed.RData")

dir.create("tables_test", showWarnings = FALSE)

# 2. Construct DESeqDataSet
metadata$condition <- relevel(metadata$condition, ref = "DMSO (NG)")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata,
                              design = ~ condition)

# 3. Process
dds <- DESeq(dds)

# 4. Extract results (Loop through doses like production)
doses <- c("0.5nM (500pM)", "1nM", "5nM")
results_list <- list()

for (dose in doses) {
  message("[TEST] Processing contrast: ", dose, " vs DMSO (NG)")
  
  res <- results(dds, contrast = c("condition", dose, "DMSO (NG)"))
  
  # Try shrinkage (standard type for test to avoid dependency issues)
  res_shrunk <- lfcShrink(dds, contrast = c("condition", dose, "DMSO (NG)"), type = "normal")
  
  results_list[[dose]] <- list(res = res, shrunk = res_shrunk)
  
  safe_name <- str_replace_all(dose, "[^a-zA-Z0-9]", "_")
  write.csv(as.data.frame(res_shrunk), 
            file = paste0("tables_test/02_test_dge_", safe_name, ".csv"))
}

# 5. Save comprehensive results (using production object names)
save(dds, results_list, file = "./.RData/02_test_deseq_results.RData")

message("[TEST] DGE analysis complete.")
