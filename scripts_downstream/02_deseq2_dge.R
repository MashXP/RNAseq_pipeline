# [[scripts_downstream/02_deseq2_dge.R]]
# Goal: Run Differential Gene Expression analysis.
# Phase 1: Multifactorial Species-level model (~ cell_line + condition).
# Phase 2: Cell-line specific models for consistency checking (UpSet plots).

library(DESeq2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 02_deseq2_dge.R <Group> [Species]")
}
group_name   <- tolower(args[1]) # Standardize to lowercase for path consistency
species_name <- if (length(args) >= 2) args[2] else str_to_title(args[1])

# Ensure output directories exist
res_dir <- paste0("../results/", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(res_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

message("--- Running DGE for Group: ", group_name, " (", species_name, ") ---")

# 1. Load data
load(paste0("./.RData/", group_name, "/01_processed_counts.RData"))

# 2. Gene Mapping Setup
extract_gtf_attribute <- function(attr_text, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  value <- sub(paste0(".*", pattern, ".*"), "\\1", attr_text)
  value[!grepl(pattern, attr_text)] <- NA_character_
  value
}

gtf_path <- if (tolower(species_name) == "human") {
  "../_data/genome/Human/Homo_sapiens.GRCh38.113.gtf"
} else {
  "../_data/genome/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
}

gene_map <- data.frame(gene_id = character(), gene_name = character())
if (file.exists(gtf_path)) {
  message("Building gene map from GTF: ", gtf_path)
  map_file <- tempfile(fileext = ".csv")
  awk_cmd <- paste0("awk -F'\\t' '$3==\"gene\" { ",
                    "match($9, /gene_id \"([^\"]+)\"/, id); ",
                    "match($9, /gene_name \"([^\"]+)\"/, name); ",
                    "print id[1] \"\\t\" name[1] }' ", gtf_path, " > ", map_file)
  system(awk_cmd)
  gene_map <- read_delim(map_file, delim = "\t", col_names = c("gene_id", "gene_name"), 
                         show_col_types = FALSE) %>%
    distinct(gene_id, .keep_all = TRUE)
  unlink(map_file)
}

# -- Comparison Engine: Match Mentor's Subset Logic --
# Logic: 
# 1. Species-level Drug vs DMSO: Subset to 12 samples (2 cell lines x 2 treatments x 3 reps).
# 2. Cell-line specific Drug vs DMSO: Subset to 6 samples (1 cell line x 2 treatments x 3 reps).

results_list <- list()

# Primary Contrasts (Species-wide)
# Each runs its own independent model (~ cell_line + condition)
drug_contrasts <- list(
  list(name = "Romidepsin_6nM_vs_DMSO_Romidepsin",       levels = c("Romidepsin_6nM", "DMSO_Romidepsin")),
  list(name = "Kromastat_6nM_vs_DMSO_Kromastat",  levels = c("Kromastat_6nM", "DMSO_Kromastat")),
  list(name = "Romidepsin_6nM_vs_Kromastat_6nM",   levels = c("Romidepsin_6nM", "Kromastat_6nM"))
)

for (drug_comp in drug_contrasts) {
  message("\n[Phase 1] Running Species-wide Model: ", drug_comp$name)
  
  # Subset to only the relevant treatments for this drug comparison
  sub_metadata <- metadata %>% filter(condition %in% drug_comp$levels)
  sub_counts <- counts_raw[, rownames(sub_metadata)]
  
  # Independent filtering and normalization
  keep_sub <- rowSums(sub_counts >= 10) >= 3
  sub_counts <- sub_counts[keep_sub, ]
  
  sub_metadata$condition <- factor(sub_metadata$condition, levels = rev(drug_comp$levels)) # Ref is last
  
  dds_sub <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_metadata, design = ~ cell_line + condition)
  dds_sub <- DESeq(dds_sub)
  
  # Extract with apeglm shrinkage (as it's species-wide and more robust)
  res <- results(dds_sub)
  res_shrunk <- lfcShrink(dds_sub, coef = resultsNames(dds_sub)[3], type = "apeglm")
  
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("Geneid") %>%
    mutate(Geneid_clean = sub("\\..*$", "", Geneid)) %>%
    left_join(gene_map, by = c("Geneid_clean" = "gene_id")) %>%
    mutate(gene_label = ifelse(is.na(gene_name) | gene_name == "", Geneid_clean, gene_name)) %>%
    select(Geneid, gene_name, gene_label, everything(), -Geneid_clean) %>%
    arrange(padj)
  
  write.csv(res_df, file = paste0(res_dir, "/tables/02_dge_", drug_comp$name, ".csv"), row.names = FALSE)
  results_list[[drug_comp$name]] <- list(res = res, shrunk = res_shrunk, df = res_df)
}

# -- Phase 2: Independent Cell-Line Specific Models --
# Runs ~ condition for every (CellLine x Drug) combination (6 samples each)
if ("cell_line" %in% colnames(metadata)) {
  cell_lines <- unique(as.character(metadata$cell_line))
  for (cl in cell_lines) {
    for (drug_comp in drug_contrasts) {
      full_name <- paste0(cl, "_", drug_comp$name)
      message("\n[Phase 2] Running Independent Subset Model: ", full_name)
      
      cl_metadata <- metadata %>% filter(cell_line == cl, condition %in% drug_comp$levels)
      cl_counts <- counts_raw[, rownames(cl_metadata)]
      
      keep_cl <- rowSums(cl_counts >= 10) >= 3
      cl_counts <- cl_counts[keep_cl, ]
      
      cl_metadata$condition <- factor(cl_metadata$condition, levels = rev(drug_comp$levels))
      
      dds_cl <- DESeqDataSetFromMatrix(countData = cl_counts, colData = cl_metadata, design = ~ condition)
      dds_cl <- DESeq(dds_cl)
      
      # Use 'normal' shrinkage for legacy cell-line subsets to match mentor
      res_cl <- results(dds_cl)
      res_cl_shrunk <- lfcShrink(dds_cl, contrast = c("condition", drug_comp$levels[1], drug_comp$levels[2]), res = res_cl, type = "normal")
      
      res_cl_df <- as.data.frame(res_cl_shrunk) %>%
        rownames_to_column("Geneid") %>%
        mutate(Geneid_clean = sub("\\..*$", "", Geneid)) %>%
        left_join(gene_map, by = c("Geneid_clean" = "gene_id")) %>%
        mutate(gene_label = ifelse(is.na(gene_name) | gene_name == "", Geneid_clean, gene_name)) %>%
        select(Geneid, gene_name, gene_label, everything(), -Geneid_clean) %>%
        arrange(padj)
      
      safe_contrast <- str_replace_all(full_name, "[^a-zA-Z0-9]", "_")
      write.csv(res_cl_df, file = paste0(res_dir, "/tables/02_dge_", safe_contrast, ".csv"), row.names = FALSE)
      results_list[[full_name]] <- list(res = res_cl, shrunk = res_cl_shrunk, df = res_cl_df)
    }
  }
}

# Save comprehensive results
# Note: dds_full is removed as results are now the composite of subsets. 
# Re-running a dummy dds_full for PCA/GSEA downstream compatibility if needed, but results_list is the primary.
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ cell_line + condition)
dds <- estimateSizeFactors(dds) # For normalized counts downstream
save(dds, results_list, species_name, group_name, metadata,
     file = paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("\n[OK] DGE analysis complete for ", group_name)
