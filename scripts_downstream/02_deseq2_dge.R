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
group_name   <- args[1]
species_name <- if (length(args) >= 2) args[2] else str_to_title(group_name)

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
  "../_data/genome/Dog/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
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

# 3. Model Definition
metadata$condition <- relevel(metadata$condition, ref = "DMSO_Romi")
results_list <- list()

# -- Phase 1: Full Multifactorial Model --
message("\n[Phase 1] Running Multifactorial Model (~ cell_line + condition)")
if ("cell_line" %in% colnames(metadata)) {
  ref_cell_line <- if (tolower(species_name) == "human") "H9" else "UL1"
  if (ref_cell_line %in% metadata$cell_line) {
    metadata$cell_line <- factor(metadata$cell_line)
    metadata$cell_line <- relevel(metadata$cell_line, ref = ref_cell_line)
  }
  
  dds_full <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ cell_line + condition)
} else {
  dds_full <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ condition)
}

dds_full <- DESeq(dds_full)

# Comparisons to extract
comparisons <- list(
  c("Romi_6nM",       "DMSO_Romi"),
  c("Kromastat_6nM",  "DMSO_Kromastat"),
  c("Romi_6nM",       "Kromastat_6nM"),
  c("DMSO_Kromastat", "DMSO_Romi")
)

process_results <- function(dds_obj, dose, ref, prefix = "") {
  contrast_name <- paste0(dose, "_vs_", ref)
  full_name <- if(prefix == "") contrast_name else paste0(prefix, "_", contrast_name)
  message("  Extracting: ", full_name)
  
  res <- results(dds_obj, contrast = c("condition", dose, ref))
  
  # Shrinkage
  all_coefs <- resultsNames(dds_obj)
  target_pattern <- paste0("condition_", make.names(dose), "_vs_", make.names(ref))
  coef_name <- all_coefs[grepl(target_pattern, all_coefs, fixed = TRUE)]
  
  if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(dds_obj, coef = coef_name, type = "apeglm")
  } else {
    res_shrunk <- lfcShrink(dds_obj, contrast = c("condition", dose, ref), type = "normal")
  }
  
  # Data Frame with Symbols
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("Geneid") %>%
    mutate(Geneid_clean = sub("\\..*$", "", Geneid)) %>%
    left_join(gene_map, by = c("Geneid_clean" = "gene_id")) %>%
    mutate(gene_label = ifelse(is.na(gene_name) | gene_name == "", Geneid_clean, gene_name)) %>%
    select(Geneid, gene_name, gene_label, everything(), -Geneid_clean) %>%
    arrange(padj)
  
  safe_contrast <- str_replace_all(full_name, "[^a-zA-Z0-9]", "_")
  write.csv(res_df, file = paste0(res_dir, "/tables/02_dge_", safe_contrast, ".csv"), row.names = FALSE)
  
  return(list(res = res, shrunk = res_shrunk, df = res_df))
}

for (comp in comparisons) {
  results_list[[paste0(comp[1], "_vs_", comp[2])]] <- process_results(dds_full, comp[1], comp[2])
}

# -- Phase 2: Subset Models for Consistency Analysis --
if ("cell_line" %in% colnames(metadata)) {
  cell_lines <- unique(as.character(metadata$cell_line))
  for (cl in cell_lines) {
    message("\n[Phase 2] Running Subset Model for Cell Line: ", cl)
    cl_metadata <- metadata[metadata$cell_line == cl, ]
    cl_counts <- counts_filtered[, rownames(cl_metadata)]
    
    # Pre-filtering subset (ensure enough reads in this specific line)
    keep_sub <- rowSums(cl_counts >= 5) >= 3
    cl_counts <- cl_counts[keep_sub, ]
    
    dds_sub <- DESeqDataSetFromMatrix(countData = cl_counts, colData = cl_metadata, design = ~ condition)
    dds_sub <- DESeq(dds_sub)
    
    # Only Romi and Krom primary contrasts for consistency
    for (comp in comparisons[1:2]) {
      results_list[[paste0(cl, "_", comp[1], "_vs_", comp[2])]] <- process_results(dds_sub, comp[1], comp[2], prefix = cl)
    }
  }
}

# Save comprehensive results (using the full dds for downstream PCA/GSEA consistency)
dds <- dds_full 
save(dds, results_list, species_name, group_name, metadata,
     file = paste0("./.RData/", group_name, "/02_deseq_results.RData"))

message("\n[OK] DGE analysis complete for ", group_name)
