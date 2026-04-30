# [[02_dge_analysis/scripts/downstream/01_data_prep.R]]
# Goal: Load raw gene counts and metadata, clean, and save for downstream DESeq2.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript 01_data_prep.R <Species>")
}
species_name <- args[1] # e.g., "Human" or "Dog"
group_name <- tolower(species_name)

# --- Dynamic Pathing Setup ---
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  # Fallback: find root by looking for 'run' script up the directory tree
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) { curr <- dirname(curr) }
  project_root <- curr
}
data_dir <- file.path(project_root, "_data")
results_root <- file.path(project_root, "results")
rdata_dir <- file.path(data_dir, "RData", "02_dge_analysis")
# -----------------------------

res_dir <- file.path(results_root, "02_dge_analysis", group_name)
dir.create(file.path(res_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load Global Metadata
metadata_raw <- read_csv(file.path(data_dir, "_metadata", "drPhuong_Sample_Data_Table.csv"), show_col_types = FALSE)

# Rename and filter for the requested species (which is in the 'Cell line' column)
metadata <- metadata_raw %>%
  rename(species_col = `Cell line`, cell_line = Group, condition = Treatment) %>%
  filter(species_col == species_name)

if (nrow(metadata) == 0) {
    stop("No samples found for species: ", species_name)
}

# Descriptive naming and formatting
metadata <- metadata %>%
  mutate(OriginalSample = str_replace(File1, "_R1_001.fastq.gz", "")) %>%
  group_by(cell_line, condition) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>%
  mutate(display_name = paste0(cell_line, "_", condition, "_", Rep)) %>%
  mutate(condition = factor(condition, levels = c("DMSO_Romidepsin", "Romidepsin_6nM", "DMSO_Kromastat", "Kromastat_6nM"))) %>%
  select(OriginalSample, display_name, condition, cell_line)

# 2. Load Gene Counts
sp_lower <- tolower(species_name)
count_file <- file.path(data_dir, "featurecounts", sp_lower, paste0(sp_lower, "_featurecounts_counts.tsv"))
meta_species_file <- file.path(data_dir, "featurecounts", sp_lower, paste0(sp_lower, "_sample_metadata.tsv"))

if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)
if (!file.exists(meta_species_file)) stop("Species metadata not found: ", meta_species_file)

raw_counts <- read_delim(count_file, delim = "\t", skip = 1, show_col_types = FALSE)

# Extract counts and GeneID
counts_matrix <- raw_counts %>%
  select(Geneid, starts_with("/mnt"), starts_with("/home")) # featureCounts headers are absolute paths

# Clean column names to just the sample ID (e.g., 27_S27)
colnames(counts_matrix) <- sub(".*_([0-9]+_S[0-9]+)_.*", "\\1", colnames(counts_matrix))
colnames(counts_matrix)[1] <- "GeneID"

# 2.1 Match Metadata with Count Matrix
# The 'OriginalSample' in metadata should match the column names in counts_matrix
# e.g., metadata$OriginalSample = "27_S27", colnames(counts_matrix) = "27_S27"

# Get sample IDs present in both
common_samples <- intersect(metadata$OriginalSample, colnames(counts_matrix))
metadata <- metadata %>% filter(OriginalSample %in% common_samples)
counts_matrix <- counts_matrix %>% select(GeneID, all_of(common_samples))

# Convert to matrix
counts_matrix <- counts_matrix %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

# Add a match_id to metadata to ensure alignment
metadata$match_id <- match(metadata$OriginalSample, colnames(counts_matrix))
metadata <- metadata %>% arrange(match_id)

# Ensure column ordering matches metadata row ordering
counts_matrix <- counts_matrix[, metadata$match_id]

# Rename columns to display_name and set metadata rownames
colnames(counts_matrix) <- metadata$display_name
metadata <- metadata %>%
  select(-OriginalSample, -match_id) %>%
  rename(Sample = display_name) %>%
  column_to_rownames("Sample")

# 3. Pre-filtering
# Keep genes with at least 10 reads in 3 or more samples
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
counts_raw <- counts_matrix # Keep unfiltered for subset-specific models
rdata_module_dir <- file.path(rdata_dir, group_name)
dir.create(rdata_module_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, counts_raw, metadata, species_name, group_name, file = file.path(rdata_module_dir, "01_processed_counts.RData"))

message("Data preparation complete for ", group_name)
