# [[test_downstream/test_01_data_prep.R]]
# Goal: Load raw counts and prepare metadata for DESeq2 analysis on specific group.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript test_01_data_prep.R <Group> <Species>")
}
group_name <- args[1]
species_name <- args[2]

message("Preparing data for Group: ", group_name, " (", species_name, ")")

# 1. Load Sample Metadata
# The CSV has some header lines we need to skip or filter
metadata_raw <- read_csv("../drPhuong_Sample_Data_Table.csv", skip = 1) %>%
  filter(!is.na(Treatment))

# Filter to only the group requested
metadata_raw <- metadata_raw %>% filter(Group == group_name)

# Create sample names that match the column headers in gene_counts.txt
# Data has File1, File2, Treatment, Cell line, Group
metadata <- metadata_raw %>%
  mutate(Sample = str_replace(File1, "_R1_001.fastq.gz", "")) %>%
  select(Sample, Treatment, Group) %>%
  mutate(condition = factor(Treatment, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM"))) %>%
  mutate(condition = relevel(condition, ref = "DMSO_Romi"))

# 2. Load Gene Counts
# Count file is now group-specific
# --- DEVIATION: Output to test-specific directories
count_file <- paste0("../_data/counts_test/gene_counts_", group_name, ".txt")
raw_counts <- read_delim(count_file, delim = "\t", skip = 1)

# Clean column names to match the Sample column in metadata
# featureCounts uses file paths, we want the sample name
counts_matrix <- raw_counts %>%
  rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col()) %>%
  select(Geneid, any_of(metadata$Sample)) %>%
  column_to_rownames("Geneid")

# Ensure ordering matches
metadata <- metadata %>% filter(Sample %in% colnames(counts_matrix))
counts_matrix <- counts_matrix[, metadata$Sample]

# 3. Pre-filtering
# Keep genes with at least 10 reads in 3 or more samples (one treatment group size)
# --- DEVIATION: Relaxed raw count filtering for small mock dataset (Chr21 subsets yield very low matrices)
keep <- rowSums(counts_matrix >= 1) >= 2
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
# --- DEVIATION: Output to test-specific directories
out_dir <- paste0("./.RData_test/", group_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, metadata, species_name, group_name, file = paste0(out_dir, "/01_processed_counts.RData"))

message("Data preparation complete. Filtered ", nrow(counts_filtered), " genes.")
print(head(metadata))
