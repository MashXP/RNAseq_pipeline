# [[scripts_downstream/01_data_prep.R]]
# Goal: Load raw counts and prepare metadata for DESeq2 analysis on specific group.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_data_prep.R <Group> <Species>")
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

# 2. Rename Samples Descriptively
# Logic: create a display_name as [Group]_[Treatment]_[Rep]
metadata <- metadata_raw %>%
  mutate(OriginalSample = str_replace(File1, "_R1_001.fastq.gz", "")) %>%
  group_by(Group, Treatment) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>%
  mutate(display_name = paste0(Group, "_", Treatment, "_", Rep)) %>%
  select(OriginalSample, display_name, Treatment, Group) %>%
  mutate(condition = factor(Treatment, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM")))

# 2. Load Gene Counts
# Count file is now group-specific
count_file <- paste0("../_data/counts/gene_counts_", group_name, ".txt")
raw_counts <- read_delim(count_file, delim = "\t", skip = 1)

# Clean column names and switch to descriptive display_names
counts_matrix <- raw_counts %>%
  rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col()) %>%
  select(Geneid, any_of(metadata$OriginalSample)) %>%
  column_to_rownames("Geneid")

# Ensure ordering and rename to display_name
counts_matrix <- counts_matrix[, metadata$OriginalSample]
colnames(counts_matrix) <- metadata$display_name

# Update metadata to use display_name as the primary key
metadata <- metadata %>%
  select(-OriginalSample) %>%
  rename(Sample = display_name) %>%
  column_to_rownames("Sample")

# 3. Pre-filtering
# Keep genes with at least 10 reads in 3 or more samples (one treatment group size)
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
out_dir <- paste0("./.RData/", group_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, metadata, species_name, group_name, file = paste0(out_dir, "/01_processed_counts.RData"))

message("Data preparation complete. Filtered ", nrow(counts_filtered), " genes.")
print(head(metadata))
