# [[scripts_downstream/01_data_prep.R]]
# Goal: Load raw species-level counts and prepare metadata for multifactorial DESeq2 analysis.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_data_prep.R <Species>")
}
species_name <- args[1]
group_name <- species_name # Use species as the master group for multifactorial analysis

message("Preparing downstream data for Species: ", species_name)

# 1. Load Metadata
# The CSV has some header lines we need to skip.
# Columns: File1, File2, Treatment, "Cell line", Group
metadata_raw <- read_csv("../drPhuong_Sample_Data_Table.csv", skip = 1) %>%
  filter(!is.na(Treatment))

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
  mutate(condition = factor(condition, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM"))) %>%
  select(OriginalSample, display_name, condition, cell_line)

# 2. Load Gene Counts
count_file <- file.path("../_data/featurecounts", species_name, paste0(species_name, "_featurecounts_counts.tsv"))
if (!file.exists(count_file)) {
    stop("Count matrix not found: ", count_file)
}
raw_counts <- read_delim(count_file, delim = "\t", skip = 1)

# Clean column names (strip BAM paths) and subset to species samples
counts_matrix <- raw_counts %>%
  rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col()) %>%
  select(Geneid, any_of(metadata$OriginalSample)) %>%
  column_to_rownames("Geneid")

# Verify sample count
if (ncol(counts_matrix) == 0) {
    stop("No columns in count matrix matched metadata OriginalSample IDs!")
}

# Ensure column ordering matches metadata row ordering
counts_matrix <- counts_matrix[, metadata$OriginalSample]

# Rename columns to display_name and set metadata rownames
colnames(counts_matrix) <- metadata$display_name
metadata <- metadata %>%
  select(-OriginalSample) %>%
  rename(Sample = display_name) %>%
  column_to_rownames("Sample")

# 3. Pre-filtering
# Keep genes with at least 10 reads in 3 or more samples
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
out_dir <- paste0("./.RData/", group_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, metadata, species_name, group_name, file = file.path(out_dir, "01_processed_counts.RData"))

message("Data preparation complete for ", group_name)
message("Filtered genes: ", nrow(counts_filtered))
message("Total samples: ", ncol(counts_filtered))
print(head(metadata))
