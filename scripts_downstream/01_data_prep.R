# [[scripts_downstream/01_data_prep.R]]
# Goal: Load raw species-level counts and prepare metadata for multifactorial DESeq2 analysis.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_data_prep.R <Species>")
}
species_name <- str_to_title(args[1]) # Force "Human" or "Dog" for CSV filtering
group_name <- tolower(args[1])      # Force "human" or "dog" for path consistency

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
sp_lower <- tolower(species_name)
count_file <- file.path("../_data/featurecounts", sp_lower, paste0(sp_lower, "_featurecounts_counts.tsv"))
meta_species_file <- file.path("../_data/featurecounts", sp_lower, paste0(sp_lower, "_sample_metadata.tsv"))

if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)
if (!file.exists(meta_species_file)) stop("Species metadata not found: ", meta_species_file)

raw_counts <- read_delim(count_file, delim = "\t", skip = 1)
meta_species <- read_delim(meta_species_file, delim = "\t") %>%
  mutate(sample_number = as.character(sample_number))

# Map OriginalSample numbers to descriptive sample_ids from species metadata
metadata <- metadata %>%
  left_join(meta_species %>% select(sample_number, sample_id), by = c("OriginalSample" = "sample_number"))

# Determine which column names to use
# If count matrix columns match BAM paths, clean them; if they match sample_id, use as is.
if (any(grepl("/", colnames(raw_counts)))) {
  message("  Detected BAM paths in count matrix columns. Cleaning...")
  counts_matrix <- raw_counts %>%
    rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col())
} else {
  message("  Detected descriptive names in count matrix columns.")
  counts_matrix <- raw_counts
}

# Subset and order
# Note: if mapping failed (sample_id is NA), fallback to OriginalSample
metadata <- metadata %>%
  mutate(match_id = ifelse(is.na(sample_id), OriginalSample, sample_id))

counts_matrix <- counts_matrix %>%
  select(Geneid, any_of(metadata$match_id)) %>%
  column_to_rownames("Geneid")

# Verify sample count
if (ncol(counts_matrix) == 0) {
    stop("No columns in count matrix matched metadata IDs!")
}

# Ensure column ordering matches metadata row ordering
counts_matrix <- counts_matrix[, metadata$match_id]

# Rename columns to display_name and set metadata rownames
colnames(counts_matrix) <- metadata$display_name
metadata <- metadata %>%
  select(-OriginalSample, -sample_id, -match_id) %>%
  rename(Sample = display_name) %>%
  column_to_rownames("Sample")

# 3. Pre-filtering
# Keep genes with at least 10 reads in 3 or more samples
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
counts_raw <- counts_matrix # Keep unfiltered for subset-specific models
out_dir <- paste0("./.RData/", group_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, counts_raw, metadata, species_name, group_name, file = file.path(out_dir, "01_processed_counts.RData"))

message("Data preparation complete for ", group_name)
message("Filtered genes: ", nrow(counts_filtered))
message("Total samples: ", ncol(counts_filtered))
print(head(metadata))
