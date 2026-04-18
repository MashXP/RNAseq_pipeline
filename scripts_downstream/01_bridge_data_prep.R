# [[scripts_downstream/00_bridge.R]]
# Goal: Bridge mentor's consolidated species-level featureCounts matrices 
#       into the .RData format used by the user's downstream pipeline.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 00_bridge.R <Species> (e.g., human or dog)")
}

species_name <- args[1]
group_name <- species_name # Using Species as the Group for multifactorial analysis
message("Bridging mentor data for Species: ", species_name)

# 1. Define Paths
base_dir <- paste0("cancer_downstream/featurecounts/", species_name)
count_file <- file.path(base_dir, paste0(species_name, "_featurecounts_counts.tsv"))
meta_file  <- file.path(base_dir, paste0(species_name, "_sample_metadata.tsv"))

if (!file.exists(count_file) || !file.exists(meta_file)) {
    stop("Input files not found in ", base_dir)
}

# 2. Load Metadata
metadata_raw <- read_delim(meta_file, delim = "\t")

# Mentor's metadata has 'treatment' and 'cell_line'.
# We need to ensure 'condition' exists for compatibility with user pipeline.
metadata <- metadata_raw %>%
  rename(condition = treatment) %>%
  mutate(condition = factor(condition, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM"))) %>%
  # set sample_id as primary key
  column_to_rownames("sample_id")

# 3. Load Gene Counts
# featureCounts file: skip 1st line, Geneid is column 1
raw_counts <- read_delim(count_file, delim = "\t", skip = 1)

# Extract matrix: use Geneid as rownames, samples are in columns 7+
counts_matrix <- raw_counts %>%
  select(Geneid, all_of(rownames(metadata))) %>%
  column_to_rownames("Geneid")

# Verify alignment
if (!all(colnames(counts_matrix) == rownames(metadata))) {
    stop("Sample mismatch between counts matrix and metadata!")
}

# 4. Pre-filtering
# User logic: Keep genes with at least 10 reads in 3 or more samples
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 5. Save processed objects
out_dir <- paste0("scripts_downstream/.RData/", group_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save(counts_filtered, metadata, species_name, group_name, file = file.path(out_dir, "01_processed_counts.RData"))

message("Bridge complete for ", group_name)
message("Filtered genes: ", nrow(counts_filtered))
message("Total samples: ", ncol(counts_filtered))
print(head(metadata))
