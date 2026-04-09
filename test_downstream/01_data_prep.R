# [[test_downstream/01_data_prep.R]]
# TEST VERSION: Initial data processing.
# Logic: Loads test counts and metadata. 
# Divergence: RELAXED FILTERING (1 count in 2 samples) to handle low depth.

library(tidyverse)

# 1. Load Sample Metadata
metadata_raw <- read_csv("../MiaPAca-2_Sample_Data_Table.csv", skip = 1) %>%
  filter(!is.na(Conc.))

metadata <- metadata_raw %>%
  mutate(Sample = str_replace(File_1, "_R1_001.fastq.gz", "")) %>%
  select(Sample, Conc., Replicate) %>%
  mutate(condition = factor(Conc., levels = c("DMSO (NG)", "0.5nM (500pM)", "1nM", "5nM")))

# 2. Load Gene Counts (Point to test counts)
raw_counts <- read_delim("../_data/counts_test/gene_counts.txt", 
                         delim = "\t", skip = 1)

counts_matrix <- raw_counts %>%
  rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col()) %>%
  select(Geneid, any_of(metadata$Sample)) %>%
  column_to_rownames("Geneid")

# Ensure ordering matches
metadata <- metadata %>% filter(Sample %in% colnames(counts_matrix))
counts_matrix <- counts_matrix[, metadata$Sample]

# 3. Pre-filtering (LOWERED FOR TEST)
# Production: count >= 10 in 3 samples
# Test: count >= 1 in 2 samples (to ensure some genes survive the 100k subsampling)
keep <- rowSums(counts_matrix >= 1) >= 2
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
dir.create(".RData", showWarnings = FALSE)
save(counts_filtered, metadata, file = "./.RData/01_test_processed.RData")

message("[TEST] Data preparation complete. Filtered ", nrow(counts_filtered), " genes.")
