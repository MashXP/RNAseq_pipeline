# [[scripts_downstream/01_data_prep.R]]
# Goal: Load raw counts and prepare metadata for DESeq2 analysis on MiaPAca-2 data.

library(tidyverse)

# 1. Load Sample Metadata
# The CSV has some header lines we need to skip or filter
metadata_raw <- read_csv("../MiaPAca-2_Sample_Data_Table.csv", skip = 1) %>%
  filter(!is.na(Conc.))
# Create sample names that match the column headers in gene_counts.txt
# (Removing the R1 suffix just like we did in the alignment script)
metadata <- metadata_raw %>%
  mutate(Sample = str_replace(File_1, "_R1_001.fastq.gz", "")) %>%
  select(Sample, Conc., Replicate) %>%
  mutate(condition = factor(Conc., levels = c("DMSO (NG)", "0.5nM (500pM)", "1nM", "5nM"))) %>%
  mutate(condition = relevel(condition, ref = "DMSO (NG)"))

# 2. Load Gene Counts
# Skip the first comment line (# program:featureCounts...)
raw_counts <- read_delim("../_data/counts/gene_counts.txt", 
                         delim = "\t", skip = 1)

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
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]

# 4. Save processed objects
dir.create(".RData", showWarnings = FALSE)
save(counts_filtered, metadata, file = "./.RData/01_processed_counts.RData")

message("Data preparation complete. Filtered ", nrow(counts_filtered), " genes.")
print(head(metadata))
