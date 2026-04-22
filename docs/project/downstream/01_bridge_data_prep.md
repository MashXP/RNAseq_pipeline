# Downstream Dissection: 01_bridge_data_prep.R

The "Bridge" script is the transition point between your mentor's workflow and your own. It takes pre-consolidated species-level matrices and reformats them so your analysis suite can recognize them.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Count Matrix**: `cancer_downstream/featurecounts/[Species]/[Species]_featurecounts_counts.tsv` (Standard `featureCounts` format).
    - **Metadata**: `cancer_downstream/featurecounts/[Species]/[Species]_sample_metadata.tsv` (Tab-delimited, mapping `sample_id` to `treatment` and `cell_line`).
- **Processing**: Renames `treatment` to `condition`, converts to factors, extracts raw counts, and filters low-signal genes (10+ reads in 3+ samples).
- **Output**: 
    - **Processed RData**: `scripts_downstream/.RData/[Species]/01_processed_counts.RData` (Contains filtered counts and metadata).

---

## 0.1 Library Rationales
This script utilizes the following libraries for data manipulation. See [**libraries.md**](libraries.md) for full technical justifications.
- `tidyverse`: Essential for piping (`%>%`) and data cleaning (`dplyr`, `readr`).

---

## 1. Parameters & Environment
```r
args <- commandArgs(trailingOnly = TRUE)
species_name <- args[1]
group_name <- species_name 
```
- **The Job**: Detects which species (Human or Dog) to process based on the command-line argument.
- **The Reasoning**: Instead of writing separate scripts for each species, we use one robust script. We set `group_name = species_name` to unify the data structure for the multifactorial analysis that follows in step 02.

---

## 2. Path Construction
```r
base_dir <- paste0("../cancer_downstream/featurecounts/", species_name)
count_file <- file.path(base_dir, paste0(species_name, "_featurecounts_counts.tsv"))
meta_file  <- file.path(base_dir, paste0(species_name, "_sample_metadata.tsv"))
```
- **The Job**: Automatically builds the path to the mentor's input data.
- **The Reasoning**: The mentor provided files in a specific folder structure (`featurecounts/human/` and `featurecounts/dog/`). This code ensures the script "finds" the right data without manual intervention.

---

## 3. Metadata Re-alignment
```r
metadata <- metadata_raw %>%
  rename(condition = treatment) %>%
  mutate(condition = factor(condition, levels = c("DMSO_Romi", "Romi_6nM", "DMSO_Kromastat", "Kromastat_6nM")))
```
- **The Job**: 
    1. Renames `treatment` to `condition` to match your pipeline's internal logic.
    2. Converts conditions into "Factors" with a fixed order.
- **Why it matters**: In R statistics, "Factors" tell the computer the order of comparison. By setting the levels explicitly, you ensure that "DMSO" (the control) is always the reference point for the drug treatments.

---

## 4. Count Matrix Extraction
```r
raw_counts <- read_delim(count_file, delim = "\t", skip = 1)
counts_matrix <- raw_counts %>%
  select(Geneid, all_of(rownames(metadata))) %>%
  column_to_rownames("Geneid")
```
- **The Job**: 
    1. Skips the first line of the featureCounts file (the tool version log).
    2. Selects only the columns that have matching entries in your metadata.
- **The Reasoning**: This creates a clean "Gene vs. Sample" matrix, discarding the extra annotation columns (Chr, Start, End) that are no longer needed for statistics.

---

## 5. Pre-filtering (Noise Removal)
```r
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]
```
- **The Job**: Keeps only genes that have at least 10 reads in 3 or more samples.
- **The Reasoning**: RNA-seq results for low-count genes are mathematically unreliable (noisy). This "pre-filtering" step removes thousands of irrelevant genes and improves the statistical power of your DESeq2 analysis.

---

## 6. Pipeline Bridge (`.RData`)
```r
out_dir <- paste0("./.RData/", group_name)
save(counts_filtered, metadata, species_name, group_name, file = file.path(out_dir, "01_processed_counts.RData"))
```
- **The Job**: Saves all the cleaned R objects into a single compressed binary file.
- **The Reasoning**: This is the "Baton Pass." The next script (`01_data_prep.R`) will load this `.RData` file instead of recalculating everything from scratch.
