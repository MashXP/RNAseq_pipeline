# Downstream Dissection: 01_data_prep.R

`01_data_prep.R` is the engine that processes the raw output from your upstream pipeline (`_data/featurecounts/`). It transforms technical filenames into human-readable sample names and prepares the data for multifactorial analysis.

---

## 0. Data Flow (I/O)
- **Input**:
    - **Count Matrix**: `_data/featurecounts/[Species]/[Species]_featurecounts_counts.tsv` (Raw featureCounts output).
    - **Metadata**: `_data/drPhuong_Sample_Data_Table.csv` (CSV, mapping groups to treatments).
- **Processing**: Filters by Group, creates `display_names`, reorders matrix columns, and filters low-signal genes.
- **Output**: 
    - **Processed RData**: `./.RData/[Group]/01_processed_counts.RData`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for data manipulation. See [**libraries.md**](libraries.md) for full technical justifications.
- `tidyverse`: Essential for string manipulation (`str_replace`, `str_remove`) and piping (`%>%`).

---

## 1. Metadata Parsing
```r
metadata_raw <- read_csv("../drPhuong_Sample_Data_Table.csv", skip = 1) %>%
  filter(!is.na(Treatment)) %>%
  filter(Group == group_name)
```
- **The Job**: Reads your experimental design table and filters it down to a single "Group" (e.g., Human-H9 or Dog-UL1).
- **The Reasoning**: Large-scale analysis is better managed one "Cell Line Group" at a time to ensure high specificity in the results.

---

## 2. Desktop-Quality Naming (Display Names)
```r
metadata <- metadata_raw %>%
  mutate(OriginalSample = str_replace(File1, "_R1_001.fastq.gz", "")) %>%
  group_by(Group, Treatment) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>%
  mutate(display_name = paste0(Group, "_", Treatment, "_", Rep))
```
- **The Job**: Automatically generates a clean name for each sample, like `Human-H9_Kromastat_6nM_1`.
- **The Reasoning**: Filenames like `S12_L001_R1_001` are impossible to read on a final graph. By creating `display_name`, you ensure that every PCA and Volcano plot is perfectly labeled and publication-ready.

---

## 3. Count Matrix Cleaning
```r
counts_matrix <- raw_counts %>%
  rename_with(~basename(.) %>% str_remove("_Aligned.sortedByCoord.out.bam"), 7:last_col()) %>%
  select(Geneid, any_of(metadata$OriginalSample))
```
- **The Job**: 
    1. Strips the long system paths from the column headers of your counts file.
    2. Identifies which BAM files belong to the current experimental group.
- **The Reasoning**: STAR output files have very long, complex headers. This cleanup step "standardizes" the matrix so it can be combined with the metadata effortlessly.

---

## 4. Primary Key Alignment
```r
counts_matrix <- counts_matrix[, metadata$OriginalSample]
colnames(counts_matrix) <- metadata$display_name
```
- **The Job**: Reorders the columns of the count matrix to match the metadata exactly and then swaps the old names for the new `display_names`.
- **Why it matters**: In R, if your matrix columns and metadata rows are out of order, the statistical model will produce total nonsense. This step is the "Safety Alignment" before analysis.

---

## 5. Pre-filtering & Serialization
```r
keep <- rowSums(counts_matrix >= 10) >= 3
counts_filtered <- counts_matrix[keep,]
save(counts_filtered, metadata, species_name, group_name, file = paste0(out_dir, "/01_processed_counts.RData"))
```
- **The Job**: 
    1. Removes low-count genes (fewer than 10 reads in 3+ samples).
    2. Packages everything into the pipeline's internal `.RData` format.
- **The Reasoning**: This ensures that regardless of whether you used the "Bridge" (mentor data) or "Data Prep" (your data), the next step in the pipeline (`02_deseq2_dge.R`) sees exactly the same format.
