# Downstream Dissection: 02_deseq2_dge.R

This is the "Brain" of the pipeline. It implements the complex statistical models required to compare Kromastat and Romidepsin across multiple cell lines while accounting for biological variation.

---

## 0. Data Flow (I/O)
- **Input**: 
    - **Processed RData**: `./.RData/[Group]/01_processed_counts.RData` (Counts and metadata).
    - **Genome Reference**: `../_data/genome/[Species]/[Files].gtf` (For gene symbology).
- **Processing**: Build gene map, run multifactorial DESeq2 (`~cell_line + condition`), run subset models for cell lines, apply APEGLM shrinkage.
- **Output**: 
    - **Results Tables**: `../results/[Group]/tables/02_dge_[Contrast].csv` (Full statistics + symbols).
    - **Results RData**: `./.RData/[Group]/02_deseq_results.RData`.

---

## 0.1 Library Rationales
This script utilizes the following libraries for statistics. See [**libraries.md**](libraries.md) for full technical justifications.
- `DESeq2`: Golden Standard for differential expression stats.
- `tidyverse`: Data wrangling and piping.

---

## 1. High-Performance Gene Mapping
```r
awk_cmd <- paste0("awk -F'\\t' '$3==\"gene\" { ",
                  "match($9, /gene_id \"([^\"]+)\"/, id); ",
                  "match($9, /gene_name \"([^\"]+)\"/, name); ",
                  "print id[1] \"\\t\" name[1] }' ", gtf_path, " > ", map_file)
system(awk_cmd)
```
- **The Job**: Uses the high-speed Linux tool `awk` to extract every Gene ID and its human-friendly Symbol (Gene Name) from the massive GTF genome file.
- **The Reasoning**: Standard R methods for reading GTF files are extremely slow. This "hybrid" approach builds a translator between Ensembl IDs (ENSG...) and Symbols (e.g., TP53) in seconds, ensuring your final tables are readable.

---

## 2. Setting the Baseline (Reference Levels)
```r
metadata$condition <- relevel(metadata$condition, ref = "DMSO_Romi")
...
ref_cell_line <- if (tolower(species_name) == "human") "H9" else "UL1"
metadata$cell_line <- relevel(metadata$cell_line, ref = ref_cell_line)
```
- **The Job**: Tells the statistical model which groups are the "Healthy Controls."
- **The Reasoning**: Without this, R picks references alphabetically. By explicitly setting `DMSO_Romi` and the healthy lines (`H9`/`UL1`) as references, your Log2FoldChange results represent the *change* caused by the drug relative to the normal state.

---

## 3. The Multifactorial Model
```r
for (drug_comp in drug_contrasts) {
  # Subset and build model
  dds_sub <- DESeqDataSetFromMatrix(..., design = ~ cell_line + condition)
  dds_sub <- DESeq(dds_sub)
  
  # Extract with apeglm shrinkage
  res_shrunk <- lfcShrink(dds_sub, coef = resultsNames(dds_sub)[3], type = "apeglm")
}
```
- **The Job**: Runs a model that looks at `condition` (the drug) while "blocking" for the differences between `cell_lines`.
- **Why it matters**: This aligns perfectly with your mentor's vision. It allows us to find drug effects that are **conserved** across both cell lines, giving the study more statistical power than looking at them separately. It uses **APEGLM** shrinkage for species-wide results to ensure maximum rigor.

---

## 4. Accurate Fold-Change (APEGLM Shrinkage)
- **The Job**: Applies a sophisticated mathematical "shrinkage" to the Fold Change results using either `apeglm` (for species-wide models) or `normal` (for cell-line subsets).
- **The Reasoning**: In RNA-seq, genes with low counts often show huge, "fake" fold changes due to noise. Shrinkage "shrinks" these noisy results toward zero while keeping high-confidence results intact, preventing "false positives" in your Volcano plots.

---

## 5. Automated Contrasts (Comparisons)
```r
comparisons <- list(
  c("Romi_6nM",       "DMSO_Romi"),
  c("Kromastat_6nM",  "DMSO_Kromastat"),
  c("Romi_6nM",       "Kromastat_6nM"),
  c("DMSO_Kromastat", "DMSO_Romi")
)
```
- **The Job**: Loops through every biological question you need to answer.
- **The Reasoning**: 
    - `Romi vs DMSO`: Does Romi work?
    - `Krom vs DMSO`: Does Krom work?
    - `Romi vs Krom`: Is your new drug (Krom) better than the old one (Romi)?

---

## 6. Phase 2: Consistency Checking (Subsets)
```r
cl_metadata <- metadata[metadata$cell_line == cl, ]
...
dds_sub <- DESeqDataSetFromMatrix(..., design = ~ condition)
```
- **The Job**: Re-runs the analysis on *each cell line individually*.
- **The Reasoning**: This is the "Rigor Check." We want to see if the drug works the same way in the aggressive cancer line as it does in the healthy line. These results are saved separately for use in **UpSet plots** later to check for conserved mechanisms.
