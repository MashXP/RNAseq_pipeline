# Utility Scripts Dissection

These utility scripts provide connective tissue for both upstream metadata parsing and downstream biological filtering.

---

## 1. `parse_samples.py`
```python
with open(csv_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    next(reader) # Skip Header
    for row in reader:
        ...
        print(f"{sample_name} {r1} {r2} {species} {group}")
```
- **The Job**: This is the "Source of Truth" for your entire pipeline. It reads your `_data/_metadata/drPhuong_Sample_Data_Table.csv` and translates each row into a single line of text.
- **Why it matters**: 
    - **Decoupling**: It keeps your Bash scripts clean. Instead of complex CSV parsing in Bash, the script just says "Hey Python, tell me what samples to run."
    - **Reliability**: By printing a clean, space-separated string, it ensures the Bash `while read` loop never skips a file or misinterprets a species name.

---

## 2. `biotype_to_multiqc.py`
```python
df = pd.read_csv(input_file, sep='\t', skiprows=1, index_col=0)
mqc_df = df.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
mqc_df = mqc_df.transpose()
...
f.write("# section_name: 'Biotype Counts'\n")
f.write("# plot_type: 'bargraph'\n")
```
- **The Job**: This script reformats the raw counts into a layout that **MultiQC** recognizes as a custom chart.
- **The Reasoning**: 
    - **Visual Audit**: Standard featureCounts output is a giant text table that MultiQC cannot convert into a chart.
    - **QC Power**: This script transposes the data and adds the formatting tags needed for MultiQC to build a "stacked bar chart." In your final report, this chart shows you instantly if a sample was "contaminated" with too much Ribosomal RNA (rRNA).

---

## 3. `generate_bam_lists.py`
```python
def parse_metadata(csv_path: Path) -> dict:
    # Extracts Sample ID, Treatment, and Cell Line (Group)
    # Maps IDs to: _data/bam/{id}/{id}_Aligned.sortedByCoord.out.bam
    # Verifies file existence before adding to list
```
- **The Job**: Specifically designed for **rMATS-turbo**, this script creates the `b1.txt` and `b2.txt` files required for each comparison.
- **The Reasoning**: 
    - **Replicate Handling**: rMATS requires replicates to be comma-separated on a single line. This script automates the grouping of the 3 biological replicates for each condition.
    - **Path Safety**: It performs a critical check (`os.path.exists`) for every BAM file. If a BAM is missing (e.g., due to a failed alignment), the script errors out early, preventing rMATS from running on incomplete data.

---

## 4. `rank_as_candidates.R`

```r
top_candidates <- df %>%
  mutate(
    avg_count1 = calc_avg_counts(IJC_SAMPLE_1, SJC_SAMPLE_1),
    avg_count2 = calc_avg_counts(IJC_SAMPLE_2, SJC_SAMPLE_2),
    min_avg_coverage = pmin(avg_count1, avg_count2),
    abs_diff = abs(IncLevelDifference)
  ) %>%
  filter(min_avg_coverage >= 20) %>%
  arrange(desc(abs_diff))
```

- **The Job**: This script filters the "master list" of thousands of significant splicing events down to a handful of top-tier candidates for experimental validation and Sashimi plotting.
- **The Reasoning**:
    - **Visual Robustness**: A significant p-value doesn't always mean a "pretty" plot. By enforcing a `min_avg_coverage` of **20 reads**, we ensure that any candidate chosen will have enough read depth to produce a clear, publication-quality Sashimi plot.
    - **Magnitude Prioritization**: By sorting by `abs(IncLevelDifference)`, we focus on the events where the drug effect is most drastic (i.e., where the isoform ratio shifts the most).
    - **Efficiency**: It automates the extraction of `geneSymbol` and `ID` for the top 50 events, saving the researcher hours of manual Excel filtering.

---

## 5. Summary
These utilities ensure that the "Engine" (the Bash and R scripts) have the right fuel and produce the right exhaust. They bridge the gap between complex raw data and clean, actionable biological insights.
