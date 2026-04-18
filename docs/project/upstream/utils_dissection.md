# Upstream Dissection: Utility Scripts

While the main `.sh` scripts do the heavy lifting, these Python utilities are the "connective tissue" that ensures your metadata is parsed correctly and your final reports are readable.

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
- **The Job**: This is the "Source of Truth" for your entire pipeline. It reads your `drPhuong_Sample_Data_Table.csv` and translates each row into a single line of text.
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

## 3. Summary
These utilities ensure that the "Engine" (the Bash scripts) has the right fuel (clean metadata) and produces the right exhaust (readable reports). Without these "small" scripts, the large-scale pipeline would be much more difficult to debug.
