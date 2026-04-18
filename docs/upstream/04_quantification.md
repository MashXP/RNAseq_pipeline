# Upstream Dissection: 04_quantification.sh

This script turns your "Biological Maps" (BAMs) into "Experimental Data" (The Gene Count Matrix). It uses `featureCounts` from the Subread package.

---

## 1. Group-Level Iteration (Species Discovery)
```bash
groups=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk '{print $5}' | sort -u)
for group in $groups; do
    species=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk -v grp="$group" '$5 == grp {print $4}' | head -n 1)
    ...
done
```
- **The Job**: This loop reads your CSV and identifies the major "Groups" (e.g., Human or Dog) and their corresponding Species.
- **The Reasoning**: 
    - **Efficiency**: It processes all samples from the same species at once.
    - **Modularity**: This allows you to generate one consolidated count table per species, which simplifies the multi-cell-line comparisons required for the "Krom vs Romi" study.

---

## 2. Dynamic GTF File Selection
```bash
if [ "$species" == "Human" ]; then
    GTF_FILE="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
else
    GTF_FILE="$GENOME_DIR/Dog/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
fi
```
- **The Job**: Automatically selects the correct mapping file for the current species.
- **The Reasoning**: Human (GRCh38) and Dog (ROS_Cfam_1.0) have entirely different genetic coordinates. Using the wrong GTF would result in exactly zero counts. This block ensures matching-model accuracy.

---

## 3. Gathering BAM Files per Group
```bash
files=""
for s in $group_samples; do
    if [ -f "$BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam" ]; then
        files="$files $BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam"
    fi
done
```
- **The Job**: This smaller internal loop scans the BAM directory and builds a long list of file paths.
- **The Reasoning**: `featureCounts` is much more efficient when it processes all samples in a group together. This line "assembles the team" of samples that belong in the final matrix.

---

## 4. The Gene-Level Quantification Call
```bash
featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_id -a "$GTF_FILE" -o "$COUNTS_DIR/gene_counts_${group}.txt" $files
```
- **The Job**: Counts exactly how many reads hit the protein-coding regions for every sample in the group.
- **Key Parameters**:
    - `-p`: Specifies that the data is "Paired-End" (R1 and R2).
    - `-s 2`: Specifies "Reversely Stranded" (Standard for modern Illumina kits).
    - `-t exon`: Only counts reads mapping to exons.
    - `-g gene_id`: Sets the row names to stable Ensembl IDs (e.g., ENSG...).

---

## 5. Biotype Classification & MultiQC Formatting
```bash
featureCounts ... -g gene_biotype -o "$COUNTS_DIR/biotype_counts_${group}.txt" ...
python3 "$UTILS_DIR/biotype_to_multiqc.py" ...
```
- **The Reasoning**: This creates the "Biotype Distribution" chart in your final report, making it immediately clear if your library is predominantly mRNA (Good) or mostly Ribosomal RNA (Bad capture).

> [!IMPORTANT]
> **MultiQC Visualization**: Raw count tables from `featureCounts` are essentially just huge text files. MultiQC cannot turn these into charts natively. We use the biotype conversion to "summarize" the data into categories that MultiQC can actually plot as a meaningful bargraph.
