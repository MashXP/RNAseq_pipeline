# Upstream Dissection: 04_quantification.sh

This script turns your "Biological Maps" (BAMs) into "Experimental Data" (The Gene Count Matrix). It uses `featureCounts` from the Subread package.

---

## 1. Species-Level Iteration
```bash
species_list=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk '{print $4}' | sort -u)
for sp in $species_list; do
    sp_lower=$(echo "$sp" | tr '[:upper:]' '[:lower:]')
    mkdir -p "$FEATURECOUNTS_DIR/$sp_lower"
    ...
done
```
- **The Job**: This loop reads your CSV and identifies the unique Species (e.g., Human or Canine).
- **The Reasoning**:
    - **Efficiency**: It processes all samples from the same species at once.
    - **Modularity**: This generates one consolidated count table per species in the `_data/featurecounts/` directory, which simplifies the multi-cell-line comparisons required for the "Krom vs Romi" study.

---

## 2. Dynamic GTF File Selection
```bash
if [ "$sp" == "Human" ]; then
    GTF_FILE="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
else
    GTF_FILE="$GENOME_DIR/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
fi
```
- **The Job**: Automatically selects the correct mapping file for the current species.
- **The Reasoning**: Human (GRCh38) and Canine (ROS_Cfam_1.0) have entirely different genetic coordinates. Using the wrong GTF would result in exactly zero counts. This block ensures matching-model accuracy.

---

## 3. Gathering BAM Files per Species
```bash
sp_samples=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk -v species="$sp" '$4 == species {print $1}')
files=""
for s in $sp_samples; do
    if [ -f "$BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam" ]; then
        files="$files $BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam"
    fi
done
```
- **The Job**: This internal loop builds a list of BAM file paths belonging to the current species.
- **The Reasoning**: `featureCounts` is most efficient when processing all samples of a species together, creating a unified matrix.

---

## 4. The Gene-Level Quantification Call
```bash
featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_id -a "$GTF_FILE" \
    -o "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_counts.tsv" $files
```
- **The Job**: Executes the core quantification at the gene level.
- **The Output**: A TSV file named `{species}_featurecounts_counts.tsv` containing the gene count matrix.

---

## 5. The Biotype-Level Quantification & Processing
```bash
featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_biotype -a "$GTF_FILE" \
    -o "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv" $files
python3 "$UTILS_DIR/biotype_to_multiqc.py" \
    "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv" \
    "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts.tsv"
```
- **The Job**: First counts reads per biotype (e.g., protein_coding, lncRNA), then processes the raw output into a cleaned format compatible with MultiQC reports.
- **The Output**: `{species}_featurecounts_biotype_counts.tsv`.
