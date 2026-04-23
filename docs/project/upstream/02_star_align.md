# Upstream Dissection: 02_star_align.sh

This script takes the clean, trimmed reads and maps them to the 3 billion letters of the Human or Canine genome. This is the core "computational" step of the upstream pipeline.

---

## 1. Resource & RAM Auto-Detection
```bash
THREADS=${SLURM_CPUS_PER_TASK:-32}
...
MEM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 128)
SORT_RAM_BYTES=$(($MEM_GB * 700000000))
```
- **The Job**: Detects the available CPUs and RAM on the system.
- **The Reasoning**: 
    - **Dynamic Threads**: Uses 32 threads by default but respects the limits of the HPC cluster if detected.
    - **STAR Sorting RAM**: STAR is notoriously memory-hungry for Human genome mapping. Calculating ~70% of your total RAM (e.g., ~$SORT_RAM_BYTES`) ensures that STAR doesn't crash from "out of memory" errors while sorting the binary data.

---

## 2. Integrity Validation
```bash
if ! gzip -t "$FASTQ_DIR/$r1" 2>/dev/null || ! gzip -t "$FASTQ_DIR/$r2" 2>/dev/null; then
    echo "CRITICAL ERROR: Trimmed FASTQ files for $sample are CORRUPTED..."
    exit 1
fi
```
- **The Job**: Uses `gzip -t` to test the integrity of the compressed data before processing.
- **The Reasoning**: If the previous trimming step crashed halfway through, the files might technically exist but be truncated. This check prevents the pipeline from "hallucinating" results from corrupted data.

---

## 3. Directory Management
```bash
OUT_DIR="$ALIGN_DIR/$sample"
mkdir -p "$OUT_DIR"
```
- **The Job**: Creates a dedicated output folder for every single sample.
- **The Reasoning**: STAR generates about 5-6 different log and data files per sample. Storing them in a per-sample directory prevents them from overwriting each other.

---

## 4. The STAR Alignment Call
```bash
STAR --genomeDir "$INDEX_DIR/$species" \
     --readFilesIn "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix "$OUT_DIR/${sample}_" \
     --outSAMtype BAM SortedByCoordinate \
     --runThreadN "$THREADS" \
     --limitBAMsortRAM "$SORT_RAM_BYTES" \
     --quantMode GeneCounts
```
- **The Job**: The main mapping engine.
- **Key Parameters**:
    - `--genomeDir`: Pulls the specific index for either Human or Canine based on your CSV metadata.
    - `--readFilesCommand gunzip -c`: Tells STAR to uncompress your data "on the fly" so you don't waste disk space.
    - `--outSAMtype BAM SortedByCoordinate`: Maps the reads AND sorts them by chromosome coordinate. This is essential for the next step (QC) and for viewing the results in a genome browser (like IGV).
    - `--quantMode GeneCounts`: Generates a "first look" at the number of reads per gene.

---

## 5. Success Metric Extraction
```bash
map_rate=$(grep "Uniquely mapped reads %" "$OUT_DIR/${sample}_Log.final.out" | cut -f2)
echo "STATUS: Uniquely Mapped Reads for $sample: $map_rate"
```
- **The Job**: Scans the STAR log file and reports the unique mapping percentage to the screen.
- **The Reasoning**: This is the first "Success Indicator" for your project.
    - **Target**: > 80% is high-quality research data.
    - **Warning**: < 50% suggests a sample or reference problem.
