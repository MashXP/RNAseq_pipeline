# Upstream Dissection: 03_alignment_qc.sh

This script is the "Health Check" for your samples. It calculates comprehensive post-alignment and RNA-seq quality metrics using Picard.

---

## 1. Resource & Concurrency Auto-detection
```bash
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_SAMPLES=4
CPUS_PER_SAMPLE=$((TOTAL_CPUS / PARALLEL_SAMPLES))
[ "$CPUS_PER_SAMPLE" -lt 1 ] && CPUS_PER_SAMPLE=1

JAVA_MEM_PER_SAMPLE=$((128 / PARALLEL_SAMPLES - 4))
[ "$JAVA_MEM_PER_SAMPLE" -lt 4 ] && JAVA_MEM_PER_SAMPLE=4
```
- **The Job**: This block calculates how many CPUs and how much RAM each Picard worker gets.
- **The Reasoning**: 
    - **Efficiency**: Picard is mostly single-threaded. Running 4 samples at once (`PARALLEL_SAMPLES=4`) uses the system's power much more effectively than running one at a time.
    - **Stability**: Picard is a Java tool and can "spike" in memory usage. We subtract 4GB from each sample's allocation to leave room for the operating system and threading overhead.

---

## 2. The `run_qc()` Function: Internal Logic

### 2.1 Smart-Skip Check
```bash
if [ -f "$OUTPUT_METRICS" ]; then
    echo "[Sample: $sample_name] Picard metrics already exist. Skipping."
    return 0
fi
```
- **The Job**: Checks if the final results already exist on disk.
- **The Reasoning**: This makes the script "idempotent"—if you add a single new sample to your 26-sample study, the script will only process the new one, saving you hours of compute time.

### 2.2 Species Detection (The "Brain")
```bash
sample_species=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | grep "^$sample_name " | awk '{print $4}')
```
- **The Job**: Uses a combination of Python, `grep`, and `awk` to extract the "Species" column for the current sample.
- **Why it matters**: In your cross-species project (Human cancer vs Dog lymphoma), the script must "know" the species of every sample *individually* as it runs through the background.

### 2.3 Dynamic refFlat Selection
```bash
if [ "$sample_species" == "Human" ]; then
    REFFLAT="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.refFlat"
else
    REFFLAT="$GENOME_DIR/Dog/Canis_lupus_familiaris.ROS_Cfam_1.0.113.refFlat"
fi
```
- **The Job**: Switches the "Map" file based on the detected species.
- **The Reasoning**: Picard's Intron/Exon counting requires the exact coordinates of every gene. This block ensures we don't accidentally check Dog reads against a Human map.

---

## 3. The Picard "CollectRnaSeqMetrics" Execution
```bash
picard -Xmx${JAVA_MEM_PER_SAMPLE}g CollectRnaSeqMetrics \
    I="$bam_file" \
    O="$SAMPLE_QC_DIR/${sample_name}.metrics.tmp" \
    REF_FLAT="$REFFLAT" \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    CHART_OUTPUT="$SAMPLE_QC_DIR/${sample_name}.coverage.pdf.tmp"
```
- **The Job**: Categorizes every mapped read into Exonic, Intronic, or UTR (Untranslated Region).
- **Key Parameters**:
    - `STRAND=SECOND_READ_TRANSCRIPTION_STRAND`: 
        - **Reasoning**: This matches the "Reverse-Stranded" library kit (Human/Dog TruSeq). If we use the wrong strand, your counts will be 0.
    - `RIBOSOMAL_INTERVALS=null`: We ignore rRNA here because we already handle biotypes in the next step (04).
    - **Atmocity (`.tmp`)**: We only write to temporary files during the run.

---

## 4. Atomic Finalization
```bash
mv "$SAMPLE_QC_DIR/${sample_name}.metrics.tmp" "$OUTPUT_METRICS"
mv "$SAMPLE_QC_DIR/${sample_name}.coverage.pdf.tmp" "$SAMPLE_QC_DIR/${sample_name}.picard_coverage.pdf"
```
- **The Job**: Atomically renames the temporary output to the final filename.
- **The Reasoning**: If the computer crashes halfway through a sample, it won't leave a "half-finished" file that looks like success. The final filename is only created if Picard finishes with zero errors.

---

## 5. Background Job Control (The Traffic Guard)
```bash
while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
    sleep 5
done
...
bash -c 'run_qc "'"$bam"'"' > "$SAMPLE_LOG" 2>&1 &
```
- **The Job**: Manages the "concurrency" of the samples.
- **The Reasoning**: 
    - The `while` loop says: "If there are already 4 jobs running, wait 5 seconds."
    - The `&` symbol at the end of the `bash -c` line sends the sample to the background so the loop can pick up the next sample immediately.

---

## 6. Final Synchronization
```bash
echo "Waiting for all Picard jobs to complete..."
wait
```
- **The Job**: Forces the script to pause until every background "worker" has finished.
- **The Reasoning**: This prevents the pipeline from moving to the next step (Quantification) before all the QC files are safely on the disk.
