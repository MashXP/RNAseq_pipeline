# Upstream Pipeline Documentation (BASH)

This document provides a technical breakdown of the automation scripts used for genome preparation, alignment, and quantification.

---

# Main Scripts (in scripts_upstream/)

## 01_genome_prep.sh

**Purpose:** Prepares the reference genome (e.g., GRCh38) and performs adapter/quality trimming on raw reads.

**Key Actions:**
- Downloads reference FASTA and GTF files (defaults to Human GRCh38).
- Generates a STAR genome index in `_data/index/`.
- Trimming with Trimmomatic: Uses adapter definitions (e.g., `TruSeq3-PE.fa`) to remove contamination and low-quality bases from all samples defined in the metadata.

**Inputs:**
- [None] (Uses URLs defined in the script for genome downloading).
- `_data/fastq/*.fastq.gz`: Raw reads for all samples.

**Outputs:**
- `_data/genome/`: Raw FASTA, GTF, and BED12 files.
- `_data/index/`: STAR binary index files.
- `_data/fastq_trimmed/`: High-quality, adapter-free reads.

---

## 02_star_align.sh

**Purpose:** Performs read alignment to the reference genome for paired-end samples.

**Key Actions:**
- Parses the provided Sample Metadata CSV via `utils/parse_samples.py`.
- Runs STAR in paired-end mode.
- Outputs coordinate-sorted BAM files for each sample.

**Inputs:**
- `_data/fastq_trimmed/*.fastq.gz`: Trimmed, high-quality reads.
- `[Sample_Data_Table].csv`: Sample metadata.
- `_data/index/`: Reference index.

**Outputs:**
- `_data/bam/[SampleName]/`: Sorted BAM file and alignment logs.

---

## 03_alignment_qc.sh

**Purpose:** Calculates comprehensive post-alignment and RNA-seq quality metrics using Picard.

**Key Actions:**
- **Picard CollectRnaSeqMetrics:** Reports detailed mapping stats (exonic, intronic, UTR) and 5'/3' bias.
- **Parallel Processing:** Handles multiple samples concurrently (standard: 4 at a time) to optimize HPC resources.
- **Log Isolation:** Captures per-sample stdout/stderr in the output directory.

**Inputs:**
- `_data/bam/*/`: Sorted BAM files.
- `_data/genome/*.refFlat`: Genome annotation in Picard RefFlat format.

**Outputs:**
- `_data/qc/[SampleName]/`: Metrics files and coverage plots (PDF).

---

## 04_quantification.sh

**Purpose:** Quantifies gene expression levels.

**Key Actions:**
- Uses `featureCounts` to count reads per gene (`gene_id`) and biotype (`gene_biotype`).
- Formats biotype counts for MultiQC integration.

**Inputs:**
- `_data/bam/*`: All aligned BAM files.
- `_data/genome/*.gtf`: Gene annotation.

**Outputs:**
- `_data/counts/gene_counts.txt`: Raw counts matrix.
- `_data/counts/biotype_counts_mqc.txt`: MultiQC-compatible biotype distribution.

---

## 05_multiqc.sh

**Purpose:** Aggregates all results into a single visual report.

**Key Actions:**
- Scans `_data/` for logs from all upstream tools.
- Generates a standalone HTML report.

**Outputs:**
- `_data/multiqc/rna_seq_pipeline_summary.html`: Final aggregated report.

---

# Utility Scripts (in scripts_upstream/utils/)

These helper scripts perform specialized data parsing and format conversions.

> [!NOTE]
> `biotype_to_multiqc.py` requires the `pandas` library. Others use standard Python libraries.

## gtf2bed.py
- **Function:** Converts Ensembl GTF files to BED12 format.
- **Why:** Useful for legacy tools that require BED12 format instead of GTF.
- **Usage:** `python3 gtf2bed.py input.gtf > output.bed`

## biotype_to_multiqc.py
- **Function:** Transforms raw `featureCounts` biotype summaries into MultiQC "Custom Content" format.
- **Why:** Allows MultiQC to display an interactive bar chart of read distributions (e.g., protein-coding vs. snoRNA).
- **Usage:** `python3 biotype_to_multiqc.py input.txt output_mqc.txt`

## parse_samples.py
- **Function:** Extracts sample names and paired-end file paths from the metadata CSV.
- **Why:** Decouples BASH logic from CSV parsing, ensuring reliable file path generation.
- **Usage:** Used internally by `02_star_align.sh`.
