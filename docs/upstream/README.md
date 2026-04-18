# Upstream Pipeline Suite: The "Genome Engine"

This directory contains the detailed in-depth dissection of the upstream automation scripts. Each file breaks down the code line-by-line, explaining the biological rationale and technical implementation of every step that turns raw sequencing data into counts.

## Step-by-Step Dissection

| Script | Biological Goal | Technical Focus |
| :--- | :--- | :--- |
| [01_genome_prep.md](01_genome_prep.md) | Indexing & Trimming | STAR, Trimmomatic |
| [02_star_align.md](02_star_align.md) | Mapping (Alignment) | STAR |
| [03_alignment_qc.md](03_alignment_qc.md) | Health Check | Picard, QC Stats |
| [04_quantification.md](04_quantification.md) | Counting Genes | featureCounts |
| [05_multiqc.md](05_multiqc.md) | Final Reporting | MultiQC |

---

## Shared Utility Logic (in scripts_upstream/utils/)

- **`parse_samples.py`**: The "Source of Truth" that maps your CSV to file paths.
- **`biotype_to_multiqc.py`**: Transforms complex counts into visual reports for MultiQC.

> [!TIP]
> **Getting Started**: Start with [01_genome_prep.md](01_genome_prep.md) to understand how the pipeline begins from a single CSV file.
