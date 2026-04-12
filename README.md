# RNA-seq Pipeline

This directory contains the automated end-to-end pipeline for analyzing cancer RNA-seq data.

## Research Purpose
The objective of this pipeline is to investigate transcriptomic responses (e.g., drug treatments across multiple doses) using high-throughput RNA-seq. The pipeline is designed to:
1.  **Identify Differentially Expressed Genes (DEGs)** for each condition relative to a control baseline.
2.  **Compare overlapping effects** across treatments using multi-way Venn Diagrams.
3.  **Perform Functional Enrichment (GSEA & ORA)** to identify activated and suppressed signaling pathways, Hallmark gene sets, and biological processes.

## Directory Structure

```text
RNA-seq_pipeline/
├── _data/                       # Centralized Data Storage (git ignored)
│   ├── fastq/                   # Raw .fastq.gz reads
│   ├── fastq_trimmed/           # Reads after Trimmomatic (Adapter/QC)
│   ├── genome/                  # Reference FASTA, GTF, and BED files
│   ├── index/                   # STAR binary genome index
│   ├── bam/                     # Aligned, Coordinate-Sorted BAM files
│   ├── rseqc/                   # Post-alignment QC logs (Read distribution, etc.)
│   ├── counts/                  # Raw gene counts (featureCounts)
│   ├── multiqc/                 # Aggregated quality reports (HTML)
│   └── logs/                    # Comprehensive pipeline run logs
├── scripts_upstream/            # Upstream Automation (BASH)
│   ├── 01_genome_prep.sh        # Build index & Trim reads
│   ├── 02_star_align.sh         # Map reads using STAR
│   ├── 03_rseqc_qc.sh           # Calculate mapping statistics
│   ├── 04_quantification.sh     # Gene-level quantification
│   ├── 05_multiqc.sh            # Generate aggregated HTML report
│   └── utils/                   # Shared utility scripts
│       ├── parse_samples.py     # Metadata CSV parser
│       ├── gtf2bed.py           # Format conversion for RSeQC
│       └── biotype_to_multiqc.py # Custom MultiQC content generator
├── scripts_downstream/          # Downstream Analysis Suite (R)
│   ├── 01_data_prep.R           # Metadata & Count filtering
│   ├── 02_deseq2_dge.R          # Diff expression (Wald test/Shrinkage)
│   ├── 03_enrichment.R          # ORA & GSEA (GO/Hallmark)
│   ├── 04_01_pca.R              # Cluster analysis (PCA)
│   ├── 04_02_volcano.R          # DEG significance (Volcano)
│   ├── 04_03_venn.R             # Multi-contrast overlap (Venn)
│   ├── 04_04_heatmap_pathway.R  # Pathway-specific expression maps
│   ├── 04_05_heatmap_variable.R # High-variance gene discovery
│   ├── 04_06_enrichment_nes.R    # Hallmark enrichment barplots
│   └── 04_07_enrichment_dotplot.R # GSEA summary dotplots
├── test_upstream/               # Verification Suite: Upstream (Chr21)
│   ├── test_01-04.sh            # Prefixed test scripts (mapped to prod)
│   └── utils/                   # Test-specific metadata parsers
├── test_downstream/             # Verification Suite: Downstream (Mock)
│   ├── test_01-03.R             # Prefixed data & enrichment tests
│   └── test_04_01-07.R          # Prefixed modular visualization tests
├── _hpc/                        # Cluster Environment Management (Slurm)
│   ├── hpc_run.sbatch           # Production run template (32 CPU, 128G)
│   ├── hpc_test_micro.sbatch    # Quick test template (4 CPU, 16G)
│   └── hpc_subsample.sbatch     # Data subsetting utility
├── Sample_Data_Table.csv        # Metadata template for the study
├── scripts_upstream_docs.md     # Upstream pipeline documentation
├── scripts_downstream_docs.md   # Downstream analysis documentation
├── environment.yml              # Conda/Mamba environment definition
├── setup_env.sh                 # Environment automation script
├── run                          # Master Production Runner
├── test                         # Master Verification Runner
└── README.md                    # Primary documentation
```

## Prerequisites & Environment Setup

It is highly recommended to use **Micromamba** or **Mamba** (or Conda) to manage the pipeline's dependencies.

1. **Initialize the Environment**:
   ```bash
   bash setup_env.sh
   ```
2. **Activate the Environment**:
   ```bash
   micromamba activate cancer_rnaseq
   ```

## Quick Start

### 1. Upstream (Alignment & Quantification)
Place raw FASTQ files in `_data/fastq/`. You can run only the upstream pipeline using:
```bash
./run up
```

### 2. Downstream (Differential Expression)
Run only the downstream R analysis using:
```bash
./run down
```

### 3. Full Pipeline
Run everything from start to finish:
```bash
./run all
```

## Verification & Testing (Local Run)

For rapid verification of the pipeline logic on local hardware, a **Chromosome 21 Test Suite** is provided.

### 1. Overview
The test suite utilizes a subsampled dataset and a single chromosome (Chr21) to minimize RAM and storage requirements.

| Mode | Target | STAR RAM | Time |
| :--- | :--- | :--- | :--- |
| **Production** | Full Genome | ~30GB | 8-10 Hours |
| **Test** | Chr21 Only | ~1.5GB | ~10 Minutes |

### 2. Usage
```bash
./test all      # Run entire verification (Upstream -> Downstream)
```

## Downstream Outputs (Module 04)

The modular visualization suite generates the following figures:

| Figure | Script | Description |
| :--- | :--- | :--- |
| PCA Plot | `04_01_pca.R` | Sample clustering |
| Volcano Plots | `04_02_volcano.R` | Individual + combined DEG plots |
| Venn Diagram | `04_03_venn.R` | DEG overlap between doses |
| Pathway Heatmap | `04_04_heatmap_pathway.R` | Top enriched pathways gene expression |
| Variable Heatmap | `04_05_heatmap_variable.R` | Top 50 most variable genes |
| Hallmark NES | `04_06_enrichment_nes.R` | Hallmark enrichment scores |
| GSEA Dotplots | `04_07_enrichment_dotplot.R` | GSEA GO enrichment dots |

## Requirements

All core dependencies are managed via **Micromamba** or **Mamba**. Refer to `environment.yml` for exact pinned versions.

### Bioinformatics Tools
- **STAR**: Ultra-fast RNA-seq aligner (requires ~32GB RAM for Human genome).
- **samtools**: Alignment indexing and processing.
- **featureCounts (subread)**: Gene-level quantification.
- **RSeQC**: Post-alignment quality metrics (Read distribution, Junctions).
- **MultiQC**: Aggregated visual reporting.
- **Trimmomatic**: Adapter and quality trimming.

### Python Environment (3.9+)
- **pandas**: Required for metadata parsing and biotype summary generation.

### R/Bioconductor Libraries
- **Differential Expression**: `DESeq2`, `apeglm`.
- **Annotation & Databases**: `org.Hs.eg.db` (Human), `org.Cf.eg.db` (Dog), `msigdbr`.
- **Enrichment Analysis**: `clusterProfiler`, `enrichplot`.
- **Visualization Core**: `ggplot2`, `ComplexHeatmap`, `circlize`, `EnhancedVolcano`.
- **Layout & Interaction**: `tidyverse`, `patchwork`, `ggVennDiagram`.

> [!TIP]
> **Automated Install**: Run `bash setup_env.sh` to automatically detect your package manager and install all requirements listed above into the `cancer_rnaseq` environment.

For technical details, see the documentation for [upstream scripts](scripts_upstream_docs.md) and [downstream analysis](scripts_downstream_docs.md).
