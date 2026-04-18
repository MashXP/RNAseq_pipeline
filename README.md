# RNA-seq Pipeline

This directory contains the automated end-to-end pipeline for analyzing cancer RNA-seq data.

## Research Purpose
The objective of this pipeline is to investigate transcriptomic responses (e.g., drug treatments across multiple doses) using high-throughput RNA-seq. The pipeline is designed to:
1.  **Identify Differentially Expressed Genes (DEGs)** for each condition relative to a control baseline.
2.  **Compare overlapping effects** across treatments using multi-way Venn Diagrams and UpSet plots.
3.  **Perform Functional Enrichment (GSEA & ORA)** to identify activated and suppressed signaling pathways, Hallmark gene sets, and biological processes.

## Directory Structure

```text
RNAseq_pipeline/
├── _data/                       # Centralized Data Storage (git ignored)
├── _hpc/                        # Cluster Environment Management (Slurm)
├── scripts_upstream/            # Upstream Automation (BASH)
├── scripts_downstream/          # Downstream Analysis Suite (R)
├── test_upstream/               # Verification Suite: Upstream (Chr21)
├── test_downstream/             # Verification Suite: Downstream (Mock)
├── docs/                        # Project Documentation
│   ├── project/                 # Project Management
│   │   ├── upstream/            # Upstream script dissections
│   │   ├── downstream/          # Downstream script dissections
│   │   ├── research_plan.md
│   │   └── TODO.md
│   └── reports/                 # Generated Analysis Reports
├── Sample_Data_Table.csv        # Metadata template for the study
├── environment.yml              # Conda/Mamba environment definition
├── setup_env.sh                 # Environment automation script
├── run                          # Master Production Runner
└── test                         # Master Verification Runner
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

| Mode | Target | STAR RAM | Time |
| :--- | :--- | :--- | :--- |
| **Production** | Full Genome | ~30GB | 8-10 Hours |
| **Test** | Chr21 Only | ~1.5GB | ~10 Minutes |

```bash
./test all      # Run entire verification (Upstream -> Downstream)
```

## Pipeline Architecture & Documentation

The pipeline is split into an upstream BASH execution engine and a downstream R statistical suite. Full, line-by-line documentation for every script can be found in the `docs/` directory.

### Upstream Pipeline Suite: The "Genome Engine"

| Script | Biological Goal | Technical Focus | Documentation |
| :--- | :--- | :--- | :--- |
| `01_genome_prep.sh` | Indexing & Trimming | STAR, Trimmomatic | [Docs](docs/project/upstream/01_genome_prep.md) |
| `02_star_align.sh` | Mapping (Alignment) | STAR | [Docs](docs/project/upstream/02_star_align.md) |
| `03_alignment_qc.sh` | Health Check | Picard, QC Stats | [Docs](docs/project/upstream/03_alignment_qc.md) |
| `04_quantification.sh` | Counting Genes | featureCounts | [Docs](docs/project/upstream/04_quantification.md) |
| `05_multiqc.sh` | Final Reporting | MultiQC | [Docs](docs/project/upstream/05_multiqc.md) |

**Shared Utilities (`scripts_upstream/utils/`)**:
- `parse_samples.py`: The "Source of Truth" that maps your CSV to file paths.
- `biotype_to_multiqc.py`: Transforms complex counts into visual reports for MultiQC.

### Downstream Pipeline Suite: R Analysis

| Script | Biological Goal | Technical Focus | Documentation |
| :--- | :--- | :--- | :--- |
| `01_bridge_data_prep.R` | Translate mentor counts | Data formatting | [Docs](docs/project/downstream/01_bridge_data_prep.md) |
| `01_data_prep.R` | Clean & organize | Metadata integration | [Docs](docs/project/downstream/01_data_prep.md) |
| `02_deseq2_dge.R` | DGE Engine | DESeq2, apeglm | [Docs](docs/project/downstream/02_deseq2_dge.md) |
| `03_enrichment.R` | Functional analysis | GSEA, GO, KEGG | [Docs](docs/project/downstream/03_enrichment.md) |

**04.x Visualization Suite**: High-resolution modular plotting.
- `04_01_pca.R`: Sample clustering / PCA. [Docs](docs/project/downstream/04_01_pca.md)
- `04_02_volcano.R`: DEG significance (Volcano). [Docs](docs/project/downstream/04_02_volcano.md)
- `04_03_venn.R`: Multi-contrast overlap (Venn). [Docs](docs/project/downstream/04_03_venn.md)
- `04_04_heatmap_pathway.R`: Top enriched pathways gene expression. [Docs](docs/project/downstream/04_04_heatmap_pathway.md)
- `04_05_heatmap_variable.R`: Top 50 most variable genes. [Docs](docs/project/downstream/04_05_heatmap_variable.md)
- `04_06_enrichment_nes.R`: Hallmark enrichment scores. [Docs](docs/project/downstream/04_06_enrichment_nes.md)
- `04_07_enrichment_dotplot.R`: GSEA Hallmark enrichment dots (Mirror Plot). [Docs](docs/project/downstream/04_07_enrichment_dotplot.md)
- `04_08_ora_dotplot.R`: ORA GO enrichment dots. [Docs](docs/project/downstream/04_08_ora_dotplot.md)
- `04_09_upset_consistency.R`: UpSet multi-way comparisons. [Docs](docs/project/downstream/04_09_upset_consistency.md)
- `04_10_correlation_plots.R`: Correlation scatters. [Docs](docs/project/downstream/04_10_correlation_plots.md)

*(Note: See [docs/project/downstream/libraries.md](docs/project/downstream/libraries.md) for full library rationales).*

## Requirements

All core dependencies are managed via **Micromamba** or **Mamba**. Refer to `environment.yml` for exact pinned versions.
Run `bash setup_env.sh` to automatically detect your package manager and install all requirements into the `cancer_rnaseq` environment.

---

## Inspiration & Related Projects

This pipeline draws inspiration from and builds upon methodological patterns established in the [YeastAnalysis](https://github.com/MashXP/YeastAnalysis) project, reflecting a shared focus on robust, automated bioinformatic workflows.
