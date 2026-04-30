# 05 — Alternative Splicing Analysis & Visualization

This stage aggregates results from rMATS-turbo, identifies top candidates for validation, and generates visualizations.

## 1. Global Summary
The script `05_01_rmats_summary.R` scans the rMATS results directory and produces:
- **`rmats_summary_counts.csv`**: A table of significant event counts per AS type for each comparison.
- **`rmats_summary_plot.png`**: A stacked bar plot visualizing the distribution of splicing events across the study.
- **`rmats_all_significant_events.csv`**: A master list of all events passing the initial filters (FDR 0.05, |DeltaPSI| 0.1).

## 2. Candidate Ranking
To narrow down thousands of events to the most biologically relevant and visually robust candidates, `utils/rank_as_candidates.R` is used. For a detailed breakdown of the ranking algorithm, see [Utility Dissection](../../../utils/utils_dissection.md#4-rank_ascandidatesr).

## 3. Sashimi Plots
Individual splicing events can be visualized using Sashimi plots via the `05_02_rmats_sashimi.sh` wrapper.

### Usage:
```bash
# bash 02_dge_analysis/scripts/downstream/05_02_rmats_sashimi.sh <cell_line> <comparison> <as_type> <event_id>
bash 02_dge_analysis/scripts/downstream/05_02_rmats_sashimi.sh CNK89 Kromastat_6nM_vs_Romidepsin_6nM SE 2318
```

### Requirements:
- Coordinate-sorted BAM files with `.bai` indexes.
- `rmats2sashimiplot` installed in the environment.

## Directory Structure
Final results are stored in `results/alt_splicing/`:
- `sashimi/`: Contains folders for each plotted event.
- `rmats_summary_counts.csv`
- `rmats_summary_plot.png`
- `top_splicing_candidates.csv`
