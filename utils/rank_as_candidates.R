# [[utils/rank_as_candidates.R]]
# Goal: Rank alternative splicing events by magnitude and visual coverage.
# Saves top candidates to results/alt_splicing/top_splicing_candidates.csv

library(tidyverse)

# --- Pathing ---
SCRIPT_DIR <- dirname(ifelse(nchar(Sys.getenv("RSTUDIO_TERM")), "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
PIPELINE_ROOT <- file.path(getwd()) # Assumes running from root
INPUT_CSV <- file.path(PIPELINE_ROOT, "results/alt_splicing/rmats_all_significant_events.csv")
OUTPUT_CSV <- file.path(PIPELINE_ROOT, "results/alt_splicing/top_splicing_candidates.csv")

if (!file.exists(INPUT_CSV)) {
  stop("Input file not found: ", INPUT_CSV)
}

message(">>> Ranking AS candidates from: ", INPUT_CSV)

# --- Load Data ---
df <- read.csv(INPUT_CSV)

# --- Helper: Calculate Average Coverage ---
calc_avg_counts <- function(ijc_str, sjc_str) {
  map2_dbl(strsplit(as.character(ijc_str), ","), 
           strsplit(as.character(sjc_str), ","), 
           ~ mean(as.numeric(.x) + as.numeric(.y)))
}

# --- Filtering and Ranking ---
top_candidates <- df %>%
  mutate(
    avg_count1 = calc_avg_counts(IJC_SAMPLE_1, SJC_SAMPLE_1),
    avg_count2 = calc_avg_counts(IJC_SAMPLE_2, SJC_SAMPLE_2),
    min_avg_coverage = pmin(avg_count1, avg_count2),
    abs_diff = abs(IncLevelDifference)
  ) %>%
  # Filter for visual robustness (min 20 reads in both groups)
  filter(min_avg_coverage >= 20) %>%
  # Sort by magnitude of change
  arrange(desc(abs_diff)) %>%
  select(cell_line, comparison, as_type, ID, geneSymbol, IncLevelDifference, FDR, min_avg_coverage)

# Save the top 50 (or all if less)
write.csv(head(top_candidates, 50), OUTPUT_CSV, row.names = FALSE)

# Print top 10 to console
message("\n>>> TOP 10 SPLICING CANDIDATES (Saved to results/alt_splicing/):")
print(head(top_candidates, 10))
