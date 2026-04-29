# [[scripts_downstream/05_01_rmats_summary.R]]
# Goal: Aggregate rMATS-turbo results and generate a global summary of splicing events.

library(tidyverse)
library(ggplot2)
library(reshape2)
library(patchwork)

# --- Pathing ---
BASE_DIR <- getwd()
PIPELINE_ROOT <- ".."
RMATS_RESULTS_DIR <- file.path(PIPELINE_ROOT, "_data/rMATS-turbo/results")
GLOBAL_RESULTS_DIR <- file.path(PIPELINE_ROOT, "results/alt_splicing")
dir.create(GLOBAL_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

message(">>> rMATS Summary: Scanning results in ", RMATS_RESULTS_DIR)

# --- Find all filtered results ---
as_types <- c("SE", "RI", "MXE", "A3SS", "A5SS")
summary_data <- list()

cell_lines <- list.dirs(RMATS_RESULTS_DIR, full.names = FALSE, recursive = FALSE)

for (cl in cell_lines) {
  comparisons <- list.dirs(file.path(RMATS_RESULTS_DIR, cl), full.names = FALSE, recursive = FALSE)
  for (comp in comparisons) {
    comp_path <- file.path(RMATS_RESULTS_DIR, cl, comp)
    
    for (as_type in as_types) {
      filtered_file <- file.path(comp_path, "filtered", paste0(as_type, ".MATS.JCEC.filtered.txt"))
      
      if (file.exists(filtered_file)) {
        # Read the file
        df <- tryCatch({
          read.table(filtered_file, header = TRUE, sep = "\t", quote = "\"", check.names = FALSE)
        }, error = function(e) return(NULL))
        
        count <- if (is.null(df)) 0 else nrow(df)
        
        summary_data[[length(summary_data) + 1]] <- data.frame(
          cell_line = cl,
          comparison = comp,
          as_type = as_type,
          significant_events = count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(summary_data) == 0) {
  stop("No rMATS results found. Check path: ", RMATS_RESULTS_DIR)
}

final_summary <- do.call(rbind, summary_data)

# --- Export CSV ---
write.csv(final_summary, file.path(GLOBAL_RESULTS_DIR, "rmats_summary_counts.csv"), row.names = FALSE)
message("[OK] Saved summary CSV: ", file.path(GLOBAL_RESULTS_DIR, "rmats_summary_counts.csv"))

# --- Visualization: Stacked Bar Plot ---

# 1. Shorten labels for the plot (Thesis-Ready Aesthetics)
# Using logic similar to 04_04_heatmap_pathway.R
final_summary$comparison_short <- final_summary$comparison %>%
  str_replace_all("DMSO_Romidepsin", "D_R") %>%
  str_replace_all("DMSO_Kromastat", "D_K") %>%
  str_replace_all("Romidepsin_6nM", "R") %>%
  str_replace_all("Kromastat_6nM", "K") %>%
  str_replace_all("_vs_", " vs ")

# Custom palette for AS types
as_palette <- c(
  "SE"   = "#E41A1C", # Red
  "RI"   = "#377EB8", # Blue
  "MXE"  = "#4DAF4A", # Green
  "A3SS" = "#984EA3", # Purple
  "A5SS" = "#FF7F00"  # Orange
)

# Calculate total counts per AS type for the legend labels
as_totals <- final_summary %>%
  group_by(as_type) %>%
  summarize(total = sum(significant_events)) %>%
  deframe()

# Updated palette with counts in labels
as_labels <- setNames(
  paste0(names(as_palette), " = ", as_totals[names(as_palette)]),
  names(as_palette)
)
as_labels[is.na(as_totals[names(as_palette)])] <- paste0(names(as_palette)[is.na(as_totals[names(as_palette)])], " = 0")

p_main <- ggplot(final_summary, aes(x = comparison_short, y = significant_events, fill = as_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~cell_line, scales = "free_x") +
  scale_fill_manual(values = as_palette, labels = as_labels) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18, margin = margin(t = 20, b = 10)),
    plot.subtitle = element_text(size = 14, margin = margin(b = 20)),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "rMATS-turbo: Significant Splicing Events",
    subtitle = "Filtering: FDR <= 0.05, |DeltaPSI| >= 0.1, Avg Counts >= 10",
    x = "Comparison",
    y = "Number of Significant Events",
    fill = "AS Type"
  )

# 2. Add Legend for Shortened names (using patchwork to combine)
total_n <- sum(final_summary$significant_events)
shorthand_text <- paste0(
  "Total Significant Events: N = ", total_n, "\n",
  "Treatments: D = DMSO (Control) | R = Romidepsin | K = Kromastat\n",
  "Splicing Types:\n",
  "  SE   = Skipped Exon\n",
  "  RI   = Retained Intron\n",
  "  MXE  = Mutually Exclusive Exons\n",
  "  A3SS = Alternative 3' Splice Site | A5SS = Alternative 5' Splice Site"
)

p_leg <- ggplot() + 
  annotate("text", x = 0, y = 0, label = shorthand_text, size = 4.2, hjust = 0, fontface = "italic", lineheight = 1.1) + 
  theme_void() +
  theme(plot.margin = margin(l = 25, b = 10, t = 10))

# Combine plots
p_final <- p_main / p_leg + plot_layout(heights = c(10, 3.5))

ggsave(file.path(GLOBAL_RESULTS_DIR, "rmats_summary_plot.png"), p_final, width = 14, height = 11, bg = "white")
message("[OK] Saved summary plot: ", file.path(GLOBAL_RESULTS_DIR, "rmats_summary_plot.png"))

# --- Export Top Candidates (Optional but useful) ---
# This part would create a big table of all significant events across all comparisons
# Use with caution if there are thousands of events.
message(">>> Generating master list of significant events...")

master_list <- list()
for (cl in cell_lines) {
  comparisons <- list.dirs(file.path(RMATS_RESULTS_DIR, cl), full.names = FALSE, recursive = FALSE)
  for (comp in comparisons) {
    comp_path <- file.path(RMATS_RESULTS_DIR, cl, comp)
    for (as_type in as_types) {
      filtered_file <- file.path(comp_path, "filtered", paste0(as_type, ".MATS.JCEC.filtered.txt"))
      if (file.exists(filtered_file)) {
        df <- tryCatch({
          read.table(filtered_file, header = TRUE, sep = "\t", quote = "\"", check.names = FALSE)
        }, error = function(e) return(NULL))
        
        if (!is.null(df) && nrow(df) > 0) {
          df$cell_line <- cl
          df$comparison <- comp
          df$as_type <- as_type
          # Reorder columns to put metadata first
          cols <- colnames(df)
          df <- df[, c("cell_line", "comparison", "as_type", cols[1:(length(cols)-3)])]
          master_list[[length(master_list) + 1]] <- df
        }
      }
    }
  }
}

if (length(master_list) > 0) {
  final_master <- bind_rows(master_list)
  write.csv(final_master, file.path(GLOBAL_RESULTS_DIR, "rmats_all_significant_events.csv"), row.names = FALSE)
  message("[OK] Saved master list: ", file.path(GLOBAL_RESULTS_DIR, "rmats_all_significant_events.csv"))
}

message("rMATS Summary Complete.")
