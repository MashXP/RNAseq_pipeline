# Downstream Dissection: Utility Scripts

These R utilities refine the raw statistical outputs into actionable biological insights and high-quality visualizations.

---

## 1. `rank_as_candidates.R`

```r
top_candidates <- df %>%
  mutate(
    avg_count1 = calc_avg_counts(IJC_SAMPLE_1, SJC_SAMPLE_1),
    avg_count2 = calc_avg_counts(IJC_SAMPLE_2, SJC_SAMPLE_2),
    min_avg_coverage = pmin(avg_count1, avg_count2),
    abs_diff = abs(IncLevelDifference)
  ) %>%
  filter(min_avg_coverage >= 20) %>%
  arrange(desc(abs_diff))
```

- **The Job**: This script filters the "master list" of thousands of significant splicing events down to a handful of top-tier candidates for experimental validation and Sashimi plotting.
- **The Reasoning**:
    - **Visual Robustness**: A significant p-value doesn't always mean a "pretty" plot. By enforcing a `min_avg_coverage` of **20 reads**, we ensure that any candidate chosen will have enough read depth to produce a clear, publication-quality Sashimi plot.
    - **Magnitude Prioritization**: By sorting by `abs(IncLevelDifference)`, we focus on the events where the drug effect is most drastic (i.e., where the isoform ratio shifts the most).
    - **Efficiency**: It automates the extraction of `geneSymbol` and `ID` for the top 50 events, saving the researcher hours of manual Excel filtering.

---

## 2. Summary
Downstream utilities are focused on **rigor** and **clarity**. They bridge the gap between "statistically significant" and "biologically compelling," ensuring that the final candidates presented in reports are both robust and visually verifiable.
