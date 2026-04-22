# Audit Principles & Shorthand

These principles ensure precise and standardized visual descriptions of gene expression patterns (blocks) in transcriptomic heatmaps.

### 1. The Expression Block (The "Wall")
When describing a block of genes for a cell line or condition, address the following three pillars:

- **Dominant Color**:
    - **Red**: Up-regulated trend.
    - **Blue**: Down-regulated trend.
- **Intensity**:
    - **Vivid**: High z-score (saturated color).
    - **Pale**: Low z-score (faded/desaturated color).
    - **White**: Neutral z-score (near zero).
- **Homogeneity**:
    - **Solid**: The entire block is uniform in color and intensity.
    - **Mixed**: The block contains a variety of intensities or alternating colors.

### 2. Outliers (The "Speckles")
After describing the dominant wall, explicitly list any genes that deviate from the trend:
- *Example*: "Vivid Red wall, but with **CCNA2** pale blue (outlier)."
- *Example*: "Solid Blue wall, with **NOP16** white."

### 3. Shorthand & Comparisons
- `//`: **Specifically means "Similar to Indolent (H9)."** Used when the Aggressive (SUPM2) cell line exhibits the same pattern as previously described for H9.
- `(D vs R)`: Comparison between DMSO and Romidepsin.
- `(D vs K)`: Comparison between DMSO and Kromastat.
