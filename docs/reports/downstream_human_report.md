---
marp: true
theme: gaia
paginate: true
backgroundColor: #fff
style: |
  section {
    font-family: 'Inter', 'Helvetica Neue', Helvetica, Arial, sans-serif;
    font-size: 28px;
    padding: 40px;
    color: #333;
  }
  img[alt~="center"] {
    position: relative;
    left: 50%;
    transform: translateX(-50%);
  }
  h1 {
    color: #004a99;
    font-size: 1.8em;
    border-bottom: 2px solid #004a99;
    margin-bottom: 0.5em;
  }
  h2 {
    color: #0066cc;
    font-size: 1.4em;
  }
  h3 {
    color: #0f1011;
    font-size: 1.4em;
    margin-bottom: -0.75em;
  }
  p {
    font-size: 1em;
    
  }
  footer {
    font-size: 0.5em;
    color: #888;
  }
  header {
    font-size: 0.5em;
    color: #004a99;
    font-weight: bold;
  }
  .columns {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 0.5rem;
  }
  .center {
    text-align: center;
  }
  .small-text {
    font-size: 0.8em;
  }
  .highlight {
    color: #d11141;
    font-weight: bold;
  }
  table {
    font-size: 0.7em;
    width: 100%;
    border-collapse: collapse;
    margin-bottom: 20px;
    font-family: 'Inter', sans-serif;
  }
  th {
    background-color: #004a99;
    color: #ffffff;
    padding: 12px 15px;
    text-align: left;
    border-bottom: 3px solid #003366;
    font-weight: 700;
  }
  td {
    padding: 10px 15px;
    border-bottom: 1px solid #eef2f7;
    color: #444;
    vertical-align: top;
  }
  tr:nth-child(even) {
    background-color: #f8faff;
  }
  tr:hover {
    background-color: #f1f5ff;
  }
  
---

<!-- _class: lead invert -->

# Cancer Treatments Project
## KROMASTAT and ROMIDEPSIN in HUMAN
### Transcriptomic Report & Analysis

---

# STUDY DESIGN

<div class="columns">

<div>

- Compare the biological effects of **Kromastat** and **Romidepsin** on lymphoma cell lines.
- **Cell Lines**: H9 and SUPM2 (Human).
- **Treatment Groups**:
  - **DMSO_Kromastat** (Control)
  - **Kromastat 6 nM**
  - **DMSO_Romidepsin** (Control)
  - **Romidepsin 6 nM**
- **Replication**: Triplicate biological replicates for statistical reliability.

</div>

<div>

![width:600px](img/human-report_0.png)

</div>

---

# METHODOLOGY: Individual Drug Response

<div class="columns">
<div>

Applied to both **Kromastat** and **Romidepsin** (H9 & SUPM2):

1. **Identify DEGs**: Significant expression changes vs. DMSO (Control).
2. **Visualize Shifts**: Mapping via Volcano and Venn diagrams.
3. **Decode Pathways**: GSEA Hallmark signature analysis.

</div>
<div>

<div style="position: relative; display: inline-block;">
  <img src="img/human-report_0.png" width="600px" />
  <svg style="position: absolute; top: 0; left: 0;" width="600" height="450" viewBox="0 0 600 450">
    <!-- H9 Kromastat 6nM Highlight -->
    <rect x="110" y="270" width="230" height="155" fill="none" stroke="#114bd1" stroke-width="8" rx="10" />
  </svg>
</div>

</div>
</div>

---

# METHODOLOGY: Individual Drug Response

<div class="columns">
<div>

Applied to both **Kromastat** and **Romidepsin** (H9 & SUPM2):

1. **Identify DEGs**: Significant expression changes vs. DMSO (Control).
2. **Visualize Shifts**: Mapping via Volcano and Venn diagrams.
3. **Decode Pathways**: GSEA Hallmark signature analysis.

</div>
<div>

<div style="position: relative; display: inline-block;">
  <img src="img/human-report_0.png" width="600px" />
  <svg style="position: absolute; top: 0; left: 0;" width="600" height="450" viewBox="0 0 600 450">
    <!-- SUP-M2 Kromastat 6nM Highlight -->
    <rect x="110" y="105" width="230" height="155" fill="none" stroke="#d11141" stroke-width="8" rx="10" />
    <!-- H9 Kromastat 6nM Highlight -->
    <rect x="110" y="270" width="230" height="155" fill="none" stroke="#114bd1" stroke-width="8" rx="10" />
  </svg>
</div>

</div>
</div>

---


# METHODOLOGY: Head-to-Head Comparison

<div class="columns">
<div>

Contrast **Romidepsin** vs **Kromastat** to define differentiation:

- **Target Selectivity**: Unique vs. shared gene sets.
- **Mechanistic Divergence**: Differentiated biological signatures.
- **Conserved Core**: Rigid "Proliferation Crash" program.

</div>
<div>

![width:600px](img/human-report_0.png)

</div>
</div>

---

# PCA: Distribution of Data

<div class="columns">
<div>

Firstly, Check the data distribution to ensure replicate clustering and group isolation.

</div>

<div class="columns">
<div>

![width:600px](../../results/human/figures/04_01_pca_combined.png)

</div>

---

# PCA: Lineage Specificity

<div class="columns">
<div class="center">

### H9 Cell Line
Clear separation of 
drug effects along PC2.

![width:500px center](../../results/human/figures/04_01_pca_H9.png)
</div>

<div class="center">

### SUPM2 Cell Line
Stronger drug-induced shift 
in the cancer background.

![width:500px center](../../results/human/figures/04_01_pca_SUPM2.png)

</div>
</div>

---

# DEGs in H9 (Healthy-like)

<div class="columns">

<div>

### Kromastat 6nM
~1,411 DEGs
![width:500px center](../../results/human/figures/04_02_volcano_H9_Kromastat_6nM_vs_DMSO_Kromastat.png)

</div>

<div>

### Romidepsin 6nM
~3,910 DEGs
![width:500px center](../../results/human/figures/04_02_volcano_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.png)

</div>

---

# Head-to-Head: H9 Lineage
### Kromastat < Romidepsin in H9

**Romidepsin** induces a **significantly broader modulation** of the H9 transcriptome at the same **molarity**.

---
<style scoped>
section {
  overflow: visible;
  position: relative;
}
img[alt~="comp_volcano"] {
  position: absolute;
  top: -50px;
  left: 50%;
  transform: translateX(-50%);
  z-index: 0; 
}
</style>
![width:850px comp_volcano](../../results/human/figures/04_02_volcano_H9_Romidepsin_6nM_vs_Kromastat_6nM.png)

---

# DEGs in SUPM2 (Cancer)

<div class="columns">

<div>

### Kromastat 6nM
~2,086 DEGs
![width:450px center](../../results/human/figures/04_02_volcano_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png)

</div>

<div>

### Romidepsin 6nM
~4,840 DEGs
![width:450px center](../../results/human/figures/04_02_volcano_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.png)

</div>

---
<style scoped>
section {
  background-position: center calc(50% - 100px);
}
</style>

# Head-to-Head: SUPM2 Lineage
### Kromastat < Romidepsin in SUPM2


Consistent with H9, Romidepsin shows a **higher potency** "overdrive" effect in the cancer background.

---
<style scoped>
section {
  overflow: visible;
  position: relative;
}
img[alt~="comp_volcano"] {
  position: absolute;
  top: -50px;
  left: 50%;
  transform: translateX(-50%);
  z-index: 0; 
}
</style>

![width:880px comp_volcano](../../results/human/figures/04_02_volcano_SUPM2_Romidepsin_6nM_vs_Kromastat_6nM.png)


---

# Global Overlap: All DEGs

Both drugs share a significant core of differentially expressed genes, though Romidepsin induces a much broader response.

<div class="columns">
<div>

- **1,218 Genes Shared** (Total)
- Significant intersection of biological targets.
- **Kromastat**: ~76% of its DEGs are shared with Romidepsin.
- **Romidepsin**: ~31% of its DEGs are shared with Kromastat.

</div>

![width:575px center](../../results/human/figures/04_03_venn_all.png)
</div>

---

# Global Overlap: Shared Activation (Up)

Both drugs are highly effective at gene activation, sharing a large core signature.

<div class="columns">
<div>

- **1,432 Genes Shared Up**
- Highly significant overlap (p < 2e-16).
- **Jaccard Index**: 0.364
- **Representation Factor**: 5.7x
- Includes core targets like `EPAS1` and `GSN`.
</div>

![width:575px center](../../results/human/figures/04_03_venn_up.png)
</div>

---

# Global Overlap: Shared Suppression (Down)

Conserved silencing program, representing an extremely rigid core response.

<div class="columns">
<div>

- **70 Genes Shared Down**
- Significant overlap (p = 2.95e-90).
- **Representation Factor**: <span class="highlight">26.9x</span>
- Striking consistency in shared silencing.
- Includes `TNFRSF8` and `PTPN7`.
</div>

![width:575px center](../../results/human/figures/04_03_venn_down.png)
</div>

---

# Pathway Analysis: H9 + Kromastat

<div class="columns">
<div class="center">

### Hallmark NES
Classic "Proliferation Crash" observed.
![width:450px](../../results/human/figures/04_06_hallmark_nes_H9_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
<div class="center">

### GSEA Dotplot
High significance in E2F and MYC targets.
![width:370px](../../results/human/figures/04_07_gsea_dotplot_H9_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
</div>

---

# Pathway Analysis: SUPM2 + Kromastat

<div class="columns">
<div class="center">

### Hallmark NES
Stronger E2F suppression in cancer.
![width:500px](../../results/human/figures/04_06_hallmark_nes_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
<div class="center">


### GSEA Dotplot
Leading edge shows broad gene engagement.
![width:370px](../../results/human/figures/04_07_gsea_dotplot_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
</div>

---
<style scoped>
table {
  font-size: 1.1em;
  margin-top: 20px;
}
td {
  padding: 15px 20px;
}
th {
  background-color: #002d5d;
}
</style>

# SUMMARY: THE RESULTS ARE GOOD!

| SUPPRESSED PATHWAYS | PROCESS | EXPECTED RESULTS |
| :--- | :--- | :--- |
| **E2F Targets** | ↓ Proliferation, ↓ Cell cycle, ↓ DNA synthesis | **Tumor growth slows** |
| **MYC Targets** | ↓ Growth, ↓ Ribosome biogenesis, ↓ Metabolism | **Tumor growth slows** |
| **G2M Checkpoint** | ↓ Mitotic division, ↑ DNA damage stress | **Mitotic arrest / Apoptosis** |

---

# Pathway Analysis: H9 + Romidepsin

<div class="columns">
<div class="center">

### Hallmark NES
Broader pathway modulation than Kromastat.
![width:450px](../../results/human/figures/04_06_hallmark_nes_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
<div class="center">

### GSEA Dotplot
Stronger immune-suppressive footprint.
![width:390px](../../results/human/figures/04_07_gsea_dotplot_H9_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
</div>

---

# Pathway Analysis: SUPM2 + Romidepsin

<div class="columns">
<div class="center">

### Hallmark NES
Maximal growth suppression signature.
![width:475px](../../results/human/figures/04_06_hallmark_nes_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
<div class="center">

### GSEA Dotplot
Apoptosis + KRAS stress activation.
![width:400px](../../results/human/figures/04_07_gsea_dotplot_SUPM2_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
</div>

---

# Comparison: SAME PATTERN

Kromastat and Romidepsin produce **broadly similar pathway-level** responses.

- **Romidepsin**: Stronger transcriptional modulation across shared pathways.
- **Kromastat**: More selective, preserving immune signaling that Romidepsin suppresses.

<div class="columns">
<div>

![width:450px](../../results/human/figures/04_06_hallmark_nes_combined.png)
</div>
<div>

![width:450px](../../results/human/figures/04_07_gsea_dotplot_combined.png)
</div>
</div>

---

# Head-to-Head: Functional Comparison

<div class="columns">

<div class="center">

### H9 Comparison
![width:530px center](../../results/human/figures/04_06_hallmark_nes_H9_Romidepsin_6nM_vs_Kromastat_6nM.png)
</div>

<div class="center">

### SUPM2 Comparison
![width:530px center](../../results/human/figures/04_06_hallmark_nes_SUPM2_Romidepsin_6nM_vs_Kromastat_6nM.png)
</div>
</div>

---

# Consistency Analysis: Shared DEGs

Intersections of DEGs across treatment groups confirm a shared core response program despite drug-specific potency variances.

<div class="columns">
<div>

| Treatment | Cell Line | Up | Down | Total |
| :--- | :--- | :---: | :---: | :---: |
| **Romidepsin** | H9 | 2,381 | 987 | 3,368 |
| **Romidepsin** | SUPM2 | 3,351 | 905 | 4,256 |
| **Kromastat** | H9 | 911 | 101 | 1,012 |
| **Kromastat** | SUPM2 | 1,323 | 140 | 1,463 |

</div>

![width:550px](../../results/human/figures/04_09_upset_consistency.png)
</div>

- **Core Signature**: A significant number of genes are universally regulated regardless of cell line or drug potency.
- **Top Conserved Targets**: `DHRS2`, `NOVA2`, `DACT3`, `TNFRSF8`, and `PTPN7`.

---

# Transcriptome Response Consistency

LFC correlation between H9 and SUPM2 confirms mechanism conservation across cell lines.

<div class="columns">
<div>

**Pearson R**: Romidepsin = 0.606, Kromastat = 0.550. High linear correlation indicates drug-induced shifts follow a strict mechanistic path independent of cancer background (Indolent vs. Aggressive).

</div>

![width:600px](../../results/human/figures/04_10_lfc_correlation.png)
</div>

---

# General Conclusions

- **Potency**: Romidepsin is a more aggressive cytotoxic agent.
- **Mechanism**: Both drugs drive a massive "Proliferation Crash."
- **Differentiation**:
  - **Romidepsin**: Metabolic crusher, broad immune silencer.
  - **Kromastat**: Immune-intelligent inhibitor, spares interferon signaling.
- **Therapeutic Potential**: Kromastat is a promising candidate for **combination with immunotherapies**.

---

# HUMAN – CANINE OVERLAP

Main matched-pair alluvial plot: visualizes conserved vs species-specific response.

---

![bg width:1200px center](../../results/comparative/figures/04_11_alluvial_gene_flow.png)

---

# Conserved Cross-Species Signal

<div class="columns">
<div class="center">

### Conserved Up
![width:500px](../../results/comparative/figures/04_11_alluvial_aggressive.png)
</div>
<div class="center">

### Conserved Down
![width:500px](../../results/comparative/figures/04_11_alluvial_indolent.png)
</div>
</div>

- **Romidepsin**: Strongest conserved cross-species signal.
- **Krom**: Smaller conserved core, more human-only genes in aggressive models.

---

# TOP CONCORDANT GENES

<div class="columns">
<div class="center">

### CONCORDANT UP
Conserved activation across species.
![width:520px](../../results/comparative/figures/04_12_ortholog_scatter_Romidepsin_6nM.png)
</div>
<div class="center">

### CONCORDANT DOWN
Conserved silencing across species.
![width:520px](../../results/comparative/figures/04_12_ortholog_scatter_Kromastat_6nM.png)
</div>
</div>

---

<!-- _class: lead invert -->

# Thank You!
## Questions?
