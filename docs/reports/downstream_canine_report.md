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
    color: #00663a;
    font-size: 1.8em;
    border-bottom: 2px solid #00663a;
    margin-bottom: 0.5em;
  }
  h2 {
    color: #009955;
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
    color: #00663a;
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
    background-color: #00663a;
    color: #ffffff;
    padding: 12px 15px;
    text-align: left;
    border-bottom: 3px solid #004d2b;
    font-weight: 700;
  }
  td {
    padding: 10px 15px;
    border-bottom: 1px solid #eef2f7;
    color: #444;
    vertical-align: top;
  }
  tr:nth-child(even) {
    background-color: #f4fff9;
  }
  tr:hover {
    background-color: #e8fff3;
  }
  
---

<!-- _class: lead invert -->

# Cancer Treatments Project
## KROMASTAT and ROMIDEPSIN in CANINE
### Transcriptomic Report & Analysis

---

# STUDY DESIGN

<div class="columns">

<div>

- Compare the biological effects of **Kromastat** and **Romidepsin** on canine lymphoma cell lines.
- **Cell Lines**: UL1 and CNK89 (Canine).
- **Treatment Groups**:
  - **DMSO_Kromastat** (Control)
  - **Kromastat 6 nM**
  - **DMSO_Romidepsin** (Control)
  - **Romidepsin 6 nM**
- **Replication**: Triplicate biological replicates for statistical reliability.

</div>

<div>

![width:600px](img/canine-report_0.png)

</div>

---

# METHODOLOGY: Individual Drug Response

<div class="columns">
<div>

Applied to both **Kromastat** and **Romidepsin** (UL1 & CNK89):

1. **Identify DEGs**: Significant expression changes vs. DMSO (Control).
2. **Visualize Shifts**: Mapping via Volcano and Venn diagrams.
3. **Decode Pathways**: GSEA Hallmark signature analysis.

</div>
<div>

<div style="position: relative; display: inline-block;">
  <img src="img/canine-report_0.png" width="600px" />
  <svg style="position: absolute; top: 0; left: 0;" width="600" height="450" viewBox="0 0 600 450">
    <!-- UL1 Kromastat 6nM Highlight -->
    <rect x="110" y="270" width="230" height="155" fill="none" stroke="#009955" stroke-width="8" rx="10" />
  </svg>
</div>

</div>
</div>

---

# METHODOLOGY: Individual Drug Response

<div class="columns">
<div>

Applied to both **Kromastat** and **Romidepsin** (UL1 & CNK89):

1. **Identify DEGs**: Significant expression changes vs. DMSO (Control).
2. **Visualize Shifts**: Mapping via Volcano and Venn diagrams.
3. **Decode Pathways**: GSEA Hallmark signature analysis.

</div>
<div>

<div style="position: relative; display: inline-block;">
  <img src="img/canine-report_0.png" width="600px" />
  <svg style="position: absolute; top: 0; left: 0;" width="600" height="450" viewBox="0 0 600 450">
    <!-- CNK89 Kromastat 6nM Highlight -->
    <rect x="110" y="105" width="230" height="155" fill="none" stroke="#d11141" stroke-width="8" rx="10" />
    <!-- UL1 Kromastat 6nM Highlight -->
    <rect x="110" y="270" width="230" height="155" fill="none" stroke="#009955" stroke-width="8" rx="10" />
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
- **Conserved Core**: Rigid \"Proliferation Crash\" program.

</div>
<div>

![width:600px](img/canine-report_0.png)

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

![width:600px](../../results/canine/figures/04_01_pca_combined.png)

</div>

---

# PCA: Lineage Specificity

<div class="columns">
<div class="center">

### UL1 Cell Line
Clear separation of 
drug effects along PC2.

![width:500px center](../../results/canine/figures/04_01_pca_UL1.png)
</div>

<div class="center">

### CNK89 Cell Line
Stronger drug-induced shift 
in the cancer background.

![width:500px center](../../results/canine/figures/04_01_pca_CNK89.png)

</div>
</div>

---

# DEGs in UL1 (Indolent)

<div class="columns">

<div>

### Kromastat 6nM
~374 DEGs
![width:500px center](../../results/canine/figures/04_02_volcano_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png)

</div>

<div>

### Romidepsin 6nM
~1,733 DEGs
![width:500px center](../../results/canine/figures/04_02_volcano_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.png)

</div>

---

# Head-to-Head: UL1 Lineage
### Kromastat < Romidepsin in UL1

**Romidepsin** induces a **significantly broader modulation** of the UL1 transcriptome at the same **molarity**.

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
![width:850px comp_volcano](../../results/canine/figures/04_02_volcano_UL1_Romidepsin_6nM_vs_Kromastat_6nM.png)

---

# DEGs in CNK89 (Aggressive)

<div class="columns">

<div>

### Kromastat 6nM
~19 DEGs
![width:450px center](../../results/canine/figures/04_02_volcano_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png)

</div>

<div>

### Romidepsin 6nM
~622 DEGs
![width:450px center](../../results/canine/figures/04_02_volcano_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.png)

</div>

---
<style scoped>
section {
  background-position: center calc(50% - 100px);
}
</style>

# Head-to-Head: CNK89 Lineage
### Kromastat << Romidepsin in CNK89


Kromastat shows **near-zero activity** in CNK89, while Romidepsin retains a robust response — a critical lineage-specific divergence.

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

![width:880px comp_volcano](../../results/canine/figures/04_02_volcano_CNK89_Romidepsin_6nM_vs_Kromastat_6nM.png)


---

# Global Overlap: All DEGs

Both drugs share a significant core of differentially expressed genes, though Romidepsin induces a much broader response.

<div class="columns">
<div>

- **367 Genes Shared** (Total)
- Significant intersection of biological targets.
- **Kromastat**: ~77% of its DEGs are shared with Romidepsin.
- **Romidepsin**: ~30% of its DEGs are shared with Kromastat.

</div>

![width:575px center](../../results/canine/figures/04_03_venn_all.png)
</div>

---

# Global Overlap: Shared Activation (Up)

Both drugs are highly effective at gene activation, sharing a large core signature.

<div class="columns">
<div>

- **365 Genes Shared Up**
- Highly significant overlap (p < 2e-16).
- **Jaccard Index**: 0.297
- **Representation Factor**: 10.75x
- Includes core targets like `ENPP2` and `HPCAL4`.
</div>

![width:575px center](../../results/canine/figures/04_03_venn_up.png)
</div>

---

# Global Overlap: Shared Suppression (Down)

Extremely narrow conserved silencing — only 2 shared down-regulated genes, with an exceptionally high representation factor.

<div class="columns">
<div>

- **2 Genes Shared Down**
- Significant overlap (p = 0.00145).
- **Representation Factor**: <span class="highlight">34.1x</span>
- Kromastat largely fails to suppress genes in canine models.
- Includes `FGF21` and `ANGPTL2`.
</div>

![width:575px center](../../results/canine/figures/04_03_venn_down.png)
</div>

---

# Pathway Analysis: UL1 + Kromastat

<div class="columns">
<div class="center">

### Hallmark NES
Classic "Proliferation Crash" observed.
![width:450px](../../results/canine/figures/04_06_hallmark_nes_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
<div class="center">

### GSEA Dotplot
High significance in E2F and MYC targets.
![width:370px](../../results/canine/figures/04_07_gsea_dotplot_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
</div>

---

# Pathway Analysis: CNK89 + Kromastat

<div class="columns">
<div class="center">

### Hallmark NES
Attenuated pathway response; minimal suppression.
![width:500px](../../results/canine/figures/04_06_hallmark_nes_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png)
</div>
<div class="center">


### GSEA Dotplot
Sparse gene engagement consistent with low DEG count.
![width:370px](../../results/canine/figures/04_07_gsea_dotplot_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png)
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
  background-color: #004d2b;
}
</style>

# SUMMARY: THE RESULTS ARE GOOD!

| SUPPRESSED PATHWAYS | PROCESS | EXPECTED RESULTS |
| :--- | :--- | :--- |
| **E2F Targets** | ↓ Proliferation, ↓ Cell cycle, ↓ DNA synthesis | **Tumor growth slows** |
| **MYC Targets** | ↓ Growth, ↓ Ribosome biogenesis, ↓ Metabolism | **Tumor growth slows** |
| **G2M Checkpoint** | ↓ Mitotic division, ↑ DNA damage stress | **Mitotic arrest / Apoptosis** |

---

# Pathway Analysis: UL1 + Romidepsin

<div class="columns">
<div class="center">

### Hallmark NES
Broader pathway modulation than Kromastat.
![width:450px](../../results/canine/figures/04_06_hallmark_nes_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
<div class="center">

### GSEA Dotplot
Stronger immune-suppressive footprint.
![width:390px](../../results/canine/figures/04_07_gsea_dotplot_UL1_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
</div>

---

# Pathway Analysis: CNK89 + Romidepsin

<div class="columns">
<div class="center">

### Hallmark NES
Retained growth suppression despite aggressive background.
![width:475px](../../results/canine/figures/04_06_hallmark_nes_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
<div class="center">

### GSEA Dotplot
Apoptosis + KRAS stress activation.
![width:400px](../../results/canine/figures/04_07_gsea_dotplot_CNK89_Romidepsin_6nM_vs_DMSO_Romidepsin.png)
</div>
</div>

---

# Comparison: SAME PATTERN

Kromastat and Romidepsin produce **broadly similar pathway-level** responses in UL1, with Kromastat nearly silent in CNK89.

- **Romidepsin**: Stronger transcriptional modulation across shared pathways.
- **Kromastat**: More selective; essentially inactive in aggressive canine background.

<div class="columns">
<div>

![width:450px](../../results/canine/figures/04_06_hallmark_nes_combined.png)
</div>
<div>

![width:450px](../../results/canine/figures/04_07_gsea_dotplot_combined.png)
</div>
</div>

---

# Head-to-Head: Functional Comparison

<div class="columns">

<div class="center">

### UL1 Comparison
![width:530px center](../../results/canine/figures/04_06_hallmark_nes_UL1_Romidepsin_6nM_vs_Kromastat_6nM.png)
</div>

<div class="center">

### CNK89 Comparison
![width:530px center](../../results/canine/figures/04_06_hallmark_nes_CNK89_Romidepsin_6nM_vs_Kromastat_6nM.png)
</div>
</div>

---

# Consistency Analysis: Shared DEGs

Intersections of DEGs across treatment groups confirm a shared core response program, though Kromastat's activity is largely confined to the UL1 lineage.

<div class="columns">
<div>

| Treatment | Cell Line | Up | Down | Total |
| :--- | :--- | :---: | :---: | :---: |
| **Romidepsin** | UL1 | 1,320 | 413 | 1,733 |
| **Romidepsin** | CNK89 | 595 | 27 | 622 |
| **Kromastat** | UL1 | 359 | 15 | 374 |
| **Kromastat** | CNK89 | 19 | 0 | 19 |

</div>

![width:550px](../../results/canine/figures/04_09_upset_consistency.png)
</div>

- **Core Signature**: Romidepsin drives the majority of the shared signature; Kromastat contributes minimally in CNK89.
- **Top Conserved Targets**: `ENPP2`, `ELAPOR1`, `HSD17B6`, `HPCAL4`, and `ANKRD29`.

---

# Transcriptome Response Consistency

LFC correlation between UL1 and CNK89 confirms mechanism conservation, though with lower concordance than the human model — reflecting Kromastat's lineage-specific inactivity.

<div class="columns">
<div>

**Pearson R**: Values computed dynamically per drug. High Romidepsin correlation indicates a conserved cytotoxic program; attenuated Kromastat correlation reflects near-absent CNK89 activity.

</div>

![width:600px](../../results/canine/figures/04_10_lfc_correlation.png)
</div>

---

# General Conclusions

- **Potency**: Romidepsin is a more aggressive cytotoxic agent in both canine cell lines.
- **Kromastat Resistance**: CNK89 shows near-complete Kromastat resistance (only 19 DEGs).
- **Mechanism**:
  - **Romidepsin**: Retained broad \"Proliferation Crash\" in both lineages.
  - **Kromastat**: Effective in Indolent (UL1), largely inactive in Aggressive (CNK89).
- **Therapeutic Implication**: Kromastat is not suitable as monotherapy for aggressive canine lymphoma; Romidepsin remains the candidate for both lineages.

---

<!-- _class: lead invert -->

# Thank You!
## Questions?
