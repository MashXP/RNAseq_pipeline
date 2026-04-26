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
    color: #5a0087;
    font-size: 1.8em;
    border-bottom: 2px solid #5a0087;
    margin-bottom: 0.5em;
  }
  h2 {
    color: #7b1fa2;
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
    color: #5a0087;
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
  .highlight-blue {
    color: #004a99;
    font-weight: bold;
  }
  .highlight-green {
    color: #00663a;
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
    background-color: #5a0087;
    color: #ffffff;
    padding: 12px 15px;
    text-align: left;
    border-bottom: 3px solid #3d0060;
    font-weight: 700;
  }
  td {
    padding: 10px 15px;
    border-bottom: 1px solid #eef2f7;
    color: #444;
    vertical-align: top;
  }
  tr:nth-child(even) {
    background-color: #faf4ff;
  }
  tr:hover {
    background-color: #f0e6ff;
  }

---

<!-- _class: lead invert -->

# Cancer Treatments Project
## Cross-Species Comparative Analysis
### KROMASTAT and ROMIDEPSIN: Human × Canine

---

# STUDY CONTEXT

This report is the **sequel** to the Human and Canine individual transcriptomic reports.

<div class="columns">
<div>

### Premise
Do Kromastat and Romidepsin induce **conserved biological programs** across mammalian species?

- **Human**: H9 (Indolent) · SUPM2 (Aggressive)
- **Canine**: UL1 (Indolent) · CNK89 (Aggressive)
- Matched pairs by disease phenotype for cross-species validity.

</div>
<div>

### Question Framework
1. **Gene-level**: Which ortholog responses are conserved?
2. **Pathway-level**: Is the Hallmark signature shared?
3. **Drug-level**: Does conservation differ between Romidepsin and Kromastat?

</div>
</div>

---

# RECAP: INDIVIDUAL DRUG RESPONSE

<div class="columns">
<div>

### <span class="highlight-blue">Human</span>
| Drug | H9 DEGs | SUPM2 DEGs |
|:---|:---:|:---:|
| Romidepsin | 3,368 | 4,256 |
| Kromastat | 1,012 | 1,463 |

- Romidepsin: stronger potency in both lines.
- Kromastat: selective, spares immune signaling.

</div>
<div>

### <span class="highlight-green">Canine</span>
| Drug | UL1 DEGs | CNK89 DEGs |
|:---|:---:|:---:|
| Romidepsin | 1,733 | 622 |
| Kromastat | 374 | **19** |

- Romidepsin: retained activity in both lines.
- Kromastat: near-silent in aggressive CNK89.

</div>
</div>

> **Key divergence**: Kromastat shows lineage-dependent resistance in the canine model absent in the human model.

---

# METHODOLOGY: Cross-Species Pipeline

<div class="columns">
<div>

1. **Gene Mapping**: Human–Canine one-to-one orthologs via Ensembl.
2. **Alluvial Flow (04_11)**: Categorize concordant, discordant, and species-specific responses per matched pair.
3. **Ortholog Scatter (04_12)**: Pearson correlation of LFC values across orthologs.
4. **NES Overlap (04_13)**: Side-by-side Hallmark pathway enrichment for matched cell lines.

</div>
<div>

### Matched Pairs
| Human | Canine | Phenotype |
|:---|:---|:---|
| H9 | UL1 | Indolent |
| SUPM2 | CNK89 | Aggressive |

</div>
</div>

---

# OVERVIEW: GENE FLOW LANDSCAPE

Main alluvial plot — all matched pairs, both drugs simultaneously.

---

![bg width:1200px center](../../results/comparative/figures/04_11_alluvial_gene_flow.png)

---

# GENE FLOW SUMMARY

<style scoped>
table { font-size: 0.7em; }
</style>

| Drug | Pair | Up ↑ | Down ↓ | Disc. | Human | Canine |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Romidepsin** | H9 / UL1 | **258** | **112** | 35 | 864 | 528 |
| **Romidepsin** | SUPM2 / CNK89 | **167** | **4** | 15 | 1103 | 237 |
| **Romidepsin** | *Global ortholog* | *331* | *22* | — | — | — |
| **Kromastat** | H9 / UL1 | **61** | **1** | 2 | 387 | 129 |
| **Kromastat** | SUPM2 / CNK89 | **4** | **0** | 0 | 454 | 10 |
| **Kromastat** | *Global ortholog* | *103* | *0* | — | — | — |

- **Romidepsin**: High conservation, particularly indolent suppression (112 genes).
- **Kromastat**: Conservation limited to indolent pair; aggressive core collapses.
- **Fidelity**: Discordance <10% — mechanism direction is highly conserved.

---

# PATHWAY FLOW: Romidepsin — Indolent Pair

H9 (Human) vs UL1 (Canine) — Romidepsin 6nM

<div class="columns">
<div>

| Pathway Status | Count |
|:---|:---:|
| Concordant Suppressed | **8** |
| Concordant Activated | **6** |
| Discordant | 1 |
| Human-only significant | 12 |
| Canine-only significant | 10 |

- **14 concordant pathways** including the full Proliferation Crash core.
- Discordance minimal: only Mitotic Spindle (reversed direction).

</div>

![width:570px](../../results/comparative/figures/04_11_alluvial_Romidepsin_H9_vs_UL1.png)
</div>

---

# PATHWAY FLOW: Romidepsin — Aggressive Pair

SUPM2 (Human) vs CNK89 (Canine) — Romidepsin 6nM

<div class="columns">
<div>

| Pathway Status | Count |
|:---|:---:|
| Concordant Suppressed | **8** |
| Concordant Activated | **11** |
| Discordant | 0 |
| Human-only significant | 20 |
| Canine-only significant | 3 |

- **19 concordant pathways** — strongest inter-species agreement of all pairs.
- Zero discordant: direction is perfectly preserved where both species respond.

</div>

![width:570px](../../results/comparative/figures/04_11_alluvial_Romidepsin_SUPM2_vs_CNK89.png)
</div>

---

# PATHWAY FLOW: Kromastat — Indolent Pair

H9 (Human) vs UL1 (Canine) — Kromastat 6nM

<div class="columns">
<div>

| Pathway Status | Count |
|:---|:---:|
| Concordant Suppressed | **8** |
| Concordant Activated | **20** |
| Discordant | 1 |
| Human-only significant | 8 |
| Canine-only significant | 5 |

- **28 concordant pathways** — Kromastat's pathway breadth is well-conserved in the indolent context.
- Single discordant: Mitotic Spindle (H9=Activated, UL1=Suppressed).

</div>

![width:570px](../../results/comparative/figures/04_11_alluvial_Kromastat_H9_vs_UL1.png)
</div>

---

# PATHWAY FLOW: Kromastat — Aggressive Pair

SUPM2 (Human) vs CNK89 (Canine) — Kromastat 6nM

<div class="columns">
<div>

| Pathway Status | Count |
|:---|:---:|
| Concordant Suppressed | **9** |
| Concordant Activated | **21** |
| Discordant | 0 |
| Human-only significant | 7 |
| Canine-only significant | 4 |

- **30 concordant pathways** despite near-absent gene-level overlap — CNK89 responds at the pathway level even with only 19 DEGs.
- Pathway-level conservation far exceeds gene-level conservation in this aggressive pair.

</div>

![width:570px](../../results/comparative/figures/04_11_alluvial_Kromastat_SUPM2_vs_CNK89.png)
</div>

---

# CONSERVED CROSS-SPECIES SIGNAL

These drug-bridge plots track each pathway across **both drugs and both species simultaneously** for each cell type pair. A pathway that is `concordant_activated` (red) or `concordant_suppressed` (blue) throughout the full 4-axis flow is stable across **treatment and biology** — the strongest conservation claim possible.

---
<div class="columns">
<div class="center">

### Indolent (H9 ↔ UL1)
![width:500px](../../results/comparative/figures/04_11_alluvial_indolent.png)
</div>
<div class="center">

### Aggressive (SUPM2 ↔ CNK89)
![width:500px](../../results/comparative/figures/04_11_alluvial_aggressive.png)
</div>
</div>

- **Concordant flows** (red/blue) represent pathways conserved in *both direction and state* regardless of which drug or which species.
- **Kromastat-specific divergence** in the Aggressive pair is expected: CNK89's near-absent gene-level response attenuates but does not eliminate pathway-level concordance.
- Flows remaining concordant in **both panels** constitute the strongest candidates for pan-species HDACi biomarkers.

---



The following Hallmark pathways are concordant across **all four cell lines** and **both drugs**.

<div class="columns">
<div>

### Universally Suppressed ↓
| Pathway |
|:---|
| MYC Targets V1 |
| MYC Targets V2 |
| E2F Targets |
| G2M Checkpoint |
| mTORC1 Signaling |
| Unfolded Protein Response |

</div>
<div>

### Universally Activated ↑
| Pathway |
|:---|
| Myogenesis |
| Heme Metabolism |
| P53 Pathway |

> **These 9 pathways define the pan-species, pan-drug anti-lymphoma core signature.**

</div>
</div>

---

# PATHWAY NES OVERLAP: H9 × CNK89

Cross-aggressiveness pair — Indolent Human (H9) vs Aggressive Canine (CNK89)

<div class="columns">
<div class="center">

### Romidepsin
![width:520px](../../results/comparative/figures/04_13_nes_bar_H9_CNK89_Romidepsin.png)
</div>
<div class="center">

### Kromastat
![width:520px](../../results/comparative/figures/04_13_nes_bar_H9_CNK89_Kromastat.png)
</div>
</div>

---

# PATHWAY NES OVERLAP: SUPM2 × UL1

Cross-aggressiveness pair — Aggressive Human (SUPM2) vs Indolent Canine (UL1)

<div class="columns">
<div class="center">

### Romidepsin
![width:520px](../../results/comparative/figures/04_13_nes_bar_SUPM2_UL1_Romidepsin.png)
</div>
<div class="center">

### Kromastat
![width:520px](../../results/comparative/figures/04_13_nes_bar_SUPM2_UL1_Kromastat.png)
</div>
</div>

---

# ORTHOLOG SCATTER: TOP CONCORDANT GENES

<div class="columns">
<div class="center">

### Romidepsin
Conserved activation across species.
![width:520px](../../results/comparative/figures/04_12_ortholog_scatter_Romidepsin_6nM.png)
</div>
<div class="center">

### Kromastat
Conserved activation — indolent-dominated.
![width:520px](../../results/comparative/figures/04_12_ortholog_scatter_Kromastat_6nM.png)
</div>
</div>

---

# ORTHOLOG HEATMAP: Conserved Targets

Heatmap of top orthologous gene pairs with concordant LFC directionality across species.

![width:900px center](../../results/comparative/figures/04_12_conserved_ortholog_heatmap.png)

---

# DRUG CONSERVATION COMPARISON

<div class="columns">
<div>

### Romidepsin
- **331 shared up**, **22 shared down** globally.
- Concordant across both Indolent and Aggressive pairs.
- Discordance confined to ~9%.
- Pathway core: **fully preserved** in both species.

</div>
<div>

### Kromastat
- **103 shared up**, **0 shared down** globally.
- Conservation **collapses** in Aggressive pair (4 genes).
- Pathway concordance high in Indolent; absent in Aggressive.
- Reveals a **species × lineage interaction** not observed with Romidepsin.

</div>
</div>

> **Conclusion**: Romidepsin is the more translationally reliable agent. Kromastat's conservation is contingent on the indolent phenotype.

---

# SPECIES-SPECIFIC DIVERGENCES

<div class="columns">
<div>

### Human-only responses (Romidepsin)
- Interferon-γ / Interferon-α suppression
- Cholesterol Homeostasis activation
- Complement activation
- Broader KRAS pathway engagement

</div>
<div>

### Canine-only responses (Romidepsin)
- Xenobiotic Metabolism activation
- Fatty Acid Metabolism activation
- Reactive Oxygen Species pathway activation
- PI3K/Akt/mTOR suppression

</div>
</div>

> Species-specific divergences likely reflect differences in basal metabolic and immune gene expression architecture rather than drug target specificity.

---

# OVERALL CONCLUSIONS

<div class="small-text">

| Finding | Romidepsin | Kromastat |
|:---|:---|:---|
| **Cross-species concordance** | High (331 up, 22 down) | Moderate (103 up, 0 down) |
| **Aggressive pair conservation** | Retained (167 up, 4 down) | Collapsed (4 up, 0 down) |
| **Universal pathway core** | Full 9-pathway signature | Full 9-pathway signature |
| **Species divergence** | Immune & metabolic modules | Immune activation absent in canine |
| **Translational reliability** | ✅ High | ⚠️ Lineage-dependent |

</div>

- The **"Proliferation Crash" core** (E2F, MYC, G2M) is unambiguously pan-mammalian.
- Kromastat's efficacy in aggressive disease is a **canine-specific gap**, not a mechanism failure.
- Canine models are valid surrogates for Romidepsin, but may **underestimate Kromastat** potency in aggressive human lymphoma.

---

<!-- _class: lead invert -->

# Thank You!
## Questions?
