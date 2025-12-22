## 🧬 Vaginal Microbiota Analysis

Important note: The data used in this repository is completely fictitious and has been generated for educational and methodological purposes. No real patient data or sensitive information is used, in compliance with data protection regulations.



## 📌 Project aim
This repository documents the complete process of statistical analysis of vaginal microbiota data, from the structure and nature of the data to the application of appropriate statistical methods to study microbial diversity and its association with different covariates.

The main objective is to show how to analytically approach microbiota data, taking into account its particular characteristics (compositional data, high dimensionality, inflation of zeros, and overdispersion), using a rigorous and reproducible statistical approach.


## Scripts folder

The `scripts/` directory includes all analysis code.

- **SimulatedMicrobiome.R**: generates simulated longitudinal microbiome data used throughout the project.
- **Models/**: contains the scripts for statistical modeling (Weighted model LM, CLMM, PERMANOVA) and visualization of selected results.



## 🧭 Analysis strategy
This repository presents a **reproducible data science workflow for microbiome analysis**, applied to simulated vaginal microbiota data. Sequencing count data are modeled as **compositional OTU abundance matrices** and analyzed using normalization, zero imputation, and **Aitchison log‑ratio transformations (CLR/ILR)**. Microbial diversity is assessed via **alpha diversity (Shannon index)** and **beta diversity (Bray–Curtis dissimilarity)**. Multivariate structure is explored using **PCoA, NMDS, and hierarchical clustering**. Group differences are tested with **permutation‑based methods (PERMANOVA, PERMDISP)**, while temporal and covariate effects are modeled using **linear and mixed‑effects models**, including multinomial extensions for taxonomic composition. The project reflects best practices in **statistical learning for microbiome data**, with an emphasis on robustness, interpretability, and reproducibility.



## 🔁 Analysis Pipeline

```text
Raw microbiome sequencing counts (simulated)
                │
                ▼
OTU / taxonomic abundance matrix
                │
                ▼
Zero handling & normalization
                │
                ▼
Log‑ratio transformation (CLR / ILR)
                │
                ▼
Exploratory analysis
                │
                ├─► Alpha diversity (Shannon index)
                │       └─ Friedman / Wilcoxon tests
                │
                ├─► Beta diversity (Bray–Curtis)
                │       ├─ NMDS / PCoA
                │       └─ PERMANOVA / PERMDISP
                │
                └─► Taxonomic composition
                        ├─ Hierarchical clustering (Ward)
                        ├─ Ordination methods
                        └─ Heatmaps
                │
                ▼
Statistical modeling
                │
                ├─ Mixed‑effects models (LME)
                ├─ Linear models
                ├─ Multinomial models (CLM / CLMM)
                └─ Forward variable selection
                │
                ▼
Interpretation, visualization & reporting
```


## 📦 Software and packages

All analyses were performed in R (v4.0.3), primarily using the following packages:

* **vegan**: `diversity`, `vegdist`, `metaMDS`, `adonis2`, `betadisper`
* **stats**: `hclust`, `cmdscale`
* **phyloseq**: `plotbar`, `plotheatmap`
* **ComplexHeatmap**: heatmaps
* **leaps**: selección de variables (`regsubsets`)
* **ordinal**: `clm`, `clmm2` models


👩‍💻 Author
Project developed for academic purposes and methodological demonstration.


