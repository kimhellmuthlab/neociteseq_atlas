# Neonatal CITE-seq atlas

Custom analysis and visualization scripts accompanying the study:

**"Human neonatal CITE-seq atlas identifies an immune transition at 32 weeks’ gestation from CD15⁺ myeloid-dominated to interferon-primed immunity"**

🔗 Preprint: https://www.biorxiv.org/content/10.64898/2026.04.01.715643v1

---

## Table of Contents
- [Overview](#overview)
- [Graphical Abstract](#graphical-abstract)
- [Abstract](#abstract)
- [Data and Materials Availability](#data-and-materials-availability)
- [Environment Setup](#environment-setup)
- [Usage](#usage)
- [Scripts](#scripts)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview
This repository provides **example analysis and visualization scripts** used in the neonatal CITE-seq atlas study.

The scripts illustrate key analysis steps and reproduce selected figures from the manuscript. 

---

## Graphical Abstract 
![Graphical Abstract](GraphicalAbstract.png)

---

## Abstract
The human neonatal immune system is developmentally specialized to balance the unique requirements of perinatal transition. Disruption of this finely tuned balance, as in preterm birth, may have profound consequences for immunity and overall health. However, the impact of prematurity on immune composition and functional responsiveness across gestational ages (GA) remains incompletely understood. Single-cell profiling has advanced our understanding of neonatal immunity, yet most studies were limited to unimodal readouts, narrow GA windows, or baseline function. Here, we present a comprehensive human neonatal CITE-seq atlas (82 samples from 25 neonates and 10 adults as controls) at the first days of life covering a wide GA range and integrating baseline and stimulated conditions. Most notably, we identify a GA-dependent immune transition point centering around 32 weeks of GA, which discriminates extremely and very preterm neonates (GA <32wks) from those of higher GA (≥32wks). In particular, early-life immunity in extremely and very preterm infants showed CD15+ granulocytic myeloid derived suppressor cell-like predominance, whereas more mature neonates exhibited interferon-primed transcriptional profiles. This was associated with divergent myeloid-to-lymphocyte signaling networks and qualitatively distinct NK- and T-cell bystander responses upon activation. Together, these findings show that intrauterine development imprints GA-specific immune programs. By defining a developmental transition around a GA of 32 weeks that regulates baseline and induced responses of neonatal immune cells, our atlas provides a framework for understanding the vulnerability of preterm infants and thus may pave the way for developing GA-adapted immunomodulatory strategies.

---

## Data and Materials Availability

All code and data associated with this study are present in the paper, the Supplementary Materials, or as indicated below:  

- **Interactive exploration and download will be possible via:** [CELLxGENE portal](https://cellxgene.cziscience.com/)  
- **Raw and processed CITE-seq data matrices will be available at:** [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXX)  
- **Custom scripts for analysis and visualization:** available in this repository

The repository contains representative scripts used for analysis and figure generation.
All other analyses were performed using publicly available software packages as described in the Methods section of the manuscript.

---

## Environment Setup

Analysis was performed in R. Below is a minimal setup:

```r
# CRAN packages
install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "forcats", "tibble",
                   "Matrix", "pheatmap", "svglite", "ggpattern", "hdf5r"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "DESeq2", "enrichR"))

# GitHub dependency
install.packages("remotes")
remotes::install_github("satijalab/MuDataSeurat")
```

---

## Usage

1. Download dataset from CELLxGENE 
2. Store locally (e.g., `data/`)
3. Adjust file paths within scripts 
4. Run scripts 

Example:
```r
source("scripts/PCAplot.R")
```

---

## Scripts

This repository contains example analysis and visualization scripts:

### `PCAplot.R`
- **Purpose:** Generate PCA plot corresponding to Fig. 1B
- **Input:** CELLxGENE dataset (“All - neonatal and adult circulating immune cells baseline and stim”)
- **Description:** Aggregates RNA expression at the sample level and performs principal component analysis

### `ndensPlot.R`
- **Purpose:** Visualize cell distributions in UMAP space
- **Input:** e.g., Fig. 2A: CELLxGENE dataset (“Myeloid cells - neonatal and adult circulating immune cells baseline”)
- **Description:** Produces normalized density plots for selected cell populations and conditions

### `StackedBarPlot_cellproportions.R`
- **Purpose:** Visualize cell type proportions across gestational age groups 
- **Input:** e.g., Fig. 2B: CELLxGENE dataset (“Myeloid cells - neonatal and adult circulating immune cells baseline”)
- **Description:** Generates stacked bar plots summarizing relative cell type abundances

### `Pheatmap.R`
- **Purpose:** Display gene expression patterns
- **Description:** Generates heatmaps of selected marker genes stratified by gestational age

### `DEGs_Analysis_Visualization.R`
- **Purpose:** Differential expression and downstream analysis
- **Input:** e.g., CELLxGENE dataset (“All - neonatal and adult circulating immune cells baseline and stim”), subset for myeloid cells. 
- **Description:** Performs age-stratified and interaction-based differential expression analysis, followed by enrichment testing and visualization

---

## Citation

If you use this code, please cite:

Rothämel P, Mattia A, Corey MI, Puzek B, Wiesel J, Michael-Kuschel P, Klein C, Sperandio M, Henneke P, Nussbaum C, Kim-Hellmuth S.  
**Human neonatal CITE-seq atlas identifies an immune transition at 32 weeks’ gestation from CD15⁺ myeloid-dominated to interferon-primed immunity.**  
bioRxiv (2026).  
https://www.biorxiv.org/content/10.64898/2026.04.01.715643v1

---

## Contact

- Name: Paula Rothämel, MD  
- Email: Paula.Rothaemel@med.uni-muenchen.de  
- Institution: Dr. von Hauner Children's Hospital, LMU Munich 
