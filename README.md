# Oocyte Isoform Transcriptome Analysis

[![R-version](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the R source code for data analysis and visualization presented in the study of human and mouse oocyte isoform landscapes.
Single oocyte full-length isoform sequencing unveils the impact of transposable elements on RNA diversity and stability during oocyte maturation (https://www.biorxiv.org/content/10.1101/2025.06.17.659919v1.full)

## 🔬 Project Summary
This project explores the transcriptomic diversity in oocytes using Long-read sequencing. Key analyses include:
- Identification of stage-specific isoforms (GV, MI, MII).
- Differential Transcript Usage (DTU) and Isoform Switching analysis.
- Characterization of Transposable Element (TE) derived transcripts.
- Multi-incoproration, ActD half-lives, Enrichment of TE-derived isoform, RBP binding prediction.

## 📂 Repository Structure

| Script Name | Main Figures | Analysis Description |
| :--- | :--- | :--- |
| `01.plot_F1.R` | Fig 1 | Expression fractions. |
| `02.plot_F2.R` | Fig 2 | SQANTI3 classification (FSM, ISM) distribution, Coding potential analysis (CPAT). |
| `03.plot_F3.R` | Fig 3 | SQANTI3 classification (FSM, ISM) distribution, Coding potential analysis (CPAT). |
| `04.plot_F4.R` | Fig 4 | Isoform switching analysis via `IsoformSwitchAnalyzeR`. |
| `05.plot_F5.R` | Fig 5 | TE-derived isoform classification of repeat class/repeat family. |
| `06.plot_F6.R` | Fig 6 | Multi-incoproration, ActD half-lives, Enrichment of TE-derived isoform, RBP binding prediction. |

## 🛠 Prerequisites

### R Environment
- **R Version**: >= 4.3.0
- **Main Packages**: `tidyverse`, `data.table`, `ggplot2`, `ggpubr`, `IsoformSwitchAnalyzeR`, `ComplexHeatmap`.

## 📊 Data Availability
The processed source data (TPM matrices, isoform annotations, etc.) required to run these scripts are provided as **Source Data** files accompanying the published version of this article in ***Nature Communications***.

### Installation
```R
install.packages(c("tidyverse", "data.table", "ggpubr", "viridis"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("IsoformSwitchAnalyzeR", "ComplexHeatmap"))

## 📄 Citation
If you use this code, please cite our study:

```bibtex
@article{author2025,
  title = {Single oocyte full-length isoform sequencing unveils the impact of transposable elements on RNA diversity and stability during oocyte maturation},
  author = {Wang et al.},
  journal = {bioRxiv},
  year = {2025},
  doi = {https://doi.org/10.1101/2025.06.17.659919},
  url = {[https://www.biorxiv.org/content/10.1101/2025.06.17.659919v1.full](https://www.biorxiv.org/content/10.1101/2025.06.17.659919v1.full)}
}
