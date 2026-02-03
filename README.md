# Oocyte-long-read-transcriptome

R scripts for the analysis and visualization of human and mouse oocyte full-length isoform landscapes.

## 🔬 Study Overview
This project focuses on:
- Single oocyte full-length isoform sequencing.
- Impact of Transposable Elements (TE) on RNA diversity.
- RNA stability and degradation during oocyte maturation (GV to MII).

## 📂 Repository Structure
- `scripts/`: Contains R scripts for generating Main and Supplementary Figures (Fig 1-6).
- `data/`: (Recommended) Place your metadata and expression matrices here.

## 🛠 Prerequisites
The following R packages are required:
- `tidyverse`, `data.table`, `ggplot2`, `ggpubr`, `RColorBrewer`, `viridis`
- `IsoformSwitchAnalyzeR`, `ComplexHeatmap` (Bioconductor)

> **Note**: Please update the working directory `setwd()` and file paths in the scripts to match your local environment before running.

## 📧 Contact
For questions, please contact **Wang et al.** at [Your Email Address].

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
