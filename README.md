# DeepClone Protocol — Figures & Analysis Code

This repository contains the code and data to reproduce the figures presented in the **protocols.io** paper:

> **DeepClone, an end-to-end protocol to study somatic mutagenesis and selection at high resolution**

## Abstract

With the introduction of DNA duplex sequencing technologies, somatic mutations that occur in one or few cells in a sample can be detected. This has opened up the study of somatic evolution in scenarios that were not previously easy to access, such as normal tissues. The general adoption of DNA duplex sequencing technologies has been hindered, among other reasons, by the lack of a complete ecosystem of computational tools that support end-to-end analysis from the sequencing raw data to the quantification of multiple aspects of mutagenesis and selection. Here, we present DeepClone, an end-to-end protocol to build DNA duplex sequencing libraries and to analyze the product of their sequencing, aimed at working with mammalian nuclear genomes. DeepClone includes one computational pipeline to efficiently call somatic mutations in a sample from DNA duplex sequencing reads, and a second pipeline that automates standard calculations carried out on those mutations across a cohort of samples. These include the deconvolution of the mutational processes active across samples, the estimation of positive selection on the mutations in genes, the influence of the exposure to exogenous agents on the clonal landscape of the samples, and the study of natural saturation mutagenesis across genes. DeepClone includes extensive quality control metrics associated with every step of the protocol that facilitate decision-making and critical assessment of its results. Any DNA duplex library preparation approach can be plugged to the computational pipelines, which guarantees the versatility of DeepClone to address different questions on the somatic evolution of tissues.

## Related repositories

| Repository | Description |
|---|---|
| [bbglab/deepUMIcaller](https://github.com/bbglab/deepUMIcaller) | Somatic mutation caller for UMI-based duplex sequencing data |
| [bbglab/deepCSA](https://github.com/bbglab/deepCSA) | Clonal Structure Analysis pipeline for duplex sequencing data |

## Repository structure

```
.
├── main/                        # Main figures
│   └── figure3/                 # Figure 3
│       ├── scripts/             # Jupyter notebooks
│       └── plots/               # Generated PDF plots
│
├── bladder/                     # Bladder urothelium case-study analyses
│
├── supplementary/               # Extended Data figures
│   ├── extended1/               # Extended Data Figure 1
│   ├── extended2/               # Extended Data Figure 2
│   ├── extended6/               # Extended Data Figure 6
│   └── extended7/               # Extended Data Figure 7
│
├── other_analysis/              # Supplementary analyses
│   ├── dupcaller/               # DupCaller diagnostic notebooks
│   ├── duplexome/               # Duplexome read count analysis
│   └── UMIcollisions/           # UMI collision frequency analysis
│
├── wetlab_metrics/              # Wet-lab metric exploration notebooks
│
├── test_datasets/               # Test data for the deepCSA pipeline
│   └── deepCSA/
│
└── data/                        # Shared reference data and QC reports
```

## Test datasets

The [`test_datasets/deepCSA/`](test_datasets/deepCSA/) directory contains a small subset (3 samples) of normal bladder urothelium data for use with the [deepCSA](https://github.com/bbglab/deepCSA) pipeline test suite. The data originates from:

> Calvet, F., Blanco Martinez-Illescas, R. et al. "Sex and smoking bias in the selection of somatic mutations in human bladder." *Nature* (2025). https://doi.org/10.1038/s41586-025-09521-x

The full dataset is available at Zenodo: [10.5281/zenodo.15836679](https://zenodo.org/records/15836679)

See [`test_datasets/deepCSA/README.txt`](test_datasets/deepCSA/README.txt) for details on how to use the test data.

## Requirements

### Python

The Jupyter notebooks are written in Python and primarily use:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `statsmodels`
- `scikit-learn`
- `matplotlib-venn`

### R

The trinucleotide-corrected mutation rate script ([`cordblood_mutrate_genome_trint_corrected.R`](supplementary/extended1/scripts/cordblood_mutrate_genome_trint_corrected.R)) requires:

- `tidyr`, `dplyr`, `stringr`
- `ggplot2`
- `Hmisc`
- `jsonlite`
- `Biostrings`
- `data.table`

