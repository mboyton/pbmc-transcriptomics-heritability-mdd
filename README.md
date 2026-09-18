# PBMC transcriptomics, heritability enrichment and Mendelian randomisation

Analysis code for the study **“NK Cell Activation as a Neuroimmune Mechanism in Depression: Insights from Single-Cell Multi-Omics.”**

This repository contains scripts used to integrate single-cell PBMC transcriptomics with genetic enrichment and Mendelian randomisation analyses of major depressive disorder and psoriasis.

## Repository structure

- `01_lawlor_preprocessing/` — preprocessing of single-cell PBMC data
- `02_gwas_preprocessing/` — preprocessing of MDD and psoriasis GWAS summary statistics
- `03_cellect/` — CELLECT heritability enrichment analyses
- `04_transcriptional_programmes/` — differential expression and WGCNA analyses
- `05_mr/` — eQTLGen instrument preparation and transcriptome-wide Mendelian randomisation
- `06_integration_and_figures/` — integration of transcriptomic and genetic results and figure generation

## Reproducibility

Analyses were conducted across local and high-performance computing environments and were not implemented as a single automated workflow. Scripts are provided to document the analytical procedures used to generate the reported results.

File paths have been generalised where appropriate, and large or externally sourced datasets are not redistributed in this repository.
