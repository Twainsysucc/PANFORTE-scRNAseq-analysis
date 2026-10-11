# PANFORTE-scRNAseq-analysis
Single-cell RNA-seq analysis and figure-generation code for the PANFORTE phase II trial.

# PANFORTE: Single-cell RNA-seq Analysis

This repository is intended to provide the analysis and figure-generation code for the exploratory single-cell RNA sequencing analyses conducted as part of the PANFORTE phase II clinical trial.

# Study overview

Pretreatment tumor samples from 48 patients with pMMR/MSS low-lying locally advanced rectal cancer were analyzed, including 32 complete responders (CR) and 16 non-complete responders (non-CR). Following quality control, 229,437 cells were retained for downstream analyses.

# Analysis workflow

The computational analyses include:

- Single-cell quality control and data integration
- Cell-type annotation and patient-level composition analysis
- MiloR differential neighborhood abundance analysis
- T-cell enrichment and functional-state analysis
- CellChat-based cell-cell communication analysis
- Clinical matching and quality-control sensitivity analyses

# Data availability

The single-cell RNA-seq datasets are deposited in the China National GeneBank DataBase (CNGBdb) under accession number CNP0009362.

The deposited datasets are currently being updated and are scheduled for release after December 30, 2026, subject to controlled-access requirements and institutional privacy restrictions.

# Code availability

The analysis and figure-generation scripts are being progressively curated and deposited in this repository.

**Repository status:** Under preparation. The complete set of analysis scripts is not yet available.

# Software

Analyses were performed using R 4.3.1, Seurat 5.4.0, Harmony 1.2.4, miloR 1.8.1, speckle 1.0.0, CellChat 1.6.1, Monocle 2.28.0, celda 1.16.1, MatchIt 4.8.0, OmicVerse 1.6.8, and Scanpy 1.10.3.
