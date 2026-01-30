# Differential Expression (RNA-seq) & Quantitative Proteomics (DIA-LFQ)  
**Analysis:** Global Z-score differences and Robust Linear Regression for Discordance Analysis

## Overview
This R Markdown pipeline performs a multi-omics integration of transcriptomic and proteomic datasets. The primary goal is to identify **discordant regulation**—instances where protein abundance diverges significantly from transcript levels, implying post-transcriptional regulation

R 4.4.2

Merging protein abundance and RNA transcript abundance dataframes and associated analyses. Uses replicate-level data and calculates foldchanges + stats.

Uniprot database must be downloaded with at least "Gene Name (primary)", "Gene Names", and the associated PTM databases. This file is also included in the example data.

The workflow merges data based on Gene Names (with synonym rescue), normalizes expression, and models the relationship using **Robust Linear Models (RLM)** to calculate discordance residuals.

## Directory Structure
 The script relies on a standard Box directory structure. You must update the `box_root` variable in the setup chunk if moving to a new machine.

```text
/R_WD
├── Reference_dbs/
│   ├── uniprot_hs_genes_PTM-DBs_17jan26.xlsx       # UniProt Mapping
│   ├── protein-half-life-Savitski-2018_high-qual.xlsx
│   └── Choudhary_ubiquitylation_TableS1.xlsx
├── AP2-327_HAEC-expression/
│   ├── 01_Input_Search_Results/
│   │   └── 327-normalized_top4-features_protein-level-data_impute.csv
│   └── SGH-JJM_RNAseq_UF-vs-DF.csv
└── POCA2/
    └── 03_graphs_plots/                            # Output Directory