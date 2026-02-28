# Total Features & Proteins Graph

## Overview
This analysis visualizes the distribution of feature-level and protein-level data across treatment conditions. It generates comparative plots showing total feature counts and protein abundance distributions, with optional stratification by protein subsets.

## Inputs
- **Feature-level data** (CSV): Individual peptide/transition measurements from DIA-NN
- **Protein-level data** (CSV): Aggregated protein intensities
- **MSstats output** (CSV): Processed feature-level data with condition assignments
- **UniProt mapping file** (XLSX): Protein ID to gene name conversion

## Methods
1. **Data Loading**: Reads feature and protein-level data, converts to Parquet for efficiency
2. **Gene Mapping**: Joins with UniProt database to add gene names
3. **Aggregation**: Counts total features per protein and condition
4. **Visualization**: Creates histograms and distribution plots
5. **Stratification**: Optional breakdown by protein subsets (e.g., cell surface, ER proteins)

## Key Parameters
- **Feature subset**: Top 4 features per protein (from DIA-NN processing)
- **Imputation**: Handles missing values in intensity data
- **Binning**: Histogram bin width for distribution visualization
- **Stratification**: Optional protein category filtering

## Outputs
- Distribution plots of feature counts across conditions
- Protein abundance histograms
- Summary statistics on feature detection rates
- Stratified comparisons by protein subset
