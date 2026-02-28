# DIA Feature Extraction and Plotting

## Overview
This analysis extracts feature-level (peptide/transition) data for proteins of interest from DIA-NN output and generates visualizations of intensity distributions across replicates and conditions.

## Inputs
- **Feature-level data** (CSV): DIA-NN output with individual peptide/transition measurements
- **Protein-level data** (CSV): Aggregated protein intensities
- **UniProt mapping file** (XLSX): Protein ID to gene name conversion
- **Protein list** (vector): Gene names of proteins to extract

## Methods
1. **Data Loading**: Reads feature-level CSV and converts to Parquet for efficient processing
2. **Gene Mapping**: Joins with UniProt database to add gene names
3. **Filtering**: Extracts features for specified proteins of interest
4. **Validation**: Checks for matches and reports unmatched genes
5. **Visualization**: Generates plots of feature-level intensity distributions

## Key Parameters
- **Feature subset**: Top 4 features per protein (from DIA-NN processing)
- **Imputation**: Handles missing values in intensity data
- **Detection threshold**: Counts valid measurements per replicate

## Outputs
- Filtered feature-level dataset for proteins of interest
- Intensity distribution plots across replicates
- Summary statistics on feature detection rates

## Dependencies
This analysis sources the following R functions from the `src/` directory:
- **src/DIA_feature_extract_barplots.R** - Generates bar plots of feature-level intensities across all replicates
- **src/DIA_feature_extract_barplots_by_group.R** - Generates bar plots of feature-level intensities stratified by treatment group
