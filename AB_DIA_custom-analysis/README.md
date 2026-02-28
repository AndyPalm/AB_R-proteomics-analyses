# AB DIA Custom Analysis

## Overview
This analysis visualizes intensity distributions from MSstats model QC output. It creates histograms comparing protein abundance distributions between treatment groups, with optional stratification by protein subsets (e.g., cell surface vs. non-cell surface proteins).

## Inputs
- **MSstats model QC output** (CSV): Fitted intensity values from MSstats data processing
- **Protein subset list** (XLSX): Optional list of proteins in specific categories (e.g., PM surface proteins)

## Methods
1. **Data Aggregation**: Groups by protein and treatment condition, calculates mean fitted intensities
2. **Subset Filtering**: Separates proteins into categories (e.g., PM surface vs. non-PM surface)
3. **Histogram Generation**: Creates inverted histograms showing distributions
4. **Visualization**: Plots subset proteins in positive direction, non-subset in negative direction
5. **Faceting**: Separates plots by treatment group for comparison

## Key Parameters
- **Grouping variables**: Protein and GROUP (treatment condition)
- **Aggregation metric**: Mean of fitted values
- **Histogram bins**: 30 bins for distribution visualization
- **Inversion**: Negative y-axis for non-subset proteins (inverted histogram)

## Outputs
- Inverted histograms showing protein abundance distributions
- Stratified comparison of protein subsets vs. total population
- Faceted plots by treatment group
- Summary statistics on protein counts by category
