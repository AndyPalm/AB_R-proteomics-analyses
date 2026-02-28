# Expression-Labeling Merge Analysis

## Overview
This comprehensive analysis merges proximity labeling data with expression data from two independent experiments. It performs ANCOVA (Analysis of Covariance) to identify proteins with altered proximity while controlling for expression changes, and generates multiple visualization types.

## Inputs
- **Experiment 1 (Labeling)**: Comparison results and protein-level data from proximity labeling
- **Experiment 2 (Expression)**: Comparison results and protein-level data from expression analysis
- **UniProt mapping file** (XLSX): Protein ID to gene name conversion
- **Target protein list** (CSV): Proteins to highlight in visualizations

## Methods
1. **Data Merging**: Joins labeling and expression fold changes on ProteinID
2. **ANCOVA**: Replicate-level analysis with expression intensity as covariate
3. **Proximity Shift Scoring**: Calculates adjusted proximity changes corrected for abundance
4. **Visualization**: Creates bubble plots, scatter plots, and waterfall plots

## Key Parameters
- **Comparison filters**: Specific condition pairs (e.g., ALOD4_DF-vs-UF, OlyA_DF-vs-UF)
- **Covariate**: Mean expression intensity
- **Multiple testing correction**: Benjamini-Hochberg method
- **Visualization thresholds**: Fold change and p-value cutoffs for highlighting

## Outputs
- Combined fold change dataset (labeling vs. expression)
- Proximity shift scores with statistical significance
- Scatter plots showing ANCOVA results
- Waterfall plots ranking proteins by proximity shift
- CSV exports of merged data and proximity scores

## Dependencies
This analysis sources the following R functions from the `src/` directory:
- **src/functions.R** - General utility functions
- **src/plot_labeling_vs_expression_colors.R** - Generates bubble plot comparing labeling vs. expression fold changes
- **src/protein_anova_analysis.R** - Performs replicate-level ANCOVA with expression as covariate
- **src/plot_ancova_scatter.R** - Creates scatter plot visualization of ANCOVA results
- **src/plot_ancova_waterfall.R** - Creates waterfall plot ranking proteins by proximity shift
