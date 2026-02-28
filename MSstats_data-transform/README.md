# MSstats Data Transform

## Overview
This analysis converts MSstats comparison results into separate dataframes for each comparison condition. It transforms the data into a format suitable for downstream visualization and analysis, handling infinite values and missing data appropriately.

## Inputs
- **MSstats comparison results** (CSV): Output from MSstats statistical analysis with fold changes and p-values

## Methods
1. **Data Extraction**: Identifies unique comparison labels from input data
2. **Subsetting**: Creates separate dataframe for each comparison condition
3. **Column Standardization**: Renames columns to consistent format (ProteinID, Label, FC, p_value, adj.pvalue, etc.)
4. **Infinite Value Handling**: Converts Inf to +10 and -Inf to -10 for visualization compatibility
5. **Data Export**: Saves each comparison to separate CSV file

## Key Parameters
- **Infinite value truncation**: Inf → +10, -Inf → -10
- **Column mapping**: log2FC → FC, pvalue → p_value, adj.pvalue → adj.pvalue
- **Comparison identifier**: Label column from MSstats output

## Outputs
- Separate CSV files for each comparison condition
- Standardized column format across all files
- Handled infinite values for downstream analysis
