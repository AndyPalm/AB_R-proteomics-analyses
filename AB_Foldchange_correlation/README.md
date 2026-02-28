# AB Foldchange Correlation

## Overview
This analysis compares fold changes between two independent experimental runs or datasets. It calculates correlation metrics and generates scatter plots to assess reproducibility and consistency of quantification across replicates.

## Inputs
- **Two comparison result files** (CSV): MSstats or FragPipe comparison results from independent runs
- **Protein labels** (optional): Gene names or custom identifiers for highlighting

## Methods
1. **Data Cleaning**: Removes rows with missing fold changes or p-values
2. **Infinite Value Handling**: Converts #NAME?, Inf, and -Inf to ±10 for visualization
3. **Data Merging**: Joins fold changes from both runs on Protein and Label
4. **Correlation Analysis**: Calculates Pearson and Spearman correlations
5. **Linear Regression**: Fits linear model to assess slope and R-squared
6. **Visualization**: Creates scatter plot with regression line and correlation statistics

## Key Parameters
- **Infinite value truncation**: #NAME? → -10, Inf → +10, -Inf → -10
- **Correlation methods**: Pearson and Spearman
- **Regression model**: log2FC_run2 ~ log2FC_run1
- **Filtering**: Removes rows with missing issues or complete missing data

## Outputs
- Correlation coefficients (Pearson and Spearman)
- Linear regression statistics (slope, intercept, R-squared)
- Scatter plot with regression line
- CSV export of paired fold changes
