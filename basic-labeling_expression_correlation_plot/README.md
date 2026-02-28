# Basic Labeling-Expression Correlation Plot

## Overview
This analysis generates correlation plots comparing fold changes between two experimental datasets. It visualizes the relationship between protein labeling intensity changes and expression level changes across multiple conditions.

## Inputs
- **Two comparison result files** (CSV): Fold change and p-value data from two independent experiments
- **Protein labels file** (CSV): List of proteins of interest for highlighting

## Methods
1. **Data Cleaning**: Removes rows with missing fold change or p-value data
2. **Categorization**: Assigns proteins to categories based on concordance of regulation between the two experiments
3. **Visualization**: Creates scatter plot with color-coded categories and highlighted proteins of interest

## Key Parameters
- **Fold change thresholds**: ±1 (log2 scale) for significance boundaries
- **Concordance categories**: 
  - Both UP (category 1)
  - Both DOWN (category 2)
  - Discordant UP/NO-DOWN (category 3)
  - Discordant DOWN/NO-UP (category 4)
  - Other (category 5)

## Outputs
- Correlation scatter plot with concordance visualization
- CSV file with categorized results and fold change comparisons
