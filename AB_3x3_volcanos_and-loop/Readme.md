# AB 3x3 Volcanos and Loop

## Overview
This analysis generates volcano plots from 3x3 experimental designs (3 control + 3 treatment replicates). It converts MSstats comparison results into a format resembling FragPipe replicate-level output and loops through multiple comparison dataframes to generate consistent volcano plots.

## Inputs
- **MSstats comparison results** (CSV): Statistical comparison output with fold changes and p-values
- **Multiple comparison dataframes** (CSV): One file per comparison condition
- **Protein labels** (XLSX): Optional custom labels for highlighting specific proteins

## Methods
1. **Data Transformation**: Converts MSstats output to FragPipe-like format
2. **Differential Expression Classification**: Categorizes proteins as UP, DOWN, or NO based on thresholds
3. **Label Assignment**: Applies custom labels to proteins of interest
4. **Volcano Plot Generation**: Creates plots with consistent styling across all comparisons
5. **Batch Processing**: Loops through all comparison files in specified directory

## Key Parameters
- **Fold change threshold**: Typically 1.5 (log2 scale)
- **P-value threshold**: 0.05 (adjusted)
- **Color scheme**: Custom colors for UP (red), DOWN (blue), NO (grey)
- **Label filtering**: Only proteins meeting FC and p-value thresholds are labeled

## Outputs
- Volcano plots for each comparison condition
- Consistent styling across all plots
- Labeled proteins of interest highlighted
- Publication-ready graphics
