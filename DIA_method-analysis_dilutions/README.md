# DIA Method Analysis - Dilution Series

## Overview
This analysis evaluates DIA (Data-Independent Acquisition) quantification accuracy across dilution series. It compares expected vs. observed fold changes across different feature subsets and processing methods (narrow vs. wide isolation windows).

## Inputs
- **Multiple DIA result files** (CSV): Feature-level, protein-level, comparison results, and model QC data
- **Dilution series metadata**: Expected fold changes based on dilution ratios

## Methods
1. **File Parsing**: Extracts subset type (HQ, top3, top4, top5, top10, top50) and table type from filenames
2. **Expected FC Calculation**: Derives expected log2 fold changes from dilution ratios (e.g., 10X vs 1X)
3. **Method Comparison**: Compares narrow vs. wide isolation window performance
4. **Accuracy Assessment**: Evaluates observed vs. expected fold changes across conditions

## Key Parameters
- **Feature subsets**: HQ, top3, top4, top5, top10, top50
- **Isolation window methods**: Narrow vs. wide
- **Imputation flag**: Tracks whether data was imputed
- **Dilution ratios**: Extracted from comparison labels (e.g., "10X_vs_1X")

## Outputs
- Comparison of quantification accuracy across methods and feature subsets
- Fold change accuracy metrics
- Method performance evaluation
