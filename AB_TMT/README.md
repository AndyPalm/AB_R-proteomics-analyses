# AB_TMT Analysis

## Overview
This analysis processes TMT (Tandem Mass Tag) quantitative proteomics data from multiple treatment conditions. The script calculates group-level statistics, fold changes, and variance metrics across biological replicates.

## Inputs
- **TMT abundance data** (TSV format): Protein-level intensity measurements from TMT-labeled samples
- **Sample groups**: Multiple treatment conditions (e.g., control, Ctx, ALO, OlyA) with 3 replicates each

## Methods
1. **Data Aggregation**: Calculates mean intensities and standard deviations for each protein within each treatment group
2. **Variance Assessment**: Computes coefficient of variation (CV) for relative variance estimation
3. **Fold Change Calculation**: Computes pairwise log2 fold changes between conditions (both linear and log2 space)
4. **Variance Categorization**: Classifies proteins as low (<0.1), medium (0.1-0.3), or high (>0.3) variance based on CV

## Key Parameters
- **CV thresholds**: Low (<0.1), Medium (0.1-0.3), High (>0.3)
- **Variance flag**: Proteins with CV > 0.3 in any condition are flagged
- **Fold change calculation**: 2^(mean_condition1 - mean_condition2)

## Outputs
- Summary table with group means, standard deviations, CVs, and pairwise fold changes
- Variance assessment with categorization flags
- CSV export of all calculated metrics
