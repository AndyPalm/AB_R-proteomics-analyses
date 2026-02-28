# MSstats Analysis Pipeline

## Overview
This analysis implements the MSstats statistical framework for quantitative proteomics data. It converts FragPipe output to MSstats format, performs data normalization and summarization, and generates statistical comparisons between conditions.

## Inputs
- **FragPipe output** (CSV): Raw peptide/protein-level data with intensity measurements
- **Sample metadata**: Condition and biological replicate assignments

## Methods
1. **Format Conversion**: Converts FragPipe output to MSstats-compatible format
2. **BioReplicate Assignment**: Maps runs to biological replicates by condition
3. **Data Processing**: 
   - Log2 transformation
   - Normalization (none, median, or quantile)
   - Feature selection (top N features per protein)
   - Summarization (TMP - Tukey's Median Polish)
4. **Statistical Comparison**: Performs t-tests between conditions with multiple testing correction

## Key Parameters
- **Log transformation**: log2
- **Normalization method**: None (or median/quantile as specified)
- **Feature subset**: Top 4 features per protein
- **Summarization method**: TMP (Tukey's Median Polish)
- **Equal feature variance**: TRUE
- **Censored intensity handling**: NA

## Outputs
- Processed protein-level data with normalized intensities
- Comparison results with fold changes and p-values
- Model QC metrics
- Statistical summaries by condition
