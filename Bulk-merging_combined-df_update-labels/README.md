# Bulk Merging & Combined DataFrame with Label Updates

## Overview
This analysis merges multiple proteomics result files into a single combined dataframe and applies gene name annotations from a UniProt reference database. Handles both homo sapiens and non-human protein entries.

## Inputs
- **UniProt lookup table** (CSV): ProteinID to gene name mapping
- **Multiple result files** (CSV): Fold change, p-value, and differential expression data from different experiments

## Methods
1. **Merge Strategy 1 (MERGE 1)**: Direct merge on ProteinID for 100% homo sapiens datasets
2. **Merge Strategy 2 (MERGE 2)**: Fallback merge that preserves original labels if UniProt match not found (for mixed-species datasets)
3. **Column Standardization**: Renames columns to consistent format across all input files
4. **Annotation**: Adds gene names from lookup table

## Key Parameters
- **Merge key**: ProteinID (UniProt accession)
- **Fallback strategy**: Preserves original labels if lookup fails
- **Column naming**: Standardized format with experiment-specific suffixes

## Outputs
- Combined dataframe with all experiments merged on ProteinID and gene name
- Handles missing annotations gracefully
- CSV export of merged results
