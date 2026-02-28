# Bubble Plot Analysis

## Overview
This analysis creates bubble plots to visualize protein abundance across multiple treatment groups. Each bubble represents a protein's mean intensity in a specific condition, with bubble size proportional to abundance and color indicating detection rate.

## Inputs
- **Protein-level data** (CSV): Log-transformed intensity measurements across replicates and groups
- **UniProt mapping file** (XLSX): Protein ID to gene name conversion table

## Methods
1. **Data Aggregation**: Calculates mean intensity and detection count for each protein-group combination
2. **Normalization**: Computes percentage detected (DetectCount / TotalReps)
3. **Visualization**: Creates bubble plots with intensity-based sizing and color gradients

## Key Parameters
- **Detection threshold**: Proteins with measurements in replicates
- **Intensity scaling**: Mean log-transformed intensities
- **Detection percentage**: Proportion of replicates with valid measurements

## Outputs
- Bubble plots showing protein abundance distribution across conditions
- Aggregated statistics table with mean intensities and detection rates
