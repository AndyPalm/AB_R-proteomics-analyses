# TMT Timecourse Analysis

## Overview
This analysis processes TMT quantitative proteomics data from a timecourse experiment. It calculates kinetic parameters including rates of change, maximum abundance, and time-to-half-maximum (T50) for each protein across multiple timepoints.

## Inputs
- **TMT abundance data** (TSV): Protein-level intensity measurements at multiple timepoints
- **Sample replicates**: Multiple replicates per timepoint (e.g., 0, 30, 60, 120, 360 minutes)

## Methods
1. **Data Aggregation**: Averages replicate measurements for each timepoint
2. **Rate Calculation**: Computes rates of change between consecutive timepoints
3. **Maximum Abundance**: Identifies peak intensity and corresponding timepoint
4. **Fold Change**: Calculates log2 fold change relative to baseline (time 0)
5. **T50 Estimation**: Interpolates time to reach half-maximum abundance using spline approximation

## Key Parameters
- **Timepoints**: 0, 30, 60, 120, 360 minutes (customizable)
- **Rate calculation**: (Intensity_t2 - Intensity_t1) / (time_t2 - time_t1)
- **T50 method**: Spline interpolation with extrapolation rule 2
- **Missing data handling**: Treats all-NA timepoints as 0

## Outputs
- Kinetic parameters table (rates, max abundance, T50, fold changes)
- Timecourse profiles for each protein
- CSV export of calculated metrics
