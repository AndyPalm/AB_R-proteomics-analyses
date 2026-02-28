# Circle Plot Heatmap

## Overview
This analysis creates circular/arc-based visualizations of protein abundance across multiple treatment conditions. Each protein is represented as a circle with arc segments showing intensity levels for different treatment groups.

## Inputs
- **Protein-level data** (CSV): Log-transformed intensity measurements across multiple conditions and replicates
- **Protein list** (vector): UniProt accession numbers of proteins to visualize

## Methods
1. **Data Aggregation**: Calculates mean intensity for each protein-condition combination
2. **Normalization**: Scales intensities to 0-1 range for consistent visualization
3. **Arc Visualization**: Creates arc bars with radius proportional to intensity
4. **Color Mapping**: Uses intensity-based color gradients to show abundance levels

## Key Parameters
- **Intensity scaling**: Power scaling (exponent 1.3) for enhanced visibility
- **Arc radius range**: 0.08 to 0.45 units based on normalized intensity
- **Color gradient**: Continuous scale from low to high intensity

## Outputs
- Circular plots showing protein abundance patterns across conditions
- Arc-based visualization with intensity-dependent sizing and coloring
