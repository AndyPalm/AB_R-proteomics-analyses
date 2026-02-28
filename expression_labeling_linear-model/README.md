# Expression-Labeling Linear Model Analysis

## Overview
This analysis performs per-protein linear regression to identify proteins with altered proximity to a specific cellular compartment or protein complex. It models labeling intensity while controlling for expression abundance to detect post-translational changes in protein localization.

## Inputs
- **Labeling data** (CSV): Proximity labeling intensities across conditions and replicates
- **Expression data** (CSV): Protein abundance measurements across same conditions and replicates

## Methods
1. **Data Reshaping**: Converts wide-format data to long format for modeling
2. **Per-Protein Linear Modeling**: Fits model: Labeling_Intensity ~ Condition + Expression_Intensity
3. **Coefficient Extraction**: Extracts condition effect (Proximity Shift Score) from model
4. **Multiple Testing Correction**: Applies Benjamini-Hochberg correction to p-values

## Key Parameters
- **Model formula**: Labeling_Intensity ~ Condition + Expression_Intensity
- **Reference level**: Parental (control condition)
- **Multiple testing correction**: Benjamini-Hochberg (BH) method
- **Effect term**: ConditionSRB1_oe (or equivalent treatment condition)

## Outputs
- Proximity Shift Scores for each protein (log2 scale, corrected for abundance)
- Adjusted p-values for statistical significance
- Ranked list of proteins with altered proximity
