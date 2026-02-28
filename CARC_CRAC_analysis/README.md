# CARC/CRAC Motif Analysis

## Overview
This analysis identifies and characterizes CARC (Cholesterol Consensus Motif) and CRAC (Cholesterol Recognition Amino acid Consensus) motifs in the human proteome. It analyzes motif positions, structural overlaps with transmembrane helices and beta-strands, and generates comprehensive statistics.

## Inputs
- **Proteome data** (XLSX): UniProt protein sequences with structural feature annotations (transmembrane helices, beta-strands)

## Methods
1. **Motif Pattern Matching**: 
   - CRAC: [LV](X0-5)Y(X0-5)[RK]
   - CARC: [RK](X0-5)[YF](X0-5)[LV]
2. **Structural Feature Parsing**: Extracts transmembrane helix and beta-strand positions from UniProt annotations
3. **Overlap Analysis**: Determines if motifs overlap with structural features
4. **Motif Characterization**: Counts motifs per protein and calculates structural overlap statistics

## Key Parameters
- **Motif patterns**: CRAC and CARC consensus sequences with flexible spacing (0-5 amino acids)
- **Structural features**: Transmembrane helices and beta-strands from UniProt
- **Overlap detection**: Boolean flag for motif-structure intersection

## Outputs
- Proteome-wide motif statistics (count, position, structural overlap)
- Protein-level annotations with motif presence/absence
- Summary statistics on motif prevalence and structural associations
