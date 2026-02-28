# Venn Diagram Visualization

## Overview
This analysis creates Venn diagrams to visualize the overlap of proteins across multiple experimental conditions or datasets. Supports 2-way, 3-way, and 4-way comparisons with customizable colors and labels.

## Inputs
- **Overlap statistics**: Counts of proteins in each region of the Venn diagram
  - Individual set sizes (area1, area2, area3, area4)
  - Pairwise intersections (n12, n13, n23, n14, n24, n34)
  - Triple intersections (n123, n124, n134, n234)
  - Quadruple intersection (n1234)

## Methods
1. **Venn Diagram Generation**: Uses VennDiagram package to create publication-quality diagrams
2. **Customization**: Applies custom colors, fonts, and sizing
3. **Export**: Saves diagrams as high-resolution PNG files

## Key Parameters
- **Diagram types**: 2-way (pairwise), 3-way (triple), 4-way (quadruple)
- **Color scheme**: Custom hex colors for each set
- **Resolution**: 1400 DPI for publication quality
- **Font**: Arial, regular style
- **Text sizing**: Customizable for intersection labels and category names

## Outputs
- High-resolution PNG files of Venn diagrams
- Visualization of protein overlap across conditions
- Publication-ready graphics
