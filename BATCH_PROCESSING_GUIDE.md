# Batch Processing Guide for .Rmd Formatting & Streamlining

## Overview
This guide provides a complete roadmap for applying formatting standardization and code streamlining to all remaining .Rmd files in the repository. Use this with a larger context window model for efficient batch processing.

## Status Summary

### ✅ COMPLETED (20/20 files)
1. **AB_TMT/AB_TMT.Rmd** - Formatting standardization applied
2. **basic-labeling_expression_correlation_plot/basic-labeling_expression_correlation_plot.Rmd** - Formatting standardization applied
3. **AB_3x3_volcanos_and-loop/AB_3x3_volcanos_and-loop.Rmd** - Formatting standardization applied
4. **AB_DIA_custom-analysis/AB_DIA_custom-analysis.Rmd** - Formatting standardization applied
5. **AB_Foldchange_correlation/AB_Foldchange_correlation.Rmd** - Formatting standardization applied
6. **Bubble_plot/Bubble_plot.Rmd** - Formatting standardization applied
7. **Bulk-merging_combined-df_update-labels/Bulk-merging_combined-df_update-labels.Rmd** - Formatting standardization applied
8. **CARC_CRAC_analysis/CARC_CRAC_analysis.Rmd** - Formatting standardization applied
9. **circle-plot_heatmap/circle-plot_heatmap.Rmd** - Formatting standardization applied
10. **DIA_feature_extract_and_plot/DIA_feature_extract_and_plot.Rmd** - Formatting standardization applied
11. **DIA_method-analysis_dilutions/DIA_method-analysis_dilutions.Rmd** - Formatting standardization applied
12. **expression_labeling_linear-model/expression_labeling_linear-model.Rmd** - Formatting standardization applied
13. **expression_labeling_merge/expression_labeling_merge.Rmd** - Formatting standardization applied
14. **mouse-human_mapping/mouse-human_mapping.Rmd** - Formatting standardization applied
15. **MSstats_AB/MSstats_AB.Rmd** - Formatting standardization applied
16. **MSstats_data-transform/MSstats_data-transform.Rmd** - Formatting standardization applied
17. **RNAseq_expression_merge/RNAseq_expression_merge.Rmd** - Formatting standardization applied
18. **TMT-timecourse-analysis/TMT-timecourse-analysis.Rmd** - Formatting standardization applied
19. **Total_features_proteins_graph/Total_features_proteins_graphs.Rmd** - Formatting standardization applied
20. **Venn_diagram/Venn_diagram.Rmd** - Formatting standardization applied

### ⏳ REMAINING (0/20 files)
All files have been processed!

---

## FORMATTING STANDARDIZATION TEMPLATE

Apply these changes to each .Rmd file in order:

### 1. Update YAML Header
**Pattern to find:**
```yaml
---
title: "..."
output: html_document
date: "..."
---
```

**Replace with:**
```yaml
---
title: "..."
author: "Andrew Becker"
output: html_document
date: "..."
R_version: "4.4.2"
description: "[Brief description of analysis from README.md]"
---
```

### 2. Consolidate Library Loading
**Pattern to find:**
- Multiple `library()` calls scattered throughout first chunk(s)
- Duplicate imports (e.g., `stringr` appears 2+ times)

**Replace with:**
```r
```{r setup, include=FALSE}
# --- SETUP: Load required libraries ---
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pbapply)
library(RColorBrewer)
library(readxl)
library(reshape2)
library(stringr)
library(tidyverse)
# [Add other libraries alphabetically]
```
```

### 3. Add Descriptive Chunk Names
**Pattern to find:**
```r
```{r}
# Some code
```
```

**Replace with:**
```r
```{r chunk_name}
# --- SECTION NAME ---
# Some code
```
```

**Standard chunk names to use:**
- `{r setup}` - Library loading and path setup
- `{r data_loading}` - Reading input files
- `{r data_cleaning}` - Data cleaning and preprocessing
- `{r data_processing}` - Main analysis calculations
- `{r analysis}` - Statistical analysis
- `{r visualization}` - Plotting and visualization
- `{r export}` - Writing output files

### 4. Add Section Headers
**Pattern to find:**
```r
# Some comment
code here
```

**Replace with:**
```r
# --- SECTION NAME ---
# Detailed comment explaining what this section does
code here
```

### 5. Remove Commented-Out Code
**Pattern to find:**
- Large blocks of `#` commented code
- Unused `#data <- read.csv(...)` examples

**Action:** Delete entirely (recoverable from git history)

---

## CODE STREAMLINING PATTERNS

### Pattern 1: Duplicate Library Imports
**Identify:** Same library imported multiple times in different chunks
**Action:** Move all to single `{r setup}` chunk, alphabetize

### Pattern 2: CSV-to-Parquet Conversion
**Identify:** This pattern repeated in 5+ files:
```r
if (!file.exists(pq_target)) {
  temp_df <- arrow::read_csv_arrow(csv_source)
  arrow::write_parquet(temp_df, pq_target)
  rm(temp_df); gc()
} else {
  message("Parquet file found...")
}
```
**Action:** Note for future extraction to `src/functions.R` as `load_or_convert_parquet()`

### Pattern 3: UniProt Lookup
**Identify:** This pattern repeated in 4+ files:
```r
df.labels <- read_excel(file.path(ref_dir, "uniprot_hs_29may25.xlsx"))
lookup <- df.labels %>% select(ProteinID, Gene)
data <- data %>% left_join(lookup, by = c("Protein" = "ProteinID"))
```
**Action:** Note for future extraction to `src/functions.R` as `add_gene_names()`

### Pattern 4: Repeated Theme Settings
**Identify:** Identical ggplot theme code in multiple files
**Action:** Note for future extraction to `src/functions.R` as `theme_custom()`

### Pattern 5: Hardcoded Box Paths
**Identify:** Paths like `"C:/Users/andbp/Box/Backus_Lab/..."`
**Action:** Note for future centralization in `src/config.R`

---

## README.md DEPENDENCY UPDATES

### Files that source .R scripts (need README updates):

1. **expression_labeling_merge/README.md** ✅ DONE
   - src/functions.R
   - src/plot_labeling_vs_expression_colors.R
   - src/protein_anova_analysis.R
   - src/plot_ancova_scatter.R
   - src/plot_ancova_waterfall.R

2. **DIA_feature_extract_and_plot/README.md** ✅ DONE
   - src/DIA_feature_extract_barplots.R
   - src/DIA_feature_extract_barplots_by_group.R

3. **Total_features_proteins_graph/README.md** ✅ DONE
   - (No .R scripts sourced)

4. **Other files** ✅ DONE
   - (No additional .R scripts sourced in remaining files)

---

## BATCH PROCESSING CHECKLIST

All 20 files have been processed with the following changes applied:

- [x] Read the .Rmd file
- [x] Update YAML header with author, R version, description
- [x] Consolidate all library() calls into single {r setup} chunk
- [x] Alphabetize library imports
- [x] Add descriptive chunk names
- [x] Add section headers with consistent formatting
- [x] Remove commented-out code blocks
- [x] Identify and note code patterns for future streamlining
- [x] Check if file sources any .R scripts
- [x] If sources .R scripts, update corresponding README.md with Dependencies section

---

## NOTES FOR NEXT SESSION

1. **Future Work - Helper Functions**: The following patterns appear in multiple files and are candidates for extraction to `src/functions.R`:
   - `load_or_convert_parquet()` - CSV-to-Parquet conversion (Bubble_plot, DIA_feature_extract_and_plot, Total_features_proteins_graph, RNAseq_expression_merge)
   - `add_gene_names()` - UniProt lookup join (Bubble_plot, DIA_feature_extract_and_plot, Total_features_proteins_graph, RNAseq_expression_merge)
   - `theme_custom()` - Repeated ggplot theme settings (Bubble_plot, circle-plot_heatmap)

2. **Future Work - Config File**: Hardcoded Box paths appear in many files. Consider creating `src/config.R` with:
   - `box_root`
   - `ref_dir`
   - Standard experiment directory structure

3. **Git History**: All changes are reversible via git if needed

---

## REFERENCE: PROPOSED_CHANGES.md

See PROPOSED_CHANGES.md for detailed explanations of:
- Why each change is needed
- What the change accomplishes
- Safety considerations

The 20 proposed changes are organized as:
- **5 Formatting Standardization changes**
- **15 Code Streamlining opportunities**
