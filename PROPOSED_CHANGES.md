# Proposed Changes for .Rmd File Standardization & Streamlining

## FORMATTING STANDARDIZATION

### 1. Add consistent YAML headers to all .Rmd files
1a. **Why**: Files have inconsistent or missing metadata (author, R version, description). Makes it hard to track which version of R was used and who created the analysis.
1b. **What**: Add standardized YAML header with title, author, date, R version, and brief description to all files.

### 2. Consolidate library loading at top of each file
1a. **Why**: Libraries are scattered throughout files (e.g., AB_TMT has libraries in first chunk, but basic-labeling_expression_correlation_plot loads tidyverse mid-file). Makes dependencies unclear.
1b. **What**: Move all `library()` calls to a single "SETUP" chunk at the beginning of each file, organized alphabetically.

### 3. Standardize section headers with consistent formatting
1a. **Why**: Section headers use inconsistent styles: `# --- SECTION ---`, `## Section`, `#Section`, etc. Reduces readability.
1b. **What**: Use consistent format: `# --- SECTION NAME ---` for major sections, `## Subsection` for subsections.

### 4. Add descriptive comments to major code blocks
1a. **Why**: Many chunks lack explanatory comments (e.g., DIA_feature_extract_and_plot has data loading without context).
1b. **What**: Add 1-2 line comments above each major code block explaining its purpose.

### 5. Standardize chunk naming conventions
1a. **Why**: Chunk names are inconsistent or missing (some use `{r}`, others use `{r setup}`, `{r plot_extraction}`). Makes navigation harder.
1b. **What**: Use consistent naming: `{r setup}`, `{r data_loading}`, `{r analysis}`, `{r visualization}`, `{r export}`.

---

## REDUNDANCIES & STREAMLINING OPPORTUNITIES

### 6. Remove duplicate library imports
1a. **Why**: Multiple files import the same libraries multiple times:
   - `stringr` appears 3+ times in AB_TMT, basic-labeling_expression_correlation_plot, Bubble_plot
   - `ggplot2`, `dplyr`, `tidyverse` duplicated in many files
1b. **What**: Consolidate all library calls into single setup chunk at top of file.

### 7. Extract common data loading pattern to helper function
1a. **Why**: CSV-to-Parquet conversion pattern is repeated identically in 5+ files (Bubble_plot, DIA_feature_extract_and_plot, Total_features_proteins_graph, etc.):
   ```r
   if (!file.exists(pq_target)) {
     temp_df <- arrow::read_csv_arrow(csv_source)
     arrow::write_parquet(temp_df, pq_target)
     rm(temp_df); gc()
   } else {
     message("Parquet file found...")
   }
   ```
1b. **What**: Create helper function `load_or_convert_parquet()` in `src/functions.R` and call it from all files.

### 8. Extract UniProt lookup pattern to helper function
1a. **Why**: Gene name mapping is repeated identically in 4+ files (Bubble_plot, DIA_feature_extract_and_plot, Total_features_proteins_graph, expression_labeling_merge):
   ```r
   df.labels <- read_excel(file.path(ref_dir, "uniprot_hs_29may25.xlsx"))
   lookup <- df.labels %>% select(ProteinID, Gene)
   data <- data %>% left_join(lookup, by = c("Protein" = "ProteinID"))
   ```
1b. **What**: Create helper function `add_gene_names()` in `src/functions.R` to standardize this operation.

### 9. Consolidate path setup into centralized configuration
1a. **Why**: Multiple files hardcode Box paths with slight variations:
   - `box_root <- "C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD"`
   - Some use `~/R_data/` instead
   - Makes it hard to update paths globally
1b. **What**: Create `config.R` in `src/` with all path definitions, source it at top of each file.

### 10. Remove commented-out code blocks
1a. **Why**: Multiple files have large blocks of commented code (e.g., basic-labeling_expression_correlation_plot has ~30 lines of commented plotting code, Venn_diagram has commented install.packages).
1b. **What**: Remove all commented-out code blocks (can be recovered from git history if needed).

### 11. Simplify redundant variance categorization logic
1a. **Why**: AB_TMT repeats variance categorization 3 times (once for Ctx, ALO, Olya) with identical logic:
   ```r
   Ctx_variance_category = case_when(
     Ctx_cv < 0.1 ~ "Low variance",
     Ctx_cv < 0.3 ~ "Medium variance",
     TRUE ~ "High variance"
   )
   ```
1b. **What**: Create helper function `categorize_variance(cv_vector, thresholds = c(0.1, 0.3))` to apply once.

### 12. Consolidate repeated ggplot theme settings
1a. **Why**: Multiple files repeat identical theme customizations:
   - `theme_minimal() + theme(axis.text.x = element_text(size = 16, color = "black", family = "Arial"), ...)`
   - Appears in basic-labeling_expression_correlation_plot, circle-plot_heatmap, and others
1b. **What**: Create `theme_custom()` function in `src/` with standard theme settings.

### 13. Simplify data cleaning patterns
1a. **Why**: Multiple files have similar data cleaning steps scattered throughout:
   - Removing rows with NA values
   - Converting columns to numeric
   - Renaming columns
1b. **What**: Create helper functions like `clean_numeric_columns()`, `remove_na_rows()` in `src/functions.R`.

### 14. Consolidate fold change calculation logic
1a. **Why**: Fold change calculations are repeated in multiple files with slight variations:
   - AB_TMT: `FC = 2^(mean1 - mean2)`
   - basic-labeling_expression_correlation_plot: Similar logic
   - AB_Foldchange_correlation: Similar logic
1b. **What**: Create `calculate_fold_changes()` helper function with consistent implementation.

### 15. Remove redundant data validation checks
1a. **Why**: Multiple files check for empty dataframes or missing columns with similar code:
   ```r
   if (nrow(data) == 0) { warning("No data...") }
   if (colnames(df)[1] == "" || is.na(colnames(df)[1])) { ... }
   ```
1b. **What**: Create `validate_data()` helper function that performs all standard checks.

### 16. Simplify repeated filtering patterns
1a. **Why**: Files repeat similar filtering logic:
   - `filter(!is.na(column))` appears in many files
   - `filter(GROUP != "control")` appears in multiple files
1b. **What**: Create helper functions like `remove_na_rows()`, `filter_by_group()` for common patterns.

### 17. Consolidate repeated color/aesthetic definitions
1a. **Why**: Color schemes are defined separately in multiple files:
   - circle-plot_heatmap: `group_colors <- c("ALOD4" = "#FF6B6B", ...)`
   - basic-labeling_expression_correlation_plot: Similar color definitions
1b. **What**: Create `color_schemes.R` in `src/` with standard color palettes.

### 18. Remove unused library imports
1a. **Why**: Some files import libraries that aren't used:
   - Venn_diagram imports `extrafont` and `Cairo` but doesn't use them
   - Several files import `pbapply` but don't use it
1b. **What**: Audit each file and remove unused imports.

### 19. Simplify repeated message/print statements
1a. **Why**: Multiple files have similar progress/status messages:
   - "Parquet file found. Loading cached data."
   - "Conversion complete."
1b. **What**: Create helper function `print_status()` for consistent messaging.

### 20. Consolidate repeated data aggregation patterns
1a. **Why**: Multiple files repeat similar aggregation logic:
   ```r
   group_by(Protein, GROUP) %>%
   summarise(mean_intensity = mean(...), count = n(), ...)
   ```
1b. **What**: Create `aggregate_by_group()` helper function with customizable metrics.

---

## SUMMARY

**Total Proposed Changes: 20**
- **Formatting Standardization: 5 changes** (headers, libraries, sections, comments, chunk names)
- **Code Streamlining: 15 changes** (remove duplicates, extract helpers, consolidate patterns)

**Safety Notes:**
- All changes preserve column names and dataframe references
- Helper functions will be non-breaking and optional
- Commented code will be removed (recoverable from git)
- No changes to analysis logic or statistical methods
