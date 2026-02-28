# =============================================================================
# config.R
# Author: Andrew Becker
# Description: Centralized configuration file for all AB_R-proteomics-analyses
#              scripts. Defines shared paths, constants, and color palettes.
#              Source this file at the top of any .Rmd that uses Box paths.
#
# Usage:
#   source(here::here("src", "config.R"))
# =============================================================================


# =============================================================================
# 1. ROOT PATHS
# =============================================================================

# Box root directory (update this if moving to a new machine)
box_root <- "C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD"

# Reference databases directory (UniProt, half-life tables, PTM databases, etc.)
ref_dir <- file.path(box_root, "Reference_dbs")


# =============================================================================
# 2. REFERENCE DATABASE FILE PATHS
# =============================================================================

# UniProt human proteome mapping (ProteinID <-> Gene Name)
uniprot_hs_path <- file.path(ref_dir, "uniprot_hs_29may25.xlsx")

# UniProt human proteome with PTM databases (for RNAseq_expression_merge)
uniprot_hs_ptm_path <- file.path(ref_dir, "uniprot_hs_genes_PTM-DBs_17jan26.xlsx")

# Protein half-life data (Savitski 2018, high-quality subset)
half_life_path <- file.path(ref_dir, "protein-half-life-Savitski-2018_high-qual.xlsx")

# Ubiquitylation sites (Choudhary 2024 Cell paper)
ubiquitylation_path <- file.path(ref_dir, "Choudhary_ubiquitylation_TableS1.xlsx")


# =============================================================================
# 3. EXPERIMENT DIRECTORY PATHS
# =============================================================================

# AP2-296: A431 3-probe experiment
ap2_296_dir        <- file.path(box_root, "AP2-296_A431-3probe")
ap2_296_input      <- file.path(ap2_296_dir, "01_Input_Search_Results")
ap2_296_parquet    <- file.path(ap2_296_dir, "02_R_Intermediates")

# AP2-327: HAEC expression experiment
ap2_327_dir        <- file.path(box_root, "AP2-327_HAEC-expression")
ap2_327_input      <- file.path(ap2_327_dir, "01_Input_Search_Results")
ap2_327_parquet    <- file.path(ap2_327_dir, "02_R_Intermediates")

# AP2-339: HAECs DF-UF ALOD4/OlyA experiment
ap2_339_dir        <- file.path(box_root, "AP2-339_HAECs_DF-UF_ALOD4_OLyA")
ap2_339_input      <- file.path(ap2_339_dir, "01_Input_Search_Results")
ap2_339_parquet    <- file.path(ap2_339_dir, "02_R_Intermediates")

# AP2-341: DIA dilution series experiment
ap2_341_dir        <- file.path(box_root, "AP2-341")
ap2_341_input      <- ap2_341_dir
ap2_341_output     <- file.path(ap2_341_dir, "analysis_outputs")

# POCA2: Combined analysis outputs
poca2_dir          <- file.path(box_root, "POCA2")
poca2_intermediates <- file.path(poca2_dir, "02_R_Intermediates")
poca2_plots        <- file.path(poca2_dir, "03_graphs_plots")


# =============================================================================
# 4. ANALYSIS THRESHOLDS & PARAMETERS
# =============================================================================

# Volcano plot thresholds
FC_THRESH    <- 1.5    # log2 fold change threshold
PVAL_THRESH  <- 0.05   # adjusted p-value threshold

# RNAseq expression thresholds
LOD_THRESH   <- 0.1    # limit of detection (CPM)
REG_THRESH   <- 0.5    # minimum CPM to be considered "regulated"
DISC_SIGMA   <- 3      # studentized residual threshold for discordance

# Imputation
TOP_N_FEATURES <- 4    # number of top features per protein (MSstats)


# =============================================================================
# 5. COLOR PALETTES
# =============================================================================

# Standard condition colors (DF = red, UF = blue)
c_df    <- "#E41A1C"   # DF condition (red)
c_uf    <- "#377EB8"   # UF condition (blue)
c_grey  <- "grey80"    # neutral/concordant

# Probe colors for POCA experiments
c_alod4   <- "#cf6a29"   # ALOD4 (orange)
c_olya    <- "#036a38"   # OlyA (green)
c_ctx     <- "#9e8599"   # Ctx (mauve)
c_control <- "#939598"   # control (grey)

# Volcano plot colors (UP, DOWN, NO)
c_up   <- "#eb6e34"    # upregulated
c_down <- "#165e3a"    # downregulated
c_no   <- "light grey" # not significant

# Proximity shift colors
c_prox_up   <- "#bd4a7f"  # upregulated proximity
c_prox_down <- "#67a5c7"  # downregulated proximity


# =============================================================================
# 6. HELPER: ENSURE OUTPUT DIRECTORIES EXIST
# =============================================================================

# Call this function to create any output directory if it doesn't exist
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }
  invisible(path)
}


# =============================================================================
# 7. CONFIRMATION MESSAGE
# =============================================================================

message("config.R loaded successfully.")
message("  box_root: ", box_root)
message("  ref_dir:  ", ref_dir)
