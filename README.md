Welcome fellow chemical biologist/mass spectrometrist/Backus lab member

This guide explains how to set up and run the R markdown files (.Rmd) for accomplishing different tasks. This requires having the R scripts (.R) accessible in the /src folder by loading the configuration file (`src/config.R`) so that all the analysis scripts in this repository work on your computer.

---

### QUICK START CHECKLIST ###
- [ ] **Open** `src/config.R` in RStudio
- [ ] **Find** the line with `box_root <- "C:/Users/andbp/..."`
- [ ] **Replace** it with YOUR data folder path
- [ ] **Save** the file
- [ ] **Open** any `.Rmd` file (e.g., `RNAseq_expression_merge/RNAseq_expression_merge.Rmd`)
- [ ] **Look for** the line `source(here::here("src", "config.R"))` in the setup chunk
- [ ] **Run** the script — it should now find your data!

---

## What is `src/config.R`?

Think of `src/config.R` as a **master settings file** that tells all the analysis scripts where to find your data files and what colors to use for plots. Instead of having the same file paths written in 20 different places, we have one central file that everyone uses.

**Why is this helpful?**
- If you move your data folder to a new location, you only need to update ONE file instead of 20
- All scripts use the same colors and settings, so your plots look consistent
- It's easier to share your work with others
- There is less bloat in the .Rmd files. Say goodbye to 1000-line markdown files where the entire text of the function is run within the .Rmd; simply call the .R file stored in /src to perform the function!

---

## Step 1: Update the Box Root Path

### What is the "Box root"?

The "Box root" is the main folder where all your data lives. In the current setup, it's:
```
C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD
```

### How to update it for your computer:

1. **Open the file** `src/config.R` in RStudio or any text editor
2. **Find this line** (it's near the top):
   ```r
   box_root <- "C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD"
   ```
3. **Replace it** with YOUR path. For example:
   ```r
   box_root <- "C:/Users/YourUsername/Box/Backus_Lab/YourName/YourFolder"
   ```
   OR if your data is on your Desktop:
   ```r
   box_root <- "C:/Users/YourUsername/Desktop/MyProteomicsData"
   ```

4. **Save the file** (Ctrl+S or Cmd+S)

**That's it!** Now all 20 analysis scripts will automatically look for data in YOUR folder.

---

## Step 2: Understand the 7 Sections of `src/config.R`

The configuration file is organized into 7 sections. Here's what each one does:

### Section 1: Root Paths
```r
box_root <- "..."  # Your main data folder
ref_dir <- "..."   # Where reference databases are stored
```
**What it does:** Defines the main folders where your data lives.

### Section 2: Reference Database File Paths
```r
uniprot_hs_path <- "..."           # Protein ID to gene name mapping
half_life_path <- "..."            # Protein stability data
ubiquitylation_path <- "..."       # Protein modification data
```
**What it does:** Points to special reference files that help identify proteins and understand their properties.

### Section 3: Experiment Directory Paths
```r
ap2_296_dir <- "..."    # Folder for experiment AP2-296
ap2_327_dir <- "..."    # Folder for experiment AP2-327
ap2_339_dir <- "..."    # Folder for experiment AP2-339
# ... and more
```
**What it does:** Defines where each experiment's data is stored. Each experiment has subfolders for:
- `01_Input_Search_Results` — Raw data from the mass spectrometer
- `02_R_Intermediates` — Processed data files (Parquet format for speed)

### Section 4: Analysis Thresholds & Parameters
```r
FC_THRESH <- 1.5       # How much a protein must change to be "significant"
PVAL_THRESH <- 0.05    # Statistical confidence level
LOD_THRESH <- 0.1      # Minimum detection level
TOP_N_FEATURES <- 4    # Number of measurements per protein
```
**What it does:** Sets the rules for what counts as a "real" change in your data. You can adjust these if you want stricter or looser criteria.

### Section 5: Color Palettes
```r
c_df <- "#E41A1C"      # Red color for DF condition
c_uf <- "#377EB8"      # Blue color for UF condition
c_alod4 <- "#cf6a29"   # Orange color for ALOD4 probe
# ... and more
```
**What it does:** Defines the colors used in all plots. This ensures all your figures look consistent.

### Section 6: Helper Function
```r
ensure_dir <- function(path) { ... }
```
**What it does:** A small helper that creates output folders if they don't exist yet.

### Section 7: Confirmation Message
When you load the config file, it prints a message confirming it worked.

---

## Step 3: How Scripts Use the Configuration File

### In the setup chunk of each .Rmd file:

Every analysis script now starts with this code:
```r
```{r setup, include=FALSE}
# --- SETUP: Load required libraries ---
library(dplyr)
library(ggplot2)
# ... other libraries ...

# --- LOAD CONFIGURATION ---
source(here::here("src", "config.R"))
```
```

**What this does:**
1. Loads all the R packages (libraries) needed for the analysis
2. Loads the configuration file, which makes all the paths and colors available

### Using the configuration in your script:

Instead of writing the full path like this:
```r
data <- read.csv("C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD/AP2-339_HAECs_DF-UF_ALOD4_OLyA/01_Input_Search_Results/339B-top4-features_comparison-results_imputed.csv")
```

You can now write:
```r
data <- read.csv(file.path(ap2_339_input, "339B-top4-features_comparison-results_imputed.csv"))
```

**Why is this better?**
- Much shorter and easier to read
- If you move your data folder, you only change ONE line in `config.R`
- All scripts automatically use the new path

---

## Step 4: Replacing Hardcoded Paths in Your Scripts

If you're updating an old script that has hardcoded paths, here's how to replace them:

### Before (hardcoded path):
```r
df.labels <- read_excel("C:/Users/andbp/Box/Backus_Lab/Andrew_Becker/ABecker_Lab_Notebook/R_WD/Reference_dbs/uniprot_hs_29may25.xlsx")
```

### After (using config):
```r
df.labels <- read_excel(uniprot_hs_path)
```

### Common replacements:

| What you want | Use this variable |
|---|---|
| UniProt mapping file | `uniprot_hs_path` |
| Protein half-life data | `half_life_path` |
| AP2-296 input folder | `ap2_296_input` |
| AP2-327 input folder | `ap2_327_input` |
| AP2-339 input folder | `ap2_339_input` |
| POCA2 output folder | `poca2_plots` |
| Reference database folder | `ref_dir` |

---

## Step 5: Quick Start Checklist

- [ ] **Open** `src/config.R` in RStudio
- [ ] **Find** the line with `box_root <- "C:/Users/andbp/..."`
- [ ] **Replace** it with YOUR data folder path
- [ ] **Save** the file
- [ ] **Open** any `.Rmd` file (e.g., `RNAseq_expression_merge/RNAseq_expression_merge.Rmd`)
- [ ] **Look for** the line `source(here::here("src", "config.R"))` in the setup chunk
- [ ] **Run** the script — it should now find your data!

---

## Troubleshooting

### Problem: "File not found" error

**Solution:** Check that your `box_root` path is correct:
1. Open File Explorer (Windows) or Finder (Mac)
2. Navigate to your data folder
3. Copy the full path from the address bar
4. Paste it into `config.R` (replace backslashes `\` with forward slashes `/`)

### Problem: "config.R not found"

**Solution:** Make sure you're running the script from the correct folder:
1. In RStudio, go to **Session → Set Working Directory → To Source File Location**
2. This tells R to look for files relative to where your `.Rmd` file is

### Problem: Colors look different in my plots

**Solution:** You can customize the colors in Section 5 of `config.R`. For example:
```r
c_df <- "#FF0000"  # Change to bright red
c_uf <- "#0000FF"  # Change to bright blue
```

---

## For Advanced Users

If you want to add new experiments or change thresholds:

1. **Add a new experiment path:**
   ```r
   ap2_400_dir <- file.path(box_root, "AP2-400_MyNewExperiment")
   ap2_400_input <- file.path(ap2_400_dir, "01_Input_Search_Results")
   ```

2. **Change analysis thresholds:**
   ```r
   FC_THRESH <- 2.0    # Stricter: only count 2-fold changes
   PVAL_THRESH <- 0.01 # Stricter: higher confidence
   ```

3. **Add new colors:**
   ```r
   c_mycolor <- "#FF00FF"  # Magenta
   ```

---

## Questions?

If something doesn't work:
1. Check that your `box_root` path is correct
2. Make sure the file exists at that path
3. Check that you saved `config.R` after making changes
4. Try running the script again
5. Load this README.md into gemini.google.com, claude.ai, or chatgpt.com along with a screenshot of your directory tree and figure it out. Don't email me.

Good luck!
