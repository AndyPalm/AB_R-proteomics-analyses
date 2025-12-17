# src/plot_labeling_vs_expression.R

# --- 0. SAFETY CHECK ---
# If variables are missing (e.g. ran chunk out of order), default to FALSE.
if (!exists("highlight_targets")) { highlight_targets <- FALSE }
if (!exists("target_file_path")) { target_file_path <- "" }

library(dplyr)
library(ggplot2)
library(ggrepel)

# --- 1. BASE PLOT CONSTRUCTION ---

# We assume 'df_combined_fc' was created in the Rmd
if (!exists("df_combined_fc")) {
  stop("Error: 'df_combined_fc' is missing. Please run the Merge Data chunk first.")
}

p_fc <- ggplot(df_combined_fc, aes(x = expr_log2FC, y = log2FC)) +
  
  # Reference Lines (Quadrants and Diagonal)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
 # geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  
  # Main Points
  geom_point(aes(color = Condition), alpha = 0.6, size = 2) +
  
  # Faceting
  facet_wrap(~Condition) +
  
  # Styling
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  labs(
    title = "Labeling vs. Expression Changes (DF vs UF)",
    subtitle = "",
    x = "Expression Log2(Fold Change)",
    y = "Labeling Log2(Fold Change)"
  )

# --- 2. CONDITIONAL HIGHLIGHTING ---

if (highlight_targets == TRUE && file.exists(target_file_path)) {
  
  message("Found target list! Applying highlights to Fold-Change plot...")
  
  # A. Read Targets
  targets <- read.csv(target_file_path, stringsAsFactors = FALSE)
  target_ids <- targets$ProteinID
  
  # B. Prepare Highlight Data
  highlight_data <- df_combined_fc %>%
    filter(Protein %in% target_ids) %>%
    # Add Gene Names
    left_join(df.labels, by = c("Protein" = "ProteinID")) %>%
    mutate(Gene = ifelse(is.na(Gene) | Gene == "", Protein, Gene))
  
  # C. Add Layers
  p_fc <- p_fc + 
    # 1. Black Circles (The "Target" indicator)
    geom_point(data = highlight_data, 
               shape = 1, color = "black", size = 2.5, stroke = 1.0) +
    
    # 2. Text Labels
    geom_text_repel(
      data = highlight_data,
      aes(label = Gene),
      size = 3.5,
      fontface = "bold",
      color = "black",
      box.padding = 0.5,
      max.overlaps = Inf,
      min.segment.length = 0
    )
}

# --- 3. OUTPUT ---
print(p_fc)