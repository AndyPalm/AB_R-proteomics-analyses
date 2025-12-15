# src/plot_ancova_scatter.R
library(dplyr)
library(ggplot2)
library(ggrepel)

# --- 1. CONFIGURATION ---
# Significance Thresholds for the base plot
alpha <- 0.05
min_shift <- 1.0 

# --- 2. DATA PREPARATION ---

# Ensure base data exists
if (!exists("master_ancova_df") || !exists("proximity_scores")) {
  stop("Error: Analysis data missing. Please run 'src/protein_anova_analysis.R' first.")
}

# Join Stats to Raw Data (for Red/Blue coloring)
plot_data <- master_ancova_df %>%
  inner_join(proximity_scores, by = c("Protein", "Treatment")) %>%
  mutate(
    Is_Significant = adj.p.value < alpha & abs(Proximity_Shift) > min_shift,
    # Define distinct categories for coloring
    Hit_Type = case_when(
      Is_Significant & Proximity_Shift > 0 ~ "Increased",
      Is_Significant & Proximity_Shift < 0 ~ "Decreased",
      TRUE ~ "Background"
    )
  )

# --- 3. BASE PLOT CONSTRUCTION ---

p_scatter <- ggplot(plot_data, aes(x = Mean_Expr, y = Label_Intensity)) +
  
  # Layer 1: Background Points (Grey, small)
  geom_point(data = subset(plot_data, Hit_Type == "Background"), 
             color = "grey85", alpha = 0.3, size = 1) +
  
  # Layer 2: Significant Hits (Colored, slightly larger)
  geom_point(data = subset(plot_data, Hit_Type != "Background"), 
             aes(color = State), alpha = 0.6, size = 1.5) +
  
  # Layer 3: Regression Lines (The ANCOVA Visual)
  geom_smooth(method = "lm", aes(color = State), se = FALSE, linewidth = 0.5) +
  
  facet_wrap(~Treatment) +
  theme_bw() +
  labs(
    title = "ANCOVA Mechanism: Labeling vs. Expression",
    subtitle = "Regression lines show global trends. Colored points are significant outliers.",
    x = "Mean Expression (Log2)",
    y = "Labeling Replicates (Log2)",
    color = "State"
  )

# --- 4. CONDITIONAL HIGHLIGHTING ---

if (highlight_targets == TRUE && file.exists(target_file_path)) {
  
  message("Found target list! Applying highlights...")
  
  # A. Read and Filter Targets
  targets <- read.csv(target_file_path, stringsAsFactors = FALSE)
  target_ids <- targets$ProteinID # Assumes column header is 'ProteinID'
  
  # B. Prepare Highlight Data (Replicates)
  highlight_points <- plot_data %>%
    filter(Protein %in% target_ids) %>%
    left_join(df.labels, by = c("Protein" = "ProteinID")) %>%
    mutate(Gene = ifelse(is.na(Gene) | Gene == "", Protein, Gene))
  
  # C. Prepare Label Data (Centroids)
  # Calculate the mean Y position so the text label sits in the middle of the vertical stack
  highlight_labels <- highlight_points %>%
    group_by(Protein, Gene, Treatment, State, Mean_Expr) %>%
    summarise(Mean_Label_Y = mean(Label_Intensity), .groups = "drop")
  
  # D. Layer onto the Plot
  p_scatter <- p_scatter +
    # Draw Black Circles around the specific replicates
    geom_point(data = highlight_points, 
               shape = 1, color = "black", size = 2, stroke = 0.8) +
    
    # Add Text Labels pointing to the centroid
    geom_text_repel(
      data = highlight_labels,
      aes(y = Mean_Label_Y, label = Gene),
      size = 3.5,
      fontface = "bold",
      color = "black",
      min.segment.length = 0, # Force connecting lines
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    labs(caption = "Black circles indicate proteins from 'target_proteins.csv'")
  
} else if (highlight_targets == TRUE) {
  warning("Highlighting was requested, but 'data/target_proteins.csv' was not found.")
}

# --- 5. OUTPUT ---
print(p_scatter)
