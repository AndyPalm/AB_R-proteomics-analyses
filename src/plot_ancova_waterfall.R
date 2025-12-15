# src/plot_ancova_waterfall.R

# Load required package for non-overlapping labels
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package 'ggrepel' is needed for this plot. Please run: install.packages('ggrepel')")
}
library(ggrepel)

# Safety Check
if (!exists("df.labels")) {
  stop("Error: 'df.labels' dataframe not found. Please load your protein-to-gene mapping file first.")
}

# 1. Prepare Data
proximity_with_genes <- proximity_scores %>%
  left_join(df.labels, by = c("Protein" = "ProteinID")) %>%
  mutate(Gene = ifelse(is.na(Gene) | Gene == "", Protein, Gene))

waterfall_data <- proximity_with_genes %>%
  group_by(Treatment) %>%
  arrange(Proximity_Shift) %>%
  mutate(
    Rank = row_number(),
    # Significance logic
    Is_Significant = adj.p.value < 0.05 & abs(Proximity_Shift) > 1.0,
    
    Point_Color = case_when(
      Is_Significant & Proximity_Shift > 0 ~ "Increased Capture",
      Is_Significant & Proximity_Shift < 0 ~ "Decreased Capture",
      TRUE ~ "NS"
    ),
    
    Point_Size = ifelse(Is_Significant, 2.5, 1.2)
  )

# 2. Filter Data for Labels (Top Hits Up and Down)
label_data <- waterfall_data %>%
  filter(Is_Significant == TRUE) %>%
  group_by(Treatment, Direction_Up = Proximity_Shift > 0) %>%
  slice_max(order_by = abs(Proximity_Shift), n = 15) %>%
  ungroup()

# 3. Build Plot
p_water <- ggplot(waterfall_data, aes(x = Rank, y = Proximity_Shift)) +
  
  geom_hline(yintercept = 0, linetype = "solid", color = "grey60", size = 0.3) +
  
  # The Points
  geom_point(aes(color = Point_Color, size = Point_Size), alpha = 0.8) +
  
  scale_color_manual(values = c("Increased Capture" = "#E41A1C", 
                                "Decreased Capture" = "#377EB8", 
                                "NS" = "grey85")) +
  scale_size_identity() +
  
  facet_wrap(~Treatment, scales = "free_x") +
  
  labs(
    title = "Proximity Shift Waterfall (Gene Level)",
    subtitle = "Labeled proteins are the top significant hits (FDR < 0.05, >2-fold shift)",
    x = "Protein Rank",
    y = "Proximity Shift (Log2)",
    color = "Status"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  ) +
  
  # 4. Add Repulsive Labels (The "Purple Arrow" Request)
  geom_text_repel(
    data = label_data, 
    aes(label = Gene),
    
    # Force horizontal text
    angle = 0, 
    size = 3.5,
    
    # Pushing logic:
    # If shift is positive (Red), push label to the LEFT (nudge_x negative)
    # If shift is negative (Blue), push label to the RIGHT (nudge_x positive)
    nudge_x = ifelse(label_data$Proximity_Shift > 0, -100, 100),
    
    # Ensure lines are drawn
    min.segment.length = 0, 
    segment.color = "black",
    segment.size = 0.3,
    
    # Spread them out vertically so they don't overlap
    direction = "y",
    
    show.legend = FALSE
  )

# Print and Save
print(p_water)