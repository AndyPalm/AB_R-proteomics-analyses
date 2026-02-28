# src/plot_labeling_vs_expression.R
#
library(dplyr)
library(ggplot2)
library(ggrepel)
#
# --- 1. DATA PREPARATION ---
# Create the unique ColorKey by combining Condition and Status
df_plot <- df_combined_fc %>%
  #
  mutate(
    #
    PointStatus = case_when(
      #
      log2FC > 1 ~ "Up",
      #
      log2FC < -1 ~ "Down",
      #
      TRUE ~ "Stable"
      #
    ),
    #
    # Ensure this exactly matches the format: "Condition_Status"
    #
    ColorKey = paste(Condition, PointStatus, sep = "_")
    #
  )
#
# --- 2. DEFINE EXPLICIT COLOR MAP ---
# Use the exact strings from your df1_clean mutate logic
color_map <- c(
  #
  # Comparison 1: ALOD4_DF-vs-UF
  #
  "ALOD4_DF-vs-UF_Up"           = "#aa571e", # Orange
  #
  "ALOD4_DF-vs-UF_Down"         = "#f08821", # Blue
  #
  "ALOD4_DF-vs-UF_Stable"       = "grey90",
  #
  # Comparison 2: OlyA_DF-vs-UF
  #
  "OlyA_DF-vs-UF_Up"            = "#09552b", # Green
  #
  "OlyA_DF-vs-UF_Down"          = "#4cb962", # Pink
  #
  "OlyA_DF-vs-UF_Stable"        = "grey90",
  #
  # Comparison 3: OlyA-DF vs ALOD4-UF
  #
  "OlyA-DF vs ALOD4-UF_Up"      = "#09552b", # Lime
  #
  "OlyA-DF vs ALOD4-UF_Down"    = "#f08821", # Yellow
  #
  "OlyA-DF vs ALOD4-UF_Stable"  = "grey90",
  #
  # Comparison 4: ALOD4-DF vs OlyA-UF
  #
  "ALOD4-DF vs OlyA-UF_Up"      = "#aa571e", # Tan
  #
  "ALOD4-DF vs OlyA-UF_Down"    = "#4cb962", # Dark Grey
  #
  "ALOD4-DF vs OlyA-UF_Stable"  = "grey90"
  #
)
#
# --- 3. PLOTTING ---
p_fc <- ggplot(df_plot, aes(x = expr_log2FC, y = log2FC)) +
  #
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  #
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  #
  geom_hline(yintercept = -1, linetype = "dashed", color = "black") +
  #
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  #
  # Map color to the specific ColorKey
  #
  geom_point(aes(color = ColorKey), alpha = 0.75, size = 1.4) +
  #
  facet_wrap(~Condition) +
  #
  scale_color_manual(values = color_map) +
  #
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"   
  ) +
  #
  labs(
    #
    title = "Labeling vs. Expression Changes (DF vs UF)",
    #
    x = "Expression Log2(Fold Change)",
    #
    y = "Labeling Log2(Fold Change)"
    #
  )
#
# --- 4. CONDITIONAL HIGHLIGHTING (Keep your existing logic) ---
if (highlight_targets == TRUE && file.exists(target_file_path)) {
  #
  targets <- read.csv(target_file_path, stringsAsFactors = FALSE)
  #
  target_ids <- targets$ProteinID
  #
  highlight_data <- df_plot %>%
    #
    filter(Protein %in% target_ids) %>%
    #
    left_join(df.labels, by = c("Protein" = "ProteinID")) %>%
    #
    mutate(Gene = ifelse(is.na(Gene) | Gene == "", Protein, Gene))
  #
  p_fc <- p_fc + 
    #
    geom_point(data = highlight_data, 
               #
               shape = 1, color = "black", size = 1.4, stroke = 0.5) +
    #
    geom_text_repel(
      #
      data = highlight_data,
      #
      aes(label = Gene),
      #
      size = 4.5, fontface = "plain", color = "black",
      #
      box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0
      #
    )
  #
}
#
print(p_fc)