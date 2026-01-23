library(ggplot2)
library(dplyr)
library(ggrepel)

generate_inf_strip <- function(inf_FC_df, comparison_label = "Inf Data", subcell_list_df = NULL, protein_labels_df = NULL) {
  
  # 1. Process Data
  plot_data <- inf_FC_df %>%
    group_by(Protein) %>%
    summarise(mean_intensity = mean(LogIntensities, na.rm = TRUE)) %>%
    arrange(mean_intensity)
  
  plot_data$Protein <- trimws(plot_data$Protein)
  
  # 2. Handle Subcellular List (Color Logic)
  # Default everyone to "Standard"
  plot_data$color_group <- "Standard"
  
  if (!is.null(subcell_list_df)) {
    # Clean whitespace for matching
    subcell_list_df$ProteinID <- trimws(subcell_list_df$ProteinID)
    
    # Identify matches and change their group
    matches_subcell <- plot_data$Protein %in% subcell_list_df$ProteinID
    plot_data$color_group[matches_subcell] <- "Subcell"
  }
  
  # 3. Handle Labeling
  plot_data$label_text <- NA
  if (!is.null(protein_labels_df)) {
    protein_labels_df$ProteinID <- trimws(protein_labels_df$ProteinID)
    matches_labels <- plot_data$Protein %in% protein_labels_df$ProteinID
    label_map <- setNames(protein_labels_df$Gene, protein_labels_df$ProteinID)
    plot_data$label_text[matches_labels] <- label_map[plot_data$Protein[matches_labels]]
  }
  
  # 4. Plotting
  p <- ggplot(data = plot_data, aes(x = 0, y = mean_intensity)) +
    
    # A. Base Points (Mapping Color to Group)
    geom_point(aes(color = color_group), size = 1.5, alpha = 0.8) + 
    
    # B. Define Colors manually
    scale_color_manual(values = c(
      "Standard" = "#eb6e34",  # Orange
      "Subcell"  = "#7570b3"   # Purple (Change this hex to whatever you like!)
    )) +
    
    # C. Themes
    theme_classic(base_size = 20) +
    labs(title = comparison_label, y = "Avg Log Intensity", x = "") +
    
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 16, family = "sans"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      legend.position = "none", # Hide the legend (optional)
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin = margin(t = 0.5, r = 3.5, b = 0.5, l = 1, unit = "cm") 
    ) +
    
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    coord_cartesian(clip = "off")
  
  # D. Black Circle Overlay
  labeled_subset <- subset(plot_data, !is.na(label_text))
  
  if (nrow(labeled_subset) > 0) {
    p <- p + geom_point(
      data = labeled_subset,
      aes(x = 0, y = mean_intensity),
      color = "black",
      size = 2.0, 
      stroke = 1.0,
      shape = 1 
    )
    
    # E. Labels
    p <- p + geom_label_repel(
      data = labeled_subset,
      aes(label = label_text),
      xlim = c(0.6, Inf), 
      direction = "y",        
      min.segment.length = 0, 
      force = 0.1,              
      max.overlaps = 50, 
      segment.size = 0.5,
      point.padding = 0.5, 
      na.rm = TRUE
    )
  }
  
  # Save
  clean_name <- gsub(" ", "_", comparison_label)
  ggsave(filename = paste0("Output_InfStrip_", clean_name, ".pdf"), 
         plot = p, width = 3.0, height = 6)
  
  return(p)
}