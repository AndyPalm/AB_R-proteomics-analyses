library(ggplot2)
library(dplyr)
library(ggrepel)

generate_volcano <- function(current_comparison_name, full_data, subcell_list_df = NULL, protein_labels_df = NULL, output_dir = ".") {
  
  # A. Filter data for JUST this comparison
  df_sub <- full_data %>% 
    filter(Comparison == current_comparison_name)
  df_sub$log2FC <- as.numeric(as.character(df_sub$log2FC))
  df_sub$adj.pvalue <- as.numeric(as.character(df_sub$adj.pvalue))
  # B. Assign UP/DOWN status (Using raw p-value for status as per your initial script logic, or adj? 
  # Your snippet used adj.pvalue for the y-axis, but your Setup chunk used p_value for diffexpressed.raw.
  # I will standardize on the logic from your snippet: log2FC and adj.pvalue for the plot axes)
  
  df_sub$diffexpressed.adj <- "NO"
  df_sub$diffexpressed.adj[df_sub$log2FC > 1 & df_sub$adj.pvalue < 0.05] <- "UP"
  df_sub$diffexpressed.adj[df_sub$log2FC < -1 & df_sub$adj.pvalue < 0.05] <- "DOWN"
  
  # C. Handle Custom Highlights (subcells)
  # Default to standard calls
  df_sub$plot_group <- df_sub$diffexpressed.adj 
  
  # Only run this block if a subcell list was provided
  if (!is.null(subcell_list_df)) {
    is_subcell <- df_sub$Protein %in% subcell_list_df$ProteinID
    
    # Mark subcells. Resulting levels: "subcell_UP", "subcell_DOWN", "subcell_NO"
    df_sub$plot_group[is_subcell] <- paste("subcell", df_sub$plot_group[is_subcell], sep = "_")
  }
  
  # Define Factor Levels explicitly to ensure color mapping stays consistent
  possible_levels <- c("NO", "UP", "DOWN", "subcell_NO", "subcell_UP", "subcell_DOWN")
  df_sub$plot_group <- factor(df_sub$plot_group, levels = possible_levels)
  
  # D. Handle Labeling
  df_sub$label_text <- NA
  
  # Only run this block if a labels file was provided
  if (!is.null(protein_labels_df)) {
    matches_selection <- df_sub$Protein %in% protein_labels_df$ProteinID
    
    # Assign the label text only to the matches
    df_sub$label_text[matches_selection] <- df_sub$GeneLabel[matches_selection]
  }
  
  # E. Define Colors (Matching your snippet's hex codes)
  my_colors <- c(
    "NO" = "light grey", 
    "UP" = "#b32d24",       # Orange
    "DOWN" = "#3d7a9c",     # Green
    "subcell_NO" = "#7570b3",     # subcells not sig
    "subcell_UP" = "#7570b3",   # subcells UP (Orange)
    "subcell_DOWN" = "#7570b3"  # subcells DOWN (Green)
  )
  
  # F. Plotting
  p <- ggplot(data = df_sub, aes(x = log2FC, y = -log10(adj.pvalue))) +
    
    # 1. Base points (colored by significance/group)
    geom_point(aes(col = plot_group), size = 1.0, alpha = 0.8) +
    
    # 2. Color Scale
    scale_color_manual(values = my_colors, drop = FALSE) + 
    
    # 3. Reference Lines
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
    
    # 4. Themes and Labels
    theme_classic(base_size = 20) +
    labs(title = current_comparison_name, y = "-log10(adj. p-value)", x = "log2 FC") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 16, family = "sans"), # 'sans' is safer than 'Arial' across OS
          axis.text.y = element_text(size = 16, family = "sans")) +
    scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1))
  
  # 5. The "Black Circle" Overlay (The requested feature)
  # We subset the data to only rows where label_text is NOT NA
  labeled_subset <- subset(df_sub, !is.na(label_text))
  
  if (nrow(labeled_subset) > 0) {
    p <- p + geom_point(
      data = labeled_subset,
      aes(x = log2FC, y = -log10(adj.pvalue)),
      color = "black",
      size = 1.4,
      stroke = 1.0,
      shape = 1 # Hollow circle
    )
  }
  
  # 6. Text Labels (Repel) - Placed on top
  p <- p + geom_label_repel(aes(label = label_text), 
                            segment.color = 'black', 
                            force = 10,
                            max.overlaps = 50, # Increased slightly to ensure labels appear
                            na.rm = TRUE)
  
  # Save the plot
  clean_name <- gsub(" ", "_", current_comparison_name)
  
  # Use the 'path' argument in ggsave to direct the file
  ggsave(
    filename = paste0("Output_", clean_name, ".pdf"), 
    plot = p, 
    path = output_dir,  # <--- This tells ggplot where to save
    width = 8, 
    height = 6
  )
  
  return(p) 
}