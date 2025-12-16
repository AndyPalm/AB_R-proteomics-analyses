# DIA_feature_extract_barplots.R

plot_top4_features <- function(data, gene_name, output_dir) {
  
  require(tidyverse)
  require(ggplot2)
  require(ggbeeswarm) 
  
  message(paste("Processing plot for:", gene_name))
  
  # 1. Filter for specific gene
  gene_data <- data %>% dplyr::filter(Gene == gene_name)
  
  # 2. Select Top 4 Features Logic
  top4_features <- gene_data %>%
    dplyr::group_by(FEATURE) %>%
    dplyr::summarize(mean_abundance = mean(ABUNDANCE, na.rm = TRUE)) %>%
    dplyr::slice_max(order_by = mean_abundance, n = 4) %>%
    dplyr::pull(FEATURE)
  
  # Filter to just these features
  plot_data <- gene_data %>%
    dplyr::filter(FEATURE %in% top4_features)
  
  # 3. Fix Factor Levels
  plot_data <- plot_data %>%
    dplyr::mutate(GROUP = fct_relevel(GROUP, "Ctx", "ALOD4", "OlyA", "control")) %>% 
    dplyr::arrange(GROUP, SUBJECT) %>%
    dplyr::mutate(SUBJECT = factor(SUBJECT, levels = unique(SUBJECT)))
  
  # 4. Generate the Plot
  p <- ggplot(plot_data, aes(x = SUBJECT, y = ABUNDANCE)) +
    
    # Bar layer
    stat_summary(geom = "bar", fun = mean, aes(fill = GROUP), 
                 alpha = 0.6, color = "black", width = 0.7) +
    
    # Error bar layer
    stat_summary(geom = "errorbar", fun.data = mean_sdl, 
                 fun.args = list(mult = 1), width = 0.2, linewidth = 0.8) +
    
    # adjust width to vary spread of individual data points
    ggbeeswarm::geom_quasirandom(
      shape = 1, 
      size = 1.5, 
      color = "black", 
      stroke = 1, 
      width = 0.35
    ) +
    
    # Zoom Y-axis if desired
    #coord_cartesian(ylim = c(19, 29)) +
    
    # Aesthetics
    theme_bw() +
    labs(
      title = paste("Top 4 Features: ", gene_name),
      subtitle = "Bar = Mean; Error Bar = +/- 1 SD; Points = Individual Features",
      y = "Log2 Abundance",
      x = "Replicate"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 5. Save as PDF
  safe_name <- gsub("[^[:alnum:]]", "_", gene_name)
  filename <- file.path(output_dir, paste0("296_Barplot_Top4_", safe_name, ".pdf"))
  
  ggsave(filename, plot = p, width = 5, height = 6)
  message(paste("Saved plot to:", filename))
  
  return(p)
}

# The batch_plot_genes function remains exactly the same
batch_plot_genes <- function(gene_list, data, output_dir) {
  require(purrr)
  message(paste("Starting batch processing for", length(gene_list), "genes..."))
  walk(gene_list, function(gene) {
    if (gene %in% data$Gene) {
      plot_top4_features(data, gene, output_dir)
    } else {
      warning(paste("Skipping", gene, "- not found in dataset."))
    }
  })
  message("Batch processing complete.")
}