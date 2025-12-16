# DIA_feature_extract_boxwhisker.R

# plot_utils.R

# -------------------------------------------------------------------
# Helper Function 1: Plot a single gene
# -------------------------------------------------------------------
plot_top4_features <- function(data, gene_name, output_dir) {
  
  require(tidyverse)
  require(ggplot2)
  
  message(paste("Processing plot for:", gene_name))
  
  # 1. Filter for specific gene
  gene_data <- data %>% filter(Gene == gene_name)
  
  # 2. Select Top 4 Features Logic
  top4_features <- gene_data %>%
    group_by(FEATURE) %>%
    summarize(mean_abundance = mean(ABUNDANCE, na.rm = TRUE)) %>%
    slice_max(order_by = mean_abundance, n = 4) %>%
    pull(FEATURE)
  
  # Filter to just these features
  plot_data <- gene_data %>%
    filter(FEATURE %in% top4_features)
  
  # 3. Fix Factor Levels (User-specified order)
  # We use the exact names and order requested: Ctx, ALOD4, OlyA, control
  plot_data <- plot_data %>%
    mutate(GROUP = fct_relevel(GROUP, "Ctx", "ALOD4", "OlyA", "control")) %>% 
    arrange(GROUP, SUBJECT) %>%
    mutate(SUBJECT = factor(SUBJECT, levels = unique(SUBJECT)))
  
  # 4. Generate the Plot (Bar + Error + Points)
  p <- ggplot(plot_data, aes(x = SUBJECT, y = ABUNDANCE)) +
    
    # Bar layer (Mean)
    stat_summary(geom = "bar", fun = mean, aes(fill = GROUP), 
                 alpha = 0.6, color = "black", width = 0.7) +
    
    # Error bar layer (Mean +/- 1 SD)
    # Note: changed 'size' to 'linewidth' to fix warning
    stat_summary(geom = "errorbar", fun.data = mean_sdl, 
                 fun.args = list(mult = 1), width = 0.2, linewidth = 0.8) +
    
    # Points layer (Individual features)
    geom_jitter(width = 0.1, shape = 1, size = 2, color = "black", stroke = 1) +
    
    # Zoom Y-axis
    coord_cartesian(ylim = c(18, 30)) +
    
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
      legend.position = "bottom"
    )
  
  # 5. Save the Plot
  safe_name <- gsub("[^[:alnum:]]", "_", gene_name)
  filename <- file.path(output_dir, paste0("Barplot_Top4_", safe_name, ".png"))
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  message(paste("Saved plot to:", filename))
  
  return(p)
}

# -------------------------------------------------------------------
# Helper Function 2: Batch process a list of genes
# -------------------------------------------------------------------
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