# ===================================================================
# === SCRIPT TO REGENERATE PLOTS FROM SAVED CSV DATA ===
# ===================================================================
# This script finds all saved data files (*_heatmap_data.csv and
# *_scatter_data.csv) in the output directory and uses them to
# recreate the heatmap and scatter plots.
#
# This allows for tweaking plot parameters without re-running the
# entire analysis pipeline.
# ===================================================================

# --- Load Required Libraries ---
library(dplyr)
library(readr)
library(svglite)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(Polychrome) # For the palette36.colors

# ===================================================================
# --- 1. SETUP PATHS ---
# ===================================================================

# IMPORTANT: This must point to the directory where your results are.
# This is both the INPUT and OUTPUT directory.
results_dir <- "/Users/wang.13246/Documents/Project/Sizun_NM_revision/MESMER_workflow_03272025/balagan_results1/out_balagan_analysis_BIDMC_phenotypic_markers_run2/BIDMC"

cat(paste("Looking for data and saving plots in:", results_dir, "\n"))

# ===================================================================
# --- 2. REGENERATE HEATMAPS ---
# ===================================================================
cat("\n--- Starting Heatmap Regeneration ---\n")

# Find all heatmap data files
heatmap_data_files <- list.files(
  path = results_dir,
  pattern = "_heatmap_data.csv$",
  full.names = TRUE # Get the full file path
)

if (length(heatmap_data_files) > 0) {
  cat(paste("Found", length(heatmap_data_files), "heatmap data files to process.\n"))
  
  for (csv_path in heatmap_data_files) {
    message(paste("Processing heatmap for:", basename(csv_path)))
    
    # Load the data
    heatmap_data <- read.csv(csv_path, row.names = 1, check.names = FALSE)
    
    # Define the output SVG filename
    svg_output_path <- gsub("_heatmap_data.csv$", "_heatmap.svg", csv_path)
    
    # Create a dynamic title
    slide_basename <- basename(gsub("_heatmap_data.csv$", "", csv_path))
    plot_title <- paste("Heatmap for", slide_basename)
    
    # Generate and save the heatmap (using parameters from your original script)
    svglite(svg_output_path, width = 10, height = 8)
    pheatmap(
      heatmap_data,
      scale = "row",
      cluster_cols = FALSE,
      angle_col = 45,
      main = plot_title
    )
    dev.off() # Close the SVG device
  }
  cat("--- Heatmap regeneration complete. ---\n")
} else {
  warning(paste("No heatmap data files ('*_heatmap_data.csv') found in:", results_dir))
}

# ===================================================================
# --- 3. REGENERATE SCATTER PLOTS ---
# ===================================================================
cat("\n--- Starting Scatter Plot Regeneration ---\n")

# Find all scatter data files
scatter_data_files <- list.files(
  path = results_dir,
  pattern = "_scatter_data.csv$",
  full.names = TRUE # Get the full file path
)

if (length(scatter_data_files) > 0) {
  cat(paste("Found", length(scatter_data_files), "scatter data files to process.\n"))
  
  for (csv_path in scatter_data_files) {
    csv_path <- "/Users/wang.13246/Documents/Project/Sizun_NM_revision/MESMER_workflow_03272025/balagan_results1/out_balagan_analysis_BIDMC_phenotypic_markers_run2/BIDMC/dataScaleSize_slide13_FOV1_scatter_data.csv"
    message(paste("Processing scatter plot for:", basename(csv_path)))
    
    # Load the data
    scatter_data <- read.csv(csv_path) %>%
      mutate(cluster = as.factor(cluster)) # Ensure cluster is a factor
    
    # Define the output SVG filename
    scatter_file <- gsub("_scatter_data.csv$", "_scatter.svg", csv_path)
    
    # Create a dynamic title
    slide_basename <- basename(gsub("_scatter_data.csv$", "", csv_path))
    plot_title <- paste("Cell Clusters for", slide_basename)
    
    # Re-create the color vector (from your original script)
    num_clusters <- length(unique(scatter_data$cluster))
    color_vector <- if (num_clusters <= 35) {
      setNames(as.character(palette36.colors(36)[-2]), sort(unique(scatter_data$cluster)))
    } else {
      scales::hue_pal()(num_clusters)
    } 
    if(slide_basename == "dataScaleSize_slide13_FOV1") {
      col_vec <- as.character(palette36.colors(36)[-2])
      color_vector <- setNames(as.character(palette36.colors(36)[-2]), sort(unique(scatter_data$cluster)))
      color_vector[3] <- col_vec[11]
      color_vector[4] <- col_vec[9]
    }
    
    # Re-build the plot
    pp1 <- ggscatter(
      scatter_data,
      x = "X_cent",
      y = "Y_cent",
      color = "cluster",
      xlab = "X Centroid",
      ylab = "Y Centroid",
      size = 0.5
    ) +
      scale_color_manual(values = color_vector) +
      labs(title = plot_title) +
      theme(legend.position = "right") +
      scale_y_reverse()
    
    # Save the plot (using parameters from your original script)
    ggsave(
      filename = scatter_file,
      plot = pp1,
      width = 12,
      height = 10
    )
  }
  cat("--- Scatter plot regeneration complete. ---\n")
} else {
  warning(paste("No scatter data files ('*_scatter_data.csv') found in:", results_dir))
}

cat("\n====================================================\n")
cat("All plot regeneration is complete.\n")
cat("====================================================\n")
