# =========================================================================================
# === SCRIPT FOR BATCH BALAGAN ANALYSIS ===
# =========================================================================================
# This script runs the balagan analysis workflow (clustering, heatmap, scatter plot, 
# subsampling) by dynamically reading and looping through a metadata file.

# --- Load Libraries ---
library(dplyr)
library(readr)
library(svglite)
library(SingleCellExperiment)
library(balagan)
library(pheatmap)
library(ggpubr)
library(Polychrome)
library(purrr)


# =========================================================================================
# === Balagan Analysis Function ===
# =========================================================================================
#' Perform Balagan Analysis on a Single Data File
#'
#' This function takes a single data frame (from one CSV file), runs the full
#' balagan analysis workflow (normalization, clustering, visualization), and saves
#' the output plots.
#'
#' @param data_df A data frame containing the cell data for one FOV.
#' @param output_prefix The base path and name for the output files (e.g., "./out_folder/dataScaleSize_slide2_FOV1").
#' @param markers_to_remove A character vector of markers to exclude from the heatmap.
#' @return The function saves plots to disk and prints messages to the console.
perform_balagan_analysis <- function(data_df, output_prefix, markers_to_remove = c()) {
  
  message(paste("Starting balagan analysis for:", basename(output_prefix)))
  
  # 1. Construct the SingleCellExperiment object
  # Columns are renamed to fit the function from the balagan package
  Expression_data <- data_df[, 5:ncol(data_df)]
  Cell_annotation <- data_df[, 1:4]
  Cell_annotation$Cell_size <- Cell_annotation$cellSize
  Cell_annotation$Location_Center_X <- Cell_annotation$X_cent
  Cell_annotation$Location_Center_Y <- Cell_annotation$Y_cent
  
  Gene_annotation <- colnames(Expression_data)
  
  sce <- SingleCellExperiment(
    assays = list(Raw_intensity = as.matrix(t(Expression_data))),
    colData = Cell_annotation,
    rowData = Gene_annotation,
    metadata = list(
      dimension = "2D",
      Bit_mode = 16,
      N_core = 6,
      Is_nuc_cyt = FALSE
    )
  )
  
  # 2. Normalize and Cluster
  sce <- Count_normalization(sce, residual_normalisation = "Pearson")
  sce <- KNN_clustering(
    sce,
    K = 15,
    clustering_method = "Louvain",
    assay_type = "Count_normalised_intensity",
    metric = "L2"
  )
  
  message(paste0("Total number of clusters found: ", length(unique(sce$label))))
  
  # 3. Create and Save the Heatmap
  expr_data <- sce@assays@data$Count_normalised_intensity
  cluster_labels <- sce$label
  df <- cbind(t(expr_data), cluster_labels)
  df <- as.data.frame(df)
  df[, 1:(ncol(df) - 1)] <- apply(df[, 1:(ncol(df) - 1)], 2, as.numeric)
  
  summary_data <- aggregate(. ~ cluster_labels, df, mean)
  summary_data$cluster_labels <- NULL
  summary_data_t <- t(summary_data)
  
  # --- DYNAMIC FEATURE SELECTION ---
  # Plot all markers except those specified for removal
  all_available_markers <- rownames(summary_data_t)
  
  # Clean up marker names from the removal list to match data
  markers_to_remove_cleaned <- gsub("[.-]", "", markers_to_remove)
  
  selected_features <- setdiff(all_available_markers, markers_to_remove_cleaned)
  
  colnames(summary_data_t) <- paste0("C", c(1:length(unique(sce$label))))
  
  # Save heatmap
  pheatmap_file <- paste0(output_prefix, "_heatmap.svg")
  svglite(pheatmap_file, width = 10, height = 8)
  pheatmap(
    summary_data_t[selected_features, ],
    scale = "row",
    cluster_cols = FALSE,
    angle_col = 45,
    main = paste("Cluster Heatmap for", basename(output_prefix))
  )
  dev.off()
  message(paste("Saved heatmap to:", pheatmap_file))
  
  # 4. Create and Save the Scatter Plot
  scatter_data <- data_df[, c("X_cent", "Y_cent")]
  scatter_data$cluster <- as.factor(sce$label)
  
  # Generate a unique color for each cluster
  num_clusters <- length(unique(scatter_data$cluster))
  color_vector <- NULL
  if (num_clusters <= 35) {
    # Use the predefined palette if we have enough colors
    color_vector <- setNames(as.character(palette36.colors(36)[-2]), sort(unique(scatter_data$cluster)))
  } else {
    # Fallback to a different color generator for more than 35 clusters
    color_vector <- scales::hue_pal()(num_clusters)
  }
  
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
    labs(title = paste("Cell Clusters for", basename(output_prefix))) +
    theme(legend.position = "right")
  
  scatter_file <- paste0(output_prefix, "_scatter.svg")
  ggsave(
    filename = scatter_file,
    plot = pp1,
    width = 12,
    height = 10
  )
  message(paste("Saved scatter plot to:", scatter_file))
  
  # 5. Run and Visualize Subsampling Analysis
  sce$ImageNumber <- 1
  Simple_sampling_analysis <- Perform_sampling_analysis(
    sce,
    Selected_image = 1,
    N_times = 20,
    N_sampling_region_vector = 1:20,
    width_FOV_vector = 400,
    height_FOV_vector = 400,
    Threshold_detection_cluster = 2
  )
  

  subsampling_file <- paste0(output_prefix, "_subsampling.svg")
  svg(paste0(output_prefix, "_subsampling.svg"),
      width = 10,
      height = 8)
  Visualize_simple_sampling(Simple_sampling_analysis)
  dev.off()
  message(paste("Saved subsampling plot to:", subsampling_file))
  
  message("---------------------------------------------------\n")
}


# =========================================================================================
# --- ANALYSIS EXECUTION ---
# =========================================================================================

### CHOOSE WHICH DATASETS TO RUN ###
# Define the specific list of "Source" names you want to analyze.
# These must match the Source column in data_mesmer/Slide_metadata.csv
sources_to_process <- c(
  # Initial Optimization datasets
  "Roche", "BIDMC", "Stanford",
  # Validation datasets - BIDMC
  "BIDMC_Tonsil", "BIDMC_DLBCL",
  # Validation datasets - Roche
  "Roche_Tonsil", "Roche_intestine",
  # Validation datasets - UKentucky
  "UKentucky_Tonsil", "UKentucky_SCC",
  # Validation datasets - ASTAR
  "ASTAR_COMET_CRC", "ASTAR_COMET_Tonsil",
  # Validation datasets - Stanford IMC
  "Stanford_IMC_Tonsil", "Stanford_IMC_OSCC",
  # Validation datasets - Stanford MIBI
  "Stanford_MIBI_Colon", "Stanford_MIBI_Liver", "Stanford_MIBI_LymphNode",
  # Validation datasets - Stanford Orion
  "Stanford_Orion_LN", "Stanford_Orion_EndometrialCancer",
  # Validation datasets - Novartis
  "Novartis_Tonsil", "Novartis_Lung_Cancer",
  # Special datasets
  "LyophilizationTest_FigS2", "Reimagedslide_FigS5", "StorageConditionsExpt"
)

# --- Load Metadata ---
# Note: Balagan analysis uses MESMER segmentation data
metadata_file <- "./data_mesmer/Slide_metadata.csv"
removal_file <- "./data_mesmer/Slide_remove_markers.csv"

if (!file.exists(metadata_file) || !file.exists(removal_file)) {
  stop("Metadata or marker removal file not found. Please check file paths.")
}
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)

cat("Starting batch analysis for", length(sources_to_process), "selected sources.\n")

# --- Main Analysis Loop ---
# This loop iterates through each source specified in 'sources_to_process'.
for(current_source in sources_to_process) {
  
  cat("\n====================================================\n")
  cat("Processing Source:", current_source, "\n")
  
  # Filter metadata for the current source
  metadata_for_source <- slide_metadata %>% 
    filter(Source == current_source, Type == "dataScaleSize")
  
  if (nrow(metadata_for_source) == 0) {
    warning(paste("No data found for source:", current_source, "- Skipping."))
    next # Skip to the next source in the list
  }
  
  # --- FOR TESTING: ONLY RUN ON THE FIRST FILE ---
  # This takes the first row from the filtered metadata for the current source.
  current_row <- metadata_for_source %>% dplyr::slice(1)
  current_filename <- current_row$Filename
  cat("TEST MODE: Running on first file only:", current_filename, "\n")
  
  # Dynamically determine data and output paths
  # Input: MESMER segmentation data from data_mesmer/
  data_folder <- file.path("./data_mesmer", current_source)
  out_folder <- file.path("./out_balagan_analysis", current_source)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Get the list of markers to remove for this specific source
  current_remove_markers <- removed_markers %>%
    filter(Source == current_source, Exclude_type == "Marker") %>%
    pull(Exclude_value)
  
  # Construct the full path to the input CSV file
  full_file_path <- file.path(data_folder, current_filename)
  
  if (file.exists(full_file_path)) {
    # Read the data for the current file
    fov_data <- read.csv(full_file_path)
    
    # Define the prefix for all output files (e.g., "out_folder/dataScaleSize_slide2_FOV1")
    output_prefix <- file.path(out_folder, sub("\\.csv$", "", current_filename))
    
    # Use a try-catch block to handle potential errors in the analysis for a single file
    tryCatch({
      # Run the entire analysis workflow on this single file
      perform_balagan_analysis(
        data_df = fov_data,
        output_prefix = output_prefix,
        markers_to_remove = current_remove_markers
      )
    }, error = function(e) {
      warning(paste("An error occurred while processing:", current_filename, "\nError:", e$message))
    })
    
  } else {
    warning(paste("File not found, skipping:", full_file_path))
  }
}

cat("\n====================================================\n")
cat("Batch analysis completed for all selected sources.\n")
cat("====================================================\n")

