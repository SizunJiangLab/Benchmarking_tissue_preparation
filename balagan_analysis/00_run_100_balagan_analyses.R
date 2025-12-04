# =========================================================================================
# === MASTER SCRIPT FOR RUNNING 100 BALAGAN ANALYSES ===
# =========================================================================================
# This script wraps the entire Balagan analysis workflow in a loop to run it multiple times.
# Each run uses a different random seed and saves its results to a unique directory,
# allowing for the assessment of result stability and the calculation of a stable
# average performance rank.

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
library(ggplot2)

# =========================================================================================
# --- MASTER CONFIGURATION ---
# =========================================================================================
NUMBER_OF_RUNS <- 100

# =========================================================================================
# --- MASTER ANALYSIS LOOP ---
# =========================================================================================

for (run_number in 1:NUMBER_OF_RUNS) {
  
  cat("\n\n#########################################################################\n")
  cat(paste("###   STARTING RUN", run_number, "OF", NUMBER_OF_RUNS, "  ###\n"))
  cat("#########################################################################\n\n")
  
  # --- Set a unique, reproducible seed for this specific run ---
  set.seed(run_number)
  
  # --- Define a unique output directory for this run ---
  output_base_dir <- paste0("./out_balagan_analysis_BIDMC_run_", run_number)
  
  # --- Setup Paths and Metadata ---
  metadata_file <- "./data_mesmer/Slide_metadata.csv"
  removal_file <- "./data_mesmer/Slide_remove_markers.csv"
  
  slide_metadata <- read_csv(metadata_file, show_col_types = FALSE)
  removed_markers <- read_csv(removal_file, show_col_types = FALSE)
  
  # --- Filter for the specific BIDMC dataset ---
  metadata_to_process <- slide_metadata %>%
    filter(Source == "BIDMC", Type == "dataScaleSize", FOV == "FOV1")
  
  data_folder <- file.path("./data_mesmer", "BIDMC")
  out_folder <- file.path(output_base_dir, "BIDMC")
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  bidmc_remove_markers <- removed_markers %>%
    filter(Source == "BIDMC", Exclude_type == "Marker") %>%
    pull(Exclude_value)
  
  # --- Initialize data frame to collect results for THIS RUN ---
  all_slides_tau_results <- data.frame()
  
  # --- Main Analysis Loop (Iterates through each slide once) ---
  cat(paste("\n--- [Run", run_number, "] Starting Complex Sampling Batch Analysis ---\n"))
  
  for (i in 1:nrow(metadata_to_process)) {
    current_row <- metadata_to_process[i, ]
    current_filename <- current_row$Filename
    full_file_path <- file.path(data_folder, current_filename)
    
    cat(paste("\n--- [Run", run_number, "] Processing Slide:", current_filename, "---\n"))
    
    if (file.exists(full_file_path)) {
      fov_data <- read.csv(full_file_path)
      output_prefix <- file.path(out_folder, sub("\\.csv$", "", current_filename))
      
      tryCatch({
        # 1. Construct SCE, Normalize, and Cluster
        phenotypic_markers <- c("CD3", "CD15", "CD8", "CD20", "CD11c", "CD68", "FoxP3", "Pax5", "CD31", "Cytokeratin")
        Expression_data_all <- fov_data[, 5:ncol(fov_data)]
        available_phenotypic_markers <- intersect(phenotypic_markers, colnames(Expression_data_all))
        Expression_data <- Expression_data_all[, available_phenotypic_markers]
        
        Cell_annotation <- fov_data[, 1:4]
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
            N_core = 12,
            Is_nuc_cyt = FALSE
          )
        )
        
        sce <- Count_normalization(sce, residual_normalisation = "Pearson")
        sce <- KNN_clustering(sce, K = 15, clustering_method = "Louvain", assay_type = "Count_normalised_intensity")
        
        # 2. Define Vectors for Complex Sampling Analysis
        fov_sizes <- seq(50, 500, 50)
        n_sampling_regions <- 1:20
        height_vector <- rep(fov_sizes, each = length(n_sampling_regions))
        width_vector <- rep(fov_sizes, each = length(n_sampling_regions))
        N_sampling_region_vector <- rep(n_sampling_regions, length(fov_sizes))
        
        # 3. Perform Complex Sampling Analysis
        sce$ImageNumber <- 1
        Complex_sampling <- Perform_sampling_analysis(
          sce, Selected_image = 1, N_times = 20,
          N_sampling_region_vector = N_sampling_region_vector,
          width_FOV_vector = width_vector,
          height_FOV_vector = height_vector,
          Threshold_detection_cluster = 2
        )
        
        write.csv(Complex_sampling, paste0(output_prefix, "_Complex_sampling.csv"))
        
        # 4. Extract Tau Values
        Parameter_table <- data.frame(Height = height_vector, Width = width_vector)
        Fitting_tau <- Visualize_complex_sampling(Complex_sampling, Parameter_table)
        
        fitting_tau_csv_path <- paste0(output_prefix, "_Fitting_tau.csv")
        write.csv(as.data.frame(Fitting_tau), file = fitting_tau_csv_path, row.names = FALSE)
        message(paste("--- Saved Fitting_tau data to:", basename(fitting_tau_csv_path)))
        
        # 5. Store Results for Final Aggregation
        slide_results <- as.data.frame(Fitting_tau) %>%
          mutate(
            Slide = current_filename,
            FoV_width = unique(width_vector)
          )
        
        all_slides_tau_results <- bind_rows(all_slides_tau_results, slide_results)
        
      }, error = function(e) {
        warning(paste("An error occurred for:", current_filename, "in run", run_number, "\nError:", e$message))
      })
    } else {
      warning(paste("File not found, skipping:", full_file_path))
    }
  }
  
  cat(paste("\n--- [Run", run_number, "] All slides processed. Aggregating results for this run. ---\n"))
  
  # --- AGGREGATE ANALYSIS FOR THIS RUN ---
  if (nrow(all_slides_tau_results) > 0) {
    # Save the combined tau data table FOR THIS RUN
    combined_csv_path <- file.path(out_folder, "AGGREGATE_all_slides_tau_per_fov.csv")
    write.csv(all_slides_tau_results, combined_csv_path, row.names = FALSE)
  }
  
  cat(paste("\n--- [Run", run_number, "] Analysis complete. Results saved to:", out_folder, "---\n"))
  
} # End of master loop

cat("\n\n#########################################################################\n")
cat(paste("###   ALL", NUMBER_OF_RUNS, "RUNS COMPLETED SUCCESSFULLY   ###\n"))
cat("#########################################################################\n")
