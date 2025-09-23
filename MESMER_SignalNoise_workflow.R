library(dplyr)
library(tidyverse)
library(matrixStats)
library(ggcorrplot)
library(ggpubr)
library(tidyr)
library(rstatix)
library(readr)
library(presto)
library(svglite)

source("helper.R")
########################################## Configuration Begin  #############################################################################

# Define paths for metadata and exclusion files
metadata_file <- "./data_03272025/Slide_metadata.csv"
removal_file <- "./data_03272025/Slide_remove_markers.csv"
exclusion_file <- "./data_03272025/Slide_exclude_markers.csv"
marker_sequence_file <- "./data_03272025/Registered_Report_marker_sequence.csv"

# Cell count files (pre-calculated by Process_Cell_Counts.R)
bidmc_cell_counts_file <- "./cell_counts/BIDMC/fov_cell_counts.csv"
roche_cell_counts_file <- "./cell_counts/Roche/fov_cell_counts.csv"
stanford_cell_counts_file <- "./cell_counts/Stanford/fov_cell_counts.csv"

# Load metadata and exclusions
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)
excluded_markers <- read_csv(exclusion_file)
marker_sequence <- read_csv(marker_sequence_file, col_names = TRUE)
marker_sequence <- marker_sequence[[1]]  # Get the first column as a vector

# Check if pre-calculated cell counts exist
if (!file.exists(bidmc_cell_counts_file) || !file.exists(roche_cell_counts_file)) {
  warning("Some pre-calculated cell counts not found. Please run Process_Cell_Counts.R for missing sources.")
}

# Load cell counts (with checks to avoid errors if files don't exist)
bidmc_cell_counts <- if (file.exists(bidmc_cell_counts_file)) read_csv(bidmc_cell_counts_file) else NULL
roche_cell_counts <- if (file.exists(roche_cell_counts_file)) read_csv(roche_cell_counts_file) else NULL
stanford_cell_counts <- if (file.exists(stanford_cell_counts_file)) read_csv(stanford_cell_counts_file) else NULL

# Filter metadata to get signal_ratios files
bidmc_snr_data <- slide_metadata %>% 
  filter(Source == "BIDMC" & Type == "signal_ratios")
roche_snr_data <- slide_metadata %>% 
  filter(Source == "Roche" & Type == "signal_ratios")
stanford_snr_data <- slide_metadata %>% 
  filter(Source == "Stanford" & Type == "signal_ratios")

# Get markers to remove for each source
roche_remove_markers <- removed_markers %>% 
  filter(Source == "Roche" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
bidmc_remove_markers <- removed_markers %>% 
  filter(Source == "BIDMC" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
stanford_remove_markers <- removed_markers %>% 
  filter(Source == "Stanford" & Exclude_type == "Marker") %>%
  pull(Exclude_value)

# Create configurations for each dataset
configurations <- list(
  BIDMC_Combined = list(
    data_folder = "./data_03272025/BIDMC/",
    out_folder = "./out_BIDMC_Combined_SNR/",
    snr_data = bidmc_snr_data,
    cell_counts = bidmc_cell_counts,
    remove_values = bidmc_remove_markers
  ),
  
  Roche_Combined = list(
    data_folder = "./data_03272025/Roche/",
    out_folder = "./out_Roche_Combined_SNR/",
    snr_data = roche_snr_data,
    cell_counts = roche_cell_counts,
    remove_values = roche_remove_markers
  ),
  
  Stanford_Combined = list(
    data_folder = "./data_03272025/Stanford/",
    out_folder = "./out_Stanford_Combined_SNR/",
    snr_data = stanford_snr_data,
    cell_counts = stanford_cell_counts,
    remove_values = stanford_remove_markers
  )
)

# Choose which configuration to use
# Options: "BIDMC_Combined", "Roche_Combined", "Stanford_Combined"
current_config_name <- "Roche_Combined"
current_config <- configurations[[current_config_name]]

########################################## Configuration End  #############################################################################
# Set up working environment
data_folder <- current_config$data_folder
out_folder <- current_config$out_folder
snr_data_meta <- current_config$snr_data
cell_counts <- current_config$cell_counts
remove_values <- current_config$remove_values

# Clean up marker names to match data processing
remove_values <- gsub("[.-]", "", remove_values)

# Create output directory if it doesn't exist
dir.create(out_folder, showWarnings = FALSE)

# Save configuration details for reference
write_csv(
  tibble(
    config_name = current_config_name,
    data_folder = data_folder,
    out_folder = out_folder,
    file_count = nrow(snr_data_meta)
  ),
  paste0(out_folder, "config_summary.csv")
)

# Main workflow execution
cat("Processing combined FOV data for", current_config_name, "...\n")

# Function to combine SNR data with cell count weights
combine_snr_with_cell_weights <- function(data_folder, snr_meta, cell_counts) {
  # Group SNR metadata by FOV
  fov1_data <- snr_meta %>% filter(FOV == "FOV1")
  fov2_data <- snr_meta %>% filter(FOV == "FOV2")
  
  # Load SNR data for each FOV
  cat("Loading FOV1 SNR data...\n")
  fov1_snr <- load_signal_ratio_data(data_folder, fov1_data$Filename, fov1_data$Name)
  
  cat("Loading FOV2 SNR data...\n")
  fov2_snr <- load_signal_ratio_data(data_folder, fov2_data$Filename, fov2_data$Name)
  
  # Get all unique markers across both datasets
  fov1_markers <- setdiff(colnames(fov1_snr), "Staining_condition")
  fov2_markers <- setdiff(colnames(fov2_snr), "Staining_condition")
  all_markers <- unique(c(fov1_markers, fov2_markers))
  
  # Get all unique conditions
  all_conditions <- unique(c(fov1_snr$Staining_condition, fov2_snr$Staining_condition))
  
  # Create output dataframe
  combined_snr <- data.frame(Staining_condition = all_conditions)
  
  # Add all marker columns
  for (marker in all_markers) {
    combined_snr[[marker]] <- NA_real_
  }
  
  # For each condition, calculate weighted average SNR
  for (i in seq_along(all_conditions)) {
    condition <- all_conditions[i]
    
    # Get corresponding cell counts and calculate weights
    condition_counts <- cell_counts %>%
      filter(Staining_condition == condition)
    
    # Skip if we don't have cell count data for this condition
    if (nrow(condition_counts) == 0) {
      cat("Warning: No cell count data for condition:", condition, "\n")
      next
    }
    
    # Get FOV-specific SNR values
    fov1_row <- which(fov1_snr$Staining_condition == condition)
    fov2_row <- which(fov2_snr$Staining_condition == condition)
    
    # For each marker, calculate weighted average
    for (marker in all_markers) {
      # Initialize values
      fov1_val <- NA
      fov2_val <- NA
      
      # Get FOV1 value if available
      if (length(fov1_row) > 0 && marker %in% colnames(fov1_snr)) {
        fov1_val <- fov1_snr[fov1_row, marker]
      }
      
      # Get FOV2 value if available
      if (length(fov2_row) > 0 && marker %in% colnames(fov2_snr)) {
        fov2_val <- fov2_snr[fov2_row, marker]
      }
      
      # Only proceed if we have at least one value
      if (is.na(fov1_val) && is.na(fov2_val)) {
        next
      }
      
      # Get cell counts and proportions
      fov1_count <- condition_counts %>% 
        filter(FOV == "FOV1") %>% 
        pull(Cell_count)
      
      fov2_count <- condition_counts %>% 
        filter(FOV == "FOV2") %>% 
        pull(Cell_count)
      
      # Handle cases where we don't have both FOVs
      if (length(fov1_count) == 0) fov1_count <- 0
      if (length(fov2_count) == 0) fov2_count <- 0
      
      total_count <- fov1_count + fov2_count
      
      # Calculate weighted average
      if (total_count > 0) {
        # If we have both values, use weighted average
        if (!is.na(fov1_val) && !is.na(fov2_val)) {
          # Calculate using proportions: (S1 * (N1/N3) + S2 * (N2/N3))
          weighted_avg <- (fov1_val * (fov1_count/total_count)) + 
            (fov2_val * (fov2_count/total_count))
        } 
        # If we only have FOV1, use that
        else if (!is.na(fov1_val)) {
          weighted_avg <- fov1_val
        } 
        # If we only have FOV2, use that
        else if (!is.na(fov2_val)) {
          weighted_avg <- fov2_val
        }
        
        combined_snr[i, marker] <- weighted_avg
      }
    }
  }
  
  return(combined_snr)
}

# Combine SNR data using cell count weights
combined_snr <- combine_snr_with_cell_weights(
  data_folder = data_folder,
  snr_meta = snr_data_meta,
  cell_counts = cell_counts
)

# Save the combined SNR data
write_csv(combined_snr, paste0(out_folder, "combined_weighted_snr.csv"))

# Create visualizations
cat("Creating signal-to-noise ratio heatmaps...\n")
heatmap_results <- create_snr_heatmaps(
  snr_data = combined_snr,
  out_folder = out_folder,
  remove_values = remove_values
)

# Print completion message
cat("\nAnalysis completed for", current_config_name, "\n")
cat("Results saved to:", out_folder, "\n")

# Function to process all configurations sequentially
process_all_configs <- function() {
  for (config_name in names(configurations)) {
    cat("\n\n==================================\n")
    cat("Processing configuration:", config_name, "\n")
    cat("==================================\n\n")
    
    # Set the current configuration
    current_config <- configurations[[config_name]]
    
    # Set up working environment
    data_folder <- current_config$data_folder
    out_folder <- current_config$out_folder
    snr_data_meta <- current_config$snr_data
    cell_counts <- current_config$cell_counts
    remove_values <- current_config$remove_values
    
    # Clean up marker names
    remove_values <- gsub("[.-]", "", remove_values)
    
    # Create output directory
    dir.create(out_folder, showWarnings = FALSE)
    
    # Save configuration details
    write_csv(
      tibble(
        config_name = config_name,
        data_folder = data_folder,
        out_folder = out_folder,
        file_count = nrow(snr_data_meta)
      ),
      paste0(out_folder, "config_summary.csv")
    )
    
    # Combine SNR data using cell count weights
    combined_snr <- combine_snr_with_cell_weights(
      data_folder = data_folder,
      snr_meta = snr_data_meta,
      cell_counts = cell_counts
    )
    
    # Save the combined SNR data
    write_csv(combined_snr, paste0(out_folder, "combined_weighted_snr.csv"))
    
    # Create visualizations
    heatmap_results <- create_snr_heatmaps(
      snr_data = combined_snr,
      out_folder = out_folder,
      remove_values = remove_values
    )
    
    cat("Analysis completed for", config_name, "\n")
    cat("Results saved to:", out_folder, "\n")
  }
  
  cat("\nAll configurations processed successfully!\n")
}

# Uncomment the line below to process all configurations
process_all_configs()
