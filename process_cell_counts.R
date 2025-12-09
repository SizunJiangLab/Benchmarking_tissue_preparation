library(dplyr)
library(tidyverse)
library(readr)
library(qs)

########################################## Configuration Begin  #############################################################################

# Define paths for metadata files
metadata_file <- "./data_mesmer/Slide_metadata.csv"

#metadata_file <- "./data_cellXpress/Slide_metadata.csv"

# Load metadata 
slide_metadata <- read_csv(metadata_file)

# Define sources to process
sources <- c("BIDMC", "Roche", "Stanford", "Stanford-scan1")

# Output directories
out_dir <- "./data_mesmer/Initial_Optimization/cell_counts"

########################################## Configuration End  #############################################################################

# Create main output directory if it doesn't exist
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source = sources[3]
# Process each source
for (source in sources) {
  cat("\n==================================\n")
  cat("Processing cell counts for", source, "\n")
  cat("==================================\n")
  
  # Create source-specific output directory
  source_dir <- paste0(out_dir, source, "/")
  dir.create(source_dir, showWarnings = FALSE)
  
  # Filter metadata to get dataScaleSize files for this source
  source_data <- slide_metadata %>% 
    filter(Source == source & Type == "dataScaleSize")
  
  # Check if we have any files
  if (nrow(source_data) == 0) {
    cat("No dataScaleSize files found for", source, "\n")
    next
  }
  
    cat("Loading dataScaleSize files for", source, "...\n")
    
    # Get data folder path - check for nested structure
    base_data_path <- "./data_mesmer/"
    possible_paths <- c(
      paste0(base_data_path, source, "/"),
      paste0(base_data_path, "Initial_Optimization/", source, "/"),
      paste0(base_data_path, "Validation/", source, "/")
    )
    
    data_folder <- possible_paths[1] # Default to first
    for (path in possible_paths) {
      if (dir.exists(path)) {
        data_folder <- path
        break
      }
    }
    
    # Initialize list to store cell counts
    cell_data_list <- list()
    
    # Process each file
    for (i in 1:nrow(source_data)) {
      filename <- source_data$Filename[i]
      condition_name <- source_data$Name[i]
      fov <- source_data$FOV[i]
      
      file_path <- paste0(data_folder, filename)
      
      # Check if file exists
      if (!file.exists(file_path)) {
        cat("Warning: File not found:", file_path, "\n")
        next
      }
      
      # Read file and count cells
      tryCatch({
        data <- read.csv(file_path, header = TRUE)
        cell_count <- nrow(data)
        
        cell_data_list[[length(cell_data_list) + 1]] <- data.frame(
          Filename = filename,
          Staining_condition = condition_name,
          FOV = fov,
          Cell_count = cell_count
        )
        
        cat("Processed:", filename, "- Found", cell_count, "cells\n")
      }, error = function(e) {
        cat("Error processing file:", filename, "- Error:", e$message, "\n")
      })
    }
    
    # Combine all cell count data
    if (length(cell_data_list) > 0) {
      cell_data <- bind_rows(cell_data_list)
    } else {
      cat("No valid files processed for", source, "\n")
      next
    }
    
    # Calculate total cells and proportions
    cell_counts <- cell_data %>%
      group_by(Staining_condition) %>%
      mutate(
        Total_cells = sum(Cell_count),
        Cell_proportion = Cell_count / Total_cells
      ) %>%
      ungroup()
    
    # Save all cell counts to CSV
    write_csv(cell_counts, paste0(source_dir, "fov_cell_counts.csv"))
    
    # Generate summary statistics by condition
    cell_summary <- cell_counts %>%
      group_by(Staining_condition) %>%
      summarise(
        Total_cells = first(Total_cells),
        FOV_count = n(),
        Min_cells = min(Cell_count),
        Max_cells = max(Cell_count),
        Mean_cells = mean(Cell_count),
        Median_cells = median(Cell_count)
      ) %>%
      ungroup() %>%
      arrange(desc(Total_cells))
    
    # Save summary to CSV
    write_csv(cell_summary, paste0(source_dir, "cell_count_summary.csv"))
    
    cat("Cell counts saved to:", source_dir, "\n")
    next
  
  
  # Process Mesmer data from qsave file
  cat("Calculating cell counts from loaded data...\n")
  
  # Group by Staining_condition and FOV
  cell_counts <- mesmer_data %>%
    group_by(Staining_condition, FOV) %>%
    summarise(
      Cell_count = n(),
      .groups = "drop"
    ) %>%
    arrange(Staining_condition, FOV)
  
  # Add proportion information
  cell_counts <- cell_counts %>%
    group_by(Staining_condition) %>%
    mutate(
      Total_cells = sum(Cell_count),
      Cell_proportion = Cell_count / Total_cells
    ) %>%
    ungroup()
  
  # Generate summary statistics by condition
  cell_summary <- cell_counts %>%
    group_by(Staining_condition) %>%
    summarise(
      Total_cells = first(Total_cells),
      FOV_count = n(),
      Min_cells = min(Cell_count),
      Max_cells = max(Cell_count),
      Mean_cells = mean(Cell_count),
      Median_cells = median(Cell_count)
    ) %>%
    ungroup() %>%
    arrange(desc(Total_cells))
  
  # Save to CSV files
  write_csv(cell_counts, paste0(source_dir, "fov_cell_counts.csv"))
  write_csv(cell_summary, paste0(source_dir, "cell_count_summary.csv"))
  
  cat("Cell counts saved to:", source_dir, "\n")
}

cat("\nCell count processing complete!\n")