library(tidyverse)

# Define possible data folders
data_folders <- c("./data_mesmer/", "./data_cellXpress/")

# Loop through each folder
for (use_folder in data_folders) {
  if (!dir.exists(use_folder)) {
    cat("Folder not found (skipping):", use_folder, "\n")
    next
  }
  
  cat("Processing folder:", use_folder, "\n")
  
  # Read the metadata file
  metadata_file <- paste0(use_folder, "Slide_metadata.csv")
  
  if (!file.exists(metadata_file)) {
    cat("Metadata file not found in:", use_folder, "\n")
    next
  }
  
  slide_metadata <- read.csv(metadata_file)
  
  # Get all unique sources that have dataScaleSize type
  if (!("Source" %in% colnames(slide_metadata)) || !("Type" %in% colnames(slide_metadata))) {
    cat("Invalid metadata columns in:", metadata_file, "\n")
    next
  }

  unique_sources <- slide_metadata %>%
    dplyr::filter(Type == "dataScaleSize") %>%
    dplyr::select(Source) %>%
    distinct() %>%
    pull(Source)
  
  # Initialize an empty dataframe to store all pairs
  all_pairs <- tibble()
  
  # Loop through each source to generate comparison pairs
  for(source in unique_sources) {
    # Get unique slide names for current source
    source_names <- slide_metadata %>%
      dplyr::filter(Source == source & Type == "dataScaleSize") %>%
      dplyr::select(Name) %>%
      distinct() %>%
      pull(Name)
    
    # Generate all possible pairs for current source
    if(length(source_names) >= 2) {
      source_pairs <- combn(source_names, 2, simplify = FALSE) %>%
        map_df(~tibble(Source = source, Compare1 = .x[1], Compare2 = .x[2]))
      
      # Add to all_pairs
      all_pairs <- bind_rows(all_pairs, source_pairs)
    }
  }
  
  # Write to CSV in the chosen folder
  output_file <- paste0(use_folder, "Slide_compare_pairs.csv")
  write.csv(all_pairs, output_file, row.names = FALSE)
  cat("Wrote pairs to:", output_file, "\n")
}