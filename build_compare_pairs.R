library(tidyverse)

# Define possible data folders
data_folders <- c("./data_03272025/", "./data_cellXpress/")

# Choose the folder that exists
use_folder <- data_folders[file.exists(data_folders)][1]

if(is.na(use_folder)) {
  stop("Neither data folder exists. Please check the paths.")
}

# Read the metadata file based on chosen folder
metadata_file <- paste0(use_folder, "Slide_metadata", 
                        ifelse(use_folder == "./data_03272025/", "", ""), 
                        ".csv")

slide_metadata <- read.csv(metadata_file)

# Get all unique sources that have dataScaleSize type
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

# Preview the result
print(all_pairs)