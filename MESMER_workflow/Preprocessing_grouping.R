library(dplyr)
library(tidyverse)
library(matrixStats)
library(ggcorrplot)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)

########################################## Configuration Begin  #############################################################################

# Set data input, output directory and the metadata table
data_folder = "./drive-download-20241104T193831Z-001/"
out_folder = "./grouping_output_11052024/"
sheet <- read_csv("grouping_metadata.csv")

########################################## Configuration End  #############################################################################
dir.create(out_folder, showWarnings = F)

# Function to read and average the data for each group
average_replicates <- function(file_list) {
  print(file_list)
  # Read all files in the group and store them in a list of data frames
  data_list <- map(file_list, ~ read_csv(file.path(data_folder, .x)))
  
  colnames_data_list <- colnames(data_list[[1]])
  # Average the rest of the columns
  averaged_data <- map_dfc(5:ncol(data_list[[1]]), ~ {
    rowMeans(do.call(cbind, lapply(data_list, function(df) df[[.x]])), na.rm = TRUE)
  })
  
  # Extract the first four columns and handle them separately
  cellLabel <- rep(NA, nrow(averaged_data))
  Y_cent <- rep(NA, nrow(averaged_data))
  X_cent <- rep(NA, nrow(averaged_data))
  cellSize <- rep(NA, nrow(averaged_data))
  
  # Combine everything into one data frame
  result <- data.frame(
    cellLabel = cellLabel,
    Y_cent = Y_cent,
    X_cent = X_cent,
    cellSize = cellSize,
    averaged_data
  )
  
  colnames(result) <- colnames_data_list
  return(result)
}

result_list <- sheet %>%
  group_by(group) %>%
  summarise(
    combined_data = list(average_replicates(file_list = filename)),
    .groups = 'drop'
  )


# Save each combined data frame to a CSV file
walk2(
  result_list$group, result_list$combined_data,
  ~ write_csv(.y, file.path(out_folder, paste0(.x, "_combined.csv")))
)
