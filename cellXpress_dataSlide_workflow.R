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
library(qs)
source("helper.R")


########################################## Configuration Begin  #############################################################################
### What data to load?
data_type <- "cellXpress"
# Define paths for metadata and exclusion files
metadata_file <- "./data_cellXpress/Slide_metadata.csv"
removal_file <- "./data_cellXpress/Slide_remove_markers.csv"
exclusion_file <- "./data_cellXpress/Slide_exclude_markers.csv"
pairs_file <- "./data_cellXpress/Slide_compare_pairs.csv"
marker_sequence_file <- "./data_cellXpress/Registered_Report_marker_sequence.csv"
cell_exclusion_file <- "./data_cellXpress/Slide_exclude_cells.csv"
# Load metadata and exclusions
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)
excluded_markers <- read_csv(exclusion_file)
comparison_pairs <- read_csv(pairs_file)
marker_sequence <- read_csv(marker_sequence_file, col_names = TRUE)
marker_sequence <- marker_sequence[[1]]  # Get the first column as a vector
# Load cell-specific exclusions if the file exists
cell_exclusions <- NULL
if (file.exists(cell_exclusion_file)) {
  cell_exclusions <- read_csv(cell_exclusion_file) %>%
    mutate(
      # Clean marker names
      Marker = gsub("[.-]", "", Marker),
      # Create a key for matching
      Key = paste0(Source, "_", Slide)
    )
}
# Filter metadata by source and FOV to create different datasets
bidmc_data <- slide_metadata %>% 
  filter(Source == "BIDMC" & Type == "dataScaleSize") %>%
  arrange(Name, FOV)
roche_data <- slide_metadata %>% 
  filter(Source == "Roche" & Type == "dataScaleSize") %>%
  arrange(Name, FOV)
stanford_data <- slide_metadata %>% 
  filter(Source == "Stanford" & Type == "dataScaleSize") %>%
  arrange(Name, FOV)
# Add Stanford-scan1 data
stanford_scan1_data <- slide_metadata %>% 
  filter(Source == "Stanford-scan1" & Type == "dataScaleSize") %>%
  arrange(Name, FOV)

# Extract comparison pairs for each source
bidmc_pairs <- comparison_pairs %>% 
  filter(Source == "BIDMC") %>% 
  select(Compare1, Compare2) %>% 
  pmap(function(Compare1, Compare2) c(Compare1, Compare2))
roche_pairs <- comparison_pairs %>% 
  filter(Source == "Roche") %>% 
  select(Compare1, Compare2) %>% 
  pmap(function(Compare1, Compare2) c(Compare1, Compare2))
stanford_pairs <- comparison_pairs %>% 
  filter(Source == "Stanford") %>% 
  select(Compare1, Compare2) %>% 
  pmap(function(Compare1, Compare2) c(Compare1, Compare2))
# Add Stanford-scan1 pairs
stanford_scan1_pairs <- comparison_pairs %>% 
  filter(Source == "Stanford-scan1") %>% 
  select(Compare1, Compare2) %>% 
  pmap(function(Compare1, Compare2) c(Compare1, Compare2))

# Process removal and exclusion markers
roche_remove_markers <- removed_markers %>% 
  filter(Source == "Roche" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
bidmc_remove_markers <- removed_markers %>% 
  filter(Source == "BIDMC" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
stanford_remove_markers <- removed_markers %>% 
  filter(Source == "Stanford" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
# Add Stanford-scan1 removal markers
stanford_scan1_remove_markers <- removed_markers %>% 
  filter(Source == "Stanford-scan1" & Exclude_type == "Marker") %>%
  pull(Exclude_value)
# If no specific markers are defined for Stanford-scan1, use Stanford's markers
if(length(stanford_scan1_remove_markers) == 0) {
  stanford_scan1_remove_markers <- stanford_remove_markers
}

# Build exclusion lists for each source
bidmc_excluded_values <- process_excluded_markers("BIDMC")
roche_excluded_values <- list()
stanford_excluded_values <- process_excluded_markers("Stanford")
# Add Stanford-scan1 excluded values
stanford_scan1_excluded_values <- process_excluded_markers("Stanford-scan1")
# If no specific exclusions are defined for Stanford-scan1, use Stanford's exclusions
if(length(stanford_scan1_excluded_values) == 0) {
  stanford_scan1_excluded_values <- stanford_excluded_values
}

# Define source-specific penalty scores
penalty_scores <- list(
  BIDMC = 30,
  Roche = 10,
  Stanford = 18,
  `Stanford-scan1` = 18  # Using the same penalty score as Stanford
)

# Create configurations for each dataset
configurations <- list(
  BIDMC_cellXpress_all = list(
    data_folder = "./data_cellXpress/BIDMC/",
    out_folder = "./out_cellXpress_BIDMC_all/",
    input_filenames = bidmc_data$Filename,
    input_note = bidmc_data$Name,
    pairs = bidmc_pairs,
    remove_values = bidmc_remove_markers,
    excluded_values = bidmc_excluded_values,
    penalty_score = penalty_scores$BIDMC
  ),
  
  Roche_cellXpress_all = list(
    data_folder = "./data_cellXpress/Roche/",
    out_folder = "./out_cellXpress_Roche_all/",
    input_filenames = roche_data$Filename,
    input_note = roche_data$Name,
    pairs = roche_pairs,
    remove_values = roche_remove_markers,
    excluded_values = roche_excluded_values,
    penalty_score = penalty_scores$Roche
  ),
  
  Stanford_cellXpress_all = list(
    data_folder = "./data_cellXpress/Stanford/",
    out_folder = "./out_cellXpress_Stanford_all/",
    input_filenames = stanford_data$Filename,
    input_note = stanford_data$Name,
    pairs = stanford_pairs,
    remove_values = stanford_remove_markers,
    excluded_values = stanford_excluded_values,
    penalty_score = penalty_scores$Stanford
  )
)

# Choose which configuration to use - can be changed to process different datasets
# Options: "BIDMC_cellXpress_all", "Roche_cellXpress_all", "Stanford_cellXpress_all"
current_config_name <- "Stanford_cellXpress_all"
current_config <- configurations[[current_config_name]]

########################################## Configuration End  #############################################################################

# Set up working environment
data_folder <- current_config$data_folder
out_folder <- current_config$out_folder
input_filenames <- current_config$input_filenames
input_note <- current_config$input_note
pairs <- current_config$pairs
excluded_values <- current_config$excluded_values
remove_values <- current_config$remove_values

dir.create(out_folder, showWarnings = FALSE)

# Save configuration details for reference
write_csv(
  tibble(
    config_name = current_config_name,
    data_folder = data_folder,
    out_folder = out_folder,
    file_count = length(input_filenames)
  ),
  paste0(out_folder, "config_summary.csv")
)

# Save list of processed files
write_csv(
  tibble(
    filename = input_filenames,
    sample_name = input_note
  ),
  paste0(out_folder, "processed_files.csv")
)

### Load CellXpress data. Try to load from qsave if available, otherwise load raw data
qsave_file <- paste0(current_config_name, "_input.qsave")
if (file.exists(qsave_file)) {
  cat("Loading data from saved qsave file:", qsave_file, "\n")
  data <- qs::qread(qsave_file)
} else {
  cat("Loading raw data files...\n")
  data <- load_cellXpress_data(data_folder, input_filenames, input_note)
  # Save data for future use
  cat("Saving loaded data to qsave file for faster future loading\n")
  qs::qsave(data, qsave_file)
}

original_marker_names <- data %>%
  rename_all(~ gsub("[.]", "-", .)) %>%
  select(5:ncol(data)) %>%
  select(-matches("Staining_condition")) %>%
  rename(
    "Na/K-ATPase" = any_of("NaKATPase"),
    "Granzyme B" = any_of("GranzymeB"),
    "Hoechst" = any_of(c("DAPI", "Nucleus","DNA"))
  ) %>%
  names()

# Remove special characters
data <- data %>%
  rename_all(~ gsub("[.-]", "", .)) 

if("DAPI" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "DAPI")
} else if("Nucleus" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "Nucleus")
} else if("DNA" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "DNA")
}



marker_names <- data %>%
  select(5:ncol(data)) %>%
  select(-Staining_condition) %>%
  names()

result = normalize_data(data)
df_norm <- result$data
marker_names <- result$marker_names

### Arcsinh (Inverse hyperbolic Sine) Transformation
df_arcsinh <- df_norm %>%
  mutate(across(all_of(marker_names), ~ asinh(.x / 0.001)))

# Universal percentile normalization
df_trans <- df_arcsinh
rng <- colQuantiles(as.matrix(df_arcsinh[, marker_names]), probs = c(0.001, 0.999))
expr <- t((t(as.matrix(df_arcsinh[, marker_names])) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr[expr < 0] <- 0
expr[expr > 1] <- 1
df_trans[, marker_names] <- expr

### Plot density plots for each staining condition and marker
p <- plot_density_plots(
  df_trans,
  marker_names,
  original_marker_names,
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.justification = "left",
  legend.rows = 6
)

ggsave(
  p,
  filename = paste0(out_folder, "Arcsinh_transformed_Hoechst_normalised_density_plots.svg"),
  width = 8,
  height = 11
)

### Perform statistical and Kruskal-Wallis tests
kruskal_results <- perform_statistical_and_kruskal_wallis_tests(df_trans, marker_names)

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, paste0(out_folder, "kruskal_pvals.csv"))


effect_size <- calc_effect_size(
  data = df_trans,
  out_folder = out_folder,
  pairs = pairs
)

# Generate heatmaps
heatmap_results <- plot_heatmaps(
  df_trans = df_trans,
  out_folder = out_folder,
  excluded_values = excluded_values,
  remove_values = remove_values,
  marker_sequence = marker_sequence,
  original_marker_names = original_marker_names
)

# Run the ranking and scoring analysis using files produced by plot_heatmaps
scoring_results <- implement_scoring_system(
  df_trans = df_trans,
  out_folder = out_folder,
  excluded_values = excluded_values,
  remove_values = remove_values,
  penalty_score = current_config$penalty_score  # Use source-specific penalty score
)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), paste0(out_folder, "session_info.txt"))

# Print completion message
cat("\nAnalysis completed for", current_config_name, "\n")
cat("Results saved to:", out_folder, "\n")