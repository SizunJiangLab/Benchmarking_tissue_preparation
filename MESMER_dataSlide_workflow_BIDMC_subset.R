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
source("helper_BIDMC_subset.R")
########################################## Configuration Begin  #############################################################################
### What data to load?
data_type <- "MESMER"
# Define paths for metadata and exclusion files
metadata_file <- "./data_03272025/Slide_metadata_03272025.csv"
removal_file <- "./data_03272025/Slide_remove_markers_03272025.csv"
#exclusion_file <- "./data_03272025/Slide_exclude_markers_03272025.csv"
exclusion_file <- "./data_03272025/Slide_exclude_markers_04212025.csv"
pairs_file <- "./data_03272025/Slide_compare_pairs_03272025.csv"
marker_sequence_file <- "./data_03272025/Registered_Report_marker_sequence.csv"
cell_exclusion_file <- "./data_03272025/Slide_exclude_cells_03272025.csv"
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

# Create subset data for slides 5, 13, 16, and 21
bidmc_subset_data <- slide_metadata %>% 
  filter(Source == "BIDMC" & Type == "dataScaleSize") %>%
  filter(grepl("slide5_|slide13_|slide16_|slide21_", Filename)) %>%
  arrange(Name, FOV)

bidmc_subset_pairs <- comparison_pairs %>% 
  filter(Source == "BIDMC") %>%
  filter(
    (Compare1 %in% c("BIDMC_5", "BIDMC_13", "BIDMC_16", "BIDMC_21") & 
       Compare2 %in% c("BIDMC_5", "BIDMC_13", "BIDMC_16", "BIDMC_21"))
  ) %>%
  select(Compare1, Compare2) %>% 
  pmap(function(Compare1, Compare2) c(Compare1, Compare2))

# Create configurations for each dataset
configurations <- list(
  BIDMC_all = list(
    data_folder = "./data_03272025/BIDMC/",
    out_folder = "./out_BIDMC_all/",
    input_filenames = bidmc_data$Filename,
    input_note = bidmc_data$Name,
    pairs = bidmc_pairs,
    remove_values = bidmc_remove_markers,
    excluded_values = bidmc_excluded_values,
    penalty_score = penalty_scores$BIDMC
  ),
  
  Roche_all = list(
    data_folder = "./data_03272025/Roche/",
    out_folder = "./out_Roche_all/",
    input_filenames = roche_data$Filename,
    input_note = roche_data$Name,
    pairs = roche_pairs,
    remove_values = roche_remove_markers,
    excluded_values = roche_excluded_values,
    penalty_score = penalty_scores$Roche
  ),
  
  Stanford_all = list(
    data_folder = "./data_03272025/Stanford/",
    out_folder = "./out_Stanford_all/",
    input_filenames = stanford_data$Filename,
    input_note = stanford_data$Name,
    pairs = stanford_pairs,
    remove_values = stanford_remove_markers,
    excluded_values = stanford_excluded_values,
    penalty_score = penalty_scores$Stanford
  ),
  
  Stanford_scan1_all = list(
    data_folder = "./data_03272025/Stanford-scan1/",
    out_folder = "./out_Stanford_scan1_all/",
    input_filenames = stanford_scan1_data$Filename,
    input_note = stanford_scan1_data$Name,
    pairs = stanford_scan1_pairs,
    remove_values = stanford_scan1_remove_markers,
    excluded_values = stanford_scan1_excluded_values,
    penalty_score = penalty_scores[["Stanford-scan1"]]
  ),
  BIDMC_subset = list(
    data_folder = "./data_03272025/BIDMC/",
    out_folder = "./out_BIDMC_subset/",
    input_filenames = bidmc_subset_data$Filename,
    input_note = bidmc_subset_data$Name,
    pairs = bidmc_subset_pairs,  # Using the same comparison pairs
    remove_values = bidmc_remove_markers,
    excluded_values = bidmc_excluded_values,
    penalty_score = 10  # Setting penalty score to 10
  )
)

# Choose which configuration to use - can be changed to process different datasets
# Options: "BIDMC_all", "Roche_all", "Stanford_all", "Stanford_scan1_all", "BIDMC_subset
current_config_name <- "BIDMC_subset"
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

### Load MESMER data. Try to load from qsave if available, otherwise load raw data
qsave_file <- paste0(current_config_name, "_input.qsave")
if (file.exists(qsave_file)) {
  cat("Loading data from saved qsave file:", qsave_file, "\n")
  data <- qs::qread(qsave_file)
} else {
  cat("Loading raw data files...\n")
  data <- load_mesmer_data(data_folder, input_filenames, input_note)
  
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
    "Hoechst" = any_of(c("DAPI", "Nucleus"))
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

# After loading data and before visualization steps, add this code:
if (current_config_name == "BIDMC_subset") {
  # Apply custom ordering for BIDMC subset
  df_trans <- reorder_bidmc_subset_conditions(df_trans)
}

# Before calling plot_density_plots function, make sure data is ordered
if (current_config_name == "BIDMC_subset") {
  # Apply ordering to ensure density plots show the correct order
  p <- plot_density_plots(
    reorder_bidmc_subset_conditions(df_trans),
    marker_names,
    original_marker_names,
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
    legend.rows = 6
  )
} else {
  # Original code for other configurations
  p <- plot_density_plots(
    df_trans,
    marker_names,
    original_marker_names,
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
    legend.rows = 6
  )
}

# Before calling plot_heatmaps function, add the custom ordering
if (current_config_name == "BIDMC_subset") {
  # Apply custom ordering for heatmap generation
  heatmap_results <- plot_heatmaps(
    df_trans = reorder_bidmc_subset_conditions(df_trans),
    out_folder = out_folder,
    excluded_values = excluded_values,
    remove_values = remove_values,
    marker_sequence = marker_sequence,
    original_marker_names = original_marker_names
  )
} else {
  # Original code for other configurations
  heatmap_results <- plot_heatmaps(
    df_trans = df_trans,
    out_folder = out_folder,
    excluded_values = excluded_values,
    remove_values = remove_values,
    marker_sequence = marker_sequence,
    original_marker_names = original_marker_names
  )
}

# Before calling implement_scoring_system, add the custom ordering
if (current_config_name == "BIDMC_subset") {
  # Apply custom ordering for scoring system
  scoring_results <- implement_scoring_system(
    df_trans = reorder_bidmc_subset_conditions(df_trans),
    out_folder = out_folder,
    excluded_values = excluded_values,
    remove_values = remove_values,
    penalty_score = current_config$penalty_score
  )
} else {
  # Original code for other configurations
  scoring_results <- implement_scoring_system(
    df_trans = df_trans,
    out_folder = out_folder,
    excluded_values = excluded_values,
    remove_values = remove_values,
    penalty_score = current_config$penalty_score
  )
}

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
