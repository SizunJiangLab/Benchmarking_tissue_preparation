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

# Create results folder if it doesn't exist
dir.create("./results", showWarnings = FALSE, recursive = TRUE)

########################################## Configuration Begin #############################################################################
### What data to load?
data_type <- "Mesmer"

# Define paths for metadata and exclusion files
metadata_file <- "./data_mesmer/Slide_metadata.csv"
removal_file <- "./data_mesmer/Slide_remove_markers.csv"
exclusion_file <- "./data_mesmer/Slide_exclude_markers.csv"
pairs_file <- "./data_mesmer/Slide_compare_pairs.csv"
marker_sequence_file <- "./data_mesmer/Registered_Report_marker_sequence.csv"


# Load all metadata and configuration files
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)
excluded_markers <- read_csv(exclusion_file)
comparison_pairs <- read_csv(pairs_file)
marker_sequence <- read_csv(marker_sequence_file, col_names = TRUE)
marker_sequence <- marker_sequence[[1]]

# Load cell-specific exclusions if the file exists
cell_exclusions <- NULL

### Define penalty scores for each source ###
### Define penalty scores for each source ###
penalty_scores <- list(
  "ASTAR_COMET_CRC" = 10,
  "ASTAR_COMET_Tonsil" = 10,
  "BIDMC" = 30,
  "BIDMC_DLBCL" = 10,
  "BIDMC_Tonsil" = 10,
  "Novartis_Lung_Cancer" = 10,
  "Novartis_Tonsil" = 10,
  "Roche" = 10,
  "Roche_intestine" = 10,
  "Roche_Tonsil" = 10,
  "Stanford" = 18,
  "Stanford_MIBI_Colon" = 10,
  "Stanford_MIBI_Liver" = 10,
  "Stanford_MIBI_LymphNode" = 10,
  "Stanford_OSCC" = 10,
  "Stanford_Orion_LN" = 10,
  "Stanford_Orion_EndometrialCancer" = 10,
  "Stanford_Tonsil" = 10,
  "UKentucky_SCC" = 10,
  "UKentucky_Tonsil" = 10,
  "StorageConditionsExpt" = 10,
  "Stanford_IMC_Tonsil" = 10,
  "Stanford_IMC_OSCC" = 10,
  "LyophilizationTest_FigS2" = 10,
  "Reimagedslide_FigS5" = 10
)

### Helper functions ###
to_pairs <- function(df) {
  select(df, Compare1, Compare2) %>% 
    pmap(function(Compare1, Compare2) c(Compare1, Compare2))
}

get_remove_markers <- function(source_name) {
  removed_markers %>% 
    filter(Source == source_name & Exclude_type == "Marker") %>% 
    pull(Exclude_value)
}

### Data-driven configuration generation ###
# Get all unique sources from metadata
all_sources <- slide_metadata %>%
  filter(Type == "dataScaleSize") %>%
  pull(Source) %>%
  unique()

# Define source-to-config-name mapping
# Define source-to-config-name mapping
source_to_config_name <- c(
  "ASTAR_COMET_CRC" = "ASTAR_COMET_CRC_all",
  "ASTAR_COMET_Tonsil" = "ASTAR_COMET_Tonsil_all",
  "BIDMC" = "BIDMC_all",
  "BIDMC_DLBCL" = "BIDMC_DLBCL_all",
  "BIDMC_Tonsil" = "BIDMC_Tonsil_all",
  "Novartis_Lung_Cancer" = "Novartis_Lung_Cancer_all",
  "Novartis_Tonsil" = "Novartis_Tonsil_all",
  "Roche" = "Roche_all",
  "Roche_intestine" = "Roche_intestine_all",
  "Roche_Tonsil" = "Roche_Tonsil_all",
  "Stanford" = "Stanford_all",
  "Stanford_MIBI_Colon" = "Stanford_MIBI_Colon_all",
  "Stanford_MIBI_Liver" = "Stanford_MIBI_Liver_all",
  "Stanford_MIBI_LymphNode" = "Stanford_MIBI_LymphNode_pooled_all",
  "Stanford_OSCC" = "Stanford_OSCC_all",
  "Stanford_Orion_LN" = "Stanford_Orion_LN_all",
  "Stanford_Orion_EndometrialCancer" = "Stanford_Orion_EndometrialCancer_all",
  "Stanford_Tonsil" = "Stanford_Tonsil_all",
  "UKentucky_SCC" = "UKentucky_SCC_all",
  "UKentucky_Tonsil" = "UKentucky_Tonsil_all",
  "StorageConditionsExpt" = "StorageConditionsExpt_all",
  "Stanford_IMC_Tonsil" = "Stanford_IMC_Tonsil_all",
  "Stanford_IMC_OSCC" = "Stanford_IMC_OSCC_all",
  "LyophilizationTest_FigS2" = "LyophilizationTest_FigS2_all",
  "Reimagedslide_FigS5" = "Reimagedslide_FigS5_all"
)

# Generate configurations dynamically
configurations <- list()

for (source in all_sources) {
  # Get config name
  config_name <- source_to_config_name[source]
  if (is.na(config_name)) {
    warning(paste("No config name mapping found for source:", source))
    next
  }
  
  # Filter metadata for this source
  source_data <- slide_metadata %>%
    filter(Source == !!source & Type == "dataScaleSize") %>%
    arrange(Name, FOV)
  
  # Get comparison pairs
  source_pairs <- comparison_pairs %>%
    filter(Source == !!source) %>%
    to_pairs()
  
  # Get markers to remove
  source_remove_markers <- get_remove_markers(source)
  
  # Get excluded markers
  source_excluded_values <- process_excluded_markers(source)
  
  # Get penalty score
  source_penalty <- penalty_scores[[source]]
  if (is.null(source_penalty)) {
    source_penalty <- 10  # default
  }
  
  # Build data folder path - checks for nested structure
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
  
  # Build output folder path
  out_folder <- paste0("./results/out_", gsub("_all$", "", config_name), "_all/")
  
  # Create configuration
  configurations[[config_name]] <- list(
    data_folder = data_folder,
    out_folder = out_folder,
    input_filenames = source_data$Filename,
    input_note = source_data$Name,
    pairs = source_pairs,
    remove_values = source_remove_markers,
    excluded_values = source_excluded_values,
    penalty_score = source_penalty
  )
}

########################################## Configuration Selection #########################################################

# Choose which configuration to use - can be changed to process different datasets
# Options: "ASTAR_COMET_CRC_all", "ASTAR_COMET_Tonsil_all",
#          "BIDMC_all", "BIDMC_DLBCL_all", "BIDMC_Tonsil_all",
#          "Novartis_Lung_Cancer_all", "Novartis_Tonsil_all", "Roche_all", 
#          "Roche_intestine_all", "Roche_Tonsil_all",
#          "Stanford_all", "Stanford_MIBI_Colon_all", "Stanford_MIBI_Liver_all",
#          "Stanford_OSCC_all", "Stanford_Orion_LN_all", "Stanford_Orion_EndometrialCancer_all",
#          "Stanford_Tonsil_all", "UKentucky_SCC_all", "UKentucky_Tonsil_all", 
#          "Stanford_MIBI_LymphNode_pooled_all", "StorageConditionsExpt_all",
#          "Stanford_IMC_Tonsil_all", "Stanford_IMC_OSCC_all",
#          "LyophilizationTest_FigS2_all", "Reimagedslide_FigS5_all"
current_config_name <- "LyophilizationTest_FigS2_all"

########################################## Configuration End ###############################################################

current_config <- configurations[[current_config_name]]

if (is.null(current_config)) {
  stop(paste("Configuration not found:", current_config_name))
}

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

### Load Mesmer data. Try to load from qsave if available, otherwise load raw data
# Create qsave directory if it doesn't exist
dir.create("./qsave_input", showWarnings = FALSE)

qsave_file <- paste0("./qsave_input/", current_config_name, "_input.qsave")
if (file.exists(qsave_file)) {
  cat("Loading data from saved qsave file:", qsave_file, "\n")
  data <- qs::qread(qsave_file)
} else {
  cat("Loading raw data files...\n")
  
  # --- UPDATED LOGIC TO CHOOSE LOADING FUNCTION ---
  if (current_config_name == "Stanford_MIBI_LymphNode_pooled_all") {
    cat("Using new pooled loading function for Tiles...\n")
    data <- load_mesmer_data_pooled(data_folder, input_filenames, input_note)
  } else {
    cat("Using standard FOV loading function...\n")
    data <- load_mesmer_data(data_folder, input_filenames, input_note)
  }
  
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
    "Hoechst" = any_of(c("DAPI", "Nucleus", "DNA1", "HH3"))
  ) %>%
  names()

# Remove special characters
data <- data %>%
  rename_all(~ gsub("[.-]", "", .))

# --- Logic to identify, rename, and position the primary nuclear marker ---

if ("DNA1" %in% colnames(data)) {
  # If DNA1 exists, prioritize it.
  # Rename it to "Hoechst", move it to column 5, and remove DNA2.
  print("Found 'DNA1'. Using it as the primary nuclear marker 'Hoechst' and removing 'DNA2'.")
  
  data <- data %>%
    rename(Hoechst = DNA1) %>%
    relocate(Hoechst, .before = 5) %>%
    select(-any_of("DNA2")) # Safely removes DNA2 only if it exists
  
} else if ("HH3" %in% colnames(data)) {
  # If DNA1 is not found, check for HH3.
  # Rename it to "Hoechst" and move it to column 5.
  print("DNA1 not found. Found 'HH3'. Using as primary nuclear marker 'Hoechst'.")
  
  data <- data %>%
    rename(Hoechst = HH3) %>%
    relocate(Hoechst, .before = 5)
  
} else if ("DAPI" %in% colnames(data)) {
  # If neither DNA1 nor HH3 is found, fall back to DAPI.
  print("DNA1 and HH3 not found. Using 'DAPI' as the nuclear marker 'Hoechst'.")
  
  data <- data %>%
    rename(Hoechst = DAPI)
  
} else if ("Nucleus" %in% colnames(data)) {
  # If none of the above are found, fall back to Nucleus.
  print("DNA1, HH3, and DAPI not found. Using 'Nucleus' as the nuclear marker 'Hoechst'.")
  
  data <- data %>%
    rename(Hoechst = Nucleus)
}


marker_names <- data %>%
  select(5:ncol(data)) %>%
  select(-Staining_condition) %>% # Ensure Staining_condition is handled or exists
  names()

# Ensure normalize_data is defined in helper.R
result = normalize_data(data)
df_norm <- result$data
marker_names <- result$marker_names

# --- Sync original_marker_names with marker_names (handling Hoechst rename, DNA2 removal, and reordering) ---
temp_clean_names <- gsub("[.-]", "", original_marker_names)
clean_to_original_map <- setNames(original_marker_names, temp_clean_names)

new_original <- c()
for (m in marker_names) {
  if (m == "Hoechst") {
    # Prioritize finding the original nuclear marker
    src <- NA
    # Check candidates in priority order (same as renaming logic)
    if ("DNA1" %in% original_marker_names) src <- "DNA1"
    else if ("HH3" %in% original_marker_names) src <- "HH3"
    else if ("DAPI" %in% original_marker_names) src <- "DAPI"
    else if ("Nucleus" %in% original_marker_names) src <- "Nucleus"
    
    if (is.na(src)) src <- "Hoechst" 
    new_original <- c(new_original, src)
  } else {
    orig <- clean_to_original_map[m]
    if (is.na(orig)) orig <- m 
    new_original <- c(new_original, orig)
  }
}
original_marker_names <- new_original

# Remove markers that are in the remove_values list
if (length(remove_values) > 0) {
  # Since order is now synced, we can just filter by matching against marker_names
  keep_mask <- ! (marker_names %in% remove_values)
  
  # Log removed
  cat("Removed markers:", paste(marker_names[!keep_mask], collapse = ", "), "\n")
  
  marker_names <- marker_names[keep_mask]
  original_marker_names <- original_marker_names[keep_mask]
  
  cat("Remaining markers (cleaned):", length(marker_names), "\n")
  cat("Remaining markers (original):", length(original_marker_names), "\n")
}

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
# Ensure plot_density_plots is defined in helper.R
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
  height = 11,
  device = svglite::svglite,
  fix_text_size = FALSE
)

### Perform statistical and Kruskal-Wallis tests
kruskal_results <- perform_statistical_and_kruskal_wallis_tests(df_trans, marker_names)

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, paste0(out_folder, "kruskal_pvals.csv"))

effect_size <- calc_effect_size(
  data = df_trans,
  out_folder = out_folder,
  pairs = pairs,
  remove_values = remove_values
)

# Generate heatmaps
heatmap_results <- plot_heatmaps(
  df_trans = df_trans,
  out_folder = out_folder,
  excluded_values = excluded_values,
  remove_values = remove_values,
  marker_sequence = marker_sequence,
  original_marker_names = original_marker_names,
  current_config_name = current_config_name
)

# Run the ranking and scoring analysis using files produced by plot_heatmaps
scoring_results <- implement_scoring_system(
  df_trans = df_trans,
  out_folder = out_folder,
  excluded_values = excluded_values,
  remove_values = remove_values,
  penalty_score = current_config$penalty_score,
  current_config_name = current_config_name
)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), paste0(out_folder, "session_info.txt"))

# Print completion message
cat("\nAnalysis completed for", current_config_name, "\n")
cat("Results saved to:", out_folder, "\n")