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
data_type <- "cellXpress"
# Define paths for metadata and exclusion files
metadata_file <- "./data_cellXpress/Slide_metadata.csv"
removal_file <- "./data_cellXpress/Slide_remove_markers.csv"
exclusion_file <- "./data_cellXpress/Slide_exclude_markers.csv"
pairs_file <- "./data_cellXpress/Slide_compare_pairs.csv"
marker_sequence_file <- "./data_cellXpress/Registered_Report_marker_sequence.csv"


# Load metadata and exclusions
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)
excluded_markers <- read_csv(exclusion_file)
comparison_pairs <- read_csv(pairs_file)
marker_sequence <- read_csv(marker_sequence_file, col_names = TRUE)
marker_sequence <- marker_sequence[[1]]  # Get the first column as a vector

# Load cell-specific exclusions if the file exists
cell_exclusions <- NULL

# --- Define all source names from metadata ---
all_sources <- unique(slide_metadata$Source)

# --- Dynamically create data filters, pairs, and exclusions for each source ---
data_filters <- list()
source_pairs <- list()
source_remove_markers <- list()
source_excluded_values <- list()

for (src in all_sources) {
  # Filter metadata for the current source
  data_filters[[src]] <- slide_metadata %>%
    filter(Source == src & Type == "dataScaleSize") %>%
    arrange(Name, FOV)
  
  # Extract comparison pairs
  source_pairs[[src]] <- comparison_pairs %>%
    filter(Source == src) %>%
    select(Compare1, Compare2) %>%
    pmap(function(Compare1, Compare2) c(Compare1, Compare2))
  
  # Process removal markers
  source_remove_markers[[src]] <- removed_markers %>%
    filter(Source == src & Exclude_type == "Marker") %>%
    pull(Exclude_value)
  
  # Process exclusion lists
  source_excluded_values[[src]] <- process_excluded_markers(src)
}


# --- Define source-specific penalty scores ---
# Initial Optimization datasets use source-specific penalty scores
# Validation datasets use standard penalty score of 10
penalty_scores <- list(
  # Initial Optimization sources
  `BIDMC` = 30,
  `Roche` = 10,
  `Stanford` = 18,
  # Validation sources - ASTAR
  `ASTAR_COMET_Tonsil` = 10,
  `ASTAR_COMET_CRC` = 10,
  # Validation sources - BIDMC
  `BIDMC_Tonsil` = 10,
  `BIDMC_DLBCL` = 10,
  # Validation sources - Roche
  `Roche_Tonsil` = 10,
  `Roche_Intestine` = 10,
  # Validation sources - UKentucky
  `UKentucky_Tonsil` = 10,
  `UKentucky_Skin` = 10,
  # Validation sources - Stanford MIBI
  `Stanford_MIBI_LymphNode_pooled` = 10,
  `Stanford_MIBI_Colon` = 10,
  `Stanford_MIBI_Liver` = 10,
  # Validation sources - Stanford Orion/Rarecyte
  `Stanford_Orion_Lymph_node` = 10,
  `Stanford_Orion_Endometrium` = 10,
  # Validation sources - Novartis
  `Novartis_LungCancer` = 10,
  `Novartis_Tonsil` = 10,
  # Validation sources - Stanford IMC
  `Stanford_IMC_OSCC` = 10,
  `Stanford_IMC_Tonsil` = 10
)

# --- Map sources to their data folders ---
# Folder structure follows README.md:
# - Initial_Optimization/BIDMC/, Initial_Optimization/Roche/, Initial_Optimization/Stanford/
# - Validation/ASTAR/, Validation/BIDMC/, Validation/Roche/, Validation/UK/
# - Validation/Stanford_MIBI/, Validation/Stanford_Orion/
# - Validation/Novartis_LungCancer/, Validation/Novartis_Tonsil/
# - Validation/Stanford_IMC_OSCC/, Validation/Stanford_IMC_Tonsil/
data_folder_map <- list(
  # Initial Optimization sources
  `BIDMC` = "./data_cellXpress/Initial_Optimization/BIDMC/",
  `Roche` = "./data_cellXpress/Initial_Optimization/Roche/",
  `Stanford` = "./data_cellXpress/Initial_Optimization/Stanford/",
  # Validation sources - ASTAR (both tissue types in same folder)
  `ASTAR_COMET_Tonsil` = "./data_cellXpress/Validation/ASTAR/",
  `ASTAR_COMET_CRC` = "./data_cellXpress/Validation/ASTAR/",
  # Validation sources - BIDMC (both tissue types in same folder)
  `BIDMC_Tonsil` = "./data_cellXpress/Validation/BIDMC/",
  `BIDMC_DLBCL` = "./data_cellXpress/Validation/BIDMC/",
  # Validation sources - Roche (both tissue types in same folder)
  `Roche_Tonsil` = "./data_cellXpress/Validation/Roche/",
  `Roche_Intestine` = "./data_cellXpress/Validation/Roche/",
  # Validation sources - UKentucky (both tissue types in same folder)
  `UKentucky_Tonsil` = "./data_cellXpress/Validation/UK/",
  `UKentucky_Skin` = "./data_cellXpress/Validation/UK/",
  # Validation sources - Stanford MIBI (all tissue types in same folder)
  `Stanford_MIBI_LymphNode_pooled` = "./data_cellXpress/Validation/Stanford_MIBI/",
  `Stanford_MIBI_Colon` = "./data_cellXpress/Validation/Stanford_MIBI/",
  `Stanford_MIBI_Liver` = "./data_cellXpress/Validation/Stanford_MIBI/",
  # Validation sources - Stanford Orion/Rarecyte (both tissue types in same folder)
  `Stanford_Orion_Lymph_node` = "./data_cellXpress/Validation/Stanford_Orion/",
  `Stanford_Orion_Endometrium` = "./data_cellXpress/Validation/Stanford_Orion/",
  # Validation sources - Novartis (separate folders per tissue type)
  `Novartis_LungCancer` = "./data_cellXpress/Validation/Novartis_LungCancer/",
  `Novartis_Tonsil` = "./data_cellXpress/Validation/Novartis_Tonsil/",
  # Validation sources - Stanford IMC (separate folders per tissue type)
  `Stanford_IMC_OSCC` = "./data_cellXpress/Validation/Stanford_IMC_OSCC/",
  `Stanford_IMC_Tonsil` = "./data_cellXpress/Validation/Stanford_IMC_Tonsil/"
)

# --- Create configurations for each dataset dynamically ---
configurations <- list()
for (src in all_sources) {
  config_name <- paste0(gsub(" ", "_", src), "_cellXpress")
  
  # Check if data exists for this source before creating a config
  if (nrow(data_filters[[src]]) > 0) {
    # Get the folder path, defaulting to a generated path if not in map
    folder_path <- data_folder_map[[src]]
    if (is.null(folder_path)) {
      # Generate a default path if not explicitly mapped
      folder_path <- paste0("./data_cellXpress/Validation/", src, "/")
    }
    
    # Get penalty score, defaulting to 10 if not defined
    penalty <- penalty_scores[[src]]
    if (is.null(penalty)) {
      penalty <- 10
    }
    
    configurations[[config_name]] <- list(
      data_folder = folder_path,
      out_folder = paste0("./results/out_", config_name, "/"),
      input_filenames = data_filters[[src]]$Filename,
      input_note = data_filters[[src]]$Name,
      pairs = source_pairs[[src]],
      remove_values = source_remove_markers[[src]],
      excluded_values = source_excluded_values[[src]],
      penalty_score = penalty
    )
  }
}


# --- Choose which configuration to use ---
# Available configurations (dynamically generated based on your metadata file):
#
# Initial Optimization:
#   "BIDMC_cellXpress", "Roche_cellXpress", "Stanford_cellXpress"
#
# Validation - ASTAR:
#   "ASTAR_COMET_Tonsil_cellXpress", "ASTAR_COMET_CRC_cellXpress"
#
# Validation - BIDMC:
#   "BIDMC_Tonsil_cellXpress", "BIDMC_DLBCL_cellXpress"
#
# Validation - Roche:
#   "Roche_Tonsil_cellXpress", "Roche_Intestine_cellXpress"
#
# Validation - UKentucky:
#   "UKentucky_Tonsil_cellXpress", "UKentucky_Skin_cellXpress"
#
# Validation - Stanford MIBI:
#   "Stanford_MIBI_LymphNode_pooled_cellXpress", "Stanford_MIBI_Colon_cellXpress", "Stanford_MIBI_Liver_cellXpress"
#
# Validation - Stanford Orion:
#   "Stanford_Orion_Lymph_node_cellXpress", "Stanford_Orion_Endometrium_cellXpress"
#
# Validation - Novartis:
#   "Novartis_LungCancer_cellXpress", "Novartis_Tonsil_cellXpress"
#
# Validation - Stanford IMC:
#   "Stanford_IMC_OSCC_cellXpress", "Stanford_IMC_Tonsil_cellXpress"

current_config_name <- "Stanford_IMC_Tonsil_cellXpress"
current_config <- configurations[[current_config_name]]

# Validate configuration exists
if (is.null(current_config)) {
  available_configs <- names(configurations)
  stop(paste0(
    "Configuration '", current_config_name, "' not found.\n",
    "Available configurations:\n  ",
    paste(available_configs, collapse = "\n  ")
  ))
}

########################################## Configuration End #############################################################################


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

### Load cellXpress data. Try to load from qsave if available, otherwise load raw data
dir.create("./qsave_input", showWarnings = FALSE)
qsave_file <- paste0("./qsave_input/", current_config_name, "_input.qsave")

if (file.exists(qsave_file)) {
  cat("Loading data from saved qsave file:", qsave_file, "\n")
  data <- qs::qread(qsave_file)
} else {
  cat("Loading raw data files...\n")
  if (grepl("_pooled", current_config_name)) {
    cat("Using POOLED loading function for cellXpress data...\n")
    data <- load_cellXpress_data_pooled(data_folder, input_filenames, input_note)
  } else {
    cat("Using standard loading function for cellXpress data...\n")
    data <- load_cellXpress_data(data_folder, input_filenames, input_note)
  }
  
  # Save data for future use
  cat("Saving loaded data to qsave file for faster future loading\n")
  qs::qsave(data, qsave_file)
}

# Handle different DNA/nuclear stain naming conventions
if("DAPI" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "DAPI")
} else if("Nucleus" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "Nucleus")
} else if("DNA" %in% colnames(data)) {
  data <- data %>%
    rename("Hoechst" = "DNA")
} else if("DNA1" %in% colnames(data)) {
  data <- data %>%
    rename(Hoechst = DNA1) %>%
    rename(PanCK = CK) %>%
    relocate(Hoechst, .before = 5) %>%
    select(-any_of("DNA2"))
}

original_marker_names <- data %>%
  rename_all(~ gsub("[.]", "-", .)) %>%
  select(5:ncol(data)) %>%
  select(-matches("Staining_condition")) %>%
  rename(
    "Na/K-ATPase" = any_of("NaKATPase"),
    "Granzyme B" = any_of("GranzymeB"),
    "Hoechst" = any_of(c("DAPI", "Nucleus", "DNA", "DNA1"))
  ) %>%
  names()

# Remove special characters
data <- data %>%
  rename_all(~ gsub("[.-]", "", .))

marker_names <- data %>%
  select(5:ncol(data)) %>%
  select(-Staining_condition) %>%
  names()

result <- normalize_data(data)
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
  penalty_score = current_config$penalty_score,
  current_config_name = current_config_name
)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), paste0(out_folder, "session_info.txt"))

# Print completion message
cat("\nAnalysis completed for", current_config_name, "\n")
cat("Results saved to:", out_folder, "\n")
