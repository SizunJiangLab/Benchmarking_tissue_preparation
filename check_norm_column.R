# Load necessary libraries for data manipulation and reading files.
# Ensure these are installed by running: install.packages(c("dplyr", "readr", "purrr"))
library(dplyr)
library(readr)
library(purrr)

# --- Configuration Setup ---
# This section mirrors your new, expanded configuration script to ensure accuracy.
# It defines all necessary file paths and metadata to build the final configurations list.

cat("--- Setting up configurations based on your script ---\n")

# Define paths for metadata and other configuration files
metadata_file <- "./data_mesmer/Slide_metadata.csv"
pairs_file <- "./data_mesmer/Slide_compare_pairs.csv"

# Check if essential files exist before proceeding
if (!file.exists(metadata_file) || !file.exists(pairs_file)) {
  stop("Essential metadata file not found. Please ensure '", metadata_file, "' and '", pairs_file, "' are in the correct directory.")
}

# Load metadata and comparison pairs
slide_metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)
comparison_pairs <- readr::read_csv(pairs_file, show_col_types = FALSE)

# --- Section 1: Filter Source Data from Metadata ---
# This mirrors the data filtering section of your script.
astar_comet_crc_data <- slide_metadata %>% filter(Source == "ASTAR_COMET_CRC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
astar_comet_tonsil_data <- slide_metadata %>% filter(Source == "ASTAR_COMET_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_data <- slide_metadata %>% filter(Source == "BIDMC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_dlbcl_data <- slide_metadata %>% filter(Source == "BIDMC_DLBCL" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_tonsil_data <- slide_metadata %>% filter(Source == "BIDMC_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
double_data <- slide_metadata %>% filter(Source == "Double" & Type == "dataScaleSize") %>% arrange(Name, FOV)
roche_data <- slide_metadata %>% filter(Source == "Roche" & Type == "dataScaleSize") %>% arrange(Name, FOV)
roche_intestine_data <- slide_metadata %>% filter(Source == "Roche_intestine" & Type == "dataScaleSize") %>% arrange(Name, FOV)
roche_tonsil_data <- slide_metadata %>% filter(Source == "Roche_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_data <- slide_metadata %>% filter(Source == "Stanford" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_colon_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_Colon" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_liver_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_Liver" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_lymphnode_tile1_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_LymphNode_Tile1" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_lymphnode_tile2_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_LymphNode_Tile2" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_lymphnode_tile3_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_LymphNode_Tile3" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_mibi_lymphnode_tile4_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_LymphNode_Tile4" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_oscc_data <- slide_metadata %>% filter(Source == "Stanford_OSCC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_rarecyte_ln_data <- slide_metadata %>% filter(Source == "Stanford_RareCyte_LN" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_rarecyte_tma_data <- slide_metadata %>% filter(Source == "Stanford_RareCyte_TMA" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_scan1_data <- slide_metadata %>% filter(Source == "Stanford-scan1" & Type == "dataScaleSize") %>% arrange(Name, FOV)
stanford_tonsil_data <- slide_metadata %>% filter(Source == "Stanford_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
ukentucky_scc_data <- slide_metadata %>% filter(Source == "UKentucky_SCC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
ukentucky_tonsil_data <- slide_metadata %>% filter(Source == "UKentucky_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)

# --- Section 2: Assemble Final Configurations List ---
# This list is taken directly from your script to ensure the check is accurate.
# Only the parameters needed for the check (data_folder, input_filenames) are used.
# Helper function to find data folder
get_data_folder <- function(source) {
  base_data_path <- "./data_mesmer/"
  possible_paths <- c(
    paste0(base_data_path, source, "/"),
    paste0(base_data_path, "Initial_Optimization/", source, "/"),
    paste0(base_data_path, "Validation/", source, "/")
  )
  
  for (path in possible_paths) {
    if (dir.exists(path)) return(path)
  }
  return(paste0(base_data_path, source, "/")) # Default fallback
}

configurations <- list(
  ASTAR_COMET_CRC_all = list(data_folder = get_data_folder("ASTAR_COMET_CRC"), input_filenames = astar_comet_crc_data$Filename),
  ASTAR_COMET_Tonsil_all = list(data_folder = get_data_folder("ASTAR_COMET_Tonsil"), input_filenames = astar_comet_tonsil_data$Filename),
  BIDMC_all = list(data_folder = get_data_folder("BIDMC"), input_filenames = bidmc_data$Filename),
  BIDMC_DLBCL_all = list(data_folder = get_data_folder("BIDMC_DLBCL"), input_filenames = bidmc_dlbcl_data$Filename),
  BIDMC_Tonsil_all = list(data_folder = get_data_folder("BIDMC_Tonsil"), input_filenames = bidmc_tonsil_data$Filename),
  Double_all = list(data_folder = get_data_folder("Double"), input_filenames = double_data$Filename),
  Roche_all = list(data_folder = get_data_folder("Roche"), input_filenames = roche_data$Filename),
  Roche_intestine_all = list(data_folder = get_data_folder("Roche_intestine"), input_filenames = roche_intestine_data$Filename),
  Roche_Tonsil_all = list(data_folder = get_data_folder("Roche_Tonsil"), input_filenames = roche_tonsil_data$Filename),
  Stanford_all = list(data_folder = get_data_folder("Stanford"), input_filenames = stanford_data$Filename),
  Stanford_MIBI_Colon_all = list(data_folder = get_data_folder("Stanford_MIBI_Colon"), input_filenames = stanford_mibi_colon_data$Filename),
  Stanford_MIBI_Liver_all = list(data_folder = get_data_folder("Stanford_MIBI_Liver"), input_filenames = stanford_mibi_liver_data$Filename),
  Stanford_MIBI_LymphNode_Tile1_all = list(data_folder = get_data_folder("Stanford_MIBI_LymphNode_Tile1"), input_filenames = stanford_mibi_lymphnode_tile1_data$Filename),
  Stanford_MIBI_LymphNode_Tile2_all = list(data_folder = get_data_folder("Stanford_MIBI_LymphNode_Tile2"), input_filenames = stanford_mibi_lymphnode_tile2_data$Filename),
  Stanford_MIBI_LymphNode_Tile3_all = list(data_folder = get_data_folder("Stanford_MIBI_LymphNode_Tile3"), input_filenames = stanford_mibi_lymphnode_tile3_data$Filename),
  Stanford_MIBI_LymphNode_Tile4_all = list(data_folder = get_data_folder("Stanford_MIBI_LymphNode_Tile4"), input_filenames = stanford_mibi_lymphnode_tile4_data$Filename),
  Stanford_OSCC_all = list(data_folder = get_data_folder("Stanford_OSCC"), input_filenames = stanford_oscc_data$Filename),
  Stanford_RareCyte_LN_all = list(data_folder = get_data_folder("Stanford_RareCyte_LN"), input_filenames = stanford_rarecyte_ln_data$Filename),
  Stanford_RareCyte_TMA_all = list(data_folder = get_data_folder("Stanford_RareCyte_TMA"), input_filenames = stanford_rarecyte_tma_data$Filename),
  Stanford_scan1_all = list(data_folder = get_data_folder("Stanford-scan1"), input_filenames = stanford_scan1_data$Filename),
  Stanford_Tonsil_all = list(data_folder = get_data_folder("Stanford_Tonsil"), input_filenames = stanford_tonsil_data$Filename),
  UKentucky_SCC_all = list(data_folder = get_data_folder("UKentucky_SCC"), input_filenames = ukentucky_scc_data$Filename),
  UKentucky_Tonsil_all = list(data_folder = get_data_folder("UKentucky_Tonsil"), input_filenames = ukentucky_tonsil_data$Filename)
)

cat("--- Verifying Original Normalization Marker for All Configurations ---\n\n")

# --- Verification Loop ---
# This loop iterates through each configuration, reads the header of the first file,
# and reports the original name of the normalization marker based on your priority list.

for (config_name in names(configurations)) {
  current_config <- configurations[[config_name]]
  
  if (length(current_config$input_filenames) == 0) {
    cat("Configuration '", config_name, "': SKIPPED (No input files).\n", sep = "")
    next
  }
  
  first_file <- current_config$input_filenames[1]
  file_path <- file.path(current_config$data_folder, first_file)
  
  if (!file.exists(file_path)) {
    cat("Configuration '", config_name, "': FAILED (File not found: '", file_path, "').\n", sep = "")
    next
  }
  
  header_info <- tryCatch({
    readr::read_csv(file_path, n_max = 1, show_col_types = FALSE)
  }, error = function(e) { NULL })
  
  if (is.null(header_info) || ncol(header_info) == 0) {
    cat("Configuration '", config_name, "': FAILED (Could not read or parse '", file_path, "').\n", sep = "")
    next
  }
  
  # Check for the normalization marker based on the prioritized list
  cols <- colnames(header_info)
  marker_found <- ""
  
  if ("DNA1" %in% cols) {
    marker_found <- "DNA1"
  } else if ("HH3" %in% cols) {
    marker_found <- "HH3"
  } else if ("DAPI" %in% cols) {
    marker_found <- "DAPI"
  } else if ("Nucleus" %in% cols) {
    marker_found <- "Nucleus"
  } else {
    marker_found <- "None of the standard markers (DNA1, HH3, DAPI, Nucleus) were found."
  }
  
  cat("Configuration '", config_name, "': Original normalization marker is '", marker_found, "'.\n", sep = "")
}

cat("\n--- Verification Complete ---\n")
