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
########################################## Configuration Begin #############################################################################
### What data to load?
data_type <- "MESMER"
# Define paths for metadata and exclusion files
metadata_file <- "./data_03272025/Slide_metadata.csv"
removal_file <- "./data_03272025/Slide_remove_markers.csv"
exclusion_file <- "./data_03272025/Slide_exclude_markers.csv"
pairs_file <- "./data_03272025/Slide_compare_pairs.csv"
marker_sequence_file <- "./data_03272025/Registered_Report_marker_sequence.csv"
cell_exclusion_file <- "./data_03272025/Slide_exclude_cells.csv"

# Load all metadata and configuration files
slide_metadata <- read_csv(metadata_file)
removed_markers <- read_csv(removal_file)
excluded_markers <- read_csv(exclusion_file)
comparison_pairs <- read_csv(pairs_file)
marker_sequence <- read_csv(marker_sequence_file, col_names = TRUE)
marker_sequence <- marker_sequence[[1]]

# Load cell-specific exclusions if the file exists
cell_exclusions <- NULL
if (file.exists(cell_exclusion_file)) {
  cell_exclusions <- read_csv(cell_exclusion_file) %>%
    mutate(
      Marker = gsub("[.-]", "", Marker),
      Key = paste0(Source, "_", Slide)
    )
}

### Define settings for each independent configuration ###

### --- Section 1: Filter Source Data from Metadata --- ###
adjusted_rarecyte_ln_data <- slide_metadata %>% filter(Source == "RareCyte_LN1_FOV1_FOV2_0.325" & Type == "dataScaleSize") %>% arrange(Name, FOV)
adjusted_rarecyte_tma_data <- slide_metadata %>% filter(Source == "RareCyte_TMA_FOV1_FOV2_0.325" & Type == "dataScaleSize") %>% arrange(Name, FOV)
astar_comet_crc_data <- slide_metadata %>% filter(Source == "ASTAR_COMET_CRC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
astar_comet_tonsil_data <- slide_metadata %>% filter(Source == "ASTAR_COMET_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_data <- slide_metadata %>% filter(Source == "BIDMC" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_dlbcl_data <- slide_metadata %>% filter(Source == "BIDMC_DLBCL" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_tonsil_data <- slide_metadata %>% filter(Source == "BIDMC_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
bidmc_tonsil_compare_data <- slide_metadata %>% filter(Source == "BIDMC_Tonsil_compare" & Type == "dataScaleSize") %>% arrange(Name, FOV)
double_data <- slide_metadata %>% filter(Source == "Double" & Type == "dataScaleSize") %>% arrange(Name, FOV)
novartis_lung_cancer_data <- slide_metadata %>% filter(Source == "Novartis_Lung_Cancer" & Type == "dataScaleSize") %>% arrange(Name, FOV)
novartis_tonsil_data <- slide_metadata %>% filter(Source == "Novartis_Tonsil" & Type == "dataScaleSize") %>% arrange(Name, FOV)
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
stanford_mibi_lymphnode_pooled_data <- slide_metadata %>% filter(Source == "Stanford_MIBI_LymphNode" & Type == "dataScaleSize") %>% arrange(Name, FOV)

# --- Section 2: Define Comparison Pairs for Each Source ---
to_pairs <- function(df) {
  select(df, Compare1, Compare2) %>% pmap(function(Compare1, Compare2) c(Compare1, Compare2))
}
adjusted_rarecyte_ln_pairs <- comparison_pairs %>% filter(Source == "RareCyte_LN1_FOV1_FOV2_0.325") %>% to_pairs()
adjusted_rarecyte_tma_pairs <- comparison_pairs %>% filter(Source == "RareCyte_TMA_FOV1_FOV2_0.325") %>% to_pairs()
astar_comet_crc_pairs <- comparison_pairs %>% filter(Source == "ASTAR_COMET_CRC") %>% to_pairs()
astar_comet_tonsil_pairs <- comparison_pairs %>% filter(Source == "ASTAR_COMET_Tonsil") %>% to_pairs()
bidmc_pairs <- comparison_pairs %>% filter(Source == "BIDMC") %>% to_pairs()
bidmc_dlbcl_pairs <- comparison_pairs %>% filter(Source == "BIDMC_DLBCL") %>% to_pairs()
bidmc_tonsil_pairs <- comparison_pairs %>% filter(Source == "BIDMC_Tonsil") %>% to_pairs()
bidmc_tonsil_compare_pairs <- comparison_pairs %>% filter(Source == "BIDMC_Tonsil_compare") %>% to_pairs()
double_pairs <- comparison_pairs %>% filter(Source == "Double") %>% to_pairs()
novartis_lung_cancer_pairs <- comparison_pairs %>% filter(Source == "Novartis_Lung_Cancer") %>% to_pairs()
novartis_tonsil_pairs <- comparison_pairs %>% filter(Source == "Novartis_Tonsil") %>% to_pairs()
roche_pairs <- comparison_pairs %>% filter(Source == "Roche") %>% to_pairs()
roche_intestine_pairs <- comparison_pairs %>% filter(Source == "Roche_intestine") %>% to_pairs()
roche_tonsil_pairs <- comparison_pairs %>% filter(Source == "Roche_Tonsil") %>% to_pairs()
stanford_pairs <- comparison_pairs %>% filter(Source == "Stanford") %>% to_pairs()
stanford_mibi_colon_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_Colon") %>% to_pairs()
stanford_mibi_liver_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_Liver") %>% to_pairs()
stanford_mibi_lymphnode_tile1_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_LymphNode_Tile1") %>% to_pairs()
stanford_mibi_lymphnode_tile2_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_LymphNode_Tile2") %>% to_pairs()
stanford_mibi_lymphnode_tile3_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_LymphNode_Tile3") %>% to_pairs()
stanford_mibi_lymphnode_tile4_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_LymphNode_Tile4") %>% to_pairs()
stanford_oscc_pairs <- comparison_pairs %>% filter(Source == "Stanford_OSCC") %>% to_pairs()
stanford_rarecyte_ln_pairs <- comparison_pairs %>% filter(Source == "Stanford_RareCyte_LN") %>% to_pairs()
stanford_rarecyte_tma_pairs <- comparison_pairs %>% filter(Source == "Stanford_RareCyte_TMA") %>% to_pairs()
stanford_scan1_pairs <- comparison_pairs %>% filter(Source == "Stanford-scan1") %>% to_pairs()
stanford_tonsil_pairs <- comparison_pairs %>% filter(Source == "Stanford_Tonsil") %>% to_pairs()
ukentucky_scc_pairs <- comparison_pairs %>% filter(Source == "UKentucky_SCC") %>% to_pairs()
ukentucky_tonsil_pairs <- comparison_pairs %>% filter(Source == "UKentucky_Tonsil") %>% to_pairs()
stanford_mibi_lymphnode_pooled_pairs <- comparison_pairs %>% filter(Source == "Stanford_MIBI_LymphNode") %>% to_pairs()

# --- Section 3: Define Markers to Remove for Each Source ---
get_remove_markers <- function(source_name) {
  removed_markers %>% filter(Source == source_name & Exclude_type == "Marker") %>% pull(Exclude_value)
}
adjusted_rarecyte_ln_remove_markers <- get_remove_markers("RareCyte_LN1_FOV1_FOV2_0.325")
adjusted_rarecyte_tma_remove_markers <- get_remove_markers("RareCyte_TMA_FOV1_FOV2_0.325")
astar_comet_crc_remove_markers <- get_remove_markers("ASTAR_COMET_CRC")
astar_comet_tonsil_remove_markers <- get_remove_markers("ASTAR_COMET_Tonsil")
bidmc_remove_markers <- get_remove_markers("BIDMC")
bidmc_dlbcl_remove_markers <- get_remove_markers("BIDMC_DLBCL")
bidmc_tonsil_remove_markers <- get_remove_markers("BIDMC_Tonsil")
bidmc_tonsil_compare_remove_markers <- get_remove_markers("BIDMC_Tonsil_compare")
double_remove_markers <- get_remove_markers("Double")
novartis_lung_cancer_remove_markers <- get_remove_markers("Novartis_Lung_Cancer")
novartis_tonsil_remove_markers <- get_remove_markers("Novartis_Tonsil")
roche_remove_markers <- get_remove_markers("Roche")
roche_intestine_remove_markers <- get_remove_markers("Roche_intestine")
roche_tonsil_remove_markers <- get_remove_markers("Roche_Tonsil")
stanford_remove_markers <- get_remove_markers("Stanford")
stanford_mibi_colon_remove_markers <- get_remove_markers("Stanford_MIBI_Colon")
stanford_mibi_liver_remove_markers <- get_remove_markers("Stanford_MIBI_Liver")
stanford_mibi_lymphnode_tile1_remove_markers <- get_remove_markers("Stanford_MIBI_LymphNode_Tile1")
stanford_mibi_lymphnode_tile2_remove_markers <- get_remove_markers("Stanford_MIBI_LymphNode_Tile2")
stanford_mibi_lymphnode_tile3_remove_markers <- get_remove_markers("Stanford_MIBI_LymphNode_Tile3")
stanford_mibi_lymphnode_tile4_remove_markers <- get_remove_markers("Stanford_MIBI_LymphNode_Tile4")
stanford_oscc_remove_markers <- get_remove_markers("Stanford_OSCC")
stanford_rarecyte_ln_remove_markers <- get_remove_markers("Stanford_RareCyte_LN")
stanford_rarecyte_tma_remove_markers <- get_remove_markers("Stanford_RareCyte_TMA")
stanford_scan1_remove_markers <- get_remove_markers("Stanford-scan1")
stanford_tonsil_remove_markers <- get_remove_markers("Stanford_Tonsil")
ukentucky_scc_remove_markers <- get_remove_markers("UKentucky_SCC")
ukentucky_tonsil_remove_markers <- get_remove_markers("UKentucky_Tonsil")
stanford_mibi_lymphnode_pooled_remove_markers <- get_remove_markers("Stanford_MIBI_LymphNode")


# --- Section 4: Define Markers to Exclude (Gray Out) for Each Source ---
adjusted_rarecyte_ln_excluded_values <- process_excluded_markers("RareCyte_LN1_FOV1_FOV2_0.325")
adjusted_rarecyte_tma_excluded_values <- process_excluded_markers("RareCyte_TMA_FOV1_FOV2_0.325")
astar_comet_crc_excluded_values <- process_excluded_markers("ASTAR_COMET_CRC")
astar_comet_tonsil_excluded_values <- process_excluded_markers("ASTAR_COMET_Tonsil")
bidmc_excluded_values <- process_excluded_markers("BIDMC")
bidmc_dlbcl_excluded_values <- process_excluded_markers("BIDMC_DLBCL")
bidmc_tonsil_excluded_values <- process_excluded_markers("BIDMC_Tonsil")
bidmc_tonsil_compare_excluded_values <- process_excluded_markers("BIDMC_Tonsil_compare")
double_excluded_values <- process_excluded_markers("Double")
novartis_lung_cancer_excluded_values <- process_excluded_markers("Novartis_Lung_Cancer")
novartis_tonsil_excluded_values <- process_excluded_markers("Novartis_Tonsil")
roche_excluded_values <- process_excluded_markers("Roche")
roche_intestine_excluded_values <- process_excluded_markers("Roche_intestine")
roche_tonsil_excluded_values <- process_excluded_markers("Roche_Tonsil")
stanford_excluded_values <- process_excluded_markers("Stanford")
stanford_mibi_colon_excluded_values <- process_excluded_markers("Stanford_MIBI_Colon")
stanford_mibi_liver_excluded_values <- process_excluded_markers("Stanford_MIBI_Liver")
stanford_mibi_lymphnode_tile1_excluded_values <- process_excluded_markers("Stanford_MIBI_LymphNode_Tile1")
stanford_mibi_lymphnode_tile2_excluded_values <- process_excluded_markers("Stanford_MIBI_LymphNode_Tile2")
stanford_mibi_lymphnode_tile3_excluded_values <- process_excluded_markers("Stanford_MIBI_LymphNode_Tile3")
stanford_mibi_lymphnode_tile4_excluded_values <- process_excluded_markers("Stanford_MIBI_LymphNode_Tile4")
stanford_oscc_excluded_values <- process_excluded_markers("Stanford_OSCC")
stanford_rarecyte_ln_excluded_values <- process_excluded_markers("Stanford_RareCyte_LN")
stanford_rarecyte_tma_excluded_values <- process_excluded_markers("Stanford_RareCyte_TMA")
stanford_scan1_excluded_values <- process_excluded_markers("Stanford-scan1")
stanford_tonsil_excluded_values <- process_excluded_markers("Stanford_Tonsil")
ukentucky_scc_excluded_values <- process_excluded_markers("UKentucky_SCC")
ukentucky_tonsil_excluded_values <- process_excluded_markers("UKentucky_Tonsil")
stanford_mibi_lymphnode_pooled_excluded_values <- process_excluded_markers("Stanford_MIBI_LymphNode")


# --- Section 5: Define Penalty Scores for Each Source ---
penalty_scores <- list(
  Adjusted_RareCyte_LN = 10,
  Adjusted_RareCyte_TMA = 10,
  ASTAR_COMET_CRC = 10,
  ASTAR_COMET_Tonsil = 10,
  BIDMC = 30,
  BIDMC_DLBCL = 10,
  BIDMC_Tonsil = 10,
  BIDMC_Tonsil_compare = 10,
  Double = 1,
  Novartis_Lung_Cancer = 10,
  Novartis_Tonsil = 10,
  Roche = 10,
  Roche_intestine = 10,
  Roche_Tonsil = 10,
  Stanford = 18,
  Stanford_MIBI_Colon = 10,
  Stanford_MIBI_Liver = 10,
  Stanford_MIBI_LymphNode_Tile1 = 10,
  Stanford_MIBI_LymphNode_Tile2 = 10,
  Stanford_MIBI_LymphNode_Tile3 = 10,
  Stanford_MIBI_LymphNode_Tile4 = 10,
  Stanford_OSCC = 10,
  Stanford_RareCyte_LN = 10,
  Stanford_RareCyte_TMA = 10,
  `Stanford-scan1` = 18,
  Stanford_Tonsil = 10,
  UKentucky_SCC = 10,
  UKentucky_Tonsil = 10,
  Stanford_MIBI_LymphNode = 10
)

# --- Section 6: Special Subset Configuration (Handled Separately) ---
bidmc_subset_data <- slide_metadata %>%
  filter(Source == "BIDMC" & Type == "dataScaleSize") %>%
  filter(grepl("slide5_|slide13_|slide16_|slide21_", Filename)) %>%
  arrange(Name, FOV)

bidmc_subset_pairs <- comparison_pairs %>%
  filter(Source == "BIDMC") %>%
  filter(
    (Compare1 %in% c("BIDMC_5", "BIDMC_13", "BIDMC_16", "BIDMC_21") &
       Compare2 %in% c("BIDMC_5", "BIDMC_13", "BIDMC_16", "BIDMC_21"))
  ) %>% to_pairs()


# --- Section 7: Assemble Final Configurations List ---
configurations <- list(
  Adjusted_RareCyte_LN_all = list(
    data_folder = "./data_03272025/RareCyte_LN1_FOV1_FOV2_0.325/", out_folder = "./out_Adjusted_RareCyte_LN_all/",
    input_filenames = adjusted_rarecyte_ln_data$Filename, input_note = adjusted_rarecyte_ln_data$Name,
    pairs = adjusted_rarecyte_ln_pairs, remove_values = adjusted_rarecyte_ln_remove_markers,
    excluded_values = adjusted_rarecyte_ln_excluded_values, penalty_score = penalty_scores$Adjusted_RareCyte_LN
  ),
  Adjusted_RareCyte_TMA_all = list(
    data_folder = "./data_03272025/RareCyte_TMA_FOV1_FOV2_0.325/", out_folder = "./out_Adjusted_RareCyte_TMA_all/",
    input_filenames = adjusted_rarecyte_tma_data$Filename, input_note = adjusted_rarecyte_tma_data$Name,
    pairs = adjusted_rarecyte_tma_pairs, remove_values = adjusted_rarecyte_tma_remove_markers,
    excluded_values = adjusted_rarecyte_tma_excluded_values, penalty_score = penalty_scores$Adjusted_RareCyte_TMA
  ),
  ASTAR_COMET_CRC_all = list(
    data_folder = "./data_03272025/ASTAR_COMET_CRC/", out_folder = "./out_ASTAR_COMET_CRC_all/",
    input_filenames = astar_comet_crc_data$Filename, input_note = astar_comet_crc_data$Name,
    pairs = astar_comet_crc_pairs, remove_values = astar_comet_crc_remove_markers,
    excluded_values = astar_comet_crc_excluded_values, penalty_score = penalty_scores$ASTAR_COMET_CRC
  ),
  ASTAR_COMET_Tonsil_all = list(
    data_folder = "./data_03272025/ASTAR_COMET_Tonsil/", out_folder = "./out_ASTAR_COMET_Tonsil_all/",
    input_filenames = astar_comet_tonsil_data$Filename, input_note = astar_comet_tonsil_data$Name,
    pairs = astar_comet_tonsil_pairs, remove_values = astar_comet_tonsil_remove_markers,
    excluded_values = astar_comet_tonsil_excluded_values, penalty_score = penalty_scores$ASTAR_COMET_Tonsil
  ),
  BIDMC_all = list(
    data_folder = "./data_03272025/BIDMC/", out_folder = "./out_BIDMC_all/",
    input_filenames = bidmc_data$Filename, input_note = bidmc_data$Name,
    pairs = bidmc_pairs, remove_values = bidmc_remove_markers,
    excluded_values = bidmc_excluded_values, penalty_score = penalty_scores$BIDMC
  ),
  BIDMC_DLBCL_all = list(
    data_folder = "./data_03272025/BIDMC_DLBCL/", out_folder = "./out_BIDMC_DLBCL_all/",
    input_filenames = bidmc_dlbcl_data$Filename, input_note = bidmc_dlbcl_data$Name,
    pairs = bidmc_dlbcl_pairs, remove_values = bidmc_dlbcl_remove_markers,
    excluded_values = bidmc_dlbcl_excluded_values, penalty_score = penalty_scores$BIDMC_DLBCL
  ),
  BIDMC_Tonsil_all = list(
    data_folder = "./data_03272025/BIDMC_Tonsil/", out_folder = "./out_BIDMC_Tonsil_all/",
    input_filenames = bidmc_tonsil_data$Filename, input_note = bidmc_tonsil_data$Name,
    pairs = bidmc_tonsil_pairs, remove_values = bidmc_tonsil_remove_markers,
    excluded_values = bidmc_tonsil_excluded_values, penalty_score = penalty_scores$BIDMC_Tonsil
  ),
  BIDMC_Tonsil_compare_all = list(
    data_folder = "./data_03272025/BIDMC_Tonsil_compare/", out_folder = "./out_BIDMC_Tonsil_compare_all/",
    input_filenames = bidmc_tonsil_compare_data$Filename, input_note = bidmc_tonsil_compare_data$Name,
    pairs = bidmc_tonsil_compare_pairs, remove_values = bidmc_tonsil_compare_remove_markers,
    excluded_values = bidmc_tonsil_compare_excluded_values, penalty_score = penalty_scores$BIDMC_Tonsil_compare
  ),
  Double_all = list(
    data_folder = "./data_03272025/Double/", out_folder = "./out_Double_all/",
    input_filenames = double_data$Filename, input_note = double_data$Name,
    pairs = double_pairs, remove_values = double_remove_markers,
    excluded_values = double_excluded_values, penalty_score = penalty_scores$Double
  ),
  Novartis_Lung_Cancer_all = list(
    data_folder = "./data_03272025/Novartis_Lung_Cancer/", out_folder = "./out_Novartis_Lung_Cancer_all/",
    input_filenames = novartis_lung_cancer_data$Filename, input_note = novartis_lung_cancer_data$Name,
    pairs = novartis_lung_cancer_pairs, remove_values = novartis_lung_cancer_remove_markers,
    excluded_values = novartis_lung_cancer_excluded_values, penalty_score = penalty_scores$Novartis_Lung_Cancer
  ),
  Novartis_Tonsil_all = list(
    data_folder = "./data_03272025/Novartis_Tonsil/", out_folder = "./out_Novartis_Tonsil_all/",
    input_filenames = novartis_tonsil_data$Filename, input_note = novartis_tonsil_data$Name,
    pairs = novartis_tonsil_pairs, remove_values = novartis_tonsil_remove_markers,
    excluded_values = novartis_tonsil_excluded_values, penalty_score = penalty_scores$Novartis_Tonsil
  ),
  Roche_all = list(
    data_folder = "./data_03272025/Roche/", out_folder = "./out_Roche_all/",
    input_filenames = roche_data$Filename, input_note = roche_data$Name,
    pairs = roche_pairs, remove_values = roche_remove_markers,
    excluded_values = roche_excluded_values, penalty_score = penalty_scores$Roche
  ),
  Roche_intestine_all = list(
    data_folder = "./data_03272025/Roche_intestine/", out_folder = "./out_Roche_intestine_all/",
    input_filenames = roche_intestine_data$Filename, input_note = roche_intestine_data$Name,
    pairs = roche_intestine_pairs, remove_values = roche_intestine_remove_markers,
    excluded_values = roche_intestine_excluded_values, penalty_score = penalty_scores$Roche_intestine
  ),
  Roche_Tonsil_all = list(
    data_folder = "./data_03272025/Roche_Tonsil/", out_folder = "./out_Roche_Tonsil_all/",
    input_filenames = roche_tonsil_data$Filename, input_note = roche_tonsil_data$Name,
    pairs = roche_tonsil_pairs, remove_values = roche_tonsil_remove_markers,
    excluded_values = roche_tonsil_excluded_values, penalty_score = penalty_scores$Roche_Tonsil
  ),
  Stanford_all = list(
    data_folder = "./data_03272025/Stanford/", out_folder = "./out_Stanford_all/",
    input_filenames = stanford_data$Filename, input_note = stanford_data$Name,
    pairs = stanford_pairs, remove_values = stanford_remove_markers,
    excluded_values = stanford_excluded_values, penalty_score = penalty_scores$Stanford
  ),
  Stanford_MIBI_Colon_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_Colon/", out_folder = "./out_Stanford_MIBI_Colon_all/",
    input_filenames = stanford_mibi_colon_data$Filename, input_note = stanford_mibi_colon_data$Name,
    pairs = stanford_mibi_colon_pairs, remove_values = stanford_mibi_colon_remove_markers,
    excluded_values = stanford_mibi_colon_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_Colon
  ),
  Stanford_MIBI_Liver_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_Liver/", out_folder = "./out_Stanford_MIBI_Liver_all/",
    input_filenames = stanford_mibi_liver_data$Filename, input_note = stanford_mibi_liver_data$Name,
    pairs = stanford_mibi_liver_pairs, remove_values = stanford_mibi_liver_remove_markers,
    excluded_values = stanford_mibi_liver_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_Liver
  ),
  Stanford_MIBI_LymphNode_Tile1_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_LymphNode_Tile1/", out_folder = "./out_Stanford_MIBI_LymphNode_Tile1_all/",
    input_filenames = stanford_mibi_lymphnode_tile1_data$Filename, input_note = stanford_mibi_lymphnode_tile1_data$Name,
    pairs = stanford_mibi_lymphnode_tile1_pairs, remove_values = stanford_mibi_lymphnode_tile1_remove_markers,
    excluded_values = stanford_mibi_lymphnode_tile1_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_LymphNode_Tile1
  ),
  Stanford_MIBI_LymphNode_Tile2_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_LymphNode_Tile2/", out_folder = "./out_Stanford_MIBI_LymphNode_Tile2_all/",
    input_filenames = stanford_mibi_lymphnode_tile2_data$Filename, input_note = stanford_mibi_lymphnode_tile2_data$Name,
    pairs = stanford_mibi_lymphnode_tile2_pairs, remove_values = stanford_mibi_lymphnode_tile2_remove_markers,
    excluded_values = stanford_mibi_lymphnode_tile2_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_LymphNode_Tile2
  ),
  Stanford_MIBI_LymphNode_Tile3_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_LymphNode_Tile3/", out_folder = "./out_Stanford_MIBI_LymphNode_Tile3_all/",
    input_filenames = stanford_mibi_lymphnode_tile3_data$Filename, input_note = stanford_mibi_lymphnode_tile3_data$Name,
    pairs = stanford_mibi_lymphnode_tile3_pairs, remove_values = stanford_mibi_lymphnode_tile3_remove_markers,
    excluded_values = stanford_mibi_lymphnode_tile3_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_LymphNode_Tile3
  ),
  Stanford_MIBI_LymphNode_Tile4_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_LymphNode_Tile4/", out_folder = "./out_Stanford_MIBI_LymphNode_Tile4_all/",
    input_filenames = stanford_mibi_lymphnode_tile4_data$Filename, input_note = stanford_mibi_lymphnode_tile4_data$Name,
    pairs = stanford_mibi_lymphnode_tile4_pairs, remove_values = stanford_mibi_lymphnode_tile4_remove_markers,
    excluded_values = stanford_mibi_lymphnode_tile4_excluded_values, penalty_score = penalty_scores$Stanford_MIBI_LymphNode_Tile4
  ),
  Stanford_OSCC_all = list(
    data_folder = "./data_03272025/Stanford_OSCC/", out_folder = "./out_Stanford_OSCC_all/",
    input_filenames = stanford_oscc_data$Filename, input_note = stanford_oscc_data$Name,
    pairs = stanford_oscc_pairs, remove_values = stanford_oscc_remove_markers,
    excluded_values = stanford_oscc_excluded_values, penalty_score = penalty_scores$Stanford_OSCC
  ),
  Stanford_RareCyte_LN_all = list(
    data_folder = "./data_03272025/Stanford_RareCyte_LN/", out_folder = "./out_Stanford_RareCyte_LN_all/",
    input_filenames = stanford_rarecyte_ln_data$Filename, input_note = stanford_rarecyte_ln_data$Name,
    pairs = stanford_rarecyte_ln_pairs, remove_values = stanford_rarecyte_ln_remove_markers,
    excluded_values = stanford_rarecyte_ln_excluded_values, penalty_score = penalty_scores$Stanford_RareCyte_LN
  ),
  Stanford_RareCyte_TMA_all = list(
    data_folder = "./data_03272025/Stanford_RareCyte_TMA/", out_folder = "./out_Stanford_RareCyte_TMA_all/",
    input_filenames = stanford_rarecyte_tma_data$Filename, input_note = stanford_rarecyte_tma_data$Name,
    pairs = stanford_rarecyte_tma_pairs, remove_values = stanford_rarecyte_tma_remove_markers,
    excluded_values = stanford_rarecyte_tma_excluded_values, penalty_score = penalty_scores$Stanford_RareCyte_TMA
  ),
  Stanford_scan1_all = list(
    data_folder = "./data_03272025/Stanford-scan1/", out_folder = "./out_Stanford_scan1_all/",
    input_filenames = stanford_scan1_data$Filename, input_note = stanford_scan1_data$Name,
    pairs = stanford_scan1_pairs, remove_values = stanford_scan1_remove_markers,
    excluded_values = stanford_scan1_excluded_values, penalty_score = penalty_scores$`Stanford-scan1`
  ),
  Stanford_Tonsil_all = list(
    data_folder = "./data_03272025/Stanford_Tonsil/", out_folder = "./out_Stanford_Tonsil_all/",
    input_filenames = stanford_tonsil_data$Filename, input_note = stanford_tonsil_data$Name,
    pairs = stanford_tonsil_pairs, remove_values = stanford_tonsil_remove_markers,
    excluded_values = stanford_tonsil_excluded_values, penalty_score = penalty_scores$Stanford_Tonsil
  ),
  UKentucky_SCC_all = list(
    data_folder = "./data_03272025/UKentucky_SCC/", out_folder = "./out_UKentucky_SCC_all/",
    input_filenames = ukentucky_scc_data$Filename, input_note = ukentucky_scc_data$Name,
    pairs = ukentucky_scc_pairs, remove_values = ukentucky_scc_remove_markers,
    excluded_values = ukentucky_scc_excluded_values, penalty_score = penalty_scores$UKentucky_SCC
  ),
  UKentucky_Tonsil_all = list(
    data_folder = "./data_03272025/UKentucky_Tonsil/", out_folder = "./out_UKentucky_Tonsil_all/",
    input_filenames = ukentucky_tonsil_data$Filename, input_note = ukentucky_tonsil_data$Name,
    pairs = ukentucky_tonsil_pairs, remove_values = ukentucky_tonsil_remove_markers,
    excluded_values = ukentucky_tonsil_excluded_values, penalty_score = penalty_scores$UKentucky_Tonsil
  ),
  BIDMC_subset = list(
    data_folder = "./data_03272025/BIDMC/", out_folder = "./out_BIDMC_subset/",
    input_filenames = bidmc_subset_data$Filename, input_note = bidmc_subset_data$Name,
    pairs = bidmc_subset_pairs, remove_values = bidmc_remove_markers,
    excluded_values = bidmc_excluded_values, penalty_score = 10
  ),
  Stanford_MIBI_LymphNode_pooled_all = list(
    data_folder = "./data_03272025/Stanford_MIBI_LymphNode/", out_folder = "./out_Stanford_MIBI_LymphNode_pooled_all/",
    input_filenames = stanford_mibi_lymphnode_pooled_data$Filename, 
    input_note = stanford_mibi_lymphnode_pooled_data$Name,
    pairs = stanford_mibi_lymphnode_pooled_pairs, 
    remove_values = stanford_mibi_lymphnode_pooled_remove_markers,
    excluded_values = stanford_mibi_lymphnode_pooled_excluded_values, 
    penalty_score = penalty_scores$Stanford_MIBI_LymphNode
  )
)

# Choose which configuration to use - can be changed to process different datasets
# Options: "Adjusted_RareCyte_LN_all", "Adjusted_RareCyte_TMA_all", "ASTAR_COMET_CRC_all", "ASTAR_COMET_Tonsil_all",
#          "BIDMC_all", "BIDMC_DLBCL_all", "BIDMC_Tonsil_all", "BIDMC_Tonsil_compare_all", "Double_all",
#          "Novartis_Lung_Cancer_all", "Novartis_Tonsil_all", "Roche_all", "Roche_intestine_all", "Roche_Tonsil_all",
#          "Stanford_all", "Stanford_MIBI_Colon_all", "Stanford_MIBI_Liver_all",
#          "Stanford_MIBI_LymphNode_Tile1_all", "Stanford_MIBI_LymphNode_Tile2_all",
#          "Stanford_MIBI_LymphNode_Tile3_all", "Stanford_MIBI_LymphNode_Tile4_all",
#          "Stanford_OSCC_all", "Stanford_RareCyte_LN_all", "Stanford_RareCyte_TMA_all",
#          "Stanford_scan1_all", "Stanford_Tonsil_all", "UKentucky_SCC_all",
#          "UKentucky_Tonsil_all", "BIDMC_subset", "Stanford_MIBI_LymphNode_pooled_all"
current_config_name <- "Stanford_all" # Set to run the new configuration
current_config <- configurations[[current_config_name]]

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

### Load MESMER data. Try to load from qsave if available, otherwise load raw data
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
  height = 11
)

### Perform statistical and Kruskal-Wallis tests
# Ensure perform_statistical_and_kruskal_wallis_tests is defined in helper.R
kruskal_results <- perform_statistical_and_kruskal_wallis_tests(df_trans, marker_names)

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, paste0(out_folder, "kruskal_pvals.csv"))

# Ensure calc_effect_size is defined in helper.R
effect_size <- calc_effect_size(
  data = df_trans,
  out_folder = out_folder,
  pairs = pairs
)

# Generate heatmaps
# Ensure plot_heatmaps is defined in helper.R
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
# Ensure implement_scoring_system is defined in helper.R
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
