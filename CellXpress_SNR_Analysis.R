library(dplyr)
library(tidyverse)
library(qs)
library(ggplot2)
library(stringr)
library(cowplot)

# Fix for font rendering issues with svglite
# Reset systemfonts cache if it exists
if (requireNamespace("systemfonts", quietly = TRUE)) {
  tryCatch({
    systemfonts::reset_font_cache()
  }, error = function(e) {
    message("Note: Could not reset font cache")
  })
}

# Set a safe default theme
theme_set(theme_minimal())

source("helper.R")

# ==============================================================================
# Custom safe ggsave function with PNG fallback
# ==============================================================================

#' Safe ggsave with PNG fallback
#' 
#' Attempts to save as SVG, falls back to PNG if font errors occur.
#' 
#' @param filename Output filename
#' @param plot ggplot object to save
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param ... Additional arguments passed to ggsave
safe_ggsave <- function(filename, plot, width = 10, height = 8, ...) {
  tryCatch({
    ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, ...)
  }, error = function(e) {
    if (grepl("vector is too large|font|locate_fonts", e$message, ignore.case = TRUE)) {
      # Fall back to PNG
      png_filename <- gsub("\\.svg$", ".png", filename)
      message("SVG save failed, falling back to PNG: ", basename(png_filename))
      ggplot2::ggsave(filename = png_filename, plot = plot, width = width, height = height, 
             dpi = 300, ...)
    } else {
      stop(e)
    }
  })
}

#' Create SNR Heatmaps (Safe Version)
#' 
#' Creates SNR heatmaps with PNG output to avoid font rendering issues.
#' This is a simplified version of create_snr_heatmaps from helper.R.
#' 
#' @param snr_data The combined, wide-format SNR data frame with pre-normalized values.
#' @param out_folder The directory to save the output heatmap files.
#' @param remove_values A vector of marker names to completely exclude from the heatmaps.
#' @param excluded_values A named list specifying marker/condition pairs to grey out (NA).
#' @return A list containing the generated ggplot objects.
create_snr_heatmaps_safe <- function(snr_data, out_folder, remove_values = c(), excluded_values = list()) {
  
  # --- 1. DATA PREPARATION ---
  snr_data_clean <- snr_data %>%
    rename_with(~ gsub("[.-]", "", .), all_of(setdiff(names(snr_data), "Staining_condition")))
  
  marker_cols <- setdiff(names(snr_data_clean), "Staining_condition")
  
  final_marker_cols <- marker_cols
  if (length(remove_values) > 0) {
    cleaned_remove_values <- gsub("[.-]", "", remove_values)
    final_marker_cols <- setdiff(marker_cols, cleaned_remove_values)
    message("Removed ", length(marker_cols) - length(final_marker_cols), " markers: ",
            paste(setdiff(marker_cols, final_marker_cols), collapse = ", "))
  }
  
  condition_order <- snr_data_clean %>%
    mutate(condition_num = as.numeric(gsub("\\D", "", Staining_condition))) %>%
    arrange(condition_num) %>%
    pull(Staining_condition) %>%
    unique()
  
  snr_long <- snr_data_clean %>%
    select(Staining_condition, all_of(final_marker_cols)) %>%
    pivot_longer(
      cols = all_of(final_marker_cols),
      names_to = "Marker",
      values_to = "SNR"
    ) %>%
    mutate(
      Staining_condition = factor(Staining_condition, levels = condition_order),
      Marker = factor(Marker, levels = final_marker_cols)  # Preserve marker order
    )
  
  # Apply exclusions (set specific marker/condition combinations to NA)
  if (length(excluded_values) > 0) {
    exclusion_count <- 0
    for (condition in names(excluded_values)) {
      markers_to_exclude <- excluded_values[[condition]]
      for (marker in markers_to_exclude) {
        # Clean marker name to match data
        cleaned_marker <- gsub("[.-]", "", marker)
        if (cleaned_marker %in% final_marker_cols) {
          snr_long <- snr_long %>%
            mutate(SNR = ifelse(Staining_condition == condition & Marker == cleaned_marker, NA, SNR))
          exclusion_count <- exclusion_count + 1
        }
      }
    }
    if (exclusion_count > 0) {
      message("Applied ", exclusion_count, " marker/slide exclusions (set to NA)")
    }
  }
  
  result_list <- list()
  
  # --- 2. CALCULATIONS ---
  
  # A. Normalize each marker's SNR to a 0-1 scale
  snr_norm_01 <- snr_long %>%
    group_by(Marker) %>%
    mutate(
      range = max(SNR, na.rm = TRUE) - min(SNR, na.rm = TRUE),
      SNR_norm_01 = if_else(range == 0, 0.5, (SNR - min(SNR, na.rm = TRUE)) / range)
    ) %>%
    ungroup()
  
  # B. Calculate the average of these 0-1 normalized values for each slide
  avg_snr_data <- snr_norm_01 %>%
    group_by(Staining_condition) %>%
    summarise(Avg_SNR = mean(SNR_norm_01, na.rm = TRUE), .groups = 'drop') %>%
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  # C. Prepare Z-score data from ORIGINAL SNR
  snr_long_zscore <- snr_long %>%
    group_by(Marker) %>%
    mutate(
      s = sd(SNR, na.rm = TRUE),
      SNR_zscore = if_else(is.na(s) | s == 0, 0, (SNR - mean(SNR, na.rm = TRUE)) / s)
    ) %>%
    ungroup()
  
  # --- 3. CREATE PLOTS ---
  
  # PLOT 1: Combined Plot (Bar chart over Z-score heatmap) - MAIN PLOT
  p_top_barplot <- ggplot(avg_snr_data, aes(x = Staining_condition, y = Avg_SNR)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Avg_SNR)), vjust = -0.4, size = 4.5) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(b = 3),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    labs(y = "Avg Score", title = "Average Normalized SNR and Z-Score Heatmap") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  p_bottom_heatmap <- ggplot(snr_long_zscore, aes(x = Staining_condition, y = Marker, fill = SNR_zscore)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = c("#330066", "#44007C", "#5B0092", "#7209B7", "#9D4EDD", "#C77DFF", 
                 "#E0AAFF", "#F3E5F5", "#FFF4E6", "#FFE082", "#FFD54F", "#FFCA28", 
                 "#FFC107", "#FFB300", "#FFA000"),
      values = scales::rescale(c(-3, -2.5, -2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3)),
      limits = c(-3, 3),
      oob = scales::squish,
      na.value = "grey80"
    ) +
    geom_text(aes(label = ifelse(is.na(SNR_zscore), "", sprintf("%.2f", SNR_zscore))), size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(0.5, 5, 5, 5)
    ) +
    labs(y = "Marker", fill = "Z-Score")
  
  combined_plot <- cowplot::plot_grid(
    p_top_barplot, p_bottom_heatmap,
    ncol = 1, align = "v",
    rel_heights = c(0.3, 0.7), axis = "lr"
  )
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_heatmap_zscore_with_avg.svg"),
    plot = combined_plot, width = 12, height = 11
  )
  result_list$combined_plot <- combined_plot
  
  # PLOT 2: Pre-normalized SNR heatmap
  heatmap_prenorm <- ggplot(snr_long, aes(x = Staining_condition, y = Marker, fill = SNR)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(low = "yellow", mid = "white", high = "#800080", 
                         midpoint = median(snr_long$SNR, na.rm = TRUE), na.value = "grey80") +
    geom_text(aes(label = ifelse(is.na(SNR), "", sprintf("%.2f", SNR))), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = "Pre-Normalized Signal-to-Noise Ratio Heatmap", x = "Staining Condition", y = "Marker", fill = "SNR")
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_heatmap_prenormalized.svg"),
    plot = heatmap_prenorm, width = 10, height = 8
  )
  result_list$heatmap_prenormalized <- heatmap_prenorm
  
  # PLOT 3: 0-1 Normalized SNR heatmap
  heatmap_01norm <- ggplot(snr_norm_01, aes(x = Staining_condition, y = Marker, fill = SNR_norm_01)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(low = "yellow", mid = "white", high = "#800080", midpoint = 0.5, na.value = "grey80") +
    geom_text(aes(label = ifelse(is.na(SNR_norm_01), "", sprintf("%.2f", SNR_norm_01))), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = "0-1 Normalized Signal-to-Noise Ratio Heatmap", x = "Staining Condition", y = "Marker", fill = "0-1 Norm. SNR")
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_heatmap_01normalized.svg"),
    plot = heatmap_01norm, width = 10, height = 8
  )
  result_list$heatmap_01normalized <- heatmap_01norm
  
  # PLOT 4: Threshold-based heatmap
  threshold_colors <- c("<0.8" = "#FF9999", "0.8-1.0" = "#FFCC99", "1.0-1.5" = "#FFFF99", 
                        "1.5-2.0" = "#CCFF99", "2.0-3.0" = "#99FF99", ">3.0" = "#009900")
  snr_long_threshold <- snr_long %>%
    mutate(SNR_category = case_when(
      SNR < 0.8 ~ "<0.8", SNR >= 0.8 & SNR < 1.0 ~ "0.8-1.0", SNR >= 1.0 & SNR < 1.5 ~ "1.0-1.5",
      SNR >= 1.5 & SNR < 2.0 ~ "1.5-2.0", SNR >= 2.0 & SNR < 3.0 ~ "2.0-3.0", SNR >= 3.0 ~ ">3.0",
      TRUE ~ NA_character_
    )) %>%
    mutate(SNR_category = factor(SNR_category, levels = names(threshold_colors)))
  
  heatmap_threshold <- ggplot(snr_long_threshold, aes(x = Staining_condition, y = Marker, fill = SNR_category)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_manual(values = threshold_colors, na.value = "grey80", drop = FALSE) +
    geom_text(aes(label = ifelse(is.na(SNR), "", sprintf("%.2f", SNR))), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = "Pre-Normalized SNR Heatmap (Threshold-Based)", x = "Staining Condition", y = "Marker", fill = "SNR Range")
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_heatmap_threshold.svg"),
    plot = heatmap_threshold, width = 10, height = 8
  )
  result_list$heatmap_threshold <- heatmap_threshold
  
  # PLOT 5: Z-Score heatmap
  heatmap_zscore <- ggplot(snr_long_zscore, aes(x = Staining_condition, y = Marker, fill = SNR_zscore)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(low = "yellow", mid = "white", high = "#800080", midpoint = 0, na.value = "grey80") +
    geom_text(aes(label = ifelse(is.na(SNR_zscore), "", sprintf("%.2f", SNR_zscore))), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    labs(title = "Z-Score of Pre-Normalized Signal-to-Noise Ratio", x = "Staining Condition", y = "Marker", fill = "Z-Score")
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_heatmap_zscore.svg"),
    plot = heatmap_zscore, width = 10, height = 8
  )
  result_list$heatmap_zscore <- heatmap_zscore
  
  # Bar plot of mean SNR per marker
  marker_means <- snr_long %>%
    group_by(Marker) %>%
    summarize(mean_SNR = mean(SNR, na.rm = TRUE)) %>%
    arrange(desc(mean_SNR))
  
  barplot_means <- ggplot(marker_means, aes(x = reorder(Marker, mean_SNR), y = mean_SNR)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.2f", mean_SNR)), hjust = -0.1, size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Mean Pre-Normalized SNR per Marker", x = "Marker", y = "Mean SNR", 
         caption = "Red dashed line indicates SNR = 1.0")
  
  safe_ggsave(
    filename = file.path(out_folder, "SNR_barplot_means.svg"),
    plot = barplot_means, width = 8, height = 8
  )
  result_list$barplot_means <- barplot_means
  
  # Save the processed data
  write_csv(snr_data_clean %>% select(Staining_condition, all_of(final_marker_cols)), 
            file.path(out_folder, "processed_snr_data.csv"))
  write_csv(marker_means, file.path(out_folder, "marker_mean_snr.csv"))
  
  message("Successfully generated and saved all plots to: ", out_folder)
  
  return(result_list)
}

# ==============================================================================
# CellXpress SNR Analysis Pipeline
# This script processes CellXpress segmentation SNR data with the same workflow
# as MESMER, including FOV cell count weighting and Z-score heatmap generation.
# ==============================================================================

# ==============================================================================
# 1. Configuration
# ==============================================================================

config <- list(
  # Path to the SNR data file (raw signal ratios from CellXpress)
  snr_file = "./data_cellXpress/SNR_BIDMC/SNR_ratios_BIDMC.csv",
  
  # Path to the metadata file mapping SNR slide names to slide identifiers
  # Format: SNR_slidename, Name
  metadata_file = "./data_cellXpress/SNR_BIDMC/metadata_BIDMC.csv",
  
  # Path to the cell counts file (pre-calculated cell counts per FOV)
  # Format: Slide, FOV, Cell count
  cell_counts_file = "./data_cellXpress/SNR_BIDMC/cell_counts_BIDMC.csv",
  
  # Paths to marker removal and exclusion files (same as MESMER)
  removal_file = "./data_mesmer/Slide_remove_markers.csv",
  exclusion_file = "./data_mesmer/Slide_exclude_markers.csv",
  
  # Source name for looking up removal/exclusion rules
  source_name = "BIDMC",
  
  # Output folder
  out_folder = "./out_CellXpress_SNR/",
  
  # Markers to remove from analysis (e.g., DAPI is used for normalization)
  # Additional markers from Slide_remove_markers.csv will be added automatically
  remove_markers = c("DAPI"),
  
  # Marker name mapping: CellXpress names -> MESMER names
  # This ensures consistent naming between pipelines
  marker_name_mapping = c(
    "Tox_Tox2" = "ToxTox2",
    "PD_1" = "PD1",
    "Na_K_ATPase" = "NAKATPase",
    "IDO_1" = "IDO1",
    "HLA_DRA" = "HLADRA",
    "Granzyme_B" = "GranzymeB",
    "FOXP3" = "FoxP3",
    "CK" = "Cytokeratin",
    "DC_SIGN" = "DCSIGN"
    # Markers with same names: pS6, Pax5, Ki67, H3K27me3, H3K27ac, CD8, CD68, 
    # CD45RA, CD45, CD4, CD31, CD3, CD20, CD15, CD138, CD11c, CD11b, CD56
  ),
  
  # Marker order to match MESMER BIDMC SNR result (reversed to match heatmap y-axis)
  # Note: CD56 and DCSIGN will be removed per Slide_remove_markers.csv for BIDMC
  marker_sequence = c(
    # MESMER markers in reverse order (so they appear correctly on y-axis from bottom to top)
    "CD11b", "CD11c", "CD138", "CD15", "CD20", "CD3", "CD31", "CD4", "CD45", "CD45RA",
    "CD68", "CD8", "Cytokeratin", "FoxP3", "GranzymeB", "H3K27ac", "H3K27me3", "HLADRA",
    "IDO1", "Ki67", "NAKATPase", "Pax5", "PD1", "pS6", "ToxTox2"
  )
)

# ==============================================================================
# 2. Helper Functions
# ==============================================================================

#' Rename Markers to Match MESMER Naming Convention
#' 
#' Applies marker name mapping to convert CellXpress names to MESMER names.
#' 
#' @param marker_names Vector of marker column names
#' @param mapping Named vector of mappings (CellXpress name -> MESMER name)
#' @return Vector of renamed marker names
rename_markers <- function(marker_names, mapping) {
  renamed <- marker_names
  for (old_name in names(mapping)) {
    renamed[renamed == old_name] <- mapping[old_name]
  }
  return(renamed)
}

#' Load Markers to Remove
#' 
#' Loads markers that should be completely removed from analysis for a given source.
#' 
#' @param removal_file Path to the Slide_remove_markers.csv file
#' @param source_name Source identifier (e.g., "BIDMC")
#' @return Vector of marker names to remove
load_removal_markers <- function(removal_file, source_name) {
  if (!file.exists(removal_file)) {
    message("Removal file not found: ", removal_file)
    return(c())
  }
  
  removed_markers <- read_csv(removal_file, show_col_types = FALSE)
  
  # Filter for the specified source and marker type
  source_markers <- removed_markers %>%
    filter(Source == source_name & Exclude_type == "Marker") %>%
    pull(Exclude_value)
  
  message("Found ", length(source_markers), " markers to remove for source '", source_name, "': ",
          paste(source_markers, collapse = ", "))
  
  return(source_markers)
}

#' Load Marker Exclusions (for specific slides)
#' 
#' Loads markers that should be excluded (greyed out) for specific slides.
#' 
#' @param exclusion_file Path to the Slide_exclude_markers.csv file
#' @param source_name Source identifier (e.g., "BIDMC")
#' @return Named list where names are slide conditions and values are vectors of markers to exclude
load_exclusion_markers <- function(exclusion_file, source_name) {
  if (!file.exists(exclusion_file)) {
    message("Exclusion file not found: ", exclusion_file)
    return(list())
  }
  
  excluded_markers <- read_csv(exclusion_file, show_col_types = FALSE)
  
  # Filter for the specified source
  source_exclusions <- excluded_markers %>%
    filter(Source == source_name)
  
  if (nrow(source_exclusions) == 0) {
    message("No marker exclusions found for source '", source_name, "'")
    return(list())
  }
  
  # Convert to list format: condition -> markers
  exclusion_list <- list()
  for (i in 1:nrow(source_exclusions)) {
    condition <- source_exclusions$Name[i]
    marker <- source_exclusions$Marker[i]
    
    # Convert BIDMC_1 to BIDMC_Slide1 format for matching
    if (grepl("^BIDMC_\\d+$", condition)) {
      slide_num <- gsub("BIDMC_", "", condition)
      condition <- paste0("BIDMC_Slide", slide_num)
    }
    
    if (!(condition %in% names(exclusion_list))) {
      exclusion_list[[condition]] <- c()
    }
    exclusion_list[[condition]] <- c(exclusion_list[[condition]], marker)
  }
  
  message("Found marker exclusions for ", length(exclusion_list), " slides")
  
  return(exclusion_list)
}

#' Load and Clean CellXpress SNR Data
#' 
#' Reads the SNR CSV and extracts signal columns based on normalization type.
#' Cleans marker names by removing the suffix and applies MESMER naming convention.
#' 
#' @param snr_file Path to the ratio_renamed.csv file
#' @param marker_mapping Named vector for renaming markers to MESMER convention
#' @param norm_type Normalization type: "Normalized" (default), "DAPInorm", or "areanorm"
#' @return A data frame with Slide_Region, Label, and marker SNR values
load_cellxpress_snr <- function(snr_file, marker_mapping = NULL, norm_type = "Normalized") {
  message("Loading SNR data from: ", snr_file)
  message("Normalization type: ", norm_type)
  
  snr_data <- read_csv(snr_file, show_col_types = FALSE)
  
  # Define column suffix based on normalization type
  suffix_pattern <- switch(norm_type,
    "Normalized" = "_Normalized_signal_invsout",
    "DAPInorm" = "_signal_invsout_DAPInorm",
    "areanorm" = "_signal_invsout_areanorm",
    stop("Unknown normalization type: ", norm_type)
  )
  
  # Select the Slide/Region, Label, and appropriate signal columns
  snr_clean <- snr_data %>%
    select(`Slide/Region`, Label, ends_with(suffix_pattern))
  
  # Clean up column names: remove the suffix
  colnames(snr_clean) <- gsub(suffix_pattern, "", colnames(snr_clean))
  
  # Apply marker name mapping if provided
  if (!is.null(marker_mapping)) {
    colnames(snr_clean) <- rename_markers(colnames(snr_clean), marker_mapping)
    message("Applied marker name mapping to match MESMER convention")
  }
  
  # Rename Slide/Region for easier handling
  snr_clean <- snr_clean %>%
    rename(Slide_Region = `Slide/Region`)
  
  message("Found ", ncol(snr_clean) - 2, " markers in SNR data")
  message("Found ", nrow(snr_clean), " total FOV/region entries")
  
  return(snr_clean)
}

#' Load Metadata Mapping
#' 
#' Reads the metadata CSV that maps SNR slide names to human-readable slide IDs.
#' 
#' @param metadata_file Path to the metadata CSV
#' @return A data frame with SNR_slidename and Name columns
load_metadata <- function(metadata_file) {
  message("Loading metadata from: ", metadata_file)
  
  meta <- read_csv(metadata_file, show_col_types = FALSE)
  
  if (!all(c("SNR_slidename", "Name") %in% colnames(meta))) {
    stop("Metadata file must contain 'SNR_slidename' and 'Name' columns.")
  }
  
  message("Found ", nrow(meta), " slide mappings in metadata")
  
  return(meta)
}

#' Load Cell Counts
#' 
#' Reads the cell counts CSV with FOV-level counts.
#' Creates a mapping from Slide + FOV to Cell count.
#' 
#' @param cell_counts_file Path to the cell counts CSV
#' @return A data frame with Slide, FOV, and Cell_count columns
load_cell_counts <- function(cell_counts_file) {
  message("Loading cell counts from: ", cell_counts_file)
  
  counts <- read_csv(cell_counts_file, show_col_types = FALSE)
  
  # Standardize column names
  counts <- counts %>%
    rename(
      Slide = Slide,
      FOV = FOV,
      Cell_count = `Cell count`
    )
  
  message("Found cell counts for ", n_distinct(counts$Slide), " slides and ", 
          nrow(counts), " total FOVs")
  
  return(counts)
}

#' Map FOV Labels to FOV Numbers
#' 
#' Determines FOV number (1 or 2) based on the Label column (e.g., A01-A12 = FOV1, A13+ = FOV2).
#' BIDMC data uses a specific pattern where the first 12 regions are FOV1.
#' 
#' @param labels A vector of Label values (e.g., "A01", "A02", etc.)
#' @return A vector of FOV numbers (1 or 2)
map_label_to_fov <- function(labels) {
  # Extract numeric part from label (e.g., "A01" -> 1, "A13" -> 13)
  label_nums <- as.numeric(gsub("[^0-9]", "", labels))
  
  # First 12 regions are FOV1, rest are FOV2
  # This may need to be adjusted based on actual data structure
  fov_nums <- ifelse(label_nums <= 12, 1, 2)
  
  return(fov_nums)
}

#' Map Region Labels to FOV Numbers
#' 
#' Maps CellXpress region labels (A01-A08) to FOV numbers (1 or 2).
#' Assumes: A01-A04 = FOV1, A05-A08 = FOV2
#' 
#' @param labels Vector of region labels (e.g., "A01", "A02", etc.)
#' @return Vector of FOV numbers (1 or 2)
map_region_to_fov <- function(labels) {
  # Extract numeric part from label (e.g., "A01" -> 1, "A05" -> 5)
  region_nums <- as.numeric(gsub("[^0-9]", "", labels))
  
  # First 4 regions (1-4) = FOV1, next 4 regions (5-8) = FOV2
  fov_nums <- ifelse(region_nums <= 4, 1, 2)
  
  return(fov_nums)
}

#' Combine SNR Data with Cell Count Weights
#' 
#' Merges SNR data with metadata and calculates weighted average SNR 
#' for each slide using CellXpress cell counts per FOV.
#' 
#' CellXpress has 8 regions per slide mapped to 2 FOVs:
#' - A01-A04 = FOV1
#' - A05-A08 = FOV2
#' 
#' @param snr_data Cleaned SNR data frame
#' @param metadata Metadata mapping data frame
#' @param cell_counts Cell counts data frame with Slide, FOV, Cell_count
#' @return A data frame with weighted average SNR per slide
combine_snr_with_weights <- function(snr_data, metadata, cell_counts = NULL) {
  
  message("Combining SNR data with cell count weights...")
  
  # Add slide name and FOV number to SNR data
  snr_with_meta <- snr_data %>%
    left_join(metadata, by = c("Slide_Region" = "SNR_slidename")) %>%
    mutate(
      # Extract slide number from Name (e.g., "BIDMC_Slide1" -> 1)
      Slide_num = as.numeric(gsub("\\D", "", Name)),
      # Map region label to FOV number (A01-A04 = FOV1, A05-A08 = FOV2)
      FOV = map_region_to_fov(Label)
    )
  
  # Check for unmapped slides
  unmapped <- snr_with_meta %>% filter(is.na(Name))
  if (nrow(unmapped) > 0) {
    warning(nrow(unmapped), " rows could not be mapped to slide names")
  }
  
  # Get marker columns (exclude metadata columns)
  marker_cols <- setdiff(colnames(snr_with_meta), 
                         c("Slide_Region", "Label", "Name", "Slide_num", "FOV"))
  
  # First, average regions within each FOV
  fov_avg <- snr_with_meta %>%
    filter(!is.na(Name)) %>%
    group_by(Name, Slide_num, FOV) %>%
    summarise(
      across(all_of(marker_cols), ~ mean(.x, na.rm = TRUE)),
      N_regions = n(),
      .groups = 'drop'
    )
  
  message("Averaged ", unique(fov_avg$N_regions), " regions per FOV")
  
  # Add cell counts if provided
  if (!is.null(cell_counts) && nrow(cell_counts) > 0) {
    message("Applying cell count weights from CellXpress...")
    
    fov_avg <- fov_avg %>%
      left_join(cell_counts, by = c("Slide_num" = "Slide", "FOV" = "FOV"))
    
    # Check for missing cell counts
    missing_counts <- fov_avg %>% filter(is.na(Cell_count))
    if (nrow(missing_counts) > 0) {
      warning(nrow(missing_counts), " FOVs missing cell counts, using equal weights")
      fov_avg <- fov_avg %>%
        mutate(Cell_count = replace_na(Cell_count, 1))
    }
    
    # Calculate weighted average SNR per slide
    weighted_snr <- fov_avg %>%
      group_by(Name) %>%
      mutate(
        Total_count = sum(Cell_count, na.rm = TRUE),
        Weight = Cell_count / Total_count
      ) %>%
      summarise(
        across(all_of(marker_cols), ~ sum(.x * Weight, na.rm = TRUE)),
        Total_cells = first(Total_count),
        .groups = 'drop'
      ) %>%
      rename(Staining_condition = Name)
    
    message("Calculated cell-count weighted SNR for ", nrow(weighted_snr), " slides")
    
  } else {
    # No cell counts - use simple average
    message("No cell counts provided, using simple average across FOVs")
    
    weighted_snr <- fov_avg %>%
      group_by(Name) %>%
      summarise(
        across(all_of(marker_cols), ~ mean(.x, na.rm = TRUE)),
        Total_cells = sum(N_regions),
        .groups = 'drop'
      ) %>%
      rename(Staining_condition = Name)
    
    message("Calculated simple average SNR for ", nrow(weighted_snr), " slides")
  }
  
  return(weighted_snr)
}

# ==============================================================================
# 3. Main Processing Function
# ==============================================================================

#' Reorder Markers According to Sequence
#' 
#' Reorders the marker columns in a data frame according to a specified sequence.
#' Markers not in the sequence are placed at the end.
#' 
#' @param data Data frame with Staining_condition and marker columns
#' @param marker_sequence Vector specifying desired marker order
#' @return Data frame with reordered columns
reorder_markers <- function(data, marker_sequence) {
  # Get current marker columns (all except Staining_condition and Total_cells)
  current_markers <- setdiff(colnames(data), c("Staining_condition", "Total_cells"))
  
  # Find common markers between data and sequence
  common_markers <- intersect(marker_sequence, current_markers)
  
  # Find markers in data but not in sequence (place at end)
  extra_markers <- setdiff(current_markers, marker_sequence)
  
  # Create final order
  final_order <- c(common_markers, extra_markers)
  
  message("Marker order applied: ", length(common_markers), " markers matched, ", 
          length(extra_markers), " extra markers at end")
  
  if (length(extra_markers) > 0) {
    message("Extra markers not in sequence: ", paste(extra_markers, collapse = ", "))
  }
  
  # Check for missing markers
  missing_markers <- setdiff(marker_sequence, current_markers)
  if (length(missing_markers) > 0) {
    message("Note: ", length(missing_markers), " markers from sequence not found in data: ",
            paste(missing_markers, collapse = ", "))
  }
  
  # Reorder data
  other_cols <- intersect(c("Staining_condition", "Total_cells"), colnames(data))
  data_reordered <- data %>%
    select(all_of(other_cols), all_of(final_order))
  
  return(data_reordered)
}

#' Process CellXpress SNR Data
#' 
#' Main function that orchestrates the entire workflow:
#' 1. Load all data files
#' 2. Combine and weight by cell counts
#' 3. Apply marker renaming and ordering
#' 4. Generate heatmaps and save results
#' 
#' @param config Configuration list with file paths and settings
#' @return List of generated plots
process_cellxpress_snr <- function(config, norm_type = "Normalized", out_suffix = "") {
  
  # Create output directory with suffix for normalization type
  out_folder <- if (out_suffix != "") {
    paste0(config$out_folder, "_", out_suffix)
  } else {
    config$out_folder
  }
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  message("\n##########################################################")
  message("# Processing normalization: ", norm_type)
  message("# Output folder: ", out_folder)
  message("##########################################################")
  
  # --- Step 1: Load marker removal and exclusion rules ---
  message("\n--- Loading marker removal/exclusion rules ---")
  removal_markers <- load_removal_markers(config$removal_file, config$source_name)
  exclusion_markers <- load_exclusion_markers(config$exclusion_file, config$source_name)
  
  # Combine removal markers with default remove_markers
  all_remove_markers <- unique(c(config$remove_markers, removal_markers))
  message("Total markers to remove: ", paste(all_remove_markers, collapse = ", "))
  
  # --- Step 2: Load all data with marker renaming ---
  message("\n--- Loading SNR data ---")
  snr_data <- load_cellxpress_snr(config$snr_file, config$marker_name_mapping, norm_type)
  metadata <- load_metadata(config$metadata_file)
  cell_counts <- load_cell_counts(config$cell_counts_file)
  
  # --- Step 3: Combine with weights ---
  weighted_snr <- combine_snr_with_weights(snr_data, metadata, cell_counts)
  
  # --- Step 4: Apply marker ordering to match MESMER ---
  message("\nApplying marker order to match MESMER BIDMC SNR result...")
  weighted_snr <- reorder_markers(weighted_snr, config$marker_sequence)
  
  # Save weighted SNR data
  write_csv(weighted_snr, file.path(out_folder, "weighted_snr_results.csv"))
  message("Saved weighted SNR results")
  
  # --- Step 5: Prepare data for heatmap functions ---
  # Remove unwanted markers and format for create_snr_heatmaps function
  snr_for_plot <- weighted_snr %>%
    select(-Total_cells)  # Remove the count column
  
  # --- Step 6: Create heatmaps ---
  message("\nGenerating SNR heatmaps...")
  
  # Use custom safe heatmap function to avoid font issues with SVG
  heatmap_results <- create_snr_heatmaps_safe(
    snr_data = snr_for_plot,
    out_folder = out_folder,
    remove_values = all_remove_markers,
    excluded_values = exclusion_markers
  )
  
  message("Done! All results saved to: ", out_folder)
  
  return(list(
    norm_type = norm_type,
    out_folder = out_folder,
    weighted_snr = weighted_snr,
    heatmap_results = heatmap_results
  ))
}

# ==============================================================================
# 4. Execute the Analysis for All Normalization Types
# ==============================================================================

# Define all normalization types with descriptions
norm_types <- list(
  list(
    type = "Normalized", 
    suffix = "Normalized",
    desc = "*_Normalized_signal_invsout - Comparable to MESMER (values 0-1)"
  ),
  list(
    type = "DAPInorm", 
    suffix = "DAPInorm",
    desc = "*_signal_invsout_DAPInorm - DAPI-normalized ratio (values can be >1)"
  ),
  list(
    type = "areanorm", 
    suffix = "AreaNorm",
    desc = "*_signal_invsout_areanorm - Area-normalized signal"
  )
)

# Process all normalization types
all_results <- list()
for (nt in norm_types) {
  message("\n")
  message("==========================================================")
  message("Processing: ", nt$desc)
  message("==========================================================")
  
  tryCatch({
    result <- process_cellxpress_snr(config, norm_type = nt$type, out_suffix = nt$suffix)
    all_results[[nt$type]] <- result
  }, error = function(e) {
    message("ERROR processing ", nt$type, ": ", e$message)
  })
}

# Print summary
message("\n========================================")
message("CellXpress SNR Analysis Complete!")
message("========================================")
message("\nProcessed ", length(all_results), " normalization types:")
for (nt in norm_types) {
  out_folder <- paste0(config$out_folder, "_", nt$suffix)
  message("\n", nt$suffix, ": ", nt$desc)
  message("  Output folder: ", out_folder)
}
message("\n\nGenerated plots in each folder:")
message("  - SNR_heatmap_zscore_with_avg.svg (Combined Z-score heatmap)")
message("  - SNR_heatmap_prenormalized.svg")
message("  - SNR_heatmap_01normalized.svg")
message("  - SNR_heatmap_threshold.svg")
message("  - SNR_heatmap_zscore.svg")
message("  - SNR_heatmap_zscore_of_01normed.svg")
message("  - SNR_barplot_means.svg")
message("  - weighted_snr_results.csv")
message("\n\n*** RECOMMENDED for comparison with MESMER: ", config$out_folder, "_Normalized ***")
