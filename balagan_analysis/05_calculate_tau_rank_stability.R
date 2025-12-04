# =========================================================================================
# === SCRIPT TO QUANTIFY TAU RANK STABILITY ACROSS FOV SIZES ===
# =========================================================================================
# This script loads the aggregated tau data from 100 runs and calculates
# a quantitative "stability score" (Standard Deviation and IQR of ranks)
# for each slide.

# --- Load Libraries ---
library(dplyr)
library(readr)
library(tidyr)

# --- 1. DEFINE PATHS ---
# INPUT: The directory where your 100-run *analysis* is saved
input_dir <- "./balagan_consistency_analysis_from_raw"
input_file <- file.path(input_dir, "AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv")

# OUTPUT: The folder where your heatmap script saves plots
output_dir <- file.path("./stable_rank_analysis_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop(paste("Required input file not found:", input_file))
}

# --- 2. LOAD AND PREP DATA ---
cat(paste("Loading aggregated tau data from:", input_file, "\n"))
all_tau_data <- read_csv(input_file, show_col_types = FALSE)

prepped_tau_data <- all_tau_data %>%
  filter(!is.na(tau), tau > 0) %>%
  mutate(Slide_Name = paste0("slide", gsub(".*slide(\\d+).*", "\\1", Slide)))

# =========================================================================================
# === PART 1: MEAN-BASED STABILITY ANALYSIS ===
# =========================================================================================
cat("\n--- Analyzing MEAN-based rank stability ---\n")

# --- Calculate stable mean tau and ranks (same as heatmap script) ---
mean_rank_data <- prepped_tau_data %>%
  group_by(Slide_Name, FoV_width) %>%
  summarise(mean_tau = mean(tau, na.rm = TRUE), .groups = 'drop') %>%
  group_by(FoV_width) %>%
  mutate(tau_rank = rank(mean_tau, ties.method = "min")) %>%
  ungroup()

# --- NEW: Calculate stability metrics ---
mean_stability_table <- mean_rank_data %>%
  group_by(Slide_Name) %>%
  summarise(
    # The average rank (used for sorting)
    mean_rank = mean(tau_rank, na.rm = TRUE),
    
    # --- QUANTITATIVE STABILITY METRICS ---
    # A low SD means the rank is consistent across FOVs
    stability_sd = sd(tau_rank, na.rm = TRUE),
    # IQR is less sensitive to one or two outlier ranks
    stability_iqr = IQR(tau_rank, na.rm = TRUE)
  ) %>%
  # Order by stability (most stable at top)
  arrange(stability_sd)

# --- Save the new table ---
output_path_mean <- file.path(output_dir, "TABLE_stable_MEAN_rank_stability.csv")
write_csv(mean_stability_table, output_path_mean)

cat(paste("Saved MEAN rank stability table to:", output_path_mean, "\n"))
print(mean_stability_table)

# =========================================================================================
# === PART 2: MEDIAN-BASED STABILITY ANALYSIS ===
# =========================================================================================
cat("\n--- Analyzing MEDIAN-based rank stability ---\n")

# --- Calculate stable median tau and ranks ---
median_rank_data <- prepped_tau_data %>%
  group_by(Slide_Name, FoV_width) %>%
  summarise(median_tau = median(tau, na.rm = TRUE), .groups = 'drop') %>%
  group_by(FoV_width) %>%
  mutate(tau_rank = rank(median_tau, ties.method = "min")) %>%
  ungroup()

# --- NEW: Calculate stability metrics ---
median_stability_table <- median_rank_data %>%
  group_by(Slide_Name) %>%
  summarise(
    mean_rank = mean(tau_rank, na.rm = TRUE),
    stability_sd = sd(tau_rank, na.rm = TRUE),
    stability_iqr = IQR(tau_rank, na.rm = TRUE)
  ) %>%
  arrange(stability_sd)

# --- Save the new table ---
output_path_median <- file.path(output_dir, "TABLE_stable_MEDIAN_rank_stability.csv")
write_csv(median_stability_table, output_path_median)

cat(paste("Saved MEDIAN rank stability table to:", output_path_median, "\n"))
print(median_stability_table)

cat("\n====================================================\n")
cat("Finished calculating all stability metrics.\n")
cat("====================================================\n")