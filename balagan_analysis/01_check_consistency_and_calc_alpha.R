# =========================================================================================
# === SCRIPT FOR ANALYZING CONSISTENCY ACROSS 100 BALAGAN RUNS (FROM RAW) ===
# =========================================================================================
# This script loads the *raw complex sampling files* from 100 separate runs,
# recalculates all tau values, and then generates plots to
# visualize the stability and consistency of the results.
#
# MODIFIED: This script now also loads CV scores to order the
# consistency boxplot by CV rank (best to worst).

# --- Load Libraries ---
library(dplyr)
library(readr)
library(purrr)   # For map_dfr
library(ggplot2)
library(svglite)
library(tidyr)   # For nest/unnest
library(viridis) # For heatmap colors
library(balagan) # For Visualize_complex_sampling

# =========================================================================================
# --- 1. CONFIGURATION ---
# =========================================================================================

# Set the base directory where all 100 run folders are located
# This should contain folders named: out_balagan_analysis_BIDMC_run_1, out_balagan_analysis_BIDMC_run_2, ...
BASE_RESULTS_DIR <- "./out_balagan_analysis"

# --- ADDED: Path to CV scores for ordering ---
# This points to the Mesmer workflow output directory with condition_summary.csv
CV_DIR <- "./out_Mesmer_BIDMC_all"

# A new directory to save this consistency analysis
OUTPUT_DIR <- "./balagan_consistency_analysis_from_raw"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define the FOV/Sampling parameters you used
# This MUST match the parameters from your master script
fov_sizes <- seq(50, 500, 50)
n_sampling_regions <- 1:20
height_vector <- rep(fov_sizes, each = length(n_sampling_regions))
width_vector <- rep(fov_sizes, each = length(n_sampling_regions))
Parameter_table <- data.frame(Height = height_vector, Width = width_vector)

cat("Starting consistency analysis for 100 runs from raw files...\n")
cat(paste("Scanning for results in:", BASE_RESULTS_DIR, "\n"))
cat(paste("Saving plots to:", OUTPUT_DIR, "\n\n"))
# =========================================================================================
# --- 2. FIND, LOAD, AND RECALCULATE TAU FROM ALL 100 RUNS (ROBUST) ---
# =========================================================================================

# Find all the raw complex sampling files across all run directories
all_sampling_files <- list.files(
  path = BASE_RESULTS_DIR,
  pattern = "_Complex_sampling\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(all_sampling_files) == 0) {
  stop(paste("No '_Complex_sampling.csv' files found in:", BASE_RESULTS_DIR))
}

cat(paste("Found", length(all_sampling_files), "total sampling files to process.\n"))


# --- NEW HELPER FUNCTION ---
# This safely fits the nls model for *one* group of data (i.e., one FOV)
safe_fit_nls <- function(data) {
  # This data is pre-filtered for one FOV size
  x <- data$N_sampling
  y <- data$Mean_number_cluster
  
  # Check for flat data (a common cause of nls failure)
  if (length(unique(y)) < 2) {
    # Return a specific reason for the failure
    return(data.frame(N = NA_real_, 
                      tau = NA_real_, 
                      R_squared = NA_real_, 
                      failure_reason = "Flat data: <2 unique cluster counts"))
  }
  
  tryCatch({
    # Use max(y) as a better starting guess for N
    expo_model <- nls(y ~ N * (1 - exp(-x/tau)), start = list(N = max(y), tau = 5))
    
    # Get results
    coeffs <- coef(expo_model)
    R_squared <- cor(predict(expo_model, newdata = x), y)^2
    
    return(data.frame(N = coeffs["N"], 
                      tau = coeffs["tau"], 
                      R_squared = R_squared, 
                      failure_reason = NA_character_)) # NA means success
    
  }, error = function(e) {
    # Return NAs and the specific error message
    return(data.frame(N = NA_real_, 
                      tau = NA_real_, 
                      R_squared = NA_real_, 
                      failure_reason = as.character(e$message)))
  })
}

# --- NEW ROBUST RECALCULATION FUNCTION ---
recalculate_tau_from_file_robust <- function(file_path, param_table) {
  
  # Extract run_id from the file path
  run_id_match <- regmatches(file_path, 
                             regexpr("out_balagan_analysis_BIDMC_run_(\\d+)", file_path))
  
  if (length(run_id_match) == 0) return(NULL) # Skip file if not in a run_X folder
  run_id <- as.integer(sub("out_balagan_analysis_BIDMC_run_", "", run_id_match[1]))
  
  # Extract slide_name
  slide_name <- basename(file_path)
  slide_name <- sub("_Complex_sampling.csv", "", slide_name)
  
  # Read data
  sampling_data <- tryCatch(read.csv(file_path), error = function(e) NULL)
  
  # Check for errors or mismatch
  if (is.null(sampling_data) || nrow(sampling_data) != nrow(param_table)) {
    warning(paste("File read error or row mismatch for:", slide_name, "run", run_id))
    return(NULL)
  }
  
  # Combine sampling data with its parameters
  full_data <- sampling_data %>%
    mutate(
      Height = param_table$Height,
      Width = param_table$Width
    )
  
  # Group by FOV width, then nest and fit the NLS curve for each group
  results_by_fov <- full_data %>%
    group_by(Width) %>%
    nest() %>%
    mutate(fit_results = map(data, safe_fit_nls)) %>%
    unnest(fit_results) %>%
    # --- MODIFIED: Added 'failure_reason' to the selection ---
    select(Width, N, tau, R_squared, failure_reason)
  
  # Rename 'Width' to 'FoV_width' and add metadata
  final_results <- results_by_fov %>%
    dplyr::rename(FoV_width = Width) %>%
    mutate(
      Slide = slide_name,
      run_id = run_id
    )
  
  return(final_results)
}

# --- MODIFIED CALL ---
# Use map_dfr to loop over all files, process them with the new robust function
cat("Recalculating tau for all files (this may take a few minutes)...\n")
all_runs_data <- map_dfr(all_sampling_files, ~ recalculate_tau_from_file_robust(.x, Parameter_table))
cat("Finished recalculating all tau values.\n")


if (nrow(all_runs_data) == 0) {
  stop("No tau data was successfully calculated. Check for errors.")
}

# Save the combined data
write_csv(all_runs_data, file.path(OUTPUT_DIR, "AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv"))
cat(paste("Successfully loaded and combined data from", length(unique(all_runs_data$run_id)), "runs.\n"))

# =========================================================================================
# --- 3. CALCULATE ALPHA SLOPE ---
# =========================================================================================
# Helper function to safely calculate slope
safe_lm <- function(data) {
  if(nrow(data) < 2) return(NA_real_)
  tryCatch({
    model <- lm(log10(tau) ~ log10(FoV_width), data = data)
    coef(model)[2]
  }, error = function(e) NA_real_)
}

cat("Calculating alpha slopes for all runs...\n")
# Calculate alpha_slope for every Slide in every run_id
alpha_slopes_all_runs <- all_runs_data %>%
  filter(!is.na(tau), tau > 0, !is.na(FoV_width), FoV_width > 0) %>%
  group_by(run_id, Slide) %>%
  nest() %>%
  mutate(alpha_slope = map_dbl(data, safe_lm)) %>%
  select(run_id, Slide, alpha_slope) %>%
  ungroup()

# Save the aggregated slopes
write_csv(alpha_slopes_all_runs, file.path(OUTPUT_DIR, "AGGREGATE_all_100_runs_alpha_slopes.csv"))

# =========================================================================================
# --- 3B. LOAD CV SCORES FOR PLOT ORDERING ---
# =========================================================================================
cat("Loading CV scores for plot ordering...\n")
cv_score_file <- file.path(CV_DIR, "condition_summary.csv")

# Load CV scores and create the 'Slide' key to match the Balagan 'Slide' column
cv_scores_for_ordering <- read_csv(cv_score_file, show_col_types = FALSE, col_select = c("Staining_condition", "Average_Score")) %>%
  mutate(Slide = paste0("slide", gsub("\\D", "", Staining_condition))) %>%
  select(Slide, Average_Score)

cat("CV scores loaded successfully.\n")

# =========================================================================================
# --- 4. GENERATE CONSISTENCY PLOTS ---
# =========================================================================================
cat("Generating consistency plots...\n")

# --- Plot 1: Success Rate Heatmap (Checks for failures) ---
# This shows how many runs (out of 100) successfully produced a tau value
fov_success_rate <- all_runs_data %>%
  group_by(Slide, FoV_width) %>%
  summarise(successful_runs = sum(!is.na(tau) & tau > 0), .groups = 'drop') %>%
  mutate(Slide = sub("dataScaleSize_", "", Slide)) # Clean up names for plot

p1 <- ggplot(fov_success_rate, aes(x = factor(FoV_width), y = Slide, fill = successful_runs)) +
  geom_tile(color = "white") +
  geom_text(aes(label = successful_runs), size = 3) +
  #scale_fill_viridis(direction = -1, limits = c(0, NA)) + # Auto-limits
  labs(
    title = "Analysis Success Rate (per Slide and FOV)",
    subtitle = "Number of successful runs that generated a valid tau value",
    x = "FOV Width (µm)",
    y = "Slide",
    fill = "Successful\nRuns"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

ggsave(file.path(OUTPUT_DIR, "plot_1_success_rate_heatmap.svg"), plot = p1, width = 10, height = 12,
       device = svglite::svglite, fix_text_size = FALSE)

# --- Plot 2: Alpha Slope Consistency (Your main consistency check) ---
# --- MODIFIED: Joined with CV scores ---
alpha_slopes_to_plot <- alpha_slopes_all_runs %>%
  filter(!is.na(alpha_slope)) %>%
  mutate(
    # First, remove "dataScaleSize_" prefix
    Slide = sub("dataScaleSize_", "", Slide), 
    # Second, remove "_FOV1" suffix to isolate the slide name
    Slide = sub("_FOV1", "", Slide), 
    # Convert alpha slope
    alpha_slope = alpha_slope * -1
  ) %>%
  # --- ADDED: Join CV scores for ordering ---
  left_join(cv_scores_for_ordering, by = "Slide") %>%
  filter(!is.na(Average_Score)) # Ensure we only plot slides that have a CV score

# --- MODIFIED: Changed 'reorder' to use Average_Score ---
p2 <- ggplot(alpha_slopes_to_plot, 
             # Order by Average_Score (ascending: best to worst)
             aes(x = reorder(Slide, Average_Score), y = alpha_slope, fill = Slide)) +
  geom_boxplot() +
  labs(
    title = "",
    subtitle = "",
    x = "Slide (Ordered by CV Score, Best-to-Worst)", # Updated X-axis label
    y = "Alpha"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "plot_2_alpha_slope_consistency_boxplot_CV_ORDERED.svg"), plot = p2, width = 14, height = 7,
       device = svglite::svglite, fix_text_size = FALSE)

# --- Plot 3: Stable Result (Your final answer) ---
# This plot will still be ordered by Alpha, which is fine.
# We join the CV scores just to ensure we are plotting the same set of slides as Plot 2.
stable_alpha_slopes <- alpha_slopes_to_plot %>%
  group_by(Slide, Average_Score) %>% # Keep Average_Score for potential use
  summarise(
    median_alpha_slope = median(alpha_slope, na.rm = TRUE),
    mean_alpha_slope = mean(alpha_slope, na.rm = TRUE),
    sd_alpha_slope = sd(alpha_slope, na.rm = TRUE),
    successful_runs = n(),
    .groups = 'drop'
  )

# Save this stable data
write_csv(stable_alpha_slopes, file.path(OUTPUT_DIR, "AGGREGATE_stable_alpha_slopes.csv"))

# Plot 3a: Ordered by Median Alpha (Original)
p3_alpha_order <- ggplot(stable_alpha_slopes, 
                         aes(x = reorder(Slide, median_alpha_slope), y = median_alpha_slope, fill = median_alpha_slope)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    aes(ymin = mean_alpha_slope - sd_alpha_slope, ymax = mean_alpha_slope + sd_alpha_slope),
    width = 0.25,
    linewidth = 0.5
  ) +
  geom_text(aes(label = successful_runs), vjust = -1, size = 3, color = "black") +
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Alpha Slope (Median of All Successful Runs)",
    subtitle = "Error bars show standard deviation. Numbers above bars = successful runs.",
    x = "Slide (Ordered by Median Alpha Slope)",
    y = "Median Alpha Slope"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "plot_3_alpha_slope_barchart_ALPHA_ORDERED.svg"), plot = p3_alpha_order, width = 14, height = 8,
       device = svglite::svglite, fix_text_size = FALSE)

# Plot 3b: Ordered by CV Score (New)
p3_cv_order <- ggplot(stable_alpha_slopes, 
                      # --- MODIFIED: Order by CV Score ---
                      aes(x = reorder(Slide, Average_Score), y = median_alpha_slope, fill = median_alpha_slope)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    aes(ymin = mean_alpha_slope - sd_alpha_slope, ymax = mean_alpha_slope + sd_alpha_slope),
    width = 0.25,
    linewidth = 0.5
  ) +
  geom_text(aes(label = round(Average_Score, 2)), vjust = 2, size = 3, color = "blue") + # Show CV score
  geom_text(aes(label = successful_runs), vjust = -1, size = 3, color = "black") + # Show run count
  scale_fill_gradient2(low = "darkgreen", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Alpha Slope (Median of All Successful Runs)",
    subtitle = "Error bars show standard deviation. Black numbers = successful runs. Blue numbers = CV score.",
    x = "Slide (Ordered by CV Score, Best-to-Worst)", # Updated X-axis label
    y = "Median Alpha Slope"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "plot_3_alpha_slope_barchart_CV_ORDERED.svg"), plot = p3_cv_order, width = 14, height = 8,
       device = svglite::svglite, fix_text_size = FALSE)


cat("\n====================================================\n")
cat("Consistency analysis complete.\n")
cat(paste("All plots saved to:", OUTPUT_DIR, "\n"))
cat("================================E====================\n")

# Display failure reasons table
cat("Summary of NLS fitting failures (if any):\n")
print(table(all_runs_data$failure_reason))