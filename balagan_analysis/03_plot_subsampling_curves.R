# =========================================================================================
# === SCRIPT TO PLOT AVERAGE SUBSAMPLING CURVES (All Slides) ===
# =========================================================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(svglite)

# --- 1. CONFIGURATION ---
# INPUT: The base directory where all 100 run folders are
BASE_RESULTS_DIR <- "/Users/wang.13246/Documents/Project/Sizun_NM_revision/MESMER_workflow_03272025/balagan_results"

# OUTPUT: The folder to save your new plot
output_dir <- "./stable_rank_analysis_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. FIND ALL RAW FILES ---
cat("Finding all complex sampling files...\n")
all_sampling_files <- list.files(
  path = BASE_RESULTS_DIR,
  pattern = "_Complex_sampling\\.csv$",  # <-- Use the REAL file name pattern
  recursive = TRUE,
  full.names = TRUE
)

# 2. Filter the full paths to keep only the ones from your 100 runs
all_sampling_files <- all_sampling_files[
  grepl("out_balagan_analysis_BIDMC_run_\\d+", all_sampling_files)
]

# Now, all_sampling_files should contain your 2400 files (24 slides * 100 runs)
if (length(all_sampling_files) == 0) {
  stop("No 'run_X' sampling files found. Check your BASE_RESULTS_DIR path.")
}

# --- 3. LOAD AND COMBINE RAW DATA ---
cat("Loading raw sampling data (this may take a moment)...\n")

# Define parameters to re-add FOV size
fov_sizes <- seq(50, 500, 50)
n_sampling_regions <- 1:20
Parameter_table <- data.frame(
  Height = rep(fov_sizes, each = length(n_sampling_regions)),
  Width = rep(fov_sizes, each = length(n_sampling_regions))
)

load_raw_data <- function(file_path, param_table) {
  # Extract slide name
  slide_name <- basename(file_path)
  slide_name <- sub("_Complex_sampling.csv", "", slide_name)
  slide_name <- sub("dataScaleSize_", "", slide_name)
  
  sampling_data <- tryCatch(read.csv(file_path), error = function(e) NULL)
  
  if (is.null(sampling_data) || nrow(sampling_data) != nrow(param_table)) {
    return(NULL)
  }
  
  # Combine with parameters
  sampling_data %>%
    mutate(
      Slide = slide_name,
      FoV_width = param_table$Width
    )
}

# Load all data
all_raw_data <- map_dfr(all_sampling_files, ~ load_raw_data(.x, Parameter_table))

# --- 4. CALCULATE AVERAGE CURVES ---
cat("Calculating average curves...\n")
avg_curves <- all_raw_data %>%
  # Group by Slide, FOV, and subsampling number
  group_by(Slide, FoV_width, N_sampling) %>%
  
  # Calculate the mean and sd of clusters found across all 100 runs
  summarise(
    mean_cluster_count = mean(Mean_number_cluster, na.rm = TRUE),
    sd_cluster_count = sd(Mean_number_cluster, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 5. PLOT THE DATA ---
cat("Generating plot...\n")
p_subsampling <- ggplot(avg_curves, 
                        aes(x = N_sampling, 
                            y = mean_cluster_count, 
                            color = Slide)) +
  geom_line(linewidth = 0.8) +
  # Add error ribbons (optional, can be noisy)
  # geom_ribbon(aes(ymin = mean_cluster_count - sd_cluster_count,
  #                 ymax = mean_cluster_count + sd_cluster_count,
  #                 fill = Slide), alpha = 0.1, linetype = 0) +
  
  # Facet by FOV size to make it readable
  facet_wrap(~ FoV_width, scales = "free_y", labeller = label_bquote(FOV: .(FoV_width) * "µm")) +
  labs(
    title = "Average Cluster Discovery Curves (Stable Mean of 100 Runs)",
    subtitle = "Comparing subsampling efficiency across all slides and FOV sizes",
    x = "Number of Subsampling Regions (N)",
    y = "Mean Number of Clusters Recovered",
    color = "Slide"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )

# --- 6. SAVE THE PLOT ---
output_path <- file.path(output_dir, "PLOT_average_subsampling_curves.svg")
ggsave(output_path, plot = p_subsampling, width = 16, height = 12)

cat(paste("\nSuccessfully saved subsampling plot to:", output_path, "\n"))