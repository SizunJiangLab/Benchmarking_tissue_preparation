# =========================================================================================
# === SCRIPT TO CREATE STABLE RANK HEATMAPS (Mean & Median, with/without CV) ===
# =========================================================================================
# This script uses the aggregated results from 100 runs to create stable,
# mean-based and median-based tau ranks and generates heatmaps:
#   - Standalone heatmaps (ordered by Tau rank)
#   - Combined heatmaps with CV rank (ordered by Tau rank)
#
# --- Load Libraries ---
library(dplyr)
library(readr)
library(ggplot2)
library(svglite)
library(tidyr)
library(patchwork) # For combining plots

# =========================================================================================
# --- 1. CONFIGURATION ---
# =========================================================================================

# --- INPUT: The directory where your 100-run *analysis* is saved ---
input_dir <- "./balagan_consistency_analysis_from_raw"
input_file <- file.path(input_dir, "AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv")

# --- INPUT: The path to your manual CV rank file (optional, for combined plots) ---
cv_rank_file <- "./data_mesmer/condition_summary.csv"

# --- OUTPUT: A new folder for these stable rank plots ---
output_dir <- file.path("./stable_rank_analysis_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Configuration: Set to TRUE to include CV rank plots ---
INCLUDE_CV_PLOTS <- file.exists(cv_rank_file)

# --- Check for input file ---
if (!file.exists(input_file)) {
  stop(paste("Required input file not found:", input_file))
}

cat(paste("Loading aggregated tau data from:", input_file, "\n"))
cat(paste("Include CV plots:", INCLUDE_CV_PLOTS, "\n"))

# =========================================================================================
# --- 2. LOAD DATA ---
# =========================================================================================

all_tau_data <- read_csv(input_file, show_col_types = FALSE)

if (nrow(all_tau_data) == 0) {
  stop("Loaded data file is empty.")
}

# Filter and prep base data once
prepped_tau_data <- all_tau_data %>%
  filter(!is.na(tau), tau > 0) %>%
  mutate(Slide_Name = paste0("slide", gsub(".*slide(\\d+).*", "\\1", Slide)))

# Load CV data if requested
cv_rank_data <- NULL
if (INCLUDE_CV_PLOTS) {
  condition_summary <- read.csv(cv_rank_file)
  staining_condition_order <- condition_summary$Staining_condition
  
  cv_rank_data <- data.frame(
    Staining_condition = staining_condition_order,
    cv_rank = 1:length(staining_condition_order)
  ) %>%
    mutate(Slide_Name = paste0("slide", gsub("\\D", "", Staining_condition)))
  
  cat("CV rank data loaded successfully.\n")
}

# =========================================================================================
# === PART 1: MEAN-BASED ANALYSIS ===
# =========================================================================================
cat("\n--- Starting MEAN Analysis ---\n")

# --- Calculate stable mean tau and ranks ---
stable_mean_tau_data <- prepped_tau_data %>%
  group_by(Slide_Name, FoV_width) %>%
  summarise(mean_tau = mean(tau, na.rm = TRUE), .groups = 'drop')

mean_rank_data <- stable_mean_tau_data %>%
  group_by(FoV_width) %>%
  mutate(tau_rank = rank(mean_tau, ties.method = "min")) %>%
  ungroup()

# --- Prepare data for plotting ---
mean_slide_order_info <- mean_rank_data %>%
  group_by(Slide_Name) %>%
  summarise(mean_rank_of_ranks = mean(tau_rank, na.rm = TRUE), .groups = 'drop')

mean_heatmap_data <- mean_rank_data %>%
  tidyr::complete(Slide_Name, FoV_width) %>%
  left_join(mean_slide_order_info, by = "Slide_Name") %>%
  mutate(Slide_Name_ordered = reorder(Slide_Name, -mean_rank_of_ranks))

# --- 1A. STANDALONE MEAN TAU RANK HEATMAP ---
p_tau_rank_mean <- ggplot(mean_heatmap_data, aes(x = factor(FoV_width), y = Slide_Name_ordered, fill = tau_rank)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient(
    low = "#fde725", high = "#0571b0", na.value = "grey85",
    name = "Stable\n(Mean) Tau Rank"
  ) +
  geom_text(aes(label = tau_rank), color = "black", size = 3.5) +
  labs(
    x = "FOV Width (µm)",
    y = "Slide (Ordered by Stable Mean Rank)"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(hjust = 1, size = 9),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

plot_mean_final <- p_tau_rank_mean +
  plot_annotation(
    title = 'Stable Mean Tau Rank Across FOVs (Mean of 100 Runs)',
    subtitle = 'Slides ordered by mean Tau Rank (best overall at top).',
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

output_path_mean <- file.path(output_dir, "AGGREGATE_stable_MEAN_tau_rank_heatmap.svg")
ggsave(output_path_mean, plot = plot_mean_final, width = 12, height = 10,
       device = svglite::svglite, fix_text_size = FALSE)
cat(paste("MEAN standalone plot saved to:", output_path_mean, "\n"))

# --- 1B. COMBINED MEAN TAU + CV RANK HEATMAP ---
if (INCLUDE_CV_PLOTS) {
  mean_heatmap_data_with_cv <- mean_rank_data %>%
    tidyr::complete(Slide_Name, FoV_width) %>%
    left_join(mean_slide_order_info, by = "Slide_Name") %>%
    left_join(cv_rank_data, by = "Slide_Name") %>%
    mutate(Slide_Name_ordered = reorder(Slide_Name, -mean_rank_of_ranks))
  
  # CV Rank Column (Left)
  p_cv_rank_mean <- ggplot(mean_heatmap_data_with_cv, aes(x = "CV Rank", y = Slide_Name_ordered, fill = cv_rank)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "#e0f3f8", high = "#4575b4", na.value = "grey85", name = "CV Rank") +
    geom_text(aes(label = cv_rank), color = "black", size = 4, fontface = "bold") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "Slide (Ordered by Stable Mean Rank)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(hjust = 1, size = 9),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # Tau Rank Heatmap (Right)
  p_tau_rank_mean_right <- ggplot(mean_heatmap_data_with_cv, aes(x = factor(FoV_width), y = Slide_Name_ordered, fill = tau_rank)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(
      low = "#fde725", high = "#0571b0", na.value = "grey85",
      name = "Stable\n(Mean) Tau Rank"
    ) +
    geom_text(aes(label = tau_rank), color = "black", size = 3.5) +
    labs(x = "FOV Width (µm)", y = NULL) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  combined_plot_mean <- p_cv_rank_mean + p_tau_rank_mean_right +
    plot_layout(widths = c(2, 10)) +
    plot_annotation(
      title = 'Stable Comparison of Slide Quality Metrics (Mean of 100 Runs)',
      subtitle = 'Manual CV Rank (left) vs. Mean Tau Rank (right). Slides ordered by mean Tau Rank.',
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  
  output_path_mean_cv <- file.path(output_dir, "AGGREGATE_stable_MEAN_cv_tau_rank_heatmap.svg")
  ggsave(output_path_mean_cv, plot = combined_plot_mean, width = 14, height = 10,
         device = svglite::svglite, fix_text_size = FALSE)
  cat(paste("MEAN combined plot (with CV) saved to:", output_path_mean_cv, "\n"))
}

# =========================================================================================
# === PART 2: MEDIAN-BASED ANALYSIS ===
# =========================================================================================
cat("\n--- Starting MEDIAN Analysis ---\n")

# --- Calculate stable median tau and ranks ---
stable_median_tau_data <- prepped_tau_data %>%
  group_by(Slide_Name, FoV_width) %>%
  summarise(median_tau = median(tau, na.rm = TRUE), .groups = 'drop')

median_rank_data <- stable_median_tau_data %>%
  group_by(FoV_width) %>%
  mutate(tau_rank = rank(median_tau, ties.method = "min")) %>%
  ungroup()

# --- Prepare data for plotting ---
median_slide_order_info <- median_rank_data %>%
  group_by(Slide_Name) %>%
  summarise(median_rank_of_ranks = mean(tau_rank, na.rm = TRUE), .groups = 'drop')

median_heatmap_data <- median_rank_data %>%
  tidyr::complete(Slide_Name, FoV_width) %>%
  left_join(median_slide_order_info, by = "Slide_Name") %>%
  mutate(Slide_Name_ordered = reorder(Slide_Name, -median_rank_of_ranks))

# --- 2A. STANDALONE MEDIAN TAU RANK HEATMAP ---
p_tau_rank_median <- ggplot(median_heatmap_data, aes(x = factor(FoV_width), y = Slide_Name_ordered, fill = tau_rank)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient(
    low = "#fde725", high = "#0571b0", na.value = "grey85",
    name = "Stable\n(Median) Tau Rank"
  ) +
  geom_text(aes(label = tau_rank), color = "black", size = 3.5) +
  labs(
    x = "FOV Width (µm)",
    y = "Slide (Ordered by Stable Median Rank)"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(hjust = 1, size = 9),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

plot_median_final <- p_tau_rank_median +
  plot_annotation(
    title = 'Stable Median Tau Rank Across FOVs (Median of 100 Runs)',
    subtitle = 'Slides ordered by median Tau Rank (best overall at top).',
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

output_path_median <- file.path(output_dir, "AGGREGATE_stable_MEDIAN_tau_rank_heatmap.svg")
ggsave(output_path_median, plot = plot_median_final, width = 12, height = 10,
       device = svglite::svglite, fix_text_size = FALSE)
cat(paste("MEDIAN standalone plot saved to:", output_path_median, "\n"))

# --- 2B. COMBINED MEDIAN TAU + CV RANK HEATMAP ---
if (INCLUDE_CV_PLOTS) {
  median_heatmap_data_with_cv <- median_rank_data %>%
    tidyr::complete(Slide_Name, FoV_width) %>%
    left_join(median_slide_order_info, by = "Slide_Name") %>%
    left_join(cv_rank_data, by = "Slide_Name") %>%
    mutate(Slide_Name_ordered = reorder(Slide_Name, -median_rank_of_ranks))
  
  # CV Rank Column (Left)
  p_cv_rank_median <- ggplot(median_heatmap_data_with_cv, aes(x = "CV Rank", y = Slide_Name_ordered, fill = cv_rank)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "#e0f3f8", high = "#4575b4", na.value = "grey85", name = "CV Rank") +
    geom_text(aes(label = cv_rank), color = "black", size = 4, fontface = "bold") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "Slide (Ordered by Median Tau Rank)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(hjust = 1, size = 9),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # Tau Rank Heatmap (Right)
  p_tau_rank_median_right <- ggplot(median_heatmap_data_with_cv, aes(x = factor(FoV_width), y = Slide_Name_ordered, fill = tau_rank)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(
      low = "#fde725", high = "#0571b0", na.value = "grey85",
      name = "Median Tau Rank"
    ) +
    geom_text(aes(label = tau_rank), color = "black", size = 3.5) +
    labs(x = "FOV Width (µm)", y = NULL) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  combined_plot_median <- p_cv_rank_median + p_tau_rank_median_right +
    plot_layout(widths = c(2, 10)) +
    plot_annotation(
      title = 'Comparison of Slide Quality Metrics (Median of 100 Runs)',
      subtitle = 'Manual CV Rank (left) vs. Median Tau Rank (right). Slides ordered by median Tau Rank.',
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  
  output_path_median_cv <- file.path(output_dir, "AGGREGATE_stable_MEDIAN_cv_tau_rank_heatmap.svg")
  ggsave(output_path_median_cv, plot = combined_plot_median, width = 14, height = 10,
         device = svglite::svglite, fix_text_size = FALSE)
  cat(paste("MEDIAN combined plot (with CV) saved to:", output_path_median_cv, "\n"))
}

cat("\n====================================================\n")
cat("Finished generating all stable rank heatmaps.\n")
cat("====================================================\n")

