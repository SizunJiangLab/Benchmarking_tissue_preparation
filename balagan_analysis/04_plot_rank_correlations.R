# =========================================================================================
# === SCRIPT TO PLOT RANK CORRELATIONS (CV vs. Alpha & Tau) ===
# =========================================================================================
# This script loads all necessary data to create a combined rank table
# and generates scatter plots to quantitatively assess the correlation
# between the CV rank and the two spatial metrics.

# --- Load Libraries ---
library(dplyr)
library(readr)
library(ggplot2)
library(svglite)
library(ggpubr) # For adding correlation stats to plots
library(tidyr)  # For nest/unnest

# =========================================================================================
# --- 1. DEFINE PATHS ---
# =========================================================================================

# --- INPUT 1: The aggregated data from all 100 runs ---
tau_data_file <- "./balagan_consistency_analysis_from_raw/AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv"

# --- INPUT 2: The stable Alpha summary ---
alpha_data_file <- "./balagan_consistency_analysis_from_raw/AGGREGATE_stable_alpha_slopes.csv"

# --- INPUT 3: The manual CV rank file ---
cv_rank_file <- "./data_mesmer/condition_summary.csv" 

# --- OUTPUT: Folder for new correlation plots ---
output_dir <- "./stable_rank_analysis_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Check for input files ---
if (!file.exists(tau_data_file)) stop(paste("File not found:", tau_data_file))
if (!file.exists(alpha_data_file)) stop(paste("File not found:", alpha_data_file))
if (!file.exists(cv_rank_file)) stop(paste("File not found:", cv_rank_file))

# =========================================================================================
# --- 2. LOAD & PREPARE ALL RANK DATA ---
# =========================================================================================

# --- Process 1: CV Rank ---
cat("Loading CV Rank data...\n")
condition_summary <- read.csv(cv_rank_file)
cv_rank_data <- data.frame(
  Staining_condition = condition_summary$Staining_condition,
  cv_rank = 1:length(condition_summary$Staining_condition)
) %>%
  mutate(Slide_Name = paste0("slide", gsub("\\D", "", Staining_condition)))

# --- Process 2: Alpha Rank ---
cat("Loading Stable Alpha data...\n")
alpha_rank_data <- read_csv(alpha_data_file, show_col_types = FALSE) %>%
  # Create the same "slideX" key for joining
  mutate(Slide_Name = paste0("slide", gsub(".*slide(\\d+).*", "\\1", Slide)), median_alpha_slope=median_alpha_slope*-1) %>%
  # Rank the slides based on their slope
  # A lower (more negative) slope is ranked first (Rank 1)
  mutate(alpha_slope_rank = rank(median_alpha_slope, ties.method = "min")) %>%
  select(Slide_Name, alpha_slope_rank, median_alpha_slope)

# --- Process 3: Stable Tau Rank (Average Tau Rank) ---
cat("Loading and processing stable tau ranks...\n")
all_tau_data <- read_csv(tau_data_file, show_col_types = FALSE)

avg_tau_rank_data <- all_tau_data %>%
  filter(!is.na(tau), tau > 0) %>%
  mutate(Slide_Name = paste0("slide", gsub(".*slide(\\d+).*", "\\1", Slide))) %>%
  group_by(Slide_Name, FoV_width) %>%
  summarise(mean_tau = mean(tau, na.rm = TRUE), .groups = 'drop') %>%
  group_by(FoV_width) %>%
  mutate(tau_rank = rank(mean_tau, ties.method = "min")) %>%
  ungroup() %>%
  group_by(Slide_Name) %>%
  summarise(avg_tau_rank = mean(tau_rank, na.rm = TRUE), .groups = 'drop')

# --- 4. COMBINE ALL RANKS INTO ONE TABLE ---
cat("Combining all ranks...\n")
combined_ranks <- cv_rank_data %>%
  left_join(alpha_rank_data, by = "Slide_Name") %>%
  left_join(avg_tau_rank_data, by = "Slide_Name") %>%
  # Remove any slides that might be missing data
  filter(!is.na(alpha_slope_rank), !is.na(avg_tau_rank))

# Save this combined table
write_csv(combined_ranks, file.path(output_dir, "TABLE_all_combined_ranks.csv"))

# =========================================================================================
# --- 5. PLOT THE CORRELATIONS ---
# =========================================================================================

# --- Plot 1: CV Rank vs. Alpha Rank ---
cat("Generating Plot 1: CV Rank vs. Alpha Rank\n")
p_corr_alpha <- ggplot(combined_ranks, aes(x = cv_rank, y = alpha_slope_rank)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", linetype = "dashed") +
  # Use ggpubr to add Spearman correlation coefficient
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Correlation: CV Rank vs. Alpha Rank",
    subtitle = "Spearman's rank correlation (Rs) shown. R near 0 = no correlation.",
    x = "CV Rank",
    y = "Alpha Rank"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "PLOT_correlation_cv_vs_alpha.svg"), 
       plot = p_corr_alpha, width = 8, height = 7,
       device = svglite::svglite, fix_text_size = FALSE)

# --- Plot 2: CV Rank vs. Stable Tau Rank ---
cat("Generating Plot 2: CV Rank vs. Stable Tau Rank\n")
p_corr_tau <- ggplot(combined_ranks, aes(x = cv_rank, y = avg_tau_rank)) +
  geom_point(size = 3, color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Correlation: CV Rank vs. Average Tau Rank",
    subtitle = "Spearman's rank correlation (Rs) shown. R near 0 = no correlation.",
    x = "CV Rank",
    y = "Average Tau Rank"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "PLOT_correlation_cv_vs_tau.svg"), 
       plot = p_corr_tau, width = 8, height = 7,
       device = svglite::svglite, fix_text_size = FALSE)

# --- NEW: Plot 3: Alpha Rank vs. Stable Tau Rank ---
cat("Generating Plot 3: Alpha Rank vs. Stable Tau Rank\n")
p_corr_alpha_vs_tau <- ggplot(combined_ranks, aes(x = avg_tau_rank, y = alpha_slope_rank)) +
  geom_point(size = 3, color = "purple") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgrey", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Correlation: Alpha Rank vs. Stable Tau Rank",
    subtitle = "Spearman's rank correlation (Rs) shown. R near 1 = strong correlation.",
    x = "Average Tau Rank",
    y = "Alpha Rank"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "PLOT_correlation_alpha_vs_tau.svg"), 
       plot = p_corr_alpha_vs_tau, width = 8, height = 7,
       device = svglite::svglite, fix_text_size = FALSE)
# --- END NEW PLOT ---

cat("\n====================================================\n")
cat("Correlation analysis complete. All 3 plots saved.\n")
cat(paste("Output directory:", output_dir, "\n"))
cat("====================================================\n")