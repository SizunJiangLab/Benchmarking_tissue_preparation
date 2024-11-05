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

source("helper.R")
########################################## Configuration Begin  #############################################################################

### What data to load?
data_type <- "MESMER"

# Step 1: Define configuration for different data sets, and add new configuration in this list:
configurations <- list(
  MESMER_09202024 = list(
    data_folder = "./data_09202024/",
    out_folder = "./out_09202024/",
    input_filenames = c("dataScaleSize_fresh.csv", "dataScaleSize_TBS.csv", "dataScaleSize_water.csv"),
    input_note = c("Fresh", "TBS", "Water"),
    pairs = list(c("TBS", "Water"), c("Fresh", "Water"), c("TBS", "Fresh")),
    excluded_values = list()
  ),
  MESMER_10032024_ROI1 = list(
    data_folder = "./data_10032024/ROI1/",
    out_folder = "./out_10032024/ROI1/",
    input_filenames = c("dataScaleSize_fresh_ROI1.csv", "dataScaleSize_lyophilized_ROI1.csv"),
    input_note = c("Fresh", "Lyophilized"),
    pairs = list(c("Fresh", "Lyophilized")),
    excluded_values = list(
      'Fresh' = c("IDO1"),
      'Lyophilized' = c("IDO1")
    )
  ),
  MESMER_10032024_ROI2 = list(
    data_folder = "./data_10032024/ROI2/",
    out_folder = "./out_10032024/ROI2/",
    input_filenames = c("dataScaleSize_fresh_ROI2.csv", "dataScaleSize_lyophilized_ROI2.csv"),
    input_note = c("Fresh", "Lyophilized"),
    pairs = list(c("Fresh", "Lyophilized")),
    excluded_values = list(
      'Fresh' = c("IDO1"),
      'Lyophilized' = c("IDO1")
    )
  ),
  MESMER_11042024 = list(
    data_folder = "./drive-download-20241104T193831Z-001/",
    out_folder = "./out_11042024/",
    input_filenames = c(
      "spleen_Slide1_ER2_3h_Core1.csv",
      "spleen_Slide1_ER2_3h_Core2.csv",
      "spleen_Slide1_ER2_3h_Core3.csv",
      "spleen_Slide2_ER1_60_Core1.csv",
      "spleen_Slide2_ER1_60_Core2.csv",
      "spleen_Slide2_ER1_60_Core3.csv",
      "spleen_Slide3_ER2_60_Core1.csv",
      "spleen_Slide3_ER2_60_Core2.csv",
      "spleen_Slide3_ER2_60_Core3.csv",
      "spleen_Slide4_ER1_3h_Core1.csv",
      "spleen_Slide4_ER1_3h_Core2.csv",
      "spleen_Slide4_ER1_3h_Core3.csv",
      "spleen_Slide5_CC1_3h_Core1.csv",
      "spleen_Slide5_CC1_3h_Core2.csv",
      "spleen_Slide5_CC1_3h_Core3.csv",
      "spleen_Slide6_CC1_60_Core1.csv",
      "spleen_Slide6_CC1_60_Core2.csv",
      "spleen_Slide6_CC1_60_Core3.csv",
      "spleen_Slide7_CC2_3h_Core1.csv",
      "spleen_Slide7_CC2_3h_Core2.csv",
      "spleen_Slide7_CC2_3h_Core3.csv",
      "spleen_Slide8_CC2_60h_Core1.csv",
      "spleen_Slide8_CC2_60h_Core2.csv",
      "spleen_Slide8_CC2_60h_Core3.csv"
    ),
    input_note = c(
      "ER2_3h_Core1",
      "ER2_3h_Core2",
      "ER2_3h_Core3",
      "ER1_60_Core1",
      "ER1_60_Core2",
      "ER1_60_Core3",
      "ER2_60_Core1",
      "ER2_60_Core2",
      "ER2_60_Core3",
      "ER1_3h_Core1",
      "ER1_3h_Core2",
      "ER1_3h_Core3",
      "CC1_3h_Core1",
      "CC1_3h_Core2",
      "CC1_3h_Core3",
      "CC1_60_Core1",
      "CC1_60_Core2",
      "CC1_60_Core3",
      "CC2_3h_Core1",
      "CC2_3h_Core2",
      "CC2_3h_Core3",
      "CC2_60h_Core1",
      "CC2_60h_Core2",
      "CC2_60h_Core3"
    ),
    pairs = list(c("CC2_60h_Core1", "CC2_60h_Core2"), c("ER2_3h_Core1", "ER2_3h_Core2")),
    excluded_values = list()
  ),
  MESMER_11052024 = list(
    data_folder = "./grouping_output_11052024/",
    out_folder = "./out_11052024/",
    input_filenames = c(
      "CC1_3h_combined.csv",
      "CC1_60_combined.csv",
      "CC2_3h_combined.csv",
      "CC2_60h_combined.csv",
      "ER1_3h_combined.csv",
      "ER1_60_combined.csv",
      "ER2_3h_combined.csv",
      "ER2_60_combined.csv"
    ),
    input_note = c(
      "CC1_3h",
      "CC1_60",
      "CC2_3h",
      "CC2_60",
      "ER1_3h",
      "ER1_60",
      "ER2_3h",
      "ER2_60"
    ),
    pairs = list(
      c("CC1_3h", "CC1_60"),
      c("CC2_3h", "CC2_60")
    ),
    excluded_values = list(
    )
  )
)

# Step 2: Choose which configuration to use
current_config <- configurations$MESMER_11052024

########################################## Configuration End  #############################################################################

data_folder <- current_config$data_folder
out_folder <- current_config$out_folder
input_filenames <- current_config$input_filenames
input_note <- current_config$input_note
pairs <- current_config$pairs
excluded_values <- current_config$excluded_values

dir.create(out_folder, showWarnings = F)

### Load MESMER data
data <- load_mesmer_data(data_folder, input_filenames, input_note)

old_marker_names <- data %>%
  rename_all(~ gsub("[.]", "-", .)) %>%
  select(5:ncol(data)) %>%
  select(-matches("Staining_condition")) %>%
  rename(
    "Na/K-ATPase" = any_of("NaKATPase"),
    "Granzyme B" = any_of("Granzyme-B"),
    "Hoechst" = any_of("DAPI")
  ) %>%
  names()

# Remove special characters
data <- data %>%
  rename_all(~ gsub("[.-]", "", .)) %>%
  rename("Hoechst" = "DAPI")

marker_names <- data %>%
  select(5:ncol(data)) %>%
  select(-Staining_condition) %>%
  names()

result = normalize_data(data)
df_norm <- result$data
marker_names <- result$marker_names

### Arcsinh (Inverse hyperbolic Sine) Transformation
df_arcsinh <- df_norm %>% 
  mutate(across(all_of(marker_names), ~asinh(.x/0.001)))

#Universal percentile normalization 
df_trans <- df_arcsinh

rng <- colQuantiles(as.matrix(df_arcsinh[,marker_names]), probs = c(0.001, 0.999))

expr <- t((t(as.matrix(df_arcsinh[,marker_names]))-rng[,1]) / (rng[,2]-rng[,1]))

expr[expr < 0] <- 0

expr[expr > 1] <- 1

df_trans[, marker_names] <- expr


### Plot density plots for each staining condition and marker
p <- plot_density_plots2(
  df_trans,
  marker_names,
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.justification = "left",
  legend.rows = 6
)

ggsave(
  p, 
  filename = 
    paste0(out_folder, "Arcsinh transformed Hoechst normalised density plots (", data_type,").svg"), width = 8, height = 11
)

### Perform statistical and Kruskal-Wallis tests
kruskal_results <- perform_statistical_and_kruskal_wallis_tests(df_trans, marker_names)

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, paste0(out_folder, "kruskal_pvals.csv"))

# plot heatmap
heatmap2 <- plot_heatmap2(data = df_trans, out_folder = out_folder, pairs = pairs, excluded_values = excluded_values)

