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



### What data to load?
data_type <- "MESMER"

## Data1: 09202024
#data_folder <- "./data_09202024/"
#out_folder <- "./out_09202024/"
#input_filenames <- c("dataScaleSize_fresh.csv", "dataScaleSize_TBS.csv", "dataScaleSize_water.csv")
#input_note <- c("Fresh", "TBS", "Water")
#pairs <- list(
#  c("TBS", "Water"),
#  c("Fresh", "Water"),
#  c("TBS", "Fresh")
#)
#excluded_values <- list() # No excluded values for 09202024

## Data2: 10032024-ROI1
#data_folder <- "./data_10032024/ROI1/"
#out_folder <- "./out_10032024/ROI1/"
#input_filenames <- c("dataScaleSize_fresh_ROI1.csv", "dataScaleSize_lyophilized_ROI1.csv")
#input_note <- c("Fresh", "Lyophilized")
#pairs <- list(
#  c("Fresh", "Lyophilized")
#)
#excluded_values <- list(
#  'Fresh' = c("IDO1"),  
#  'Lyophilized' = c("IDO1")
#)

# Data3: 10032024-ROI2
data_folder <- "./data_10032024/ROI2/"
out_folder <- "./out_10032024/ROI2/"
input_filenames <- c("dataScaleSize_fresh_ROI2.csv", "dataScaleSize_lyophilized_ROI2.csv")
input_note <- c("Fresh", "Lyophilized")
pairs <- list(
  c("Fresh", "Lyophilized")
)
excluded_values <- list(
  'Fresh' = c("IDO1"),  
  'Lyophilized' = c("IDO1") 
)

dir.create(out_folder, showWarnings = F)
### Use MESMER data
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
p <- plot_density_plots2(df_trans, marker_names)

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

#Save the plot as svg

