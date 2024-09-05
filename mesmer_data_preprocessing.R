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

source("R-scripts/helper.R")

### Remove all special characters
data_colnames <- c(
  "cellLabel", "Y_cent", "X_cent", "cellSize",
  "Hoechst", "CD3",  "aSMA", "CD15", "CD4", "CD8",
  "CD11b", "CD11c", "CD20", "CD21", "H3K27me3",
  "Ki67", "HLADRA", "HistoneH3", "CD68", "DCSIGN",
  "Foxp3", "PD1", "CD163", "H3K27ac", "GranzymeB",
  "CD31", "CD206", "CD138", "NaKATPase", "CD45RA",
  "CD45", "Cytokeratin"
)

data_folder <- "./data/"
out_folder <- "./out/"

### What data to load?
data_type <- "MESMER" ### MESMER or cX2

### Use MESMER data
data <- load_mesmer_data(data_colnames, data_folder)
marker_names <- colnames(data[, 6:33])

result = normalize_data(data)
df_norm <- result$data
marker_names <- result$marker_names

### Arcsinh (Inverse hyperbolic Sine) Transformation
df_arcsinh <- df_norm %>% 
  mutate(across(all_of(marker_names), ~asinh(.x/0.001)))

#Universal percentile normalization 
df_trans <- df_arcsinh

rng <- colQuantiles(as.matrix(df_arcsinh[,5:28]), probs = c(0.001, 0.999))

expr <- t((t(as.matrix(df_arcsinh[,5:28]))-rng[,1]) / (rng[,2]-rng[,1]))

expr[expr < 0] <- 0

expr[expr > 1] <- 1

df_trans[,5:28] <- expr


### Plot density plots for each staining condition and marker
p <- plot_density_plots(df_trans, marker_names)

ggsave(
  p, 
  filename = 
    paste0(out_folder, "Arcsinh transformed Hoechst normalised density plots (", data_type,").svg"), width = 8, height = 10
)

### Perform statistical and Kruskal-Wallis tests
kruskal_results <- perform_statistical_and_kruskal_wallis_tests(df_trans, marker_names)

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, paste0(out_folder, "kruskal_pvals.csv"))

heatmap <- plot_heatmap(df_trans)

#Save the plot as svg
ggsave(
  filename = paste0(out_folder, "Heatmap of arcsinh transform normalised data Z scores of CV_with excluded values.svg"), 
  width = 10, 
  height = 8
)