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

data_folder <- "/home/ubuntu/project/temp/Benchmarking_tissue_preparation_data"
out_folder <- "./"

### What data to load?
data_type <- "cX2" ### MESMER or cX2

### Use cX2 data
marker_names <- c(
  "DAPI", "CD3", "aSMA", "CD15", "CD4", "CD8",
  "CD11b", "CD11c", "CD20", "CD21", "H3K27me3",
  "Ki67", "HLADRA", "HistoneH3", "CD68", "DCSIGN",
  "Foxp3", "PD1", "CD163", "H3K27ac", "GranzymeB",
  "CD31", "CD206", "CD138", "NaKATPase", "CD45RA",
  "CD45", "Cytokeratin"
)

feature_names <- paste0(
  "mean_intensity:", marker_names, ":cell_region"
)

std_features <- c(
  "label:mask:cell_region",
  "pos_y:mask:cell_region",
  "pos_x:mask:cell_region",
  "area:mask:cell_region"
)

data <- load_cX2_data(all_of(c(std_features, feature_names)), 'data')

colnames(data) <- c(
  "cellLabel", "Y_cent", "X_cent", "cellSize",
  marker_names,
  "Staining_condition"
)


df_norm = normalize_data(data)
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
