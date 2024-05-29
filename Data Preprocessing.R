library(dplyr)
library(tidyverse)
library(matrixStats)
library(ggcorrplot)
library(ggpubr)
library(tidyr)
library(rstatix)
library(readr)

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

### What data to load?
data_type <- "MESMER" ### MESMER or cX2

if (data_type == "MESMER") {
    ### Use MESMER data
    data_slide1_FOV1 <- read.csv("dataScaleSize_slide1_FOV1.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide1_FOV2 <- read.csv("dataScaleSize_slide1_FOV2.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide1 <- rbind(data_slide1_FOV1, data_slide1_FOV2)
    data_slide1 <- mutate(data_slide1, Staining_condition = "20 min HIER, 1h 37C stain")
    head(data_slide1)

    data_slide2_FOV1 <- read.csv("dataScaleSize_slide2_FOV1.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide2_FOV2 <- read.csv("dataScaleSize_slide2_FOV2.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide2 <- rbind(data_slide2_FOV1, data_slide2_FOV2)
    data_slide2 <- mutate(data_slide2, Staining_condition = "20 min HIER, ON 4C stain")
    head(data_slide2)

    data_slide3_FOV1 <- read.csv("dataScaleSize_slide3_FOV1.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide3_FOV2 <- read.csv("dataScaleSize_slide3_FOV2.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide3 <- rbind(data_slide3_FOV1, data_slide3_FOV2)
    data_slide3 <- mutate(data_slide3, Staining_condition = "10 min HIER, ON 4C stain")
    head(data_slide3)

    data_slide4_FOV1 <- read.csv("dataScaleSize_slide4_FOV1.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide4_FOV2 <- read.csv("dataScaleSize_slide4_FOV2.csv",
        col.names = data_colnames, header = TRUE
    )
    data_slide4 <- rbind(data_slide4_FOV1, data_slide4_FOV2)
    data_slide4 <- mutate(data_slide4, Staining_condition = "40 min HIER, ON 4C stain")
    head(data_slide4)


    data <- rbind(data_slide1, data_slide2, data_slide3, data_slide4)
    str(data)

    marker_names <- colnames(data[, 6:33])

} else {
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

    data_slide1 <- read_csv(
        file = "data/cX2_Slide1_rawdata.csv",
        col_sel = all_of(c(std_features, feature_names)),
        show_col_types = FALSE
    ) |>
    mutate(
          Staining_condition = "20 min HIER, 1h 37C stain"
    )

    data_slide2 <- read_csv(
        file = "data/cX2_Slide2_rawdata.csv",
        col_sel = all_of(c(std_features, feature_names)),
        show_col_types = FALSE
    ) |>
    mutate(
          Staining_condition = "20 min HIER, ON 4C stain"
    )

    data_slide3 <- read_csv(
        file = "data/cX2_Slide3_rawdata.csv",
        col_sel = all_of(c(std_features, feature_names)),
        show_col_types = FALSE
    ) |>
    mutate(
          Staining_condition = "10 min HIER, ON 4C stain"
    )

    data_slide4 <- read_csv(
        file = "data/cX2_Slide4_rawdata.csv",
        col_sel = all_of(c(std_features, feature_names)),
        show_col_types = FALSE
    ) |>
    mutate(
          Staining_condition = "40 min HIER, ON 4C stain"
    )

    data <- bind_rows(data_slide1, data_slide2, data_slide3, data_slide4)

    colnames(data) <- c(
        "cellLabel", "Y_cent", "X_cent", "cellSize",
        marker_names,
        "Staining_condition"
    )
}

# reorder columns
data$Staining_condition <- factor(
    data$Staining_condition,
    levels = c(
    "20 min HIER, 1h 37C stain",
    "20 min HIER, ON 4C stain",
    "10 min HIER, ON 4C stain",
    "40 min HIER, ON 4C stain"
    )
)
head(data)

# Remove markers shown to be not working for all conditions
data <- data %>% 
  dplyr::select(-'aSMA',
                -'HistoneH3',
                -'CD163',
                -'CD45')
str(data)
marker_names <- colnames(data[,5:28])

# Quantify how many observations have 0 nuclear marker signal

data %>% 
  filter('Hoechst' == 0) %>% 
  dplyr::count()

# Filter those cells out

data <- data %>%
  filter('Hoechst' > 0)

#Visualize distribution of Hoechst signal
data %>% 
  select(Staining_condition, `Hoechst`) %>% 
  ggplot(aes(x = as.factor(Staining_condition), y = `Hoechst`, fill = as.factor(Staining_condition))) + 
  geom_boxplot(outlier.alpha = 0.3) + 
  theme_bw() + 
  theme(legend.position = 'none') + 
  labs(x = 'Staining Condition')

#Filter out cells with nuclear signal below Q1 - 1.5IQR, where Q1 is the first quantile and IQR is the interquartile range.

data %>%
  group_by(Staining_condition) %>%
  mutate(min = quantile(Hoechst, 0.25) - 1.5 * IQR(Hoechst)) %>%
  filter(Hoechst < min) %>%
  ungroup() %>%
  dplyr::count()

data <- data %>% 
  group_by(Staining_condition) %>% 
  mutate(min = quantile(Hoechst, 0.25) - 1.5 * IQR(Hoechst)) %>% 
  filter(Hoechst >= min) %>% 
  ungroup()

#Mean Nuclear Marker Signal Normalization
mean_marker_data <- data %>% 
  dplyr::select(Staining_condition, all_of(marker_names)) %>% 
  pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value') %>% 
  group_by(Staining_condition, marker) %>% 
  summarise(marker_mean = mean(value)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Staining_condition, names_from = 'marker', values_from = 'marker_mean')

str(mean_marker_data)

df_norm <- data
for (i in 1:length(unique(data$Staining_condition))){
  row_index <- which(data$Staining_condition == unique(data$Staining_condition)[i])
  df_norm[row_index, 5:28] <- data[row_index, 5:28]/median(data$Hoechst[row_index])
}

str(df_norm)

#Arcsinh (Inverse hyperbolic Sine) Transformation

df_arcsinh <- df_norm %>% 
  mutate(across(all_of(marker_names), ~asinh(.x/0.001)))

#Universal percentile normalization 
df_trans <- df_arcsinh

rng <- colQuantiles(as.matrix(df_arcsinh[,5:28]), probs = c(0.001, 0.999))

expr <- t((t(as.matrix(df_arcsinh[,5:28]))-rng[,1]) / (rng[,2]-rng[,1]))

expr[expr < 0] <- 0

expr[expr > 1] <- 1

df_trans[,5:28] <- expr


# Plot density plots for each staining condition and marker
p <- df_trans %>%
    select(Staining_condition, all_of(marker_names)[1:24]) %>%
    pivot_longer(
        cols = -c(Staining_condition), names_to = "marker", values_to = "value"
    ) %>%
    mutate(
        marker = factor(
            marker,
            levels = c(
                "CD11b", "CD11c", "CD138", "CD15", "CD20",
                "CD206", "CD21", "CD3", "CD31", "CD4",
                "CD45RA", "CD68", "CD8", "Cytokeratin",
                "DCSIGN", "Foxp3", "GranzymeB", "H3K27ac", "H3K27me3",
                "HLADRA", "Ki67", "NaKATPase", "PD1", "DAPI"
            ),
            labels = c(
                "CD11b", "CD11c", "CD138", "CD15", "CD20",
                "CD206", "CD21", "CD3", "CD31", "CD4",
                "CD45RA", "CD68", "CD8", "Cytokeratin",
                "DC-SIGN", "Foxp3", "Granzyme B", "H3K27ac", "H3K27me3",
                "HLA-DRA", "Ki67", "Na/K-ATPase", "PD-1", "Hoechst"
            )
        )
    ) %>%
    # filter(marker == "CD3", Staining_condition == "20 min HIER, 1h 37C stain") %>%
    ggplot(aes(x = value, fill = Staining_condition)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~marker, ncol = 3, scales = "free") +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal"
    ) +
    labs(x = "Signal", y = "Density", fill = "Staining Condition") +
    guides(
        fill = guide_legend(
            nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5
        )
    )
    
ggsave(p, filename = 
paste0("Arcsinh transformed Hoechst normalised density plots (", data_type,").svg"), width = 8, height = 10)

# Perform statistical tests
df_long <- df_trans %>%
  select(Staining_condition, all_of(marker_names)[1:24]) %>%
  pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value')

# Perform Kruskal-Wallis test for each marker
kruskal_results <- df_long %>%
  group_by(marker) %>%
  kruskal_test(value ~ Staining_condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Save Kruskal-Wallis p-values to CSV
write_csv(kruskal_results, "kruskal_pvals.csv")

# Filter markers with significant p-values after adjustment
significant_markers <- kruskal_results %>%
  filter(p.adj < 0.05) %>%
  pull(marker)

# Initialize list to store Dunn's test results
dunn_results <- list()

# Perform Dunn's test for significant markers
for (marker in significant_markers) {
  dunn_test <- df_long %>%
    filter(marker == marker) %>%
    dunn_test(value ~ Staining_condition, p.adjust.method = "BH")
  dunn_test$marker <- marker
  dunn_results[[marker]] <- dunn_test
}

# Combine all Dunn's test results into a single dataframe
dunn_results_combined <- bind_rows(dunn_results)

# Save Dunn's test p-values to CSV
write_csv(dunn_results_combined, "dunn_pvals.csv")

# Calculate mean of each marker for each staining condition
mean_df <- df_trans %>%
  group_by(Staining_condition) %>%
  summarise(across(.cols = Hoechst:Cytokeratin, .fns = \(x) mean(x, na.rm = TRUE)))

# Calculate sd of each marker for each staining condition
sd_df <- df_trans %>%
  group_by(Staining_condition) %>%
  summarise(across(.cols = Hoechst:Cytokeratin, .fns = \(x) sd(x, na.rm = TRUE)))

# Calculate the number of cells for each staining condition
n_df <- df_trans %>%
  group_by(Staining_condition) %>%
  summarise(n = n()) %>%
  ungroup()

# Combine mean, sd, and sample size data frames
combined_df <- mean_df %>%
  rename_with(~ paste0(., "_mean"), -Staining_condition) %>%
  left_join(sd_df %>% rename_with(~ paste0(., "_sd"), -Staining_condition), by = "Staining_condition") %>%
  left_join(n_df, by = "Staining_condition")

# Define Cohen's d calculation function
calculate_cohens_d <- function(mean1, mean2, sd1, sd2, n1, n2) {
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  cohens_d <- (mean1 - mean2) / pooled_sd
  return(cohens_d)
}

# Define function to categorize effect size
categorize_effect_size <- function(d) {
  if (d < 0) {
    return("")
  } else if (d >= 0 & d <= 0.2) {
    return("small")
  } else if (d > 0.2 & d < 0.8) {
    return("medium")
  } else if (d >= 0.8) {
    return("large")
  }
}

# Define pairs of staining conditions
pairs <- list(
  c("10 min HIER, ON 4C stain", "20 min HIER, 1h 37C stain"),
  c("20 min HIER, 1h 37C stain", "10 min HIER, ON 4C stain"),
  c("10 min HIER, ON 4C stain", "20 min HIER, ON 4C stain"),
  c("20 min HIER, ON 4C stain", "10 min HIER, ON 4C stain"),
  c("10 min HIER, ON 4C stain", "40 min HIER, ON 4C stain"),
  c("40 min HIER, ON 4C stain", "10 min HIER, ON 4C stain"),
  c("20 min HIER, 1h 37C stain", "20 min HIER, ON 4C stain"),
  c("20 min HIER, ON 4C stain", "20 min HIER, 1h 37C stain"),
  c("20 min HIER, 1h 37C stain", "40 min HIER, ON 4C stain"),
  c("40 min HIER, ON 4C stain", "20 min HIER, 1h 37C stain"),
  c("20 min HIER, ON 4C stain", "40 min HIER, ON 4C stain"),
  c("40 min HIER, ON 4C stain", "20 min HIER, ON 4C stain")
)


# Define function to calculate the standard error of Cohen's d
calculate_se_cohens_d <- function(d, n1, n2) {
  se_d <- sqrt(((n1 + n2) / (n1 * n2)) + (d^2 / (2 * (n1 + n2))))
  return(se_d)
}

# Calculate Cohen's d for each marker and pair of staining conditions, including 95% CI
cohens_d_results <- data.frame()

for (pair in pairs) {
  condition1 <- pair[1]
  condition2 <- pair[2]
  
  data1 <- combined_df %>% filter(Staining_condition == condition1)
  data2 <- combined_df %>% filter(Staining_condition == condition2)
  
  for (marker in colnames(mean_df)[-1]) {
    mean1 <- as.numeric(data1[[paste0(marker, "_mean")]])
    mean2 <- as.numeric(data2[[paste0(marker, "_mean")]])
    sd1 <- as.numeric(data1[[paste0(marker, "_sd")]])
    sd2 <- as.numeric(data2[[paste0(marker, "_sd")]])
    n1 <- as.numeric(data1$n)
    n2 <- as.numeric(data2$n)
    
    cohens_d <- calculate_cohens_d(mean1, mean2, sd1, sd2, n1, n2)
    
    size_of_effect <- categorize_effect_size(cohens_d)
    
    se_d <- calculate_se_cohens_d(cohens_d, n1, n2)
    lower_ci <- cohens_d - 1.96 * se_d
    upper_ci <- cohens_d + 1.96 * se_d
    
    cohens_d_results <- rbind(cohens_d_results,
                              data.frame(Marker = marker,
                                         Condition1 = condition1,
                                         Condition2 = condition2,
                                         Cohens_d = cohens_d,
                                         Size_of_effect = size_of_effect,
                                         Lower_CI = lower_ci,
                                         Upper_CI = upper_ci)) %>%
      arrange(Marker)
  }
}

# View the Cohen's d results with 95% confidence intervals
print(cohens_d_results)

# Export the results to a CSV file
write.csv(cohens_d_results, "cohens_d_results.csv", row.names = FALSE)

# Filter rows where the Size_of_effect is "large"
large_effect_results <- cohens_d_results %>%
  filter(Size_of_effect == "large") %>%
  arrange(Marker)

# View the filtered results
print(large_effect_results)

# Write the filtered results to a new CSV file
write.csv(large_effect_results, "large_effect_results.csv", row.names = FALSE)

######Calculate Z scores and plot heatmap

# Calculate mean of each marker for each staining condition
mean_df_2 <- df_trans %>%
  group_by(Staining_condition) %>%
  summarise(across(.cols = CD3:Cytokeratin, .fns = \(x) mean(x, na.rm = TRUE)))

excluded_values <- list(
  '20 min HIER, 1h 37C stain' = c("DCSIGN"),  # Exclude Marker1 values for Condition1
  '10 min HIER, ON 4C stain' = c("DCSIGN", "Cytokeratin")   # Exclude Marker2 values for Condition3
)

#Replace values to exclude with NA
for (condition in names(excluded_values)) {
  for (marker in excluded_values[[condition]]) {
    mean_df_2[mean_df$Staining_condition == condition, marker] <- NA
  }
}
# Step 3: Calculate Z-scores, excluding NA values
means <- mean_df_2 %>%
  summarise(across(-Staining_condition, mean, na.rm = TRUE))

sds <- mean_df_2 %>%
  summarise(across(-Staining_condition, sd, na.rm = TRUE))

z_scores <- mean_df_2 %>%
  mutate(across(-Staining_condition, ~ (.-means[[cur_column()]]) / sds[[cur_column()]]))


# Step 5: Plot the heatmap
z_scores_long <- z_scores %>%
  pivot_longer(cols = -Staining_condition, names_to = "Marker", values_to = "Z_score")


heatmap <- ggplot(z_scores_long, aes(x = Staining_condition, y = Marker, fill = Z_score)) +
  geom_tile(color = "black", lwd = 0.4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Staining Condition", y = "Marker", fill = "Z-Score")

heatmap


#Save the plot as svg
ggsave(filename = "Heatmap of arcsinh transform normalised data Z scores_with excluded values.svg", width = 10, height = 8)