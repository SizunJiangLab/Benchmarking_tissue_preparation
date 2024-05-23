library(dplyr)
library(tidyverse)
library(matrixStats)
library(ggcorrplot)
library(ggpubr)
library(tidyr)
install.packages("rstatix")
library(rstatix)
library(readr)


data_slide1_FOV1 <- read.csv("dataScaleSize_slide1_FOV1.csv", header = TRUE)
data_slide1_FOV2 <- read.csv("dataScaleSize_slide1_FOV2.csv", header = TRUE)
names(data_slide1_FOV1)
names(data_slide1_FOV2)
col_names <- c("a.SMA", "H3K27me3", "Ki67") 
names(data_slide1_FOV2)[c(7, 15, 16)] <- col_names
names(data_slide1_FOV2)
data_slide1 <- rbind(data_slide1_FOV1, data_slide1_FOV2) 
data_slide1 <- mutate(data_slide1, Staining_condition = "20 min HIER, 1h 37C stain") 
head(data_slide1)

data_slide2_FOV1 <- read.csv("dataScaleSize_slide2_FOV1.csv", header = TRUE)
data_slide2_FOV2 <- read.csv("dataScaleSize_slide2_FOV2.csv", header = TRUE)
names(data_slide2_FOV1)
names(data_slide2_FOV2)
col_names <- c("a.SMA", "H3K27me3", "Ki67") 
names(data_slide2_FOV2)[c(7, 15, 16)] <- col_names
names(data_slide2_FOV2)
data_slide2 <- rbind(data_slide2_FOV1, data_slide2_FOV2) 
data_slide2 <- mutate(data_slide2, Staining_condition = "20 min HIER, ON 4C stain") 
head(data_slide2)

data_slide3_FOV1 <- read.csv("dataScaleSize_slide3_FOV1.csv", header = TRUE)
data_slide3_FOV2 <- read.csv("dataScaleSize_slide3_FOV2.csv", header = TRUE)
names(data_slide3_FOV1)
names(data_slide3_FOV2)
col_names <- c("a.SMA", "H3K27me3", "Ki67") 
names(data_slide3_FOV2)[c(7, 15, 16)] <- col_names
names(data_slide3_FOV2)
data_slide3 <- rbind(data_slide3_FOV1, data_slide3_FOV2) 
data_slide3 <- mutate(data_slide3, Staining_condition = "10 min HIER, ON 4C stain") 
head(data_slide3)

data_slide4_FOV1 <- read.csv("dataScaleSize_slide4_FOV1.csv", header = TRUE)
data_slide4_FOV2 <- read.csv("dataScaleSize_slide4_FOV2.csv", header = TRUE)
names(data_slide4_FOV1)
names(data_slide4_FOV2)
col_names <- c("a.SMA", "H3K27me3", "Ki67") 
names(data_slide4_FOV2)[c(7, 15, 16)] <- col_names
names(data_slide4_FOV2)
data_slide4 <- rbind(data_slide4_FOV1, data_slide4_FOV2) 
data_slide4 <- mutate(data_slide4, Staining_condition = "40 min HIER, ON 4C stain") 
head(data_slide4)


data <- rbind(data_slide1, data_slide2, data_slide3, data_slide4)
str(data)

marker_names <- colnames(data[,6:33])

#reorder columns
data$Staining_condition <- factor(data$Staining_condition,levels = c("20 min HIER, 1h 37C stain", 
                                                                     "20 min HIER, ON 4C stain", 
                                                                     "10 min HIER, ON 4C stain", 
                                                                     "40 min HIER, ON 4C stain"))
head(data)

data %>% 
  select(Staining_condition,  all_of(marker_names)[1:27]) %>% 
  # transform the dataframe longer for plotting
  pivot_longer(cols = -c(Staining_condition), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, fill = as.factor(Staining_condition))) + 
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, ncol = 3, nrow = 9, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'top')

data <- data %>% 
  dplyr::select(-'a.SMA',
                -'Histone.H3',
                -'CD163',
                -'CD45')
str(data)
marker_names <- colnames(data[,5:28])



df_arcsinh <- data %>% 
  mutate(across(all_of(marker_names), ~asinh(.x/0.001)))

#Universal percentile normalization 
df_trans <- df_arcsinh

rng <- colQuantiles(as.matrix(df_arcsinh[,5:28]), probs = c(0.001, 0.999))

expr <- t((t(as.matrix(df_arcsinh[,5:28]))-rng[,1]) / (rng[,2]-rng[,1]))

expr[expr < 0] <- 0

expr[expr > 1] <- 1

df_trans[,5:28] <- expr

df_trans %>% 
  select(Staining_condition,  all_of(marker_names)[1:24]) %>% 
  # transform the dataframe longer for plotting
  pivot_longer(cols = -c(Staining_condition), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = Staining_condition, y = marker, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Staining Condition", y = "Marker", fill = "Signal")

#Save the plot as svg
ggsave(filename = "Heatmap of signals without DAPI normalisation.svg", width = 10, height = 8)


# Quantify how many observations have 0 nuclear marker singal

data %>% 
  filter('DAPI' == 0) %>% 
  dplyr::count()

# Filter those cells out

data <- data %>%
  filter('DAPI' > 0)

data %>%
  select(Staining_condition, DAPI) %>%
  ggplot(aes(x = DAPI, fill = as.factor(Staining_condition))) + 
  geom_histogram(bins = 100, alpha = 0.5) + 
  theme_classic() + 
  theme(legend.position = 'right')

data %>% 
  select(Staining_condition, DAPI) %>% 
  ggplot(aes(x = as.factor(Staining_condition), y = DAPI , fill = as.factor(Staining_condition))) + 
  geom_boxplot(outlier.alpha = NULL) + 
  theme(legend.position = 'none') + 
  theme_classic() + 
  labs(x = 'Staining Condition') +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill = FALSE)

data %>%
  group_by(Staining_condition) %>%
  mutate(min = quantile(DAPI, 0.25) - 1.5 * IQR(DAPI)) %>%
  filter(DAPI < min) %>%
  ungroup() %>%
  dplyr::count()

data <- data %>% 
  group_by(Staining_condition) %>% 
  mutate(min = quantile(DAPI, 0.25) - 1.5 * IQR(DAPI)) %>% 
  filter(DAPI >= min) %>% 
  ungroup()

#Median Nuclear Marker Signal Normalization
mean_marker_data <- data %>% 
  dplyr::select(Staining_condition, all_of(marker_names)) %>% 
  pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value') %>% 
  group_by(Staining_condition, marker) %>% 
  summarise(marker_mean = mean(value)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = Staining_condition, names_from = 'marker', values_from = 'marker_mean')

str(mean_marker_data)

ggcorrplot(as.matrix(cor(mean_marker_data[,2:25], method = 'spearman')['DAPI', ])) + 
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_blank())

df_norm <- data

for (i in 1:length(unique(data$Staining_condition))){
  row_index <- which(data$Staining_condition == unique(data$Staining_condition)[i])
  df_norm[row_index, 5:28] <- data[row_index, 5:28]/median(data$DAPI[row_index])
}

str(df_norm)

df_norm %>% 
  select(Staining_condition,  all_of(marker_names)[1:24]) %>% 
  # transform the dataframe longer for plotting
  pivot_longer(cols = -c(Staining_condition), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, fill = as.factor(Staining_condition))) + 
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, ncol = 3, nrow = 9, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'right')


#Arcsinh (Inverse Hyperbolic Sine) Transformation
#Visualizing transformation of some markers to determine the optimal factors. 

markers_to_trans <- c('CD3')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('CD8')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('CD11b')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('CD206')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('CD31')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = '')

markers_to_trans <- c('DC.SIGN')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'right')

markers_to_trans <- c('CD20')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('CD4')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

markers_to_trans <- c('Granzyme.B')

df_norm %>%
  dplyr::select(Staining_condition, all_of(markers_to_trans)) %>% 
  mutate(
    across(
      all_of(markers_to_trans),
      list(
        "1000" = ~ asinh(.x / 1000),
        "100" = ~ asinh(.x / 100),
        "10" = ~ asinh(.x / 10),
        "0.1" = ~ asinh(.x / 0.1),
        "0.01" = ~ asinh(.x / 0.01),
        "0.001" = ~ asinh(.x / 0.001)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>% 
  pivot_longer(cols = -c('Staining_condition'), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = value, color = as.factor(Staining_condition))) +
  geom_histogram(bins = 100, alpha = 0.3) + 
  facet_wrap(~marker, scales = 'free') + 
  theme_classic() + 
  theme(legend.position = 'none')

#Visualize and choose optimal cofactor 

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
df_trans %>%
  select(Staining_condition, all_of(marker_names)[1:24]) %>%
  pivot_longer(cols = -c(Staining_condition), names_to = 'marker', values_to = 'value') %>%
  ggplot(aes(x = value, fill = Staining_condition)) +
  geom_density(alpha = 0.5) +  
  facet_wrap(~ marker, ncol = 3, scales = 'free') +
  theme_classic() +
  labs(x = "Signal", y = "Density", fill = "Staining Condition")
ggsave(filename = "Arcsinh transformed DAPI normalised density plots.svg", width = 10, height = 8)

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

#Plot Violin plots with p-values
# Define the comparisons for pairwise tests
my_comparisons <- list(
  c("10 min HIER, ON 4C stain", "20 min HIER, 1h 37C stain"),
  c("10 min HIER, ON 4C stain", "20 min HIER, ON 4C stain"),
  c("10 min HIER, ON 4C stain", "40 min HIER, ON 4C stain"),
  c("20 min HIER, 1h 37C stain", "20 min HIER, ON 4C stain"),
  c("20 min HIER, 1h 37C stain", "40 min HIER, ON 4C stain"),
  c("20 min HIER, ON 4C stain", "40 min HIER, ON 4C stain")
)

# Transform data to long format and create the plot
df_trans %>%
  select(Staining_condition, all_of(marker_names)[1:24]) %>%
  pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value') %>%
  ggplot(aes(x = Staining_condition, y = value)) +
  geom_violin(aes(fill = Staining_condition), trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  scale_y_continuous(name = "signal", labels = scales::comma) + 
  scale_fill_manual(values = c("#fbb4ae", "#ccebc5", "#b3cde3", "#FFCCFF")) + # specify your color values
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~marker, ncol = 3, scales = "free_y") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(p.adjust.method = "BH"), label.y = seq(1, by = 0.2, length.out = length(my_comparisons)), size = 2) +
  stat_compare_means(label.y = 2, size = 2) 

#Heatmap
df_trans %>% 
  select(Staining_condition,  all_of(marker_names)[2:24]) %>% 
  # transform the dataframe longer for plotting
  pivot_longer(cols = -c(Staining_condition), names_to = 'marker', values_to = 'value') %>% 
  ggplot(aes(x = Staining_condition, y = marker, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Staining Condition", y = "Marker", fill = "Signal")


# Calculate median of each marker for each staining condition
mean_df <- df_trans %>%
  group_by(Staining_condition) %>%
  summarise(across(.cols = DAPI:Cytokeratin, .fns = \(x) mean(x, na.rm = TRUE)))

excluded_values <- list(
  '20 min HIER, 1h 37C stain' = c("DC.SIGN"),  # Exclude Marker1 values for Condition1
  '10 min HIER, ON 4C stain' = c("DC.SIGN", "Cytokeratin")   # Exclude Marker2 values for Condition3
)

#Replace values to exclude with NA
for (condition in names(excluded_values)) {
  for (marker in excluded_values[[condition]]) {
    mean_df[mean_df$Staining_condition == condition, marker] <- NA
  }
}
# Step 3: Calculate Z-scores, excluding NA values
means <- mean_df %>%
  summarise(across(-Staining_condition, mean, na.rm = TRUE))

sds <- mean_df %>%
  summarise(across(-Staining_condition, sd, na.rm = TRUE))

z_scores <- mean_df %>%
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

