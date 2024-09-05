
# load mesmer data
load_mesmer_data <- function(data_colnames, data_folder) {
  data_slide <- function(slide_num, staining_condition) {
    data_slide_FOV1 <- read.csv(
      paste0(data_folder, "dataScaleSize_slide", slide_num, "_FOV1.csv"),
      col.names = data_colnames, header = TRUE
    )
    data_slide_FOV2 <- read.csv(
      paste0(data_folder, "dataScaleSize_slide", slide_num, "_FOV2.csv"),
      col.names = data_colnames, header = TRUE
    )
    data_slide <- rbind(data_slide_FOV1, data_slide_FOV2)
    data_slide <- mutate(data_slide, Staining_condition = staining_condition)
    return(data_slide)
  }

  data_slide1 <- data_slide(1, "20 min HIER, 1h 37C stain")
  data_slide2 <- data_slide(2, "20 min HIER, ON 4C stain")
  data_slide3 <- data_slide(3, "10 min HIER, ON 4C stain")
  data_slide4 <- data_slide(4, "40 min HIER, ON 4C stain")

  data <- bind_rows(data_slide1, data_slide2, data_slide3, data_slide4)
  return(data)
}

# load cX2 data
load_cX2_data <- function(data_colnames, data_folder) {
  
  load_data <- function(slide_num, staining_condition) {
    data_slide <- read_csv(
      file = paste0(data_folder, "cX2_Slide", slide_num, "_rawdata.csv"),
      col_sel = data_colnames,
      show_col_types = FALSE
      ) |>
      mutate(Staining_condition = staining_condition)
      
      return(data_slide)
    }

    data_slide1 <- load_data(1, "20 min HIER, 1h 37C stain")
    data_slide2 <- load_data(2, "20 min HIER, ON 4C stain")
    data_slide3 <- load_data(3, "10 min HIER, ON 4C stain")
    data_slide4 <- load_data(4, "40 min HIER, ON 4C stain")

    data <- bind_rows(data_slide1, data_slide2, data_slide3, data_slide4)
    return(data)
}

normalize_data <- function (data) {
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
   
  # Filter out the cells have 0 nuclear marker signal
  data <- data %>%
    filter('Hoechst' > 0)
  
  # Check if "Hoechst" column exists in the dataset
  if ("Hoechst" %in% colnames(data)) {
    ### Filter out cells with nuclear signal below Q1 - 1.5IQR, where Q1 is the first quantile and IQR is the interquartile range.
    data <- data %>% 
      group_by(Staining_condition) %>% 
      mutate(min = quantile(Hoechst, 0.25) - 1.5 * IQR(Hoechst)) %>% 
      filter(Hoechst >= min) %>% 
      ungroup()
    df_norm <- data
    for (i in 1:length(unique(data$Staining_condition))){
      row_index <- which(data$Staining_condition == unique(data$Staining_condition)[i])
      df_norm[row_index, 5:28] <- data[row_index, 5:28]/median(data$Hoechst[row_index])
    }
    str(df_norm)

    return(list(data=df_norm, marker_names=marker_names))
  } else {
    print("Column 'Hoechst' not found in the dataset.")
    return(list(data=data, marker_names=marker_names))
  }  
}

# Function to plot density plots for each staining condition and marker
plot_density_plots <- function(data, marker_names) {
  p <- data %>%
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
  
  return(p)
}

perform_statistical_and_kruskal_wallis_tests <- function(data, marker_names) {
  ### Perform statistical tests
  df_long <- data %>%
    select(Staining_condition, all_of(marker_names)[1:24]) %>%
    pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value')

  # Perform Kruskal-Wallis test for each marker
  kruskal_results <- df_long %>%
    group_by(marker) %>%
    kruskal_test(value ~ Staining_condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()

  return(kruskal_results)
}

plot_heatmap <- function(data, out_folder) {
  data <- data %>% mutate(id = str_c("c", 1:nrow(.)))
  counts <- data %>%
    column_to_rownames("id") %>%
    select(Hoechst:Cytokeratin) %>%
    t()
  group.by <- data$Staining_condition
  comb_2 <- combn(unique(data$Staining_condition), 2)
  res <- vector("list", ncol(comb_2))
  for (i in seq_along(res)) {
    res[[i]] <- wilcoxauc(counts, group.by, comb_2[, 1]) %>%
      mutate(ident.1 = comb_2[1, i], ident.2 = comb_2[2, i], )
  }
  res %>% bind_rows() -> wilcox_results
  
  #Prepare data for Cohen's d effect size test
  # Calculate mean of each marker for each staining condition
  mean_df <- data %>%
    group_by(Staining_condition) %>%
    summarise(across(.cols = Hoechst:Cytokeratin, .fns = \(x) mean(x, na.rm = TRUE)))
  
  # Calculate sd of each marker for each staining condition
  sd_df <- data %>%
    group_by(Staining_condition) %>%
    summarise(across(.cols = Hoechst:Cytokeratin, .fns = \(x) sd(x, na.rm = TRUE)))
  
  # Calculate the number of cells for each staining condition
  n_df <- data %>%
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
  write.csv(cohens_d_results, paste0(out_folder, "cohens_d_results.csv"), row.names = FALSE)
  
  # Filter rows where the Size_of_effect is "large"
  large_effect_results <- cohens_d_results %>%
    filter(Size_of_effect == "large") %>%
    arrange(Marker)
  
  # View the filtered results
  print(large_effect_results)
  
  # Write the filtered results to a new CSV file
  write.csv(large_effect_results, paste0(out_folder, "large_effect_results.csv"), row.names = FALSE)
  
  ####Calculate Coefficient of Variation (CV) and plot heatmap
  
  # Calculate mean of each marker for each staining condition
  cv_df_1 <- df_trans %>%
    group_by(Staining_condition) %>%
    summarise(across(.cols = CD3:Cytokeratin, .fns = function(x) {
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    }))
  excluded_values <- list(
    '20 min HIER, 1h 37C stain' = c("DCSIGN"),  # Exclude Marker1 values for Condition1
    '10 min HIER, ON 4C stain' = c("DCSIGN", "Cytokeratin")   # Exclude Marker2 values for Condition3
  )

  #Replace values to exclude with NA
  for (condition in names(excluded_values)) {
    for (marker in excluded_values[[condition]]) {
      cv_df_1[mean_df$Staining_condition == condition, marker] <- NA
    }
  }

  #Calculate Z-scores, excluding NA values
  means <- cv_df_1 %>%
    summarise(across(-Staining_condition, mean, na.rm = TRUE))
  
  sds <- cv_df_1 %>%
    summarise(across(-Staining_condition, sd, na.rm = TRUE))
  
  z_scores <- cv_df_1 %>%
    mutate(across(-Staining_condition, ~ (.-means[[cur_column()]]) / sds[[cur_column()]]))
  
  #Plot the heatmap
  cv_z_scores_long <- z_scores %>%
    pivot_longer(cols = -Staining_condition, names_to = "Marker", values_to = "Z_score")

  heatmap <- ggplot(cv_z_scores_long, aes(x = Staining_condition, y = Marker, fill = Z_score)) +
    geom_tile(color = "black", lwd = 0.4) +
    scale_fill_distiller(palette = 'PiYG', direction = -1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Staining Condition", y = "Marker", fill = "Z-Score")
  
  return(heatmap)
}

