load_mesmer_data <- function(data_folder, input_filenames, input_note) {
  # Extract slide identifiers from filenames (without FOV part)
  slide_ids <- gsub("_FOV\\d+\\.csv$", "", input_filenames)
  
  # Create a grouping dataframe
  file_groups <- data.frame(
    filename = input_filenames,
    staining_condition = input_note,
    slide_id = slide_ids
  )
  
  # Group by slide_id and staining_condition
  grouped_files <- file_groups %>%
    group_by(slide_id, staining_condition) %>%
    summarize(filenames = list(filename)) %>%
    ungroup()
  
  # Load and combine data for each group
  data_list <- vector("list", nrow(grouped_files))
  
  for (i in seq_along(data_list)) {
    group_filenames <- unlist(grouped_files$filenames[i])
    group_condition <- grouped_files$staining_condition[i]
    
    # Load and combine all FOVs for this slide
    fov_data_list <- vector("list", length(group_filenames))
    
    for (j in seq_along(group_filenames)) {
      file_path <- paste0(data_folder, group_filenames[j])
      fov_data <- read.csv(file_path, header = TRUE)
      fov_data_list[[j]] <- fov_data
    }
    
    # Combine all FOVs and add staining condition
    combined_data <- bind_rows(fov_data_list)
    combined_data <- mutate(combined_data, Staining_condition = group_condition)
    data_list[[i]] <- combined_data
  }
  
  # Combine all slides
  data <- bind_rows(data_list)
  return(data)
}

load_cellXpress_data <- function(data_folder, input_filenames, input_note) {
  # Create a grouping dataframe
  file_groups <- data.frame(
    slide_id = input_filenames,
    staining_condition = input_note,
    stringsAsFactors = FALSE
  )
  
  # Load and combine data for each group
  data_list <- vector("list", nrow(file_groups))
  
  for (i in seq_along(data_list)) {
    slide_id <- file_groups$slide_id[i]
    condition <- file_groups$staining_condition[i]
    
    # Find files matching the pattern for this slide
    pattern <- paste0(slide_id, ".*-raw_data\\.qs$")
    matching_files <- list.files(
      path = data_folder, 
      pattern = pattern, 
      full.names = TRUE
    )
    
    if (length(matching_files) == 0) {
      message("No files found matching pattern: ", pattern, " in folder: ", data_folder)
      next
    }
    
    message("Found ", length(matching_files), " file(s) for slide: ", slide_id)
    
    # Load and process each matching file
    file_data_list <- vector("list", length(matching_files))
    
    for (j in seq_along(matching_files)) {
      file_path <- matching_files[j]
      message("Loading file: ", basename(file_path))
      
      # Load the data using qs package
      slide_data <- qs::qread(file_path)
      
      # Process: keep first 4 columns and only whole-cell markers
      processed_data <- slide_data %>%
        # Step 1: Select only first 4 columns and those containing "whole-cell"
        select(c(1,3,2,4), contains("whole-cell")) %>%
        # Step 2: Rename columns to remove "(whole-cell)" suffix
        rename_with(~ str_replace(.x, " \\(whole-cell\\)", ""), 
                    contains("whole-cell")) %>%
        # Step 3: Add staining condition
        mutate(Staining_condition = condition)
      
      file_data_list[[j]] <- processed_data
    }
    
    # Combine all files for this slide
    combined_data <- bind_rows(file_data_list)
    data_list[[i]] <- combined_data
  }
  
  # Combine all slides
  data <- bind_rows(data_list)
  
  # Log the columns that were kept
  message("Final columns: ", paste(head(colnames(data), 10), collapse=", "), 
          "... (", length(colnames(data)), " total columns)")
  
  return(data)
}
normalize_data <- function (data) {
  marker_names <- data %>%
    select(5:ncol(data)) %>%
    select(-Staining_condition) %>%
    names()
  
  # Filter out the cells have 0 nuclear marker signal
  data <- data %>%
    filter('Hoechst' > 0)
  
  # Check if "Hoechst" column exists in the dataset
  if ("Hoechst" %in% colnames(data)) {
    print("Column 'Hoechst' found, applied extra normalization")
    ### Filter out cells with nuclear signal below Q1 - 1.5IQR, where Q1 is the first quantile and IQR is the interquartile range.
    data <- data %>% 
      group_by(Staining_condition) %>% 
      mutate(min = quantile(Hoechst, 0.25) - 1.5 * IQR(Hoechst)) %>% 
      filter(Hoechst >= min) %>% 
      ungroup()
    df_norm <- data
    for (i in 1:length(unique(data$Staining_condition))){
      row_index <- which(data$Staining_condition == unique(data$Staining_condition)[i])
      df_norm[row_index, marker_names] <- data[row_index, marker_names]/median(data$Hoechst[row_index])
    }
    str(df_norm)
    
    return(list(data=df_norm, marker_names=marker_names))
  } else {
    print("Column 'Hoechst' not found in the dataset.")
    return(list(data=data, marker_names=marker_names))
  } 
}

# Function to plot density plots for each staining condition and marker
plot_density_plots <- function(data,
                                marker_names,
                                original_marker_names,
                                legend.position = "bottom",
                                legend.direction = "horizontal",
                                legend.justification = "center",
                                legend.rows = 2) {
  p <- data %>%
    select(Staining_condition, all_of(marker_names)) %>%
    pivot_longer(
      cols = -c(Staining_condition), names_to = "marker", values_to = "value"
    ) %>%
    mutate(
      marker = factor(
        marker,
        levels = marker_names,
        labels = original_marker_names
      )
    ) %>%
    ggplot(aes(x = value, fill = Staining_condition)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~marker, ncol = 3, scales = "free") +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      legend.position = legend.position,
      legend.justification = legend.justification,
      legend.direction = legend.direction
    ) +
    labs(x = "Signal", y = "Density", fill = "Staining Condition") +
    guides(
      fill = guide_legend(
        nrow = legend.rows, byrow = TRUE, title.position = "top", title.hjust = 0.5
      )
    )
  
  return(p)
}

perform_statistical_and_kruskal_wallis_tests <- function(data, marker_names) {
  ### Perform statistical tests
  df_long <- data %>%
    select(Staining_condition, all_of(marker_names)) %>%
    pivot_longer(cols = -Staining_condition, names_to = 'marker', values_to = 'value')
  
  # Perform Kruskal-Wallis test for each marker
  kruskal_results <- df_long %>%
    group_by(marker) %>%
    kruskal_test(value ~ Staining_condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  return(kruskal_results)
}

calc_effect_size <- function(data, out_folder, pairs) {
  data <- data %>% mutate(id = str_c("c", 1:nrow(.)))
  counts <- data %>%
    column_to_rownames("id") %>%
    select(5:(ncol(data)-3)) %>%
    t()
  group.by <- data$Staining_condition
  comb_2 <- combn(unique(data$Staining_condition), 2)
  res <- vector("list", ncol(comb_2))
  for (i in seq_along(res)) {
    res[[i]] <- wilcoxauc(counts, group.by, comb_2[, 1]) %>%
      mutate(ident.1 = comb_2[1, i], ident.2 = comb_2[2, i], )
  }
  wilcox_results <- res %>% bind_rows()
  # Export the wilcox_results to a CSV file
  write.csv(wilcox_results, paste0(out_folder, "wilcox_results.csv"), row.names = FALSE)
  
  #Prepare data for Cohen's d effect size test
  # Calculate mean of each marker for each staining condition
  cols <- names(data)[5:(ncol(data) - 3)]
  
  mean_df <- data %>%
    group_by(Staining_condition) %>%
    summarise(across(.cols = all_of(cols), .fns = \(x) mean(x, na.rm = TRUE)))
  # Calculate sd of each marker for each staining condition
  sd_df <- data %>%
    group_by(Staining_condition) %>%
    summarise(across(.cols = all_of(cols), .fns = \(x) sd(x, na.rm = TRUE)))
  
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
  
  return(large_effect_results)
}


# Function to generate multiple heatmaps with different color schemes
plot_heatmaps <- function(df_trans, out_folder, excluded_values, remove_values, marker_sequence = NULL, original_marker_names = NULL) {
  # Select column names (markers)
  cols <- names(df_trans)[6:(ncol(df_trans) - 2)]
  
  # Remove specified markers if needed
  if (length(remove_values) > 0) {
    cols <- setdiff(cols, remove_values)
  }
  
  # Create mapping between processed names and original names
  name_mapping <- NULL
  if (!is.null(original_marker_names)) {
    # Create a mapping using exact matches only
    name_mapping <- cols
    names(name_mapping) <- cols
    
    # For each column, find its corresponding original name
    for (i in seq_along(cols)) {
      col_name <- cols[i]
      # Find exact match in original_marker_names
      if (col_name %in% original_marker_names) {
        name_mapping[col_name] <- col_name
      }
    }
  }
  
  # Order markers according to sequence file if provided
  if (!is.null(marker_sequence)) {
    # Find intersections between cols and marker_sequence
    common_markers <- intersect(cols, marker_sequence)
    
    # Order the common markers according to sequence
    ordered_common <- common_markers[order(match(common_markers, marker_sequence))]
    
    # Add any remaining markers
    remaining_cols <- setdiff(cols, ordered_common)
    cols <- c(ordered_common, remaining_cols)
  }
  
  
  # Calculate mean of each marker for each staining condition
  mean_df <- df_trans %>%
    group_by(Staining_condition) %>%
    summarise(across(
      .cols = all_of(cols),
      .fns = \(x) mean(x, na.rm = TRUE)
    ))
  
  # Get numeric ordering of staining conditions
  condition_order <- mean_df %>%
    mutate(
      condition_num = as.numeric(gsub("\\D", "", Staining_condition))
    ) %>%
    arrange(condition_num) %>%
    pull(Staining_condition)
  
  # Calculate CV of each marker for each staining condition
  cv_df <- df_trans %>%
    group_by(Staining_condition) %>%
    summarise(across(
      .cols = all_of(cols),
      .fns = function(x) {
        sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      }
    ))
  
  # Handle excluded markers
  for (condition in names(excluded_values)) {
    for (marker in excluded_values[[condition]]) {
      if (condition %in% mean_df$Staining_condition &&
          marker %in% names(mean_df)) {
        mean_df[mean_df$Staining_condition == condition, marker] <- NA
        cv_df[cv_df$Staining_condition == condition, marker] <- NA
      }
    }
  }
  
  # Calculate Z-scores for means
  mean_means <- mean_df %>%
    summarise(across(-Staining_condition, mean, na.rm = TRUE))
  
  mean_sds <- mean_df %>%
    summarise(across(-Staining_condition, sd, na.rm = TRUE))
  
  mean_z_scores <- mean_df %>%
    mutate(across(-Staining_condition, ~ (. - mean_means[[cur_column()]]) / mean_sds[[cur_column()]]))
  
  # Calculate Z-scores for CVs
  cv_means <- cv_df %>%
    summarise(across(-Staining_condition, mean, na.rm = TRUE))
  
  cv_sds <- cv_df %>%
    summarise(across(-Staining_condition, sd, na.rm = TRUE))
  
  cv_z_scores <- cv_df %>%
    mutate(across(-Staining_condition, ~ (. - cv_means[[cur_column()]]) / cv_sds[[cur_column()]]))
  
  # Convert to long format for plotting
  mean_long <- mean_df %>%
    pivot_longer(
      cols = -Staining_condition,
      names_to = "Marker",
      values_to = "Mean"
    ) %>%
    # Set factor levels to ensure correct ordering
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  mean_z_long <- mean_z_scores %>%
    pivot_longer(
      cols = -Staining_condition,
      names_to = "Marker",
      values_to = "Z_score"
    ) %>%
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  
  cv_long <- cv_df %>%
    pivot_longer(
      cols = -Staining_condition,
      names_to = "Marker",
      values_to = "CV"
    ) %>%
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  cv_z_long <- cv_z_scores %>%
    pivot_longer(
      cols = -Staining_condition,
      names_to = "Marker",
      values_to = "Z_score"
    ) %>%
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  # Join CV long with Z-score data to have both in one dataframe
  combined_cv_data <- cv_z_long %>%
    left_join(cv_long, by = c("Staining_condition", "Marker"))
  
  # 1. Heatmap of mean without Z-score transformation (white to red)
  mean_heatmap <- ggplot(mean_long, aes(x = Staining_condition, y = Marker, fill = Mean)) +
    geom_tile(color = "black", lwd = 0.4) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      na.value = "grey60"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Mean fluorescence intensity\n for each marker (without Z score)",
      x = "Condition",
      y = "Marker",
      fill = "Expression"
    )
  
  # Save mean heatmap
  ggsave(
    filename = paste0(out_folder, "Heatmap_mean_white_to_red.svg"),
    plot = mean_heatmap,
    width = 10,
    height = 8
  )
  
  # 2. Heatmap of mean with Z-score transformation (blue, white, red)
  mean_z_heatmap <- ggplot(mean_z_long,
                           aes(x = Staining_condition, y = Marker, fill = Z_score)) +
    geom_tile(color = "black", lwd = 0.4) +
    scale_fill_gradient2(
      low = "#5657A4",     # Custom blue
      mid = "white",       # White for 0 values
      high = "#EB2024",    # Custom red
      midpoint = 0,        # Ensure 0 is white
      na.value = "grey60"
    )  +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Z-scores of mean fluorescence\n intensity per marker",
      x = "Condition",
      y = "Marker",
      fill = "Z-Score"
    )
  
  # Save mean Z-score heatmap
  ggsave(
    filename = paste0(out_folder, "Heatmap_mean_zscore_blue_white_red.svg"),
    plot = mean_z_heatmap,
    width = 10,
    height = 8
  )
  
  # 3. Heatmap of CV with Z-score transformation (purple and green)
  cv_z_heatmap <- ggplot(combined_cv_data,
                         aes(x = Staining_condition, y = Marker, fill = Z_score)) +
    geom_tile(color = "black", lwd = 0.4) +
    scale_fill_gradientn(
      colors = c("#4D9221", "#A1D76A", "#E6F5D0", 
                 "#FFFFFF",  # Pure white at zero (instead of #F7F7F7)
                 "#FDE0EF", "#E9A3C9", "#C51B7D"),
      values = scales::rescale(c(min(combined_cv_data$Z_score, na.rm = TRUE), 
                                 -0.5, -0.35, 
                                 0,  # Zero point
                                 0.35, 0.5, 
                                 max(combined_cv_data$Z_score, na.rm = TRUE))),
      na.value = "grey60",
      guide = "colorbar"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Z-scores of coefficient of variation per marker",
      x = "Condition",
      y = "Marker",
      fill = "Z-Score of CV"
    )
  
  # Save CV Z-score heatmap
  ggsave(
    filename = paste0(out_folder, "Heatmap_CV_zscore_purple_green.svg"),
    plot = cv_z_heatmap,
    width = 12,
    height = 8
  )
  
  # 4. Heatmap of raw CV values (without Z-score transformation)
  cv_min <- min(cv_long$CV, na.rm = TRUE)
  cv_max <- max(cv_long$CV, na.rm = TRUE)
  cv_heatmap <- ggplot(cv_long, aes(x = Staining_condition, y = Marker, fill = CV)) +
    geom_tile(color = "black", lwd = 0.4) +
    geom_text(aes(label = sprintf("%.2f", CV)), size = 2.5) +
    scale_fill_gradientn(
      colors = c("white", "#FBEDF6", "#F5D0E6", "#EAB2D5", "#D378B4", "#C11C84"),
      values = c(0, 0.4, 0.6, 0.75, 0.85, 1),
      na.value = "grey60"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Coefficient of variation\nfor each marker (without Z score)",
      x = "Condition",
      y = "Marker",
      fill = "CV"
    )
  
  # Save raw CV heatmap
  ggsave(
    filename = paste0(out_folder, "Heatmap_CV_raw_values.svg"),
    plot = cv_heatmap,
    width = 12,
    height = 8
  )
  
  # Save the data frames for further analysis if needed
  write.csv(mean_df, paste0(out_folder, "mean_values.csv"), row.names = FALSE)
  write.csv(cv_df, paste0(out_folder, "cv_values.csv"), row.names = FALSE)
  write.csv(mean_z_scores,
            paste0(out_folder, "mean_z_scores.csv"),
            row.names = FALSE)
  write.csv(cv_z_scores,
            paste0(out_folder, "cv_z_scores.csv"),
            row.names = FALSE)
  
  # Return a list of the four heatmap objects
  return(
    list(
      mean_heatmap = mean_heatmap,
      mean_z_heatmap = mean_z_heatmap,
      cv_z_heatmap = cv_z_heatmap,
      cv_heatmap = cv_heatmap,
      cv_df = cv_df,
      cv_z_scores = cv_z_scores
    )
  )
}

# Function to load and process signal ratio data
load_signal_ratio_data <- function(data_folder, input_filenames, input_note) {
  data_list <- vector("list", length(input_filenames))
  
  for (i in seq_along(input_filenames)) {
    file_path <- paste0(data_folder, input_filenames[i])
    
    # Read the signal ratio CSV file
    data_file <- read_csv(file_path, show_col_types = FALSE)
    
    # Check if we have the expected columns
    if("Filename" %in% colnames(data_file) && 
       "Signal_withinvsoutsidecellmask" %in% colnames(data_file)) {
      
      # Clean up filenames to remove extensions
      data_file <- data_file %>%
        mutate(Filename = gsub("\\.tiff$|\\.tif$|\\.jpg$|\\.png$", "", Filename))
      
      # Pivot to wide format
      data_wide <- data_file %>%
        pivot_wider(
          names_from = Filename,
          values_from = Signal_withinvsoutsidecellmask
        )
      
      # Add staining condition
      data_wide$Staining_condition <- input_note[i]
      
      data_list[[i]] <- data_wide
    } else {
      warning(paste("Unexpected column format in file:", input_filenames[i]))
      # Create an empty dataframe with the staining condition
      data_list[[i]] <- tibble(Staining_condition = input_note[i])
    }
  }
  
  # Combine all dataframes
  if(length(data_list) > 0) {
    # Find common columns across all dataframes
    all_cols <- unique(unlist(lapply(data_list, names)))
    
    # Fill missing columns with NA and bind rows
    for(i in seq_along(data_list)) {
      missing_cols <- setdiff(all_cols, names(data_list[[i]]))
      for(col in missing_cols) {
        data_list[[i]][[col]] <- NA
      }
    }
    
    combined_data <- bind_rows(data_list)
    return(combined_data)
  } else {
    return(tibble())
  }
}

# Function to create multiple heatmap visualizations for signal-to-noise ratio data
create_snr_heatmaps <- function(snr_data, out_folder, remove_values = c()) {
  # First, identify marker columns (all except Staining_condition)
  marker_cols <- setdiff(names(snr_data), "Staining_condition")
  
  # Remove specified markers if requested
  if(length(remove_values) > 0) {
    marker_cols <- setdiff(marker_cols, remove_values)
  }
  
  # Clean marker names (remove special characters)
  snr_data <- snr_data %>%
    rename_with(~ gsub("[.-]", "", .), all_of(marker_cols))
  
  snr_data <- snr_data %>%
    select(-any_of(remove_values))
  
  # Update marker_cols with cleaned names
  marker_cols <- setdiff(names(snr_data), "Staining_condition")
  
  # Get original marker names for display
  old_marker_names <- marker_cols
  
  # Get numeric ordering of staining conditions
  condition_order <- snr_data %>%
    mutate(
      condition_num = as.numeric(gsub("\\D", "", Staining_condition))
    ) %>%
    arrange(condition_num) %>%
    pull(Staining_condition) %>%
    unique()
  
  # Convert data to long format for ggplot
  snr_long <- snr_data %>%
    pivot_longer(
      cols = all_of(marker_cols),
      names_to = "Marker",
      values_to = "SNR"
    ) %>%
    # Set factor levels to ensure correct ordering
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  # Check if nuclear marker (DAPI/Hoechst/Nucleus) is available for normalization
  has_nuclear_marker <- "DAPI" %in% marker_cols || "Hoechst" %in% marker_cols || "Nucleus" %in% marker_cols
  
  # Determine which nuclear marker column to use
  nuclear_marker_col <- if ("Hoechst" %in% marker_cols) {
    "Hoechst"
  } else if ("DAPI" %in% marker_cols) {
    "DAPI"
  } else if ("Nucleus" %in% marker_cols) {
    "Nucleus"
  } else {
    NULL
  }
  
  # Create normalized versions of all visualizations if DAPI/Hoechst is available
  if (!is.null(nuclear_marker_col)) {
    # Create a wide format dataset by marker for normalization
    snr_wide_by_marker <- snr_long %>%
      pivot_wider(
        names_from = Marker,
        values_from = SNR
      )
    
    # Normalize each marker by DAPI
    snr_dapi_normalized <- snr_wide_by_marker %>%
      mutate(across(all_of(setdiff(marker_cols, nuclear_marker_col)), 
                    ~ . / .data[[nuclear_marker_col]], 
                    .names = "{.col}_norm"))
    
    # Extract normalized columns
    norm_cols <- names(snr_dapi_normalized)[grepl("_norm$", names(snr_dapi_normalized))]
    
    # Convert back to long format for visualization
    snr_norm_long <- snr_dapi_normalized %>%
      select(Staining_condition, all_of(norm_cols)) %>%
      pivot_longer(
        cols = -Staining_condition,
        names_to = "Marker",
        values_to = "Normalized_SNR"
      ) %>%
      # Clean up marker names
      mutate(Marker = gsub("_norm$", "", Marker))
    
    # Calculate min/max for consistent scales
    min_norm_snr <- min(snr_norm_long$Normalized_SNR, na.rm = TRUE)
    max_norm_snr <- max(snr_norm_long$Normalized_SNR, na.rm = TRUE)
    
    # 1. Normalized heatmap with blue-white-red color scheme
    heatmap_dapi_norm <- ggplot(snr_norm_long, 
                                aes(x = Staining_condition, y = Marker, fill = Normalized_SNR)) +
      geom_tile(color = "black", linewidth = 0.4) +
      scale_fill_gradient2(
        low = "blue", 
        mid = "white", 
        high = "red", 
        midpoint = median(snr_norm_long$Normalized_SNR, na.rm = TRUE),
        na.value = "grey80",
        limits = c(min_norm_snr, max_norm_snr)
      ) +
      geom_text(aes(label = sprintf("%.2f", Normalized_SNR)), size = 3) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        title = paste0("Signal-to-Noise Ratio Normalized to ", nuclear_marker_col),
        x = "Staining Condition",
        y = "Marker",
        fill = "Normalized SNR"
      )
    
    # 2. Normalized threshold-based heatmap
    # Define thresholds for normalized values
    norm_threshold_colors <- c(
      "<0.5" = "#FF9999",   # Light red
      "0.5-0.8" = "#FFCC99", # Light orange
      "0.8-1.0" = "#FFFF99", # Light yellow
      "1.0-1.5" = "#CCFF99", # Light green
      "1.5-2.0" = "#99FF99", # Medium green
      ">2.0" = "#009900"     # Dark green
    )
    
    # Add normalized threshold category
    snr_norm_long <- snr_norm_long %>%
      mutate(SNR_category = case_when(
        Normalized_SNR < 0.5 ~ "<0.5",
        Normalized_SNR >= 0.5 & Normalized_SNR < 0.8 ~ "0.5-0.8",
        Normalized_SNR >= 0.8 & Normalized_SNR < 1.0 ~ "0.8-1.0",
        Normalized_SNR >= 1.0 & Normalized_SNR < 1.5 ~ "1.0-1.5",
        Normalized_SNR >= 1.5 & Normalized_SNR < 2.0 ~ "1.5-2.0",
        Normalized_SNR >= 2.0 ~ ">2.0",
        TRUE ~ NA_character_
      ))
    
    # Convert to ordered factor for proper legend ordering
    snr_norm_long$SNR_category <- factor(
      snr_norm_long$SNR_category,
      levels = names(norm_threshold_colors)
    )
    
    heatmap_norm_threshold <- ggplot(snr_norm_long, 
                                     aes(x = Staining_condition, y = Marker, fill = SNR_category)) +
      geom_tile(color = "black", linewidth = 0.4) +
      scale_fill_manual(
        values = norm_threshold_colors,
        na.value = "grey80",
        drop = FALSE
      ) +
      geom_text(aes(label = sprintf("%.2f", Normalized_SNR)), size = 3) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        title = paste0(nuclear_marker_col, "-Normalized SNR Heatmap (Threshold-Based)"),
        subtitle = "Values > 1.0 indicate better SNR than nuclear marker",
        x = "Staining Condition",
        y = "Marker",
        fill = "Normalized SNR Range"
      )
    
    # 3. Normalized bar plot of mean SNR per marker
    norm_marker_means <- snr_norm_long %>%
      group_by(Marker) %>%
      summarize(mean_norm_SNR = mean(Normalized_SNR, na.rm = TRUE)) %>%
      arrange(desc(mean_norm_SNR))
    
    barplot_norm_means <- ggplot(norm_marker_means, 
                                 aes(x = reorder(Marker, mean_norm_SNR), y = mean_norm_SNR)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = sprintf("%.2f", mean_norm_SNR)), hjust = -0.1, size = 3) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste0("Mean Signal-to-Noise Ratio Normalized to ", nuclear_marker_col),
        subtitle = "Values > 1.0 indicate better SNR than nuclear marker",
        x = "Marker",
        y = paste0("Mean ", nuclear_marker_col, "-Normalized SNR"),
        caption = "Red line indicates equal performance to nuclear marker"
      )
    
    # Save the normalized visualizations
    ggsave(
      filename = paste0(out_folder, "SNR_heatmap_DAPI_normalized.svg"),
      plot = heatmap_dapi_norm,
      width = 10,
      height = 8
    )
    
    ggsave(
      filename = paste0(out_folder, "SNR_heatmap_DAPI_norm_threshold.svg"),
      plot = heatmap_norm_threshold,
      width = 10,
      height = 8
    )
    
    ggsave(
      filename = paste0(out_folder, "SNR_barplot_DAPI_normalized.svg"),
      plot = barplot_norm_means,
      width = 8,
      height = 8
    )
    
    # Save the normalized data
    write_csv(snr_dapi_normalized, paste0(out_folder, "snr_dapi_normalized.csv"))
    write_csv(norm_marker_means, paste0(out_folder, "marker_mean_norm_snr.csv"))
  }
  
  # Calculate min/max for consistent scales for original SNR
  min_snr <- min(snr_long$SNR, na.rm = TRUE)
  max_snr <- max(snr_long$SNR, na.rm = TRUE)
  
  # 1. Basic heatmap with Purple-Yellow color scheme
  heatmap_purpleyellow <- ggplot(snr_long, aes(x = Staining_condition, y = Marker, fill = SNR)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "purple", 
      mid = "white", 
      high = "yellow", 
      midpoint = median(snr_long$SNR, na.rm = TRUE),
      na.value = "grey80",
      limits = c(min_snr, max_snr)
    ) +
    geom_text(aes(label = sprintf("%.2f", SNR)), size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Signal-to-Noise Ratio Heatmap (Purple-Yellow)",
      x = "Staining Condition",
      y = "Marker",
      fill = "SNR"
    )
  
  # Save the first heatmap
  ggsave(
    filename = paste0(out_folder, "SNR_heatmap_purpleyellow.svg"),
    plot = heatmap_purpleyellow,
    width = 10,
    height = 8
  )
  
  # 2. Threshold-based heatmap (clear visual for SNR > 1)
  # Define thresholds for coloring
  threshold_colors <- c(
    "<0.8" = "#FF9999",  # Light red
    "0.8-1.0" = "#FFCC99", # Light orange
    "1.0-1.5" = "#FFFF99", # Light yellow
    "1.5-2.0" = "#CCFF99", # Light green
    "2.0-3.0" = "#99FF99", # Medium green
    ">3.0" = "#009900"     # Dark green
  )
  
  # Add threshold category
  snr_long <- snr_long %>%
    mutate(SNR_category = case_when(
      SNR < 0.8 ~ "<0.8",
      SNR >= 0.8 & SNR < 1.0 ~ "0.8-1.0",
      SNR >= 1.0 & SNR < 1.5 ~ "1.0-1.5",
      SNR >= 1.5 & SNR < 2.0 ~ "1.5-2.0",
      SNR >= 2.0 & SNR < 3.0 ~ "2.0-3.0",
      SNR >= 3.0 ~ ">3.0",
      TRUE ~ NA_character_
    ))
  
  # Convert to ordered factor for proper legend ordering
  snr_long$SNR_category <- factor(
    snr_long$SNR_category,
    levels = names(threshold_colors)
  )
  
  heatmap_threshold <- ggplot(snr_long, aes(x = Staining_condition, y = Marker, fill = SNR_category)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_manual(
      values = threshold_colors,
      na.value = "grey80",
      drop = FALSE
    ) +
    geom_text(aes(label = sprintf("%.2f", SNR)), size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Signal-to-Noise Ratio Heatmap (Threshold-Based)",
      x = "Staining Condition",
      y = "Marker",
      fill = "SNR Range"
    )
  
  # Save the heatmap
  ggsave(
    filename = paste0(out_folder, "SNR_heatmap_threshold.svg"),
    plot = heatmap_threshold,
    width = 10,
    height = 8
  )
  
  # 3. Bar plot of mean SNR per marker
  # Calculate mean SNR for each marker
  marker_means <- snr_long %>%
    group_by(Marker) %>%
    summarize(mean_SNR = mean(SNR, na.rm = TRUE)) %>%
    arrange(desc(mean_SNR))
  barplot_means <- ggplot(marker_means, aes(x = reorder(Marker, mean_SNR), y = mean_SNR)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.2f", mean_SNR)), hjust = -0.1, size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Mean Signal-to-Noise Ratio per Marker",
      x = "Marker",
      y = "Mean SNR",
      caption = "Red dashed line indicates SNR = 1.0"
    )
  
  ggsave(
    filename = paste0(out_folder, "SNR_barplot_means.svg"),
    plot = barplot_means,
    width = 8,
    height = 8
  )
  
  # Save the processed data for reference
  write_csv(snr_data, paste0(out_folder, "processed_snr_data.csv"))
  
  # Save marker statistics
  write_csv(marker_means, paste0(out_folder, "marker_mean_snr.csv"))
  
  # Return a list of plot objects for potential further use
  result_list <- list(
    heatmap_purpleyellow = heatmap_purpleyellow,
    heatmap_threshold = heatmap_threshold,
    barplot_means = barplot_means
  )
  
  # Add DAPI-normalized visualizations to results if created
  if (!is.null(nuclear_marker_col)) {
    result_list$heatmap_dapi_norm <- heatmap_dapi_norm
    result_list$heatmap_norm_threshold <- heatmap_norm_threshold
    result_list$barplot_norm_means <- barplot_norm_means
  }
  
  return(result_list)
}
# Function to implement the scoring system directly as described
implement_scoring_system <- function(df_trans, out_folder, excluded_values = list(), remove_values = c(), penalty_score = 10, marker_sequence = NULL) {
  # Select column names (markers)
  cols <- names(df_trans)[6:(ncol(df_trans) - 2)]
  
  # Remove specified markers if needed
  if (length(remove_values) > 0) {
    cols <- setdiff(cols, remove_values)
  }
  
  # Order markers according to sequence file if provided
  if (!is.null(marker_sequence)) {
    # Find intersections between cols and marker_sequence
    common_markers <- intersect(cols, marker_sequence)
    
    # Order the common markers according to sequence
    ordered_common <- common_markers[order(match(common_markers, marker_sequence))]
    
    # Add any remaining markers
    remaining_cols <- setdiff(cols, ordered_common)
    cols <- c(ordered_common, remaining_cols)
  }
  
  # Calculate CV of each marker for each staining condition
  cv_df <- df_trans %>%
    group_by(Staining_condition) %>%
    summarise(across(
      .cols = all_of(cols),
      .fns = function(x) {
        sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      }
    ))
  
  # Get numeric ordering of staining conditions
  condition_order <- cv_df %>%
    mutate(
      condition_num = as.numeric(gsub("\\D", "", Staining_condition))
    ) %>%
    arrange(condition_num) %>%
    pull(Staining_condition)
  
  # Create a tracking matrix for excluded markers
  excluded_matrix <- matrix(FALSE, 
                            nrow = nrow(cv_df), 
                            ncol = length(cols), 
                            dimnames = list(cv_df$Staining_condition, cols))
  
  # Mark excluded markers in the tracking matrix and set to NA in cv_df
  for (condition in names(excluded_values)) {
    for (marker in excluded_values[[condition]]) {
      if (condition %in% cv_df$Staining_condition && marker %in% names(cv_df)) {
        # Mark as excluded
        excluded_matrix[condition, marker] <- TRUE
        # Set to NA in the data
        cv_df[cv_df$Staining_condition == condition, marker] <- NA
      }
    }
  }
  
  # Convert to long format for ranking
  cv_long <- cv_df %>%
    pivot_longer(
      cols = -Staining_condition,
      names_to = "Marker",
      values_to = "CV"
    ) %>%
    # Set factor levels to ensure correct ordering
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  # Add excluded status
  cv_long <- cv_long %>%
    mutate(
      Is_Excluded = mapply(function(cond, mrk) excluded_matrix[cond, mrk], 
                           as.character(Staining_condition), Marker)
    )
  
  # Calculate Z-scores for CV values
  cv_z_long <- cv_long %>%
    group_by(Marker) %>%
    mutate(Z_score = (CV - mean(CV, na.rm = TRUE)) / sd(CV, na.rm = TRUE)) %>%
    ungroup()
  
  # Rank CVs for each marker (1 = lowest/best, 4 = highest/worst)
  cv_ranked <- cv_long %>%
    group_by(Marker) %>%
    mutate(
      # Simplified rank calculation - only rank non-NA values
      Rank = if (all(is.na(CV))) {
        rep(NA_real_, length(CV))
      } else {
        # Create a vector for ranking that excludes NAs
        cv_no_na <- CV
        ranks <- rep(NA_real_, length(CV))
        non_na_indices <- which(!is.na(CV))
        
        if (length(non_na_indices) > 0) {
          # Only calculate ranks for non-NA values
          ranks[non_na_indices] <- rank(CV[non_na_indices], ties.method = "min")
        }
        
        ranks
      },
      
      # Apply penalty for non-working markers (both NA and explicitly excluded)
      Score = case_when(
        Is_Excluded ~ penalty_score,  # Excluded markers get penalty
        is.na(CV) ~ penalty_score,    # NA values get penalty
        TRUE ~ Rank                   # Otherwise use rank
      )
    ) %>%
    ungroup()
  
  # Convert back to wide format for easier viewing
  cv_ranks_wide <- cv_ranked %>%
    select(Staining_condition, Marker, Rank) %>%
    pivot_wider(
      names_from = Marker,
      values_from = Rank
    )
  
  cv_scores_wide <- cv_ranked %>%
    select(Staining_condition, Marker, Score) %>%
    pivot_wider(
      names_from = Marker,
      values_from = Score
    )
  
  # Calculate summary statistics for each condition
  condition_summary <- cv_ranked %>%
    group_by(Staining_condition) %>%
    summarise(
      Total_Markers = n_distinct(Marker),
      Working_Markers = sum(!is.na(CV) & !Is_Excluded),
      Excluded_Markers = sum(Is_Excluded),
      NonWorking_Markers = sum(is.na(CV) & !Is_Excluded),
      Total_Score = sum(Score, na.rm = TRUE),
      Average_Score = Total_Score / Total_Markers,
      Average_Rank = mean(Rank, na.rm = TRUE),
      Penalty_Count = sum(Score == penalty_score)
    ) %>%
    arrange(Average_Score) # Sort from best to worst
  
  # Calculate summary statistics for each marker
  marker_summary <- cv_ranked %>%
    group_by(Marker) %>%
    summarise(
      Best_Condition = as.character(Staining_condition[which.min(Score)]),
      Worst_Condition = as.character(Staining_condition[which.max(Score)]),
      Min_CV = min(CV, na.rm = TRUE),
      Max_CV = max(CV, na.rm = TRUE),
      CV_Range = Max_CV - Min_CV,
      Failed_Conditions = sum(is.na(CV)),
      Excluded_Conditions = sum(Is_Excluded)
    ) %>%
    arrange(desc(CV_Range))
  
  # Create total Z-score ranks
  cv_z_ranks <- cv_z_long %>%
    group_by(Marker) %>%
    mutate(
      # Robust Z-rank calculation
      Z_Rank = if (all(is.na(Z_score))) {
        rep(NA_real_, length(Z_score))
      } else {
        # Create a vector for ranking that excludes NAs
        ranks <- rep(NA_real_, length(Z_score))
        non_na_indices <- which(!is.na(Z_score))
        
        if (length(non_na_indices) > 0) {
          # Only calculate ranks for non-NA values
          ranks[non_na_indices] <- rank(Z_score[non_na_indices], ties.method = "min")
        }
        
        ranks
      }
    ) %>%
    ungroup() %>%
    select(Staining_condition, Marker, Z_Rank) %>%
    pivot_wider(
      names_from = Marker,
      values_from = Z_Rank
    )
  # Calculate total rank for each condition
  total_ranks <- cv_z_ranks %>%
    rowwise() %>%
    mutate(Total_Z_Rank = sum(c_across(-Staining_condition), na.rm = TRUE)) %>%
    arrange(Total_Z_Rank)
  
  # Save all results to CSV files
  write.csv(cv_df, paste0(out_folder, "cv_values.csv"), row.names = FALSE)
  write.csv(cv_long, paste0(out_folder, "cv_values_long.csv"), row.names = FALSE)
  write.csv(cv_z_long, paste0(out_folder, "cv_z_scores_long.csv"), row.names = FALSE)
  write.csv(cv_ranked, paste0(out_folder, "cv_ranks_and_scores.csv"), row.names = FALSE)
  write.csv(cv_ranks_wide, paste0(out_folder, "cv_ranks_wide.csv"), row.names = FALSE)
  write.csv(cv_scores_wide, paste0(out_folder, "cv_scores_wide.csv"), row.names = FALSE)
  write.csv(condition_summary, paste0(out_folder, "condition_summary.csv"), row.names = FALSE)
  write.csv(marker_summary, paste0(out_folder, "marker_summary.csv"), row.names = FALSE)
  write.csv(total_ranks, paste0(out_folder, "total_z_ranks.csv"), row.names = FALSE)
  
  # Create visualizations
  
  # 1. Score heatmap
  score_heatmap <- ggplot(cv_ranked, aes(x = Staining_condition, y = Marker, fill = Score)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "white", 
      high = "red",
      na.value = "grey60"
    ) +
    # Add markers for excluded cells
    geom_text(aes(label = ifelse(Is_Excluded, "X", round(Score, 1))), 
              size = 3, fontface = ifelse(cv_ranked$Is_Excluded, "bold", "plain")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Marker Score Heatmap (Lower is Better)",
      subtitle = paste("Non-working markers penalized with score of", penalty_score),
      x = "Staining Condition",
      y = "Marker",
      fill = "Score"
    )
  
  # Save the score heatmap
  ggsave(
    filename = paste0(out_folder, "marker_score_heatmap.svg"),
    plot = score_heatmap,
    width = 12,
    height = 8
  )
  
  # 2. Average score bar chart ordered by score value
  avg_score_plot <- ggplot(condition_summary,
                           aes(x = reorder(Staining_condition, Average_Score), y = Average_Score)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f", Average_Score)), vjust = -0.5, size = 3.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Average Score by Staining Condition",
      subtitle = "Lower scores indicate better performance (lower CV)",
      x = "Staining Condition",
      y = "Average Score"
    )
  
  # Save the average score plot
  ggsave(
    filename = paste0(out_folder, "average_score_by_condition.svg"),
    plot = avg_score_plot,
    width = 10,
    height = 6
  )
  
  # 3. CV values heatmap with rankings
  cv_heatmap <- ggplot(cv_ranked, aes(x = Staining_condition, y = Marker, fill = CV)) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = c("white", "#FBEDF6", "#F5D0E6", "#EAB2D5", "#D378B4", "#C11C84"),
      values = c(0, 0.4, 0.6, 0.75, 0.85, 1),
      na.value = "grey60"
    ) +
    # Modify text display for excluded cells
    geom_text(aes(label = ifelse(Is_Excluded, 
                                 "X", 
                                 sprintf("%.2f\n(R:%d)", CV, Rank))), 
              size = 2.5, 
              fontface = ifelse(cv_ranked$Is_Excluded, "bold", "plain")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "CV Values with Rankings",
      subtitle = "R:1 = best (lowest CV), R:4 = worst (highest CV)",
      x = "Staining Condition",
      y = "Marker",
      fill = "CV"
    )
  
  # Save the CV heatmap with rankings
  ggsave(
    filename = paste0(out_folder, "cv_heatmap_with_ranks.svg"),
    plot = cv_heatmap,
    width = 12,
    height = 8
  )
  
  #4. # Create the combined CV visualization
  combined_cv_vis <- create_combined_cv_plot(
    condition_summary = condition_summary,
    cv_z_long = cv_z_long,
    condition_order = condition_order
  )
  
  # Save the combined plot
  ggsave(
    filename = paste0(out_folder, "combined_cv_ranks_heatmap.svg"),
    plot = combined_cv_vis$combined_plot,
    width = 12,
    height = 10
  )
  
  # Return the results
  return(list(
    cv_df = cv_df,
    cv_long = cv_long,
    cv_ranked = cv_ranked,
    condition_summary = condition_summary,
    marker_summary = marker_summary,
    total_ranks = total_ranks,
    score_heatmap = score_heatmap,
    avg_score_plot = avg_score_plot,
    cv_heatmap = cv_heatmap
  ))
}

# Process excluded markers (markers to gray out in specific conditions)
process_excluded_markers <- function(source) {
  # Create a default empty result list
  result <- list()
  
  # Filter the exclusions for the given source
  exclusions <- excluded_markers %>%
    filter(Source == source)
  
  # Check if we have any exclusions
  if (nrow(exclusions) > 0) {
    # Process each exclusion
    for (i in 1:nrow(exclusions)) {
      condition <- exclusions$Name[i]  # Using Name instead of Condition
      marker <- exclusions$Marker[i]
      
      # Initialize the list element if it doesn't exist
      if (!(condition %in% names(result))) {
        result[[condition]] <- c()
      }
      
      # Add the marker to the exclusion list for this condition
      result[[condition]] <- c(result[[condition]], marker)
    }
  }
  
  return(result)
}
# Function to create a combined visualization with average scores above CV heatmap
create_combined_cv_plot <- function(condition_summary, cv_z_long, condition_order) {
  # Prepare the condition summary data for plotting
  rank_data <- condition_summary %>%
    mutate(Staining_condition = factor(Staining_condition, levels = condition_order))
  
  # Create the top barplot for average scores
  p_ranks <- ggplot(rank_data, aes(x = Staining_condition, y = Average_Score)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", Average_Score)), vjust = -0.3, size = 4.5) + # Increased font size, one decimal place
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(b = 0),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 12), # Larger font for y-axis values
      axis.title.y = element_text(size = 14)  # Larger font for y-axis title
    ) +
    labs(y = "Avg Score") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  # Create the bottom CV z-score heatmap
  p_heatmap <- ggplot(cv_z_long, aes(x = Staining_condition, y = Marker, fill = Z_score)) +
    geom_tile(color = "black", lwd = 0.4) +
    scale_fill_gradientn(
      colors = c("#4D9221", "#A1D76A", "#E6F5D0", 
                 "#FFFFFF",  # Pure white at zero
                 "#FDE0EF", "#E9A3C9", "#C51B7D"),
      values = scales::rescale(c(min(cv_z_long$Z_score, na.rm = TRUE), 
                                 -0.5, -0.35, 
                                 0,  # Zero point
                                 0.35, 0.5, 
                                 max(cv_z_long$Z_score, na.rm = TRUE))),
      na.value = "grey60",
      guide = "colorbar"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Larger condition names
      axis.text.y = element_text(size = 16), # Larger marker names
      axis.title = element_text(size = 14),  # Larger axis titles
      plot.title = element_text(hjust = 0.5, size = 16), # Larger title
      legend.title = element_text(size = 12), # Larger legend title
      legend.text = element_text(size = 10)   # Larger legend text
    ) +
    labs(
      x = "Condition",
      y = "Marker",
      fill = "Z-Score of CV"
    )
  
  # Combine the plots
  combined_plot <- cowplot::plot_grid(
    p_ranks, p_heatmap,
    ncol = 1,
    align = "v",
    rel_heights = c(0.2, 0.8),
    axis = "lr"
  )
  
  return(list(
    combined_plot = combined_plot,
    rank_plot = p_ranks,
    heatmap = p_heatmap
  ))
}


# Custom function to reorder staining conditions for BIDMC subset
reorder_bidmc_subset_conditions <- function(data) {
  # Define the desired order
  custom_order <- c("BIDMC_16", "BIDMC_13", "BIDMC_5", "BIDMC_21")
  
  # Check if the data contains Staining_condition column
  if ("Staining_condition" %in% colnames(data)) {
    # Convert Staining_condition to factor with custom order
    data <- data %>%
      mutate(Staining_condition = factor(Staining_condition, levels = custom_order))
    
    # Sort the dataframe by the new factor levels
    data <- data %>%
      arrange(Staining_condition)
  }
  
  return(data)
}