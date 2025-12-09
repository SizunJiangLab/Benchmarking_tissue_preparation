# Balagan Spatial Analysis Workflow

This directory contains R scripts for comprehensive spatial analysis of tissue samples using the [Balagan](https://github.com/PierreBSC/Balagan) framework. The workflow runs 100 independent Balagan analyses to generate stable, reproducible metrics for comparing tissue preparation conditions.

> **Note**: This workflow analyzes **Mesmer segmentation data** (from `./data_mesmer/`). It does NOT use cellXpress data.

## Overview

The Balagan framework assesses spatial heterogeneity and subsampling efficiency in multiplexed imaging data through:

- **Tau (τ)**: Sampling efficiency metric - indicates how quickly subsampling recovers all cell clusters
- **Alpha (α)**: Spatial heterogeneity metric - quantifies the relationship between field-of-view size and cluster discovery

Lower τ and higher α indicate more heterogeneous tissue organization.

## Directory Structure

### Input Directories

```
../data_mesmer/                      # Mesmer segmentation data (INPUT)
  ├── Slide_metadata.csv             # Maps files to sources and conditions
  ├── Slide_remove_markers.csv       # Markers to exclude from analysis
  └── {Source}/                      # Source-specific data folders (e.g., BIDMC/)
      └── *.csv                      # Single-cell feature data
```

### Output Directories

```
../out_balagan_analysis/             # Base output from Step 0
  └── out_balagan_analysis_BIDMC_run_{1-100}/  # Results per run
      └── BIDMC/
          ├── *_Complex_sampling.csv # Raw subsampling data
          └── *_Fitting_tau.csv      # Fitted tau values

./balagan_consistency_analysis_from_raw/  # Output from Step 1
  ├── AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv
  ├── AGGREGATE_stable_alpha_slopes.csv
  └── plot_*.svg

./stable_rank_analysis_plots/        # Output from Steps 2-5
  ├── TABLE_*.csv                    # Combined rank tables
  └── PLOT_*.svg                     # Rank comparison plots

./advanced_correlation_analysis/     # Output from Step 6
  ├── TABLE_combined_core_metrics.csv
  └── PLOT_*.svg
```

## Workflow

### Step 0: Run 100 Balagan Analyses

**Script**: `00_run_100_balagan_analyses.R`

- Runs Balagan analysis 100 times with different random seeds
- Each run performs complex sampling analysis across multiple FOV sizes (50-500 µm)
- Outputs tau values for each slide and FOV combination
- **Runtime**: ~20-30 hours for 24 slides × 100 runs

**Input**:

- `./data_mesmer/Slide_metadata.csv`
- `./data_mesmer/Slide_remove_markers.csv`
- `./data_mesmer/{Source}/*.csv` (single-cell data)

**Output**:

- `./out_balagan_analysis_BIDMC_run_{N}/BIDMC/`
  - `*_Complex_sampling.csv` - Raw subsampling data
  - `*_Fitting_tau.csv` - Fitted tau values per FOV
  - `AGGREGATE_all_slides_tau_per_fov.csv` - Combined results per run

**Configuration**:

```r
NUMBER_OF_RUNS <- 100
output_base_dir <- "./out_balagan_analysis_BIDMC_run_"
```

---

### Step 1: Check Consistency and Calculate Alpha

**Script**: `01_check_consistency_and_calc_alpha.R`

- Loads raw complex sampling files from all 100 runs
- Recalculates tau values using non-linear least squares fitting
- Calculates alpha slope (exponent of power-law relationship)
- Generates consistency plots to assess result stability

**Input**:

- `./out_balagan_analysis/out_balagan_analysis_BIDMC_run_{1-100}/`
- `./out_MESMER_BIDMC_all/condition_summary.csv` (optional, for CV ordering)

**Output**:

- `./balagan_consistency_analysis_from_raw/`
  - `AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv` - All tau values
  - `AGGREGATE_all_100_runs_alpha_slopes.csv` - Alpha slopes per slide
  - `AGGREGATE_stable_alpha_slopes.csv` - Median alpha across runs
  - `plot_1_success_rate_heatmap.svg` - Success rate per slide/FOV
  - `plot_2_alpha_slope_consistency_boxplot_CV_ORDERED.svg` - Alpha consistency
  - `plot_3_alpha_slope_barchart_CV_ORDERED.svg` - Stable alpha ranks

**Key Metrics**:

- **Alpha (positive)**: Spatial heterogeneity exponent (higher = more heterogeneous)
- **Success rate**: Percentage of runs that generated valid tau values
- **Stability (SD)**: Standard deviation of alpha across 100 runs

---

### Step 2: Plot Stable Rank Heatmaps

**Script**: `02_plot_stable_rank_heatmaps.R`

- Calculates stable mean and median tau ranks across all 100 runs
- Generates standalone and CV-combined heatmaps
- Orders slides by stable performance rank

**Input**:

- `./balagan_consistency_analysis_from_raw/AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv`
- `./data_mesmer/condition_summary.csv` (optional, for CV comparison)

**Output**:

- `./stable_rank_analysis_plots/`
  - `AGGREGATE_stable_MEAN_tau_rank_heatmap.svg` - Mean-based ranks
  - `AGGREGATE_stable_MEDIAN_tau_rank_heatmap.svg` - Median-based ranks
  - `AGGREGATE_stable_MEAN_cv_tau_rank_heatmap.svg` - Mean ranks with CV comparison
  - `AGGREGATE_stable_MEDIAN_cv_tau_rank_heatmap.svg` - Median ranks with CV comparison

---

### Step 3: Plot Subsampling Curves

**Script**: `03_plot_subsampling_curves.R`

- Visualizes average cluster discovery curves across all 100 runs
- Faceted by FOV size to compare subsampling efficiency
- Shows how many regions are needed to recover all clusters

**Input**:

- `./out_balagan_analysis/out_balagan_analysis_BIDMC_run_{1-100}/`

**Output**:

- `./stable_rank_analysis_plots/PLOT_average_subsampling_curves.svg`

**Interpretation**:

- Steeper curves = more efficient sampling (clusters discovered quickly)
- Plateaus = all clusters recovered
- Slide-to-slide differences indicate heterogeneity variation

---

### Step 4: Plot Rank Correlations

**Script**: `04_plot_rank_correlations.R`

- Tests correlations between different ranking systems:
  - CV Rank vs. Alpha Rank
  - CV Rank vs. Tau Rank
  - Alpha Rank vs. Tau Rank

**Input**:

- `./balagan_consistency_analysis_from_raw/AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv`
- `./balagan_consistency_analysis_from_raw/AGGREGATE_stable_alpha_slopes.csv`
- `./data_mesmer/condition_summary.csv`

**Output**:

- `./stable_rank_analysis_plots/`
  - `TABLE_all_combined_ranks.csv` - All ranking metrics
  - `PLOT_correlation_cv_vs_alpha.svg`
  - `PLOT_correlation_cv_vs_tau.svg`
  - `PLOT_correlation_alpha_vs_tau.svg`

---

### Step 5: Calculate Tau Rank Stability

**Script**: `05_calculate_tau_rank_stability.R`

- Quantifies rank stability across FOV sizes for each slide
- Calculates standard deviation (SD) and interquartile range (IQR) of ranks
- Identifies slides with consistent vs. variable performance

**Input**:

- `./balagan_consistency_analysis_from_raw/AGGREGATE_all_100_runs_RECALCULATED_tau_data.csv`

**Output**:

- `./stable_rank_analysis_plots/`
  - `TABLE_stable_MEAN_rank_stability.csv`
  - `TABLE_stable_MEDIAN_rank_stability.csv`

---

### Step 6: Comprehensive Correlation Analysis

**Script**: `06_correlation_analysis.R`

- Performs multiple advanced correlation analyses
- Tests relationships between technical quality, spatial metrics, and biological features

**Analyses**:

1. Mean Intensity vs. CV Rank/Score
2. CV Score vs. Alpha/Tau Rank
3. Unknown/Mixed Cell Fraction vs. Alpha (optional)
4. All Markers vs. 10 Markers (optional)
5. Quadrant Classification (Efficiency vs. Heterogeneity)

**Input**:

- `./balagan_consistency_analysis_from_raw/`
- `./out_MESMER_BIDMC_all/` (Mesmer workflow output)
- `./stable_rank_analysis_plots/`
- `./data_mesmer/`

**Output**:

- `./advanced_correlation_analysis/`
  - `TABLE_combined_core_metrics.csv`
  - `PLOT_01a_intensity_vs_cv_rank.svg` through `PLOT_05_quadrant_*`

---

### Step 7: Regenerate Scatter/Heatmap Plots

**Script**: `07_regenerate_scatter_heatmap_plots.R`

- Utility script to regenerate plots from saved CSV data
- Allows tweaking plot parameters without re-running analysis

**Input**:

- `./out_balagan_analysis/{Source}/` (with `*_heatmap_data.csv` and `*_scatter_data.csv`)

**Output**:

- `*_scatter.svg` - Cell cluster spatial distribution
- `*_heatmap.svg` - Marker expression per cluster

---

### Standalone Analysis Script

**Script**: `balagan_analysis.R`

- Single-run Balagan analysis for testing/exploration
- Processes one file from each source at a time
- Uses Mesmer data from `./data_mesmer/`

---

## Prerequisites

### R Packages

```r
# Core packages
install.packages(c(
  "dplyr", "readr", "tidyr", "purrr",
  "ggplot2", "svglite", "ggpubr", "ggrepel",
  "pheatmap", "patchwork", "Polychrome"
))

# Balagan (for spatial analysis)
devtools::install_github("PierreBSC/Balagan")

# SingleCellExperiment (dependency)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

### Input Data Requirements

1. **Single-cell feature data** (from Mesmer segmentation): CSV files with columns:

   - `cellLabel`, `Y_cent`, `X_cent`, `cellSize`
   - Marker expression columns (e.g., `CD3`, `CD20`, etc.)

2. **Metadata files**:

   - `./data_mesmer/Slide_metadata.csv`
   - `./data_mesmer/Slide_remove_markers.csv`

3. **Optional files**:
   - `./data_mesmer/condition_summary.csv` (from Mesmer workflow, for CV comparison)
   - `./out_MESMER_{Source}_all/` (MESMER workflow outputs)

---

## Expected Runtime

- **Step 0**: ~20-30 hours (for 24 slides × 100 runs)
- **Step 1**: ~10-15 minutes (loading and recalculating 2400 files)
- **Steps 2-7**: < 5 minutes each

**Total**: ~24-30 hours for complete workflow

---

## Interpretation Guide

### Tau Rank (Sampling Efficiency)

- **Lower rank = Better efficiency** (fewer regions needed to find all clusters)
- Ranks are calculated per FOV size, then averaged across FOVs
- Stable across 100 runs indicates robust metric

### Alpha (Spatial Heterogeneity)

- **Positive alpha = Power-law scaling** (higher = more heterogeneous)
- **Alpha > 0.5**: Highly heterogeneous spatial organization
- **Alpha < 0.3**: Relatively homogeneous tissue
- Calculated as slope of log(tau) vs. log(FOV size)

### CV Score (Technical Quality)

- **Lower score = Better technical quality**
- Based on coefficient of variation across markers
- Independent metric from main Mesmer workflow

### Quadrant Classification

- **Top-left**: Inefficient & Heterogeneous (challenging samples)
- **Top-right**: Efficient & Heterogeneous (good heterogeneity, easy to sample)
- **Bottom-left**: Inefficient & Homogeneous (unexpected, investigate)
- **Bottom-right**: Efficient & Homogeneous (ideal for focused studies)
