# CellXpress Segmentation Data

This folder contains cell segmentation data processed using the CellXpress pipeline. The data includes single-cell measurements from multiplexed imaging experiments, organized by batch and source institution.

**Note**: The actual data files are not included in this GitHub repository. Data files will be made available on Zenodo, and download links will be provided in the main [README.md](../README.md) once available. Please download the data from Zenodo and place it in this directory to use the analysis workflows.

## Folder Structure

```
data_cellXpress/
├── README.md                          # This file
├── Slide_metadata.csv                 # Metadata mapping files to sources and conditions
├── Slide_compare_pairs.csv            # Comparison pairs for statistical analysis
├── Slide_remove_markers.csv           # Markers to exclude from analysis
├── Slide_exclude_markers.csv          # Marker/slide-specific exclusions
├── Registered_Report_marker_sequence.csv  # Marker ordering for visualization
│
├── B1_ASTAR/                          # Batch 1: ASTAR data
├── B1_BIDMC/                          # Batch 1: BIDMC data
├── B1_Roche/                          # Batch 1: Roche data
├── B1_UK/                             # Batch 1: University of Kentucky data
│
├── B2_Stanford_IMC_Tonsil/           # Batch 2: Stanford IMC tonsil
├── B2_Stanford_MIBI/                  # Batch 2: Stanford MIBI data
├── B2_Stanford_Rarecyte/              # Batch 2: Stanford RareCyte data
│
├── B3_Novartis_LuCa/                  # Batch 3: Novartis lung cancer
├── B3_Novartis_tonsil/                # Batch 3: Novartis tonsil
├── B3_Stanford_OSCC/                   # Batch 3: Stanford OSCC
│
├── BIDMC/                             # BIDMC main dataset
├── Roche/                              # Roche main dataset
├── Stanford/                          # Stanford main dataset
│
└── SNR_BIDMC/                         # Signal-to-noise ratio data for BIDMC
    ├── SNR_ratios_BIDMC.csv           # SNR ratio data
    ├── metadata_BIDMC.csv             # Metadata mapping
    └── cell_counts_BIDMC.csv          # Cell counts per FOV
```

## Data Format

Each source folder contains `.qs` (quick serialization) files with the following naming convention:
- `*_raw_data.qs` - Main segmentation data (cell-level measurements)
- `*_ROI_info.qs` - Region of Interest (ROI) metadata

### Data File Structure

The `.qs` files contain R data structures:
- **raw_data.qs**: Data frame with columns:
  - Cell identifiers and spatial coordinates
  - Marker expression values (normalized fluorescence intensity)
  - Staining condition information
- **ROI_info.qs**: Metadata about regions of interest

To load these files in R:
```r
library(qs)
data <- qread("path/to/file_raw_data.qs")
roi_info <- qread("path/to/file_ROI_info.qs")
```

## Batch Organization

Data is organized into batches (B1, B2, B3) representing different experimental runs or time periods:
- **Batch 1 (B1)**: Early experimental batches
- **Batch 2 (B2)**: Stanford platform-specific data (MIBI, IMC, RareCyte)
- **Batch 3 (B3)**: Later experimental batches

## SNR Data

The `SNR_BIDMC/` folder contains pre-calculated signal-to-noise ratio data for BIDMC samples:
- **SNR_ratios_BIDMC.csv**: Signal-to-noise ratios per marker and slide
- **metadata_BIDMC.csv**: Maps SNR slide names to slide identifiers
- **cell_counts_BIDMC.csv**: Cell counts per slide and FOV for weighted averaging

See `SNR_BIDMC/README.md` for more details on SNR analysis.

## Metadata Files

- **Slide_metadata.csv**: Maps filenames to sources, types, FOVs, and sample names
- **Slide_compare_pairs.csv**: Defines comparison pairs for statistical tests
- **Slide_remove_markers.csv**: Lists markers to completely exclude from analysis
- **Slide_exclude_markers.csv**: Lists marker/slide combinations to exclude (set to NA)
- **Registered_Report_marker_sequence.csv**: Defines marker order for heatmap visualization

## Usage

To use this data with the analysis workflows:

1. Ensure all data files are in their respective source folders
2. Run `cellXpress_dataSlide_workflow.R` for main analysis
3. Run `CellXpress_SNR_Analysis.R` for SNR analysis (BIDMC only)

The workflows will automatically:
- Load data based on `Slide_metadata.csv`
- Apply marker exclusions and removals
- Perform normalization and transformation
- Generate statistical comparisons and visualizations

## Data Sources

This dataset includes contributions from:
- **ASTAR**: ASTAR COMET platform data
- **BIDMC**: Beth Israel Deaconess Medical Center
- **Novartis**: Novartis research data
- **Roche**: Roche research data
- **Stanford**: Stanford University data (MIBI, IMC, RareCyte platforms)
- **UKentucky**: University of Kentucky data

## Notes

- Data files use `.qs` format for efficient serialization (faster than CSV for large datasets)
- Each slide typically has multiple regions/ROIs
- Marker names follow CellXpress conventions (see marker mapping in workflows)
- The `SNR_BIDMC` folder contains processed SNR data, not raw segmentation data

