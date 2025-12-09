# CellXpress2 Segmentation Data

This folder contains single-cell marker signal intensity data following cellXpress2 cell segmentation. The data includes single-cell measurements from multiplexed imaging experiments across multiple institutions and tissue types.

## Dataset Organization

The data is organized into two main categories:

### 1. Initial Optimization

Primary dataset for benchmarking tissue preparation protocols.

| Source   | Slides | Tissue     | Description                                      |
| -------- | ------ | ---------- | ------------------------------------------------ |
| BIDMC    | 24     | Tonsil     | Beth Israel - 24 HIER and Ab staining conditions |
| Roche    | 6      | Tonsil     | Roche - Selected conditions                      |
| Stanford | 12     | Lymph Node | Stanford - Selected conditions                   |

### 2. Cross-site Validation of Top Three and Bottom One Conditions

Independent validation across different institutions, platforms, and tissue types.

| Source                         | Institution | Platform | Tissue Type                 |
| ------------------------------ | ----------- | -------- | --------------------------- |
| ASTAR_COMET_Tonsil             | ASTAR       | COMET    | Tonsil                      |
| ASTAR_COMET_CRC                | ASTAR       | COMET    | Colorectal Cancer           |
| BIDMC_Tonsil                   | BIDMC       | Fusion   | Tonsil                      |
| BIDMC_DLBCL                    | BIDMC       | Fusion   | DLBCL (Lymph Node)          |
| Roche_Tonsil                   | Roche       | Fusion   | Tonsil                      |
| Roche_Intestine                | Roche       | Fusion   | Intestine                   |
| UKentucky_Tonsil               | U. Kentucky | COMET    | Tonsil                      |
| UKentucky_Skin                 | U. Kentucky | COMET    | Skin (SCC)                  |
| Novartis_Tonsil                | Novartis    | Fusion   | Tonsil                      |
| Novartis_LungCancer            | Novartis    | Fusion   | Lung Cancer                 |
| Stanford_IMC_Tonsil            | Stanford    | IMC      | Tonsil                      |
| Stanford_IMC_OSCC              | Stanford    | IMC      | Oral SCC                    |
| Stanford_MIBI_LymphNode_pooled | Stanford    | MIBI     | Lymph Node (4 tiles pooled) |
| Stanford_MIBI_Colon            | Stanford    | MIBI     | Colon                       |
| Stanford_MIBI_Liver            | Stanford    | MIBI     | Liver                       |
| Stanford_Orion_Lymph_node      | Stanford    | Orion    | Lymph Node                  |
| Stanford_Orion_Endometrium     | Stanford    | Orion    | Endometrial Cancer          |

## Folder Structure

```
data_cellXpress/
├── README.md
├── Slide_metadata.csv
├── Slide_compare_pairs.csv
├── Slide_remove_markers.csv
├── Slide_exclude_markers.csv
├── Registered_Report_marker_sequence.csv
│
├── Initial_Optimization/
│   ├── BIDMC/                    # 24 slides
│   ├── Roche/                    # 6 slides
│   └── Stanford/                 # 12 slides
│
└── Validation/
    ├── ASTAR/                    # ASTAR COMET (Tonsil + CRC)
    ├── BIDMC/                    # BIDMC Fusion (Tonsil + DLBCL)
    ├── Roche/                    # Roche Fusion (Tonsil + Intestine)
    ├── UK/                       # U. Kentucky COMET (Tonsil + Skin)
    ├── Stanford_MIBI/            # Stanford MIBI (LymphNode + Colon + Liver)
    ├── Stanford_Orion/           # Stanford Orion (LymphNode + Endometrium)
    ├── Stanford_IMC_OSCC/        # Stanford IMC OSCC
    ├── Stanford_IMC_Tonsil/      # Stanford IMC Tonsil
    ├── Novartis_LungCancer/      # Novartis Lung Cancer
    └── Novartis_Tonsil/          # Novartis Tonsil
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

The workflows will automatically:

- Load data based on `Slide_metadata.csv`
- Apply marker exclusions and removals
- Perform normalization and transformation
- Generate statistical comparisons and visualizations

## Data Sources

This dataset includes contributions from:

- **ASTAR**: Agency for Science, Technology and Research (Singapore)
- **BIDMC**: Beth Israel Deaconess Medical Center
- **Novartis**: Novartis Institutes for BioMedical Research
- **Roche**: Roche Tissue Diagnostics
- **Stanford**: Stanford University (MIBI, IMC, Orion platforms)
- **UKentucky**: University of Kentucky

## Notes

- Data files use `.qs` format for efficient serialization (faster than CSV for large datasets)
- Each slide typically has one or more ROIs
- Marker names follow CellXpress conventions (see marker mapping in workflows)
- Stanford_MIBI_LymphNode_pooled has 4 tiles per slide that are pooled together for analysis
