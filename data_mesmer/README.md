# Mesmer Segmentation Data

This folder contains single-cell marker signal intensity data following Mesmer cell segmentation. The data includes single-cell measurements from multiplexed imaging experiments across multiple institutions and tissue types.

## Dataset Organization

The data is organized into three main categories:

### 1. Initial Optimization

Primary dataset for benchmarking tissue preparation protocols with signal-to-noise ratio analysis.

| Source   | Slides | Tissue     | Description                           |
| -------- | ------ | ---------- | ------------------------------------- |
| BIDMC    | 24     | Tonsil     | Beth Israel - Full optimization panel |
| Roche    | 6      | Tonsil     | Roche - Selected conditions           |
| Stanford | 12     | Lymph Node | Stanford - Selected conditions        |

### 2. Cross-site Validation of Top Three and Bottom One Conditions

Independent validation across different institutions, platforms, and tissue types.

| Source                           | Institution | Platform | Tissue Type             |
| -------------------------------- | ----------- | -------- | ----------------------- |
| ASTAR_COMET_Tonsil               | ASTAR       | COMET    | Tonsil                  |
| ASTAR_COMET_CRC                  | ASTAR       | COMET    | Colorectal Cancer       |
| BIDMC_Tonsil                     | BIDMC       | Fusion   | Tonsil                  |
| BIDMC_DLBCL                      | BIDMC       | Fusion   | DLBCL (Lymph Node)      |
| Roche_Tonsil                     | Roche       | Fusion   | Tonsil                  |
| Roche_intestine                  | Roche       | Fusion   | Intestine               |
| UKentucky_Tonsil                 | U. Kentucky | COMET    | Tonsil                  |
| UKentucky_SCC                    | U. Kentucky | COMET    | Squamous Cell Carcinoma |
| Novartis_Tonsil                  | Novartis    | Fusion   | Tonsil                  |
| Novartis_Lung_Cancer             | Novartis    | Fusion   | Lung Cancer             |
| Stanford_IMC_Tonsil              | Stanford    | IMC      | Tonsil                  |
| Stanford_IMC_OSCC                | Stanford    | IMC      | Oral SCC                |
| Stanford_MIBI_LymphNode          | Stanford    | MIBI     | Lymph Node (4 tiles)    |
| Stanford_MIBI_Colon              | Stanford    | MIBI     | Colon                   |
| Stanford_MIBI_Liver              | Stanford    | MIBI     | Liver                   |
| Stanford_Orion_LN                | Stanford    | Orion    | Lymph Node              |
| Stanford_Orion_EndometrialCancer | Stanford    | Orion    | Endometrial Cancer      |

### 3. Supplementary Experiments

Experiments for supplementary figures.

| Source                       | Description                                          | Figure  |
| ---------------------------- | ---------------------------------------------------- | ------- |
| LyophilizationTest_FigS2     | Fresh vs Lyophilized manual vs Lyophilized automated | Fig S2  |
| Reimagedslide_FigS5          | BIDMC Tonsil Run 1 vs Run 2 comparison               | Fig S5  |
| StorageConditionsExpt_FigS16 | Comparing different slide storage conditions         | Fig S16 |

## Folder Structure

```
data_mesmer/
├── README.md
├── Slide_metadata.csv
├── Slide_compare_pairs.csv
├── Slide_remove_markers.csv
├── Slide_exclude_markers.csv
├── Registered_Report_marker_sequence.csv
│
├── Initial_Optimization/
│   ├── BIDMC/                    # 24 slides, 2 FOVs each
│   ├── Roche/                    # 6 slides, 2 FOVs each
│   ├── Stanford/                 # 12 slides, 2 FOVs each
│   └── cell_counts/              # Pre-calculated cell counts (from process_cell_counts.R)
│
├── Validation/
│   ├── ASTAR_COMET_CRC/
│   ├── ASTAR_COMET_Tonsil/
│   ├── BIDMC_DLBCL/
│   ├── BIDMC_Tonsil/
│   ├── Novartis_Lung_Cancer/
│   ├── Novartis_Tonsil/
│   ├── Roche_Tonsil/
│   ├── Roche_intestine/
│   ├── Stanford_IMC_OSCC/
│   ├── Stanford_IMC_Tonsil/
│   ├── Stanford_MIBI_Colon/
│   ├── Stanford_MIBI_Liver/
│   ├── Stanford_MIBI_LymphNode/
│   ├── Stanford_Orion_EndometrialCancer/
│   ├── Stanford_Orion_LN/
│   ├── UKentucky_SCC/
│   └── UKentucky_Tonsil/
│
├── LyophilizationTest_FigS2/
├── Reimagedslide_FigS5/
└── StorageConditionsExpt_FigS16/
```

## Data Format

Each source folder contains CSV files with the following naming convention:

- `dataScaleSize_slide[#]_FOV[1|2].csv` - Main segmentation data
- `signal_ratios_slide[#]_FOV[1|2].csv` - Signal-to-noise ratio data (Initial Optimization only)

The conditions corresponding to each slide number can be found in Supplementary Table 18.

### CSV File Structure

Each CSV file contains:

- **Columns 1-4**: Cell identifiers and spatial coordinates
  - `CellID`: Unique cell identifier
  - `X`, `Y`: Spatial coordinates
  - `Area`: Cell area
- **Column 5+**: Marker expression values
  - Marker names (e.g., `CD3`, `CD8`, `CD20`, `Cytokeratin`, etc.)
  - Values represent normalized fluorescence intensity
- **Last column**: `Staining_condition` - Experimental condition identifier

## Metadata Files

- **Slide_metadata.csv**: Maps filenames to sources, types, FOVs, and sample names
- **Slide_compare_pairs.csv**: Defines comparison pairs for statistical tests
- **Slide_remove_markers.csv**: Lists markers to completely exclude from analysis
- **Slide_exclude_markers.csv**: Lists marker/slide combinations to exclude (set to NA)
- **Registered_Report_marker_sequence.csv**: Defines marker order for heatmap visualization

## Cell Counts

The `Initial_Optimization/cell_counts/` folder contains pre-calculated cell counts per FOV, generated by `process_cell_counts.R`. This data is used by `Mesmer_SignalNoise_workflow.R` to calculate weighted average signal-to-noise ratios across FOVs.

Each source subfolder (BIDMC/, Roche/, Stanford/) contains:

- **fov_cell_counts.csv**: Cell counts per FOV with columns: `Filename`, `Staining_condition`, `FOV`, `Cell_count`, `Total_cells`, `Cell_proportion`
- **cell_count_summary.csv**: Summary statistics per condition including total cells, FOV count, min/max/mean/median cells

## Usage

To use this data with the analysis workflows:

1. Ensure all data files are in their respective source folders
2. Run `Mesmer_dataSlide_workflow.R` for main analysis
3. Run `Mesmer_SignalNoise_workflow.R` for SNR analysis (requires signal_ratios files)

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

- Each row in `dataScaleSize` is a cell and each column a marker.
- Signal intensities have been normalized by cell size.
- Each slide typically has 2 FOVs (Fields of View), except Stanford_MIBI_LymphNode which has 4 tiles.
- Marker names may vary slightly between sources (see marker mapping in workflows)
- Initial Optimization sources include both `dataScaleSize` and `signal_ratios` data types
