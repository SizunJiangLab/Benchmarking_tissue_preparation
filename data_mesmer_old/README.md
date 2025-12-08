# MESMER Segmentation Data

This folder contains cell segmentation data processed using the MESMER (Mesmerizing Segmentation) pipeline. The data includes single-cell measurements from multiplexed imaging experiments across multiple institutions and tissue types.

**Note**: The actual data files are not included in this GitHub repository. Data files will be made available on Zenodo, and download links will be provided in the main [README.md](../README.md) once available. Please download the data from Zenodo and place it in this directory to use the analysis workflows.

## Folder Structure

```
data_mesmer/
├── README.md                          # This file
├── Slide_metadata.csv                 # Metadata mapping files to sources and conditions
├── Slide_compare_pairs.csv            # Comparison pairs for statistical analysis
├── Slide_remove_markers.csv           # Markers to exclude from analysis
├── Slide_exclude_markers.csv          # Marker/slide-specific exclusions
├── Registered_Report_marker_sequence.csv  # Marker ordering for visualization
├── condition_summary.csv              # Summary of staining conditions
│
├── ASTAR_COMET_CRC/                   # ASTAR COMET CRC data
├── ASTAR_COMET_Tonsil/                # ASTAR COMET Tonsil data
├── BIDMC/                             # BIDMC main dataset (24 slides)
├── BIDMC_DLBCL/                       # BIDMC DLBCL dataset
├── BIDMC_Tonsil/                      # BIDMC Tonsil dataset
├── BIDMC_Tonsil_compare/              # BIDMC Tonsil comparison dataset
├── Double/                            # Double staining condition data
├── Novartis_Lung_Cancer/              # Novartis lung cancer data
├── Novartis_Tonsil/                   # Novartis tonsil data
├── RareCyte_LN1_FOV1_FOV2_0.325/     # RareCyte lymph node data
├── RareCyte_TMA_FOV1_FOV2_0.325/     # RareCyte TMA data
├── Roche/                             # Roche main dataset
├── Roche_Tonsil/                      # Roche tonsil data
├── Roche_intestine/                   # Roche intestine data
├── Stanford/                          # Stanford main dataset
├── Stanford-scan1/                    # Stanford scan 1 data
├── Stanford_MIBI_Colon/               # Stanford MIBI colon data
├── Stanford_MIBI_Liver/               # Stanford MIBI liver data
├── Stanford_MIBI_LymphNode/           # Stanford MIBI lymph node (pooled)
├── Stanford_MIBI_LymphNode_Tile1-4/   # Stanford MIBI lymph node tiles
├── Stanford_OSCC/                     # Stanford OSCC data
├── Stanford_RareCyte_LN/              # Stanford RareCyte lymph node
├── Stanford_RareCyte_TMA/             # Stanford RareCyte TMA
├── Stanford_Tonsil/                   # Stanford tonsil data
├── StorageConditionsExpt/             # Storage conditions experiment
├── UKentucky_SCC/                     # University of Kentucky SCC data
└── UKentucky_Tonsil/                  # University of Kentucky tonsil data
```

## Data Format

Each source folder contains CSV files with the following naming convention:

- `dataScaleSize_slide[#]_FOV[1|2].csv` - Main segmentation data
- `signal_ratios_slide[#]_FOV[1|2].csv` - Signal-to-noise ratio data (for some sources)

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

## Usage

To use this data with the analysis workflows:

1. Ensure all data files are in their respective source folders
2. Run `MESMER_dataSlide_workflow.R` for main analysis
3. Run `MESMER_SignalNoise_workflow.R` for SNR analysis (requires signal_ratios files)

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
- **RareCyte**: RareCyte platform data

## Notes

- Data has been pre-processed and normalized
- Each slide typically has 2 FOVs (Fields of View)
- Marker names may vary slightly between sources (see marker mapping in workflows)
- Some sources include both `dataScaleSize` and `signal_ratios` data types
