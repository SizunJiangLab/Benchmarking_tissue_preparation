# Mesmer Segmentation Data

This folder contains single-cell marker signal intensity data following Mesmer cell segmentation. The data includes single-cell measurements from multiplexed imaging experiments across multiple institutions and tissue types.

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
├── Initial_Optimization/BIDMC/                   # BIDMC Initial Optimization data (24 slides)
├── Initial_Optimization/Roche/                   # Roche Initial Optimization data (6 slides)
├── Initial_Optimization/Stanford/                # Stanford Initial Optimization data (12 slides)
├── LyophilizationTest_FigS2/                     # Comparison of Fresh vs Lyophilized manual vs Automated 
├── Reimagedslide_FigS5/                          # BIDMC Tonsil Run 1 vs Run 2 data
├── StorageConditionsExpt/                        # Storage conditions experiment
├── Validation/ASTAR_COMET_CRC/                   # ASTAR COMET CRC data (Validation)
├── Validation/ASTAR_COMET_Tonsil/                # ASTAR COMET Tonsil data
├── Validation/BIDMC_DLBCL/                       # BIDMC DLBCL data
├── Validation/BIDMC_Tonsil/                      # BIDMC Tonsil data
├── Validation/Novartis_Lung_Cancer/              # Novartis lung cancer data
├── Validation/Novartis_Tonsil/                   # Novartis tonsil data
├── Validation/Roche_Tonsil/                      # Roche tonsil data
├── Validation/Roche_Intestine/                   # Roche intestine data
├── Validation/Stanford_IMC_OSCC/                 # Stanford IMC OSCC data
├── Validation/Stanford_IMC_Tonsil/               # Stanford IMC tonsil data
├── Validation/Stanford_MIBI_Colon/               # Stanford MIBI colon data
├── Validation/Stanford_MIBI_Liver/               # Stanford MIBI liver data
├── Validation/Stanford_MIBI_LymphNode/           # Stanford MIBI lymph node data
├── Validation/Stanford_Orion_LN/                 # Stanford Orion lymph node data
├── Validation/Stanford_Orion_EndometrialCancer/  # Stanford Orion endometrial cancer data
├── Validation/UKentucky_SCC/                     # University of Kentucky SCC data
├── Validation/UKentucky_Tonsil/                  # University of Kentucky tonsil data
```

## Data Format

Each source folder contains CSV files with the following naming convention:

- `dataScaleSize_slide[#]_FOV[1|2].csv` - Main segmentation data
- `signal_ratios_slide[#]_FOV[1|2].csv` - Signal-to-noise ratio data (for Initial Optimization dataset)

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

- **ASTAR**: Agency for Science, Technology and Research
- **BIDMC**: Beth Israel Deaconess Medical Center
- **Novartis**: Novartis 
- **Roche**: Roche 
- **Stanford**: Stanford University (MIBI, IMC, Orion platforms)
- **UKentucky**: University of Kentucky

## Notes
- Each row in `dataScaleSize` is a cell and each column a marker. 
- Signal intensities have been normalized by cell size.
- Each slide typically has 2 FOVs (Fields of View)
- Marker names may vary slightly between sources (see marker mapping in workflows)
- Some sources include both `dataScaleSize` and `signal_ratios` data types
