# MESMER Data

This folder holds input CSVs and metadata for the MESMER workflow. Raw data are large and are NOT tracked in git. Only small metadata CSVs are versioned.

Data availability (raw): see Zenodo `https://doi.org/10.5281/zenodo.11391050`. Download and place files locally as described below.

## Folder Structure

Place source-specific CSVs under subfolders and keep metadata CSVs at the root of `data_mesmer/`.

```
data_mesmer/
  BIDMC/                             # per-source CSVs (e.g., dataScaleSize_*_FOV1.csv)
  Roche/
  Stanford/
  ...
  Slide_metadata.csv                 # maps filenames to Source/Type/FOV/Name
  Slide_compare_pairs.csv            # auto-generated (see below)
  Slide_exclude_markers.csv          # markers to gray out in plots
  Slide_remove_markers.csv           # markers to drop before analysis
  Registered_Report_marker_sequence.csv  # desired marker order (one column)
  Slide_exclude_cells.csv            # optional; per-slide cell exclusions
```

Notes:

- FOV files typically follow the pattern `dataScaleSize_<slide>_FOV1.csv` and `..._FOV2.csv`.
- Some datasets have tiles instead of FOVs (e.g., Stanford MIBI tiles). The workflow supports both.

## Metadata Files

- Slide_metadata.csv

  - Columns: `Filename, Source, Type, FOV, Name`
  - `Type` must be `dataScaleSize` for inclusion in MESMER analysis.

- Slide_compare_pairs.csv

  - Columns: `Source, Compare1, Compare2`
  - Auto-generated from unique `Name` per `Source` where `Type == dataScaleSize`.

- Slide_exclude_markers.csv and Slide_remove_markers.csv

  - Columns: `Source, Exclude_type, Exclude_value`

- Registered_Report_marker_sequence.csv

  - One column listing marker order to display in outputs.

- Slide_exclude_cells.csv (optional)
  - Columns: `Source, Slide, Marker`. Marker names are normalized (remove `.` and `-`).

## Setup Steps

1. Download raw data from Zenodo: `https://doi.org/10.5281/zenodo.11391050`.
2. Extract and place source CSVs into the corresponding subfolders under `data_mesmer/` (e.g., `data_mesmer/BIDMC/`).
3. Ensure the metadata CSVs listed above are present and up to date.
4. Generate comparison pairs:
   ```bash
   Rscript build_compare_pairs.R
   ```
   This writes `data_mesmer/Slide_compare_pairs.csv`.

## Example Slide_metadata.csv

```
Filename,Source,Type,FOV,Name
dataScaleSize_slide10_FOV1.csv,Roche,dataScaleSize,FOV1,Roche_10
dataScaleSize_slide10_FOV2.csv,Roche,dataScaleSize,FOV2,Roche_10
```

## Next Steps

- Choose `current_config_name` in `MESMER_dataSlide_workflow.R` and run:
  ```bash
  Rscript MESMER_dataSlide_workflow.R
  ```

## Git Tracking

The repository `.gitignore` excludes raw data in `data_mesmer/`. The following are tracked: `README.md`, `Slide_metadata.csv`, `Slide_exclude_markers.csv`, `Slide_remove_markers.csv`, `Registered_Report_marker_sequence.csv`. The generated `Slide_compare_pairs.csv` remains ignored by default.
