# cellXpress Data

This folder mirrors `data_mesmer/` but is intended for the cellXpress workflow. Raw input CSVs live in subfolders by Source and are NOT tracked in git. Only small metadata CSVs are versioned.

## Folder Structure

```
data_cellXpress/
  BIDMC/
  Roche/
  Stanford/
  ...
  Slide_metadata.csv
  Slide_exclude_markers.csv
  Slide_remove_markers.csv
  Registered_Report_marker_sequence.csv
```

Notes:

- Filenames in `Slide_metadata.csv` should correspond to CSVs in the per-Source folders.
- The schema parallels MESMER; `Type` can still be `dataScaleSize` for consistency.

## Metadata Files

- Slide_metadata.csv

  - Columns: `Filename, Source, Type, FOV, Name`

- Slide_exclude_markers.csv and Slide_remove_markers.csv

  - Columns: `Source, Exclude_type, Exclude_value`

- Registered_Report_marker_sequence.csv
  - One column listing desired marker order for plotting.

## Usage

1. Populate per-Source folders with your cellXpress CSVs.
2. Fill in `Slide_metadata.csv` with entries for your files.
3. Configure and run the workflow:
   ```bash
   Rscript cellXpress_dataSlide_workflow.R
   ```

## Git Tracking

The repository `.gitignore` excludes raw data in `data_cellXpress/`. The following are tracked: `README.md`, `Slide_metadata.csv`, `Slide_exclude_markers.csv`, `Slide_remove_markers.csv`, `Registered_Report_marker_sequence.csv`.
