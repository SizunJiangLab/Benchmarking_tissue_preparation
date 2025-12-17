# Image Preprocessing (Stage 1)

This folder contains Python scripts for processing raw QPTIFF images: cropping FOV regions, running Mesmer cell segmentation, extracting single-cell features, and calculating signal intensity ratios.

## Overview

The preprocessing pipeline takes raw multiplexed imaging data and produces single-cell CSV files that can be analyzed by the R workflows (Stage 2).

```
Raw QPTIFF Images → Crop FOVs → Mesmer Segmentation → Single-cell CSVs
```

## Requirements

```bash
pip install -r ../requirements.txt
pip install deepcell  # For Mesmer segmentation
```

## Scripts

| Script                                            | Purpose                                                              |
| ------------------------------------------------- | -------------------------------------------------------------------- |
| `crop_mesmer_featureextraction_signaltonoise.py`  | Main preprocessing pipeline                                          |
| `pylibs/tissue_preparation.py`                    | Helper functions for image reading, segmentation, feature extraction |
| `segmentation_scFeature_extraction_4slides.ipynb` | Jupyter notebook example                                             |

## Inputs

- **Raw QPTIFF images**: Download from [BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/) (accession TBD)
- **FOV coordinates**: Defined in the script's `crop_coords_dict` or reference [`Master_metadata_Mesmer.csv`](../Master_metadata_Mesmer.csv)

## Configuration

Edit these variables in `crop_mesmer_featureextraction_signaltonoise.py`:

```python
# Path to folder containing .qptiff files
data_folder = "/path/to/qptiff/files"

# Output directory for results
output_folder = "/path/to/output"

# FOV crop coordinates: {slide_key: {"FOV1": (x_min, x_max, y_min, y_max), "FOV2": ...}}
crop_coords_dict = {
    "1": {
        "FOV1": (10000, 17000, 5500, 11000),
        "FOV2": (14000, 21000, 11000, 16500),
    },
    # ... more slides
}

# Marker names in channel order
markers = ['DAPI', 'CD3', 'CD45', ...]
```

## Usage

```bash
cd preprocessing
python crop_mesmer_featureextraction_signaltonoise.py
```

To process specific FOVs only, modify the script:

```python
# Process only FOV1 for all slides
process_and_segment(data_folder, output_folder, crop_coords_dict, markers, fovs=["FOV1"])
```

## Outputs

For each slide and FOV, the script generates:

| Output                                                                   | Description                            |
| ------------------------------------------------------------------------ | -------------------------------------- |
| `MESMER_outputs/{slide}_{FOV}/MESMER_mask.tiff`                          | Cell segmentation mask                 |
| `MESMER_outputs/{slide}_{FOV}/seg_overlay.tiff`                          | Segmentation overlay visualization     |
| `MESMER_outputs/{slide}_{FOV}/seg_outline.tiff`                          | Cell outlines visualization            |
| `data_slide{slide}_{FOV}.csv`                                            | Single-cell features (raw)             |
| `dataScaleSize_slide{slide}_{FOV}.csv`                                   | Single-cell features (size-normalized) |
| `Individualtiff_slide{slide}_{FOV}/{marker}.tiff`                        | Per-marker cropped images              |
| `Individualtiff_slide{slide}_{FOV}/signal_ratios_slide{slide}_{FOV}.csv` | Signal intensity ratios                |

## Signal Intensity Ratios

The `signal_ratios` CSV contains per-marker ratios comparing signal inside vs outside cell masks:

| Column                      | Description                                                                           |
| --------------------------- | ------------------------------------------------------------------------------------- |
| `Marker`                    | Marker name                                                                           |
| `Normalized_signal_invsout` | Ratio of (signal inside cells / DAPI inside) to (signal outside cells / DAPI outside) |

## Mesmer Parameters

Default segmentation parameters (can be adjusted in script):

```python
image_mpp = 0.50           # Microns per pixel
maxima_threshold = 0.075   # Cell detection sensitivity
interior_threshold = 0.2   # Cell boundary threshold
```

## Note

Pre-generated CSV outputs are available on [Zenodo](https://zenodo.org/) (TBD). Most users can skip this preprocessing step and download the CSVs directly for use with Stage 2 R workflows.
