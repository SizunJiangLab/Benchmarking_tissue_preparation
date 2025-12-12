import os
import re
import tifffile
import numpy as np
import pandas as pd
import skimage.io
from skimage import img_as_ubyte
from pathlib import Path
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay
from pylibs.tissue_preparation import (
    read_qtiff_image,
    generate_nuclear_and_membrane,
    stack_nuclear_and_membrane,
    run_segmentation,
    extract_sc_features,
)

# Mesmer PARAMETERS
image_mpp = 0.50
maxima_threshold = 0.075
interior_threshold = 0.2


def process_and_segment(data_folder, output_folder, crop_coords_dict, markers, fovs=None):
    """
    Process QPTIFF images: crop, segment with Mesmer, extract features, and calculate SNR.
    
    Args:
        data_folder: Path to folder containing .qptiff files
        output_folder: Path to save outputs
        crop_coords_dict: Nested dict {slide_key: {"FOV1": (x_min, x_max, y_min, y_max), "FOV2": ...}}
        markers: List of marker names in channel order
        fovs: List of FOVs to process, e.g. ["FOV1", "FOV2"]. If None, process all FOVs in dict.
    """
    os.makedirs(output_folder, exist_ok=True)

    data_folder = Path(data_folder)
    all_qptiffs = list(data_folder.glob("*.qptiff"))

    for slide_key, fov_coords in crop_coords_dict.items():
        # Find matching QPTIFF file
        qptiff_path = None
        for qptiff_f in all_qptiffs:
            match = re.search(rf"(^.*Registered_Report_Slide{slide_key}.qptiff)", str(qptiff_f))
            if match:
                qptiff_path = os.path.join(data_folder, f"{match.group(0)}")
                break
        
        if qptiff_path is None:
            print(f"⚠️ No matching qptiff file for slide {slide_key}, skipping.")
            continue

        if not os.path.exists(qptiff_path):
            print(f"⚠️ Missing qptiff file for slide {slide_key}, skipping.")
            continue

        print(f"\n{'='*60}")
        print(f"Processing slide {slide_key}: {qptiff_path}")
        print(f"{'='*60}")

        # Read qptiff image (once per slide)
        slide = read_qtiff_image(qptiff_path)
        print(f"  - Slide shape: {slide.shape}")

        # Generate nuclear and membrane images (once per slide)
        nuclear, membrane = generate_nuclear_and_membrane(slide)
        stack = stack_nuclear_and_membrane(nuclear=nuclear, membrane=membrane)
        print(f"  - Stacked shape: {stack.shape}")

        # Transpose for feature extraction (once per slide)
        array_list = np.transpose(slide, (1, 2, 0))

        # Determine which FOVs to process
        fovs_to_process = fovs if fovs else list(fov_coords.keys())

        for fov in fovs_to_process:
            if fov not in fov_coords:
                print(f"  ⚠️ {fov} not in coordinates for slide {slide_key}, skipping.")
                continue

            coords = fov_coords[fov]
            x_min, x_max, y_min, y_max = coords

            print(f"\n  Processing {fov} with coords: {coords}")

            # Validate cropping boundaries
            _, height, width, _ = stack.shape
            if x_min < 0 or y_min < 0 or x_max > width or y_max > height:
                print(f"  ⚠️ Invalid crop coordinates for slide {slide_key} {fov}: {coords}")
                continue

            cropped_stack = stack[:, y_min:y_max, x_min:x_max, :]

            if cropped_stack.size == 0:
                print(f"  ⚠️ Cropped stack for slide {slide_key} {fov} is empty! Skipping.")
                continue

            # Run Mesmer segmentation
            predictions = run_segmentation(cropped_stack)

            # Create RGB image and overlay
            rgb_image = create_rgb_image(cropped_stack, channel_colors=["green", "blue"])
            overlay = make_outline_overlay(rgb_data=rgb_image, predictions=predictions)

            # Save Mesmer outputs
            output_dir = os.path.join(
                output_folder, 
                f'MESMER_outputs/{maxima_threshold}maxima_{slide_key}_{interior_threshold}interior_{slide_key}_{fov}/'
            )
            os.makedirs(output_dir, exist_ok=True)

            skimage.io.imsave(os.path.join(output_dir, "seg_outline.tiff"), img_as_ubyte(overlay[0, ..., 0]), check_contrast=False)
            skimage.io.imsave(os.path.join(output_dir, "seg_overlay.tiff"), img_as_ubyte(overlay[0, ...]), check_contrast=False)
            skimage.io.imsave(os.path.join(output_dir, "MESMER_mask.tiff"), predictions[0, ..., 0], check_contrast=False)

            print(f"  ✅ Saved Mesmer outputs for slide {slide_key} {fov}")

            # Single-cell feature extraction
            mask = skimage.io.imread(os.path.join(output_dir, "MESMER_mask.tiff"))
            cropped_array_list = array_list[y_min:y_max, x_min:x_max, :]

            data, dataScaleSize, cell_props, cellSizes = extract_sc_features(cropped_array_list, mask=mask, markers=markers)

            data_df = pd.DataFrame(data)
            data_df.columns = markers
            data_full = pd.concat((
                pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
                pd.DataFrame(cellSizes, columns=["cellSize"]),
                data_df
            ), axis=1)

            dataScaleSize_df = pd.DataFrame(dataScaleSize)
            dataScaleSize_df.columns = markers
            dataScaleSize_full = pd.concat((
                pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
                pd.DataFrame(cellSizes, columns=["cellSize"]),
                dataScaleSize_df
            ), axis=1)

            # Save feature CSVs
            save_dir = os.path.join(output_folder, 'extracted_features')
            os.makedirs(save_dir, exist_ok=True)

            data_full.to_csv(os.path.join(save_dir, f'data_slide{slide_key}_{fov}.csv'), index=False)
            dataScaleSize_full.to_csv(os.path.join(save_dir, f'dataScaleSize_slide{slide_key}_{fov}.csv'), index=False)

            print(f"  ✅ Saved single-cell features for slide {slide_key} {fov}")

            # Export individual TIFF files and calculate SNR
            individual_tiff_dir = os.path.join(output_folder, f'Individualtiff_slide{slide_key}_{fov}')
            os.makedirs(individual_tiff_dir, exist_ok=True)

            raw_signals = {}
            area_in = np.sum(mask != 0)
            area_out = np.sum(mask == 0)

            for i in range(cropped_array_list.shape[2]):
                marker = markers[i]
                im_marker = cropped_array_list[:, :, i]

                # Save TIFF
                path_marker = os.path.join(individual_tiff_dir, f"{marker}.tiff")
                tifffile.imwrite(path_marker, im_marker)

                # Calculate signals
                signal_in = np.sum(im_marker[mask != 0])
                signal_out = np.sum(im_marker[mask == 0])

                raw_signals[marker] = {
                    "signal_in": signal_in,
                    "signal_out": signal_out,
                    "area_in": area_in,
                    "area_out": area_out
                }

            # Calculate SNR ratios
            results = []
            dapi_in = raw_signals["DAPI"]["signal_in"]
            dapi_out = raw_signals["DAPI"]["signal_out"]

            for marker, vals in raw_signals.items():
                norm_ratio = np.nan

                if dapi_in > 0 and dapi_out > 0 and vals["signal_in"] > 0 and vals["signal_out"] > 0:
                    # DAPI + Area normalized ratio
                    norm_ratio = (vals["signal_in"] / dapi_in / area_in) / (vals["signal_out"] / dapi_out / area_out)

                # Handle DAPI itself
                if marker == "DAPI":
                    norm_ratio = 1.0

                results.append({
                    "Marker": marker,
                    "Normalized_signal_invsout": norm_ratio
                })

            # Save SNR CSV
            df = pd.DataFrame(results)
            output_csv_path = os.path.join(individual_tiff_dir, f"signal_ratios_slide{slide_key}_{fov}.csv")
            df.to_csv(output_csv_path, index=False)

            print(f"  ✅ Saved SNR ratios for slide {slide_key} {fov}")

    print(f"\n{'='*60}")
    print("Processing complete!")
    print(f"{'='*60}")


# =============================================================================
# CONFIGURATION
# =============================================================================

data_folder = "/mnt/nfs/storage/Fusion_Registered_Report/Initial_Optimization_Stage2/BIDMC"
output_folder = "/mnt/nfs/home/jialelee/Registered_Report/Initial Optimization/BIDMC/CheckSNR_Normalizedbyarea/Output"

# Crop coordinates: {slide_key: {"FOV1": (x_min, x_max, y_min, y_max), "FOV2": (...)}}
crop_coords_dict = {
    "1": {
        "FOV1": (10000, 17000, 5500, 11000),
        "FOV2": (14000, 21000, 11000, 16500),
    },
    "2": {
        "FOV1": (12500, 19500, 6250, 11750),
        "FOV2": (15600, 22600, 12350, 17850),
    },
    "3": {
        "FOV1": (9000, 16000, 6000, 11500),
        "FOV2": (12300, 19300, 13000, 18500),
    },
    "4": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "5": {
        "FOV1": (9500, 16500, 5500, 11000),
        "FOV2": (12000, 19000, 12400, 17900),
    },
    "6": {
        "FOV1": (9000, 16000, 4500, 10000),
        "FOV2": (12500, 19500, 11000, 16500),
    },
    "7": {
        "FOV1": (9500, 16500, 6000, 11500),
        "FOV2": (12000, 19000, 13000, 18500),
    },
    "8": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "9": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "10": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "11": {
        "FOV1": (9000, 16000, 4000, 9500),
        "FOV2": (13000, 20700, 10200, 15700),
    },
    "12": {
        "FOV1": (12000, 19000, 6250, 11750),
        "FOV2": (15500, 22500, 13000, 18500),
    },
    "13": {
        "FOV1": (9000, 16000, 5500, 11000),
        "FOV2": (12500, 19500, 12000, 17500),
    },
    "14": {
        "FOV1": (9000, 16000, 8000, 13500),
        "FOV2": (12600, 19600, 15000, 20500),
    },
    "15": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "16": {
        "FOV1": (10000, 17000, 5500, 11000),
        "FOV2": (12500, 19500, 13000, 18500),
    },
    "17": {
        "FOV1": (9000, 16000, 4500, 10000),
        "FOV2": (12500, 19500, 11800, 17300),
    },
    "18": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 12500, 18000),
    },
    "19": {
        "FOV1": (10000, 17000, 5500, 11000),
        "FOV2": (13600, 20600, 12200, 17700),
    },
    "20": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "21": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 12000, 17500),
    },
    "22": {
        "FOV1": (9000, 16000, 4500, 10000),
        "FOV2": (12500, 19500, 11500, 17000),
    },
    "23": {
        "FOV1": (9000, 16000, 5000, 10500),
        "FOV2": (12500, 19500, 12000, 17500),
    },
    "24": {
        "FOV1": (9500, 16500, 6000, 11500),
        "FOV2": (12000, 19000, 13000, 18500),
    },
}

markers = [
    'DAPI', 'CD3', 'CD45', 'CD15', 'CD4', 'CD8', 'CD11b', 'Tox.Tox2', 
    'CD20', 'CD11c', 'H3K27me3', 'Ki-67', 'HLA-DRA', 'pS6', 'CD68', 
    'DC-SIGN', 'FoxP3', 'PD-1', 'Pax5', 'H3K27ac', 'GranzymeB', 'CD31', 
    'IDO-1', 'CD56', 'CD138', 'NAK-ATPase', 'CD45RA', 'Cytokeratin'
]

# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    # Process all FOVs for all slides
    process_and_segment(data_folder, output_folder, crop_coords_dict, markers)
    
    # Or process specific FOVs only:
    # process_and_segment(data_folder, output_folder, crop_coords_dict, markers, fovs=["FOV1"])
    # process_and_segment(data_folder, output_folder, crop_coords_dict, markers, fovs=["FOV2"])
