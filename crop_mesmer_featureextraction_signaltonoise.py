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

def process_and_segment(data_folder, output_folder, crop_coords_dict, markers):
    os.makedirs(output_folder, exist_ok=True)

    data_folder = Path(data_folder)
    all_qptiffs = list(data_folder.glob("*.qptiff"))

    for key in crop_coords_dict.keys():
        for qptiff_f in all_qptiffs:
            match = re.search(rf"(^.*Registered_Report_Slide{key}.qptiff)", str(qptiff_f))
            if match:
                break
        
        if match:
        # If the match is successful, construct the file path with the correct Slide{key} part
            qptiff_path = os.path.join(data_folder, f"{match.group(0)}")  # match.group(0) gives the full match
            print(f"Path: {qptiff_path}")
    
        else:
            print(f"No match for key: {key}")


        #if key in ["18", "22"]:
            #qptiff_path = os.path.join(data_folder, f"slide{key}_Scan2.qptiff")
        #else:
            #qptiff_path = os.path.join(data_folder, f"Slide{key}_Scan1.qptiff")


        
        output_crop_path = os.path.join(output_folder, f"cropped_stack{key}.tiff")
        output_segmentation_path = os.path.join(output_folder, f"segmentation_{key}.tiff")

        if not os.path.exists(qptiff_path):
            print(f"⚠️ Missing qptiff file for {key}, skipping.")
            continue

        print(f"Processing {key}...")

        # Read qptiff image
        slide = read_qtiff_image(qptiff_path)
        print(f"  - Slide shape: {slide.shape}")

        # Generate nuclear and membrane images
        nuclear, membrane = generate_nuclear_and_membrane(slide)

        # Use imported function to create the stack
        stack = stack_nuclear_and_membrane(nuclear=nuclear, membrane=membrane)
        print(f"  - Stacked shape: {stack.shape}")  

        # Extract crop coordinates
        x_min, x_max, y_min, y_max = crop_coords_dict[key]

        # Validate cropping boundaries
        _, height, width, _ = stack.shape
        if x_min < 0 or y_min < 0 or x_max > width or y_max > height:
            print(f"⚠️ Invalid crop coordinates for {key}: {crop_coords_dict[key]}")
            continue

        cropped_stack = stack[:, y_min:y_max, x_min:x_max, :]

        # Check if cropped stack is empty
        if cropped_stack.size == 0:
            print(f"⚠️ Cropped stack for {key} is empty! Skipping.")
            continue

        # Run Mesmer segmentation
        predictions = run_segmentation(cropped_stack)

        # Create RGB image and overlay
        rgb_image = create_rgb_image(cropped_stack, channel_colors=["green", "blue"])
        overlay = make_outline_overlay(rgb_data=rgb_image, predictions=predictions)

        # Save Mesmer outputs
        output_dir = os.path.join(output_folder, f'Mesmer_outputs/{maxima_threshold}maxima_{key}_{interior_threshold}interior_{key}_FOV2/')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        skimage.io.imsave(os.path.join(output_dir, "seg_outline.tiff"), img_as_ubyte(overlay[0, ..., 0]), check_contrast=False)  # segmentation outline
        skimage.io.imsave(os.path.join(output_dir, "seg_overlay.tiff"), img_as_ubyte(overlay[0, ...]), check_contrast=False)  # segmentation overlay (nuc + membrane + outline)
        skimage.io.imsave(os.path.join(output_dir, "Mesmer_mask.tiff"), predictions[0, ..., 0], check_contrast=False)  # Mesmer mask

        print(f"✅ Processed, cropped, segmented, and saved results for {key}.")

        # Single-cell feature extraction
        mask = skimage.io.imread(os.path.join(output_dir, "Mesmer_mask.tiff"))

        # Transpose the stack dimensions 
        array_list = np.transpose(slide, (1, 2, 0))

        # Crop the array to match cropping done earlier
        cropped_array_list = array_list[y_min:y_max, x_min:x_max, :]

        # Extract single-cell features
        data, dataScaleSize, cell_props, cellSizes = extract_sc_features(cropped_array_list, mask=mask, markers=markers)

        data_df = pd.DataFrame(data)
        data_df.columns = markers
        data_full = pd.concat((pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]), pd.DataFrame(cellSizes, columns=["cellSize"]), data_df), axis=1)

        dataScaleSize_df = pd.DataFrame(dataScaleSize)
        dataScaleSize_df.columns = markers
        dataScaleSize_full = pd.concat((pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]), pd.DataFrame(cellSizes, columns=["cellSize"]), dataScaleSize_df), axis=1)

        # Save the dataframes
        save_dir = os.path.join(output_folder, 'extracted_features')
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        data_full.to_csv(os.path.join(save_dir, f'data_slide{key}_FOV1.csv'), index=False)
        dataScaleSize_full.to_csv(os.path.join(save_dir, f'dataScaleSize_slide{key}_FOV2.csv'), index=False)

        print(f"✅ Extracted and saved single-cell features for {key}.")

        # Export individual TIFF files for each channel for the FOV
        individual_tiff_dir = os.path.join(output_folder, f'Individualtiff_slide{key}_FOV1')
        os.makedirs(individual_tiff_dir, exist_ok=True)

        # Save individual TIFFs and collect raw signals
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

        # Initialize results
        results = []

        # Extract DAPI reference
        dapi_in = raw_signals["DAPI"]["signal_in"]
        dapi_out = raw_signals["DAPI"]["signal_out"]

        # Loop through markers and compute ratios
        for marker, vals in raw_signals.items():

            # Default NaNs in case of invalid divisions
            signal_invsout_DAPInorm = np.nan
            signal_invsout_areanorm = np.nan
            norm_ratio = np.nan

            if dapi_in > 0 and dapi_out > 0 and vals["signal_in"] > 0 and vals["signal_out"] > 0:
              # 1️⃣ DAPI-normalized ratio (no area correction)
              signal_invsout_DAPInorm = (vals["signal_in"] / dapi_in) / (vals["signal_out"] / dapi_out)

             # 2️⃣ Area-normalized ratio (no DAPI correction)
              signal_invsout_areanorm = (vals["signal_in"] / area_in) / (vals["signal_out"] / area_out)

             # 3️⃣ DAPI + Area normalized ratio (original definition)
              norm_ratio = (vals["signal_in"] / dapi_in / area_in) / (vals["signal_out"] / dapi_out / area_out)

            # 4️⃣ Handle DAPI itself (avoid dividing by itself)
            if marker == "DAPI":
                signal_invsout_DAPInorm = 1.0
                norm_ratio = 1.0

            # Append results
            results.append({
                "Marker": marker,
                "signal_invsout_DAPInorm": signal_invsout_DAPInorm,
                "signal_invsout_areanorm": signal_invsout_areanorm,
                "Normalized_signal_invsout": norm_ratio
             })

        # Convert results to a DataFrame
        df = pd.DataFrame(results)

        # Save to CSV
        output_csv_path = os.path.join(individual_tiff_dir, f"signal_ratios_slide{key}_FOV2.csv")
        df.to_csv(output_csv_path, index=False)

        print(f"✅ Ratios have been calculated and saved to {output_csv_path}")

# Example usage
data_folder = "/mnt/nfs/storage/Fusion_Registered_Report/Initial_Optimization_Stage2/BIDMC"
output_folder = "/mnt/nfs/home/jialelee/Registered_Report/Initial Optimization/BIDMC/CheckSNR_Normalizedbyarea/Original_FOV2_251031"


crop_coords_dict = {
    #"1": (10000, 17000, 5500, 11000), 
    #"2": (12500, 19500, 6250, 11750),  
    #"3": (9000, 16000, 6000, 11500), 
    #"4": (9000, 16000, 5000, 10500),
    #"5": (9500, 16500, 5500, 11000), 
    #"6": (9000, 16000, 4500, 10000), 
    #"7": (9500, 16500, 6000, 11500), 
    #"8": (9000, 16000, 5000, 10500), 
    #"9": (9000, 16000, 5000, 10500), 
    #"10": (9000, 16000, 5000, 10500), 
   #"11": (9000, 16000, 4000, 9500), 
    #"12": (12000, 19000, 6250, 11750), 
   #"13": (9000, 16000, 5500, 11000), 
   #"14": (9000, 16000, 8000, 13500), 
   #"15": (9000, 16000, 5000, 10500), 
   #"16": (10000, 17000, 5500, 11000), 
   #"17": (9000, 16000, 4500, 10000), 
   #"18": (9000, 16000, 5000, 10500), 
   #"19": (10000, 17000, 5500, 11000), 
    #"20": (9000, 16000, 5000, 10500), 
   #"21": (9000, 16000, 5000, 10500), 
    #"22": (9000, 16000, 4500, 10000), 
    #"23": (9000, 16000, 5000, 10500), 
    #"24": (9500, 16500, 6000, 11500) 
    "1": (14000, 21000, 11000, 16500),#FOV2
    "2": (15600, 22600, 12350, 17850),#FOV2
    "3": (12300, 19300, 13000, 18500),#FOV2
    "4": (12500, 19500, 11500, 17000),#FOV2
    "5": (12000, 19000, 12400, 17900),#FOV2
    "6": (12500, 19500, 11000, 16500),#FOV2
    "7": (12000, 19000, 13000, 18500), #FOV2
    "8": (12500, 19500, 11500, 17000),#FOV2
    "9": (12500, 19500, 11500, 17000), #FOV2
    "10": (12500, 19500, 11500, 17000), #FOV2
    "11": (13000, 20700, 10200, 15700),#FOV2
    "12": (15500, 22500, 13000, 18500), #FOV2
    "13": (12500, 19500, 12000, 17500), #FOV2
    "14": (12600, 19600, 15000, 20500), #FOV2
    "15": (12500, 19500, 11500, 17000),#FOV2
    "16": (12500, 19500, 13000, 18500), #FOV2
    "17": (12500, 19500, 11800, 17300), #FOV2
    "18": (12500, 19500, 12500, 18000), #FOV2
    "19": (13600, 20600, 12200, 17700), #FOV2
    "20": (12500, 19500, 11500, 17000), #FOV2
    "21": (12500, 19500, 12000, 17500), #FOV2
   "22": (12500, 19500, 11500, 17000), #FOV2
   "23": (12500, 19500, 12000, 17500), #FOV2
    "24": (12000, 19000, 13000, 18500) #FOV2
}

markers = ['DAPI', 'CD3', 'CD45', 'CD15', 'CD4', 'CD8', 'CD11b', 'Tox.Tox2', 'CD20', 'CD11c', 'H3K27me3', 'Ki-67', 'HLA-DRA', 'pS6', 'CD68', 'DC-SIGN', 'FoxP3', 'PD-1', 'Pax5', 'H3K27ac', 'GranzymeB', 'CD31', 'IDO-1', 'CD56', 'CD138', 'NAK-ATPase', 'CD45RA','Cytokeratin']

process_and_segment(data_folder, output_folder, crop_coords_dict, markers)
