from typing import Tuple, Optional, List
import matplotlib.pyplot as plt # for plotting images
from os import path
from tifffile import tifffile
import imagecodecs
import numpy as np
import pandas as pd
import glob
from skimage.io import imsave, imread
import skimage
import skimage.io
import skimage.measure
import skimage.morphology
from tqdm import tqdm

from deepcell.applications import Mesmer # Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay

nuclear_markers = [0] # these are indices of the channels we want to use as the nuclear signal (only one here)
membrane_markers = [2, 6, 12, 27] # CD45, HLA-DR, Cytokeratin, CD11b these are the indices we want to use as the membrane signal (referenced above)


def read_qtiff_image(fn: str)->np.ndarray: 
    """
    Read qtiff image
    Args:
        fn (str): image file path
    Returns:
        ndarray: image content in
    """
    if not path.exists(fn):
        raise ValueError(f"File {fn} not found")
    
    img = glob.glob(fn)
    return tifffile.imread(img[0])

def generate_nuclear_and_membrane(img: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate nuclear array and membrane array
    Args:
        img (ndarray): image content
    Returns:
        ndarray: nuclear array and membrane array
    """
    shp = (img.shape[1], img.shape[2])
    nuclear = np.zeros(shp)
    membrane = np.zeros(shp)

    for chn_index in range(len(img)):
        arr = img[chn_index, :, :] # slice the corresponding index in the first dimension
        if chn_index in nuclear_markers: # if the chn_index is in our list of nuclear markers, stack it to the nuclear array
            nuclear = np.add(nuclear, arr)
        elif chn_index in membrane_markers: # if the chn_index is a membrane_marker, stack it to the membrane array
            membrane = np.add(membrane, arr)
        else:
            pass

    return nuclear, membrane

def stack_nuclear_and_membrane(nuclear: np.ndarray, membrane: np.ndarray) -> np.ndarray:
    """
    Stack the nuclear and membrane arrays into one array
    Args:
        nuclear (ndarray): nuclear array
        membrane (ndarray): membrane array
    Returns:
        ndarray: Stacked array
    """
    # stack
    stack = np.stack((nuclear, membrane), axis=-1)
    # also expand to 4 dimensions
    return np.expand_dims(stack, 0)

def run_segmentation(
    img: np.ndarray, 
    image_mpp: Optional[float] = 0.50,
    maxima_threshold: Optional[float] = 0.075, 
    interior_threshold: Optional[float] = 0.2,
):
    """
    Identify cell segmentation with MESMER
    Args:
        image (ndarray): image content
        image_mpp (Optional, float): microns per pixel resolution
        maxima_threshold (Optional, float): controls what is considered a 
            unique cell (lower values = more separate cells, higher values
             = fewer cells)
        interior_threshold (Optional, float): determines what is considered 
            background/not part of a cell (lower value = larger cells)
    Returns:
        ndarray: Image content with cell segmentation
    """
    maxima_threshold = maxima_threshold
    interior_threshold = interior_threshold
    app = Mesmer()
    return app.predict(
        img, 
        image_mpp=image_mpp,
        postprocess_kwargs_whole_cell={
            "maxima_threshold": maxima_threshold,
            "interior_threshold": interior_threshold,
        },
    )

def extract_sc_features(
    img: np.ndarray,
    mask: np.ndarray,
    markers: List[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract sc features from image
    Args:
        img (ndarray): image content
        mask (ndarray): MESMER mask
        markers (List[str]): markers
    
    Returns:
        Tuple[ndarray, ndarray, ndarray, ndarray]: data - array of marker expression sum for cells
                                                   dataScaleSize: array of scaled marker expression sum for cells
                                                   cell_props: array of cell props
                                                   cellSizes: array of cell sizes
    """
    stats = skimage.measure.regionprops(mask)
    cell_count = len(stats) # number of actual cells not always equal to np.max(mask) 
    marker_count = len(markers)

    # empty containers of zeros
    data = np.zeros((cell_count, marker_count))
    dataScaleSize = np.zeros((cell_count, marker_count))
    cellSizes = np.zeros((cell_count, 1))
    cell_props = np.zeros((cell_count, 3))

    # extract info
    for i in tqdm(range(cell_count)): # tqdm creates the progress bar
        cellLabel = stats[i].label
        label_counts = [img[coord[0], coord[1], :] for coord in stats[i].coords] # all markers for this cell
        data[i, 0:marker_count] = np.sum(label_counts, axis = 0) # sum the marker expression for this cell
        dataScaleSize[i, 0:marker_count] = np.sum(label_counts, axis = 0) / stats[i].area # scale the sum by size
        cellSizes[i] = stats[i].area # cell size
        cell_props[i, 0] = cellLabel
        cell_props[i, 1] = stats[i].centroid[0] # Y centroid
        cell_props[i, 2] = stats[i].centroid[1] # X centroid
    
    return data, dataScaleSize, cell_props, cellSizes

if __name__ == "__main__":
    img = read_qtiff_image("/home/ubuntu/project/temp/Benchmarking_tissue_preparation_data/Slide 1_20 min HIER 1h RT stain_Scan1.qptiff")
    n, m = generate_nuclear_and_membrane(img)
    print(n)
    print(m)
