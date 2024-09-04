# Benchmark Tissue Preparation

### Table of Contents
- [Abstract](#abstract)
- [Environment](#environment)
- [Reproducing manuscript figure](#reproducing)
- [Directory structure](#directory)

<a name="abstract"></a>
# Abstract
 
The advent of multiplexed protein imaging has revolutionized the study of biological processes within the native tissue context by enabling the simultaneous detection of over 50 protein markers on single cells, an approach often termed “spatial proteomics”. These approaches, now commercially available with widespread adoption, have advanced fields ranging from developmental biology to viral pathogenesis and cancer biology. With the increased adoption of multiplexed protein imaging, particularly in clinical archival tissues, a robust assessment of various tissue preparation steps, including antigen retrieval and antibody staining conditions, on antibody performance is direly needed to guide the development of standardized and robust protocols for maximal performance and reduced batch-to-batch variability. Here, we propose to perform a large-scale assessment of 24 tissue preparation and staining conditions initially on the commercially available co-detection by indexing (CODEX) fluorescence multiplexed imaging platform, on adjacent formalin-fixed paraffin-embedded sections of human tonsil tissue stained using a commercial-grade, immune-focused 23-plex antibody panel. We will then perform a quantitative image analysis to evaluate the impact of these conditions on marker performance, including specificity and signal-to-noise ratio. To assess robustness and reproducibility, we will then repeat the top and bottom 2 conditions on other commercial platforms (including the Ionpath MIBI and Lunaphore Comet) across multiple academic and pharmaceutical institutions. This resource will be invaluable in guiding adopters of multiplexed imaging in the design of experiments, selection of reagents and experimental conditions, and quantitative evaluation of results. 


The following is an overview of the code in this repository. Original data to reproduce the figures can be downloaded from Zenodo (Link to pilot data: https://doi.org/10.5281/zenodo.11391050). To reproduce the analysis, file paths have to be changed accordingly.  



| Name | Description |
| -------------------------------- | -------------------------------- | 
| segmentation_scFeature_extraction_4slides.ipynb | This script uses stitched and background subtracted images in qptiff format and is designed for segmentation and single cell feature extraction of four slides simultaneously. It includes cropping of ROIs, segmentation using Mesmer, extraction of single cell features, and export as .csv files. In this study’s pipeline, these files are subjected to preprocessing and feature comparison using “Data Preprocessing.R”. Produces cell masks for Fig. 1E and input for Fig. 1C, Fig. 1D, Supp Fig. 1F |
| Data Preprocessing.R | This script uses extracted per-cell signals (in .csv format) from either Mesmer or cellXpress as input for subsequent outlier removal, mean nuclear signal normalization, arcsinh transformation, and universal percentile normalization. This code also performs all statistical tests including Kruskal-Wallis with Dunn’s post-hoc tests and Benjamin-Hochberg correction. Visualizations include density plots of signal distribution for each marker and heatmaps of Z-scores of mean marker expression across all conditions. Produces Fig. 1C, Fig. 1D, Supp Fig. 1F, Supp Fig. 1G. R version 4.3.2 was used in this study. |

<a name="environment"></a>
# Environment
1. Python: should be compatible with Python ^3.9.12
2. R: should be equal to or greater than 4.3.2
3. Pre-install:
   
 - Python environment:
```
pip install -r ./requirements.txt
```
 - R environment:
```
install.packages(c("dplyr", "tidyverse", "matrixStats", "ggcorrplot", "ggpubr", "tidyr", "rstatix", "readr", "svglite", "devtools"))
```
```
devtools::install_github("immunogenomics/presto")
```

<a name="reproducing"></a>
# Reproducing manuscript figure

To reproduce the figures in manuscript, please first download original data from [Zenodo](https://doi.org/10.5281/zenodo.11391050) to `./data` folder.

 1. Run segmentation_scFeature_extraction_4slides.ipynb:

This script generate cell masks for Fig. 1E and input for Fig. 1C, Fig. 1D, Supp Fig. 1F.

It processes stitched and background subtracted images in qptiff format and is designed for segmentation and single cell feature extraction of four slides simultaneously. It processes stitched and background-subtracted images in qptiff format for segmentation and single-cell feature extraction across four slides simultaneously.

The results are written to `./out/extracted_features/dataScaleSize_slide{1,2,3,4}.csv` and `./out/extracted_features/data_slide{1,2,3,4}.csv`. 

**Note**: the result files `dataScaleSize_slide{1,2,3,4}.csv` have been included in the original data downloaded from [Zenodo](https://doi.org/10.5281/zenodo.11391050).

This notebook is compatible with python 3.9.12

2. Run mesmer_data_preprocessing.R:

This script produces Fig. 1C, Fig. 1D, Supp Fig. 1F. It processes extracted per-cell signals (in .csv format) from Mesmer for outlier removal, normalization, transformation, and statistical analysis. Visualizations include density plots and heatmaps. 

This script is compatible with R 4.3.2.

3. Run cX2_data_preprocessing.R:

This script produces Supp Fig. 1G. It processes extracted per-cell signals (in .csv format) from cellXpress2 for outlier removal, normalization, transformation. Visualizations include density plots. 

<a name="directory"></a>
# Directory structure

```
Benchmark_tissue_preparation
|---data: default data directory
|---out: default output directory
|---pylibs
| |---tissue_preparation.py: python module containing reusable functions
|---R-scripts
| |---helper.R: R scripts containing reusable functions
|---requirements.txt: dependent python packages
|---segmentation_scFeature_extraction_4slides: python notebook
|---mesmer_data_preprocessing.R: R script to generate the figures for MESMER cell segmentations
|---cX2_data_preprocessing.R: R script to generate the figures for cellXpress2 cell segmentations
|---README.md: this file
```
The default data and output directory can be overriden in python notebook "segmentation_scFeature_extraction_4slides":

```
# To set data folder and out folder, un-comment the following code
# os.environ["DATA_FOLDER"] = "/project/temp/Benchmarking_tissue_preparation_data/"
# os.environ["OUT_FOLDER"] = "/project/temp/Benchmarking_tissue_preparation_out/"
```
and in mesmer_data_preprocessing.R and cX2_data_preprocessing.R:
```
data_folder = "/project/temp/Benchmarking_tissue_preparation_data/"
out_folder = "/project/temp/Benchmarking_tissue_preparation_out/"
```

