# Benchmark Tissue Preparation

### Table of Contents
- [Environment](#environment)
- [Reproducing Manuscript Figures](#reproducing)
- [Directory Structure](#directory)
- [Contributor](#contributor)

<a name="environment"></a>
# Environment
1. Python: should be compatible with Python ^3.9.12
2. R: should be equal to or greater than 4.3.2
3. Pre-installation:
   
 - Python environment:
```
pip install -r ./requirements.txt
```
 - R environment:
```
install.packages(c("dplyr", "tidyverse", "matrixStats", "ggcorrplot", "ggpubr", "tidyr", "rstatix", "readr", "svglite", "devtools"))

devtools::install_github("immunogenomics/presto")
```

<a name="reproducing"></a>
# Reproducing Manuscript Figures

To reproduce the figures in manuscript, please first download original data from [Zenodo](https://doi.org/10.5281/zenodo.11391050) to `./data` folder.

 1. Run the `segmentation_scFeature_extraction_4slides.ipynb` notebook:

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
# Directory Structure

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
The default data and output directory can be overriden in python notebook `segmentation_scFeature_extraction_4slides.ipynb`:

```
# To set data folder and out folder, un-comment the following code
# os.environ["DATA_FOLDER"] = "/project/temp/Benchmarking_tissue_preparation_data/"
# os.environ["OUT_FOLDER"] = "/project/temp/Benchmarking_tissue_preparation_out/"
```
and in `mesmer_data_preprocessing.R` and `cX2_data_preprocessing.R`:
```
data_folder = "/project/temp/Benchmarking_tissue_preparation_data/"
out_folder = "/project/temp/Benchmarking_tissue_preparation_out/"
```

<a name="contributor"></a>
# Contributor
* [Johanna Schaffenrath](https://github.com/johannaschaffenrath)
* [Shaohong Feng](https://github.com/fengsh27)
* [Cankun Wang](https://github.com/Wang-Cankun)
