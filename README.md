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

# To run balagan_analysis.R:
devtools::install_github("PierreBSC/Balagan")
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

2. Run `mesmer_data_preprocessing.R`:

This script produces Fig. 1C, Fig. 1D, Supp Fig. 1F. It processes extracted per-cell signals (in .csv format) from Mesmer for outlier removal, normalization, transformation, and statistical analysis. Visualizations include density plots and heatmaps. 

This script is compatible with R 4.3.2.

3. Run `cX2_data_preprocessing.R`:

This script produces Supp Fig. 1G. It processes extracted per-cell signals (in .csv format) from cellXpress2 for outlier removal, normalization, transformation. Visualizations include density plots. 

4. Run `balagan_analysis.R`:

This script produces Supp Fig. 1I. It processes "dataScaleSize_slide2_FOV1.csv" using the balagan R package to run spatial clustering and subsampling analysis. The heatmap, clustering scatter plot, and subsampling plot were generated. 

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
|---balagan_analysis.R: R script to generate the supplementary fig 1I using the balagan R package
|---Random.seed.txt: Random seed used to reproduce balagan_analysis.R
|---README.md: this file
```

<a name="contributor"></a>
# Contributor
* [Johanna Schaffenrath](https://github.com/johannaschaffenrath)
* [Cankun Wang](https://github.com/Wang-Cankun)
* [Shaohong Feng](https://github.com/fengsh27)

If you have any questions or feedback regarding this repository, please contact Sizun Jiang via sjiang3@bidmc.harvard.edu.
