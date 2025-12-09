# How to run:

## 1️⃣ 01_clustering.py

**Purpose:**  
Performs unsupervised clustering (PhenoGraph, k=200) per slide, generating heatmaps and GeoJSON files.

**Inputs:**  
- `input/h5ad/slide_{XX}_adata.h5ad`  
- `input/geojson/{slide_prefix}.geojson`  

**Outputs:**  
- `output/slide_{XX}_adata_k=200/`  
  - `clustering_results/*.csv` — cluster membership per cell  
  - `heatmap_{slide}_dual_raw.png` — raw marker heatmap  
  - `heatmap_{slide}_dual_vmin_vmax.png` — z-score heatmap  
  - `geojson/` — updated GeoJSONs with cluster assignments  

**Run manually:**  
```bash
python 01_clustering.py <slide_number>
```

**Run on SLURM:**  
```bash
sbatch 01_submit_clustering.slurm
```

---

## 2️⃣ 02_annotation.py

**Purpose:**  
Loads `cluster_annotations.json` and applies cell type labels to each cluster.

**Inputs:**  
- `cluster_annotations.json`  
- Clustering results from Step 1 (`clustering_results/`)  
- `input/h5ad/*.h5ad`  

**Outputs:**  
- Annotated clustering heatmaps:  
  - `heatmap_*_full_annotations.png/svg`  
  - `heatmap_*_combined_annotations.png/svg`  
- Updated clustering objects with annotations  

**Run manually:**  
```bash
python 02_annotation.py
```

**Run on SLURM:**  
```bash
sbatch 02_submit_annotation.slurm
```

---

## 3️⃣ 03_phenotype_map.py

**Purpose:**  
Renders per-FOV phenotype maps using segmentation masks and annotations.

**Inputs:**  
- `output/slide_{XX}_adata_k=200/clustering_results/*.csv`  
- `Mesmer_overlay/*.tiff` and `Mesmer_mask/*.tiff`  

**Outputs:**  
- `output/slide_{XX}_adata_k=200/*.png, *.svg` — per-FOV colored phenotype maps  
  (each cell colored by its annotation)

**Run manually:**  
```bash
python 03_phenotype_map.py --slide 7
```

**Run on SLURM:**  
```bash
sbatch 03_submit_phenotype_map.slurm
```

---

## 4️⃣ 04_stacked_bar_plots.py

**Purpose:**  
Generates stacked bar charts showing the fractional composition of annotations per slide.

**Inputs:**  
- `cluster_annotations.json`  
- `input/h5ad/*.h5ad`  
- Annotated clustering results from Step 2  

**Outputs:**  
- `annotation_mixed_vs_single_vs_unknown.svg`  
  - Stacked bar plot with annotation fractions across slides  
  - Categories: cell type, Mixed, Unknown  

**Run manually:**  
```bash
python 04_stacked_bar_plots.py
```

**Run on SLURM:**  
```bash
sbatch 04_submit_stacked_bar_plot.slurm
```

---

## 5️⃣ 05_enrichment_plots.py

**Purpose:**  
Computes cell-type enrichment and produces comparative plots across conditions.

**Inputs:**  
- Annotated clustering results (`output/*/`)  
- Cell counts per annotation  

**Outputs:**  
- Enrichment heatmaps (`*.png`, `*.svg`)  
- CSV tables with enrichment values per condition  

**Run manually:**  
```bash
python 05_enrichment_plots.py
```

**Run on SLURM:**  
```bash
sbatch 05_submit_enrichment.slurm
```
