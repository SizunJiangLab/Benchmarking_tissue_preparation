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
  - `clustering_stash/*.pickle` — clustering objects for future reloading   
  - `clustering_sequence.txt` — log of the clustering order
  - `heatmap_{slide}_dual_raw.png` — raw marker heatmap  
  - `heatmap_{slide}_dual_vmin_vmax.png` — z-score heatmap  
  - `geojson/` — updated GeoJSONs with cluster assignments  
  - `per_class_geojson/` - per each cluster per FOV separate geojson file

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
- `MESMER_overlay/*.tiff` and `MESMER_mask/*.tiff`  

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
- Categories: cell type, Mixed, Other  

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

**Outputs:**  
- log2 ratio per slide of marker enrichment (`*.png`, `*.svg`) 
- CSV tables with enrichment values per condition (`*.csv`)

**Run manually:**  
```bash
python 05_enrichment_plots.py
```

**Run on SLURM:**  
```bash
sbatch 05_submit_enrichment.slurm
```

---

## 6️⃣ 06_umap_all_slides.py

**Purpose:**  
Concatenates all 24 slides, attaches PhenoGraph cluster IDs and cell-type
annotations from the per-slide clustering CSVs, and runs UMAP. Provides the
resolved cell-type annotations on dimensionality-reduction plots requested by
the reviewers. Two analysis modes are produced: `all_cells` and
`filtered_annotations` (keeping only cells with a definite annotation:
CD8+ T cells, CD8- T cells, Tregs, B cells, Epithelial, Macrophages).
SVGs are exported with editable text for Illustrator.

**Inputs:**  
- `input/h5ad/slide_{XX}_adata.h5ad`  
- `output/slide_{XX}_adata_k=200/clustering_results/*.csv`  

**Outputs:**  
- `output/umap/all_cells/` and `output/umap/filtered_annotations/`, each with:  
  - `combined_umap.h5ad` — combined AnnData object for the analysis mode
  - `umap_by_cluster.svg` — UMAP colored by PhenoGraph cluster ID
  - `umap_by_annotation.svg` — UMAP colored by cell-type annotation
  - `umap_by_slide.svg` — UMAP colored by slide ID for batch checking
  - `umap_marker_expression_grid.svg` — grid of UMAPs colored by marker intensity

**Run manually:**  
Paths default to this folder's `input/` and `output/`. To point at data stored
elsewhere, set `REGISTERED_REPORT_DIR` (no paths are hard-coded in the script);
column names can still be adjusted in the `CONFIG` block. Then run:
```bash
python 06_umap_all_slides.py
```

---

## 7️⃣ 07_umap_all_slides_individual_stitched.py

**Purpose:**  
Generates UMAP visualizations for each slide individually and stitches them
into A4 pages for side-by-side comparison across conditions. Each per-slide
panel pairs a cell-type annotation UMAP with a marker-expression grid.

**Inputs:**  
- `input/h5ad/slide_{XX}_adata.h5ad`  
- `output/slide_{XX}_adata_k=200/clustering_results/*.csv`  

**Outputs:**  
- `output/umap/slide_{XX}/`  
  - `slide_{XX}_umap.h5ad` — slide-specific AnnData with UMAP coordinates and metadata
  - `slide_{XX}_umap_annotation_and_markers.svg` — annotation UMAP (left) plus marker-expression panels (right)
- `output/umap/slides_01_04_A4.svg` … `slides_21_24_A4.svg` — A4 portrait pages stitching four slides each

**Run manually:**  
Paths default to this folder's `input/` and `output/`. To point at data stored
elsewhere, set `REGISTERED_REPORT_DIR` (no paths are hard-coded in the script);
column names can still be adjusted in the `CONFIG` block. Then run:
```bash
python 07_umap_all_slides_individual_stitched.py
```
