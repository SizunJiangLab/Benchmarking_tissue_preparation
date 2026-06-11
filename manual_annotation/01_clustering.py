import pandas as pd
import importlib
from pathlib import Path
import numpy as np
import random
import anndata as ad
from IPython.core.interactiveshell import InteractiveShell
from clustering import subcluster
import re
import sys

import json
from pathlib import Path
import geopandas as gpd
import glob

importlib.reload(subcluster)

if True:
    from clustering.subcluster import (
        ClusteringResultManager,
        plot_clustering_heatmap_2,
        run_clustering,
        update_geojson_from_clustering_result,
    )

InteractiveShell.ast_node_interactivity = "all"

random.seed(0)
np.random.seed(0)

markers_all = [
    "CD3",
    "CD15",
    "CD8",
    "CD20",
    "CD11c",
    "CD68",
    "FoxP3",
    "Pax5",
    "CD31",
    "Cytokeratin",
]

# Get slide number from command-line arguments
slide_index = int(sys.argv[1])  # e.g., 1, 2, ..., 24
slide_filename = f"slide_{slide_index:02d}_adata.h5ad"

input_dir = Path("/registered_report/input/h5ad")
output_root = Path("/registered_report/output")

h5ad_path = input_dir / slide_filename

if not h5ad_path.exists():
    print(f"File {h5ad_path} not found. Exiting.")
    sys.exit(1)

slide_name = h5ad_path.stem  # e.g., slide_01_adata
print(f"Processing {slide_name}")

adata = ad.read_h5ad(h5ad_path)

output_dir = output_root / f"{slide_name}_k=200"
output_dir.mkdir(parents=True, exist_ok=True)

present_markers = [m for m in markers_all if m in adata.var_names]
missing_markers = [m for m in markers_all if m not in adata.var_names]
if missing_markers:
    print(f"Warning: markers not found and will be skipped: {missing_markers}")
if not present_markers:
    print("Error: none of the requested markers were found in adata.var_names. Exiting.")
    sys.exit(1)

features = present_markers
unit_ids = adata.obs.index
manager = ClusteringResultManager(output_dir=output_dir, unit_ids=unit_ids)

# Run clustering
clustering_result = run_clustering(
    adata,
    unit_ids,
    features,
    method="phenograph",
    method_params={"k": 200, "n_jobs": 8},
    output_dir=output_dir,
)

unit_ids = np.unique(clustering_result.unit_ids)
matches = np.array(["_".join(s.split("_")[:2]) for s in unit_ids])
matches = np.unique(matches)

for unit_id_prefix in matches:
    n_match = len(
        [
            unit_id
            for unit_id in clustering_result.unit_ids
            if unit_id.startswith(unit_id_prefix)
        ]
    )
    if n_match > 0:
        print(f"Found {n_match} units starting with {unit_id_prefix} in this clustering")
    else:
        print(f"No unit found starting with '{unit_id_prefix} in this clustering")

heatmap_raw = plot_clustering_heatmap_2(
    adata,
    clustering_result,
    features,
    col_cellsize="cellSize",
    figsize=(20, 8),
    col_gap=30,
    legend_hpad=60,
)
heatmap_raw.savefig(f"{output_dir}/heatmap_{slide_name}_dual_raw.png", dpi=300, bbox_inches="tight")


non_col_cluster = {"col_cluster": False, "col_dendrogram": False}
heatmap_vmin_vmax = plot_clustering_heatmap_2(
    adata,
    clustering_result,
    features,
    col_cellsize="cellSize",
    figsize=(30, 8),
    col_gap=30,
    legend_hpad=60,
    kwargs_zscore={"vmin": -3, "center": 0, "vmax": 3} | non_col_cluster,
    kwargs_mean={"vmin": 0, "vmax": 0.5} | non_col_cluster,
)
heatmap_vmin_vmax.savefig(f"{output_dir}/heatmap_{slide_name}_dual_vmin_vmax.png", dpi=300, bbox_inches="tight")

for prefix in matches:
    print(f"/registered_report/input/geojson/{prefix}.geojson")
    geojson_file = Path(
        f"/registered_report/input/geojson/{prefix}.geojson"
    )
    update_geojson_from_clustering_result(
        geojson_file,
        clustering_result,
        output_dir,
        unit_id_prefix=f"{prefix}_",
    )

    matches = glob.glob(f"{output_dir}/geojson/*/{prefix}.geojson")
    if not matches:
        raise FileNotFoundError(f"No geojson found for prefix {prefix}")

    in_path = matches[0]
    gdf = gpd.read_file(in_path)
    out_dir = Path(f"/registered_report/output/{slide_name}_k=200/per_class_geojson")
    out_dir.mkdir(parents=True, exist_ok=True)


    def to_class_dict(val):
        if isinstance(val, dict):
            return val
        if isinstance(val, str):
            val = val.strip()
            if val.startswith("{"):
                try:
                    return json.loads(val)
                except json.JSONDecodeError:
                    pass
        return None  # unparseable / missing

    gdf["classification_dict"] = gdf["classification"].apply(to_class_dict)

    # Drop rows without a valid classification dict
    gdf = gdf[gdf["classification_dict"].notna()].copy()

    # Normalize fields for grouping
    def to_key(d):
        # d like {"name": "11", "color": [242, 72, 72]}
        name = str(d.get("name", "unknown"))
        color = tuple(d.get("color") or [])  # make hashable
        return (name, color)

    gdf["class_key"] = gdf["classification_dict"].apply(to_key)

    summary = []
    for (name, color), sub in gdf.groupby("class_key"):
        # build a safe filename: e.g., name-11_color-242-72-72.geojson
        color_str = "-".join(map(str, color)) if color else "none"
        fname = f"name-{name}_color-{color_str}_prefix-{prefix}.geojson"
        out_path = out_dir / fname

        # write just these features
        sub.drop(columns=[c for c in ["classification_dict", "class_key"] if c in sub.columns]) \
        .to_file(out_path, driver="GeoJSON")

        summary.append((name, color, len(sub), str(out_path)))

    # --- print summary ---
    print("Written files:")
    for name, color, n, p in summary:
        print(f"  name={name:>4}  color={list(color)}  count={n:>6}  -> {p}")

# Stash result
clustering_result.stash(output_dir=output_dir)

print(f"Finished: {slide_name}\n")
