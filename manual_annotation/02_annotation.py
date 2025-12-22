import pandas as pd
import importlib
from pathlib import Path
import numpy as np
import random
import anndata as ad
import json
import re
import sys

from clustering import subcluster
importlib.reload(subcluster)

from clustering.subcluster import (
    ClusteringResult,
    ClusteringResultManager,
    plot_clustering_heatmap_2,
)

# Set seeds for reproducibility
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

def reduce_annotation_for_category(label: str) -> str:
    if label == "Other":
        return "Other"
    if "and" in label:
    #if "and" in label or label == "T cells":
        return "Mixed"
    return label  # keep actual annotation name (e.g. "B cells", "Tregs", etc.)


# Paths
annotation_json_path = Path("/registered_report/cluster_annotations.json")
output_root = Path("/registered_report/output")
input_dir = Path("/registered_report/input/h5ad")

# Load annotations
with open(annotation_json_path) as f:
    annotations_all = json.load(f)

for slide, cluster_dict in annotations_all.items():
    for k in list(cluster_dict.keys()):
        if cluster_dict[k] == "Artifacts":
            annotations_all[slide][k] = ""

per_slide_annotation_counts = {}
# Process each slide in annotations
for slide_key, annotations in annotations_all.items():
    print(f"Processing {slide_key}...")
    print(annotations)
    print(slide_key)

    slide_number = str(int(slide_key.split("_")[1]))  # '01' -> '1'

    # Setup paths
    slide_filename = f"{slide_key}_adata.h5ad"
    h5ad_path = input_dir / slide_filename
    output_dir = output_root / f"{slide_key}_adata_k=200"

    if not h5ad_path.exists():
        print(f"H5AD file {h5ad_path} not found. Skipping.")
        continue

    # Load data
    adata = ad.read_h5ad(h5ad_path)
    # derive features from the fixed list, intersected with adata
    present_markers = [m for m in markers_all if m in adata.var_names]
    missing_markers = [m for m in markers_all if m not in adata.var_names]
    if missing_markers:
        print(f"Warning for {slide_key}: markers not found and skipped: {missing_markers}")
    if not present_markers:
        print(f"No requested markers present in {slide_key}; skipping.")
        continue

    features = present_markers

    # Determine clustering_id by checking geojson folder
    geojson_dir = output_dir / "geojson"
    clustering_id = None

    if geojson_dir.is_dir():
        clustering_dirs = [d for d in geojson_dir.iterdir() if d.is_dir()]
        if clustering_dirs:
            newest_dir = max(clustering_dirs, key=lambda d: d.stat().st_mtime)
            clustering_id = newest_dir.name

    if not clustering_id:
        print(f"No clustering ID folder found in {geojson_dir}, skipping.")
        continue


    # Load clustering result
    clustering_result = ClusteringResult.pop(
        clustering_id=clustering_id,
        output_dir=output_dir
    )

    df = clustering_result.cluster_df
    print("Empty cluster_ids:", df["cluster_ids"].isna().sum())
    print("Blank string cluster_ids:", (df["cluster_ids"] == "").sum())
    print("Cluster ID counts:\n", df["cluster_ids"].value_counts(dropna=False))

    print("Unique cluster IDs in result:")
    print(sorted(clustering_result.cluster_df["cluster_ids"].unique()))
    print("Type of one cluster ID:", type(clustering_result.cluster_df["cluster_ids"].iloc[0]))

    print("Keys in annotations:", list(annotations.keys()))
    print("Types in annotations:", {type(k) for k in annotations.keys()})

    # Normalize both to clean string keys
    cluster_ids = clustering_result.cluster_df["cluster_ids"].astype(str).str.strip().unique()
    normalized_annotations = {str(k).strip(): v for k, v in annotations.items()}

    # Check coverage
    missing = set(cluster_ids) - set(normalized_annotations)
    if missing:
        print(f"Missing annotations for cluster IDs: {missing}")


    print("N clusters in result:", clustering_result.cluster_df["cluster_ids"].nunique())
    print("N unique annotations passed:", len(annotations))

    # Apply annotation
    clustering_result.add_annotation(normalized_annotations)

    non_col_cluster = {"col_cluster": False, "col_dendrogram": False}
    heatmap_full_annotations = plot_clustering_heatmap_2(
        adata,
        clustering_result,
        features,
        figsize=(30, 8),
        x_label="tag",
        col_cellsize="cellSize",
        col_gap=30,
        legend_hpad=60,
        kwargs_zscore={"vmin": -3, "center": 0, "vmax": 3} | non_col_cluster,
        kwargs_mean={"vmin": 0, "vmax": 0.5} | non_col_cluster,
    )
    heatmap_full_annotations.savefig(output_dir / f"heatmap_{slide_filename}_full_annotations.png", dpi=300, bbox_inches="tight")

    # Save result
    clustering_result.save(output_dir)

    # Create summary manager and filtered result
    manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)
    print(f"Unannotated units: {len(manager.non_explicit_df)}")
    print(manager.summary_df.annotation.value_counts())

    # Add per-annotation breakdown (excluding unannotated)
    annotation_counts_clean = manager.summary_df.annotation[manager.summary_df.annotation != ""].value_counts()

    # Normalize to string for consistent JSON keys (optional)
    annotation_breakdown = {k: int(v) for k, v in annotation_counts_clean.items()}

    per_slide_annotation_counts[slide_key] = annotation_breakdown

    reduced_annotation_counts = {}

    for slide, counts in per_slide_annotation_counts.items():
        reduced = {}
        for label, count in counts.items():
            reduced_label = reduce_annotation_for_category(label)
            if reduced_label == "":  # ignore empty/unannotated
                continue
            reduced[reduced_label] = reduced.get(reduced_label, 0) + count
        reduced_annotation_counts[slide] = reduced

    idx_annotated = manager.summary_df.annotation != ""
    result = ClusteringResult(
        clustering_id="annotated_summary",
        method="manual_annot",
        unit_ids=manager.summary_df.index[idx_annotated],
        cluster_ids=manager.summary_df.annotation[idx_annotated],
    )

    # Plot heatmap of combined annotations
    heatmap_combined_annotations = plot_clustering_heatmap_2(
        adata,
        result,
        features,
        figsize=(20, 8),
        col_cellsize="cellSize",
        x_label="annotation",
        kwargs_zscore={"vmin": -3, "center": 0, "vmax": 3} | non_col_cluster,
        kwargs_mean={"vmin": 0, "vmax": 0.5} | non_col_cluster,
    )
    heatmap_combined_annotations.savefig(output_dir / f"heatmap_{slide_filename}_combined_annotations.png", dpi=300, bbox_inches="tight")

    # Full annotation

    heatmap_full_annotations.savefig(output_dir / f"heatmap_{slide_filename}_full_annotations.svg", bbox_inches="tight")

    heatmap_combined_annotations.savefig(output_dir / f"heatmap_{slide_filename}_combined_annotations.svg", bbox_inches="tight")

print("Done.")