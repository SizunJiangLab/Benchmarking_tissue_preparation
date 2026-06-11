"""
06_umap_all_slides.py
-----------------------
Concatenates all per-slide h5ad files, attaches cluster IDs and annotations
from the per-slide clustering CSVs, then plots UMAP coloured by:
  1. Cluster ID
  2. Cell-type annotation
  3. Slide ID  (batch check)

Directory layout expected:
  h5ad:  <H5AD_DIR>/slide_NN_adata.h5ad
  csvs:  <CSV_ROOT>/slide_NN_adata_k=200/clustering_results/<uuid>.csv
         Each CSV must contain a 'unit_ids' column that matches
         adata.obs['slide_FOV_cell'] (or adata.obs.index).

Edit the CONFIG block at the top then run:
    python 06_umap_all_slides.py
"""

# =============================================================================
# CONFIG - edit these if paths or column names differ
# =============================================================================

import os
from pathlib import Path

# Base directory holding the Annotation `input/` and `output/` folders.
# Defaults to this script's folder so the repository carries no absolute paths.
# Override with the REGISTERED_REPORT_DIR environment variable when the data lives
# elsewhere, e.g.  REGISTERED_REPORT_DIR=/path/to/data python 06_umap_all_slides.py
BASE_DIR = Path(os.environ.get("REGISTERED_REPORT_DIR", Path(__file__).resolve().parent))
H5AD_DIR = str(BASE_DIR / "input" / "h5ad")
CSV_ROOT = str(BASE_DIR / "output")
OUT_DIR  = str(BASE_DIR / "output" / "umap")

N_SLIDES       = 24
ONLY_SLIDE     = ""          # e.g. "slide_07"; leave blank to keep all slides
SUBSAMPLE_N    = 240000     # total cells to keep (stratified per slide); 0 = keep all
N_NEIGHBORS    = 30
UMAP_MIN_DIST  = 0.5
USE_HARMONY    = False        # batch-correct PCA by slide_id before neighbours/UMAP
MARKER_CMAP    = "viridis"   # colormap for marker-expression UMAPs
PNG_DPI        = 300
SVG_DPI        = 600          # raises embedded raster resolution inside SVG files

# Column names
OBS_CELL_ID_COL = "slide_FOV_cell"  # column in adata.obs matching CSV unit_ids
                                     # set "" to fall back to adata.obs.index
CSV_UNIT_COL    = "unit_ids"
CSV_CLUSTER_COL = "cluster_ids"     # adjust if your column is named differently
CSV_ANNOT_COL   = "annotation"      # adjust if your column is named differently

MARKERS_ALL = [
    "CD3", "CD15", "CD8", "CD20", "CD11c",
    "CD68", "FoxP3", "Pax5", "CD31", "Cytokeratin",
]

ANNOTATIONS_TO_KEEP = [
    "CD8+ T cells",
    "CD8- T cells",
    "Tregs",
    "B cells",
    "Epithelial",
    "Macrophages",
]

ANNOTATION_COLOR_DICT = {
    "CD8+ T cells": "#F781BF",
    "CD8- T cells": "#4DAF4A",
    "B cells": "#FEB24C",
    "Tregs": "#2196F0",
    "Epithelial": "#A65628",
    "Macrophages": "#E41A1C",
}
# =============================================================================

import glob
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"   # keep text editable in Illustrator
matplotlib.rcParams["image.interpolation"] = "nearest"
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def load_cluster_csv(slide_tag: str) -> "pd.DataFrame | None":
    """
    Read all CSVs under <CSV_ROOT>/slide_NN_adata_k=200/clustering_results/*.csv,
    concatenate them, and return a DataFrame indexed by unit_ids.
    """
    pattern = str(
        Path(CSV_ROOT) / f"{slide_tag}_adata_k=200" / "clustering_results" / "*.csv"
    )
    files = glob.glob(pattern)
    if not files:
        print(f"    [warn] no CSVs matched: {pattern}")
        return None

    parts = []
    for f in files:
        try:
            parts.append(pd.read_csv(f))
        except Exception as e:
            print(f"    [warn] could not read {f}: {e}")

    if not parts:
        return None

    df = pd.concat(parts, ignore_index=True)

    if CSV_UNIT_COL not in df.columns:
        print(f"    [warn] column '{CSV_UNIT_COL}' not found in CSVs for {slide_tag}")
        print(f"           available columns: {list(df.columns)}")
        return None

    df = df.drop_duplicates(subset=CSV_UNIT_COL).set_index(CSV_UNIT_COL)
    print(f"    CSV rows: {len(df):,}  (from {len(files)} file(s))")
    return df


def attach_metadata(adata: ad.AnnData, csv_df: "pd.DataFrame | None") -> ad.AnnData:
    """Map cluster_id / annotation from csv_df onto adata.obs via unit_ids."""
    if OBS_CELL_ID_COL and OBS_CELL_ID_COL in adata.obs.columns:
        keys = adata.obs[OBS_CELL_ID_COL]
    else:
        keys = pd.Series(adata.obs.index, index=adata.obs.index)

    for col, default in [(CSV_CLUSTER_COL, "unknown"), (CSV_ANNOT_COL, "unknown")]:
        if csv_df is not None and col in csv_df.columns:
            adata.obs[col] = keys.map(csv_df[col]).fillna(default).values
        else:
            adata.obs[col] = default

    if CSV_CLUSTER_COL in adata.obs.columns and CSV_CLUSTER_COL != "cluster_id":
        adata.obs.rename(columns={CSV_CLUSTER_COL: "cluster_id"}, inplace=True)
    if CSV_ANNOT_COL in adata.obs.columns and CSV_ANNOT_COL != "annotation":
        adata.obs.rename(columns={CSV_ANNOT_COL: "annotation"}, inplace=True)

    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)
    adata.obs["annotation"] = adata.obs["annotation"].astype(str)
    return adata


def stratified_subsample(
    adata: ad.AnnData, n_total: int, key: str = "slide_id"
) -> ad.AnnData:
    """Sample n_total cells proportionally from each group."""
    counts = adata.obs[key].value_counts()
    per_grp = max(1, n_total // len(counts))
    rng = np.random.default_rng(0)
    keep = []
    for grp, cnt in counts.items():
        idx = np.where(adata.obs[key] == grp)[0]
        keep.extend(rng.choice(idx, size=min(per_grp, cnt), replace=False).tolist())
    return adata[np.sort(keep)].copy()


def filter_definite_annotations(adata: ad.AnnData, annotations_to_keep: list[str]) -> ad.AnnData:
    """Keep only cells whose annotation is one of the definite cell types."""
    if not annotations_to_keep:
        return adata

    before_n = adata.n_obs
    adata.obs["annotation"] = adata.obs["annotation"].astype(str).str.strip()
    keep_mask = adata.obs["annotation"].isin(annotations_to_keep)
    filtered = adata[keep_mask].copy()

    kept_by_annotation = filtered.obs["annotation"].value_counts().reindex(
        annotations_to_keep, fill_value=0
    )
    print(
        "Filtered to definite annotations: "
        f"{filtered.n_obs:,} / {before_n:,} cells kept "
        f"({filtered.n_obs / before_n:.1%})"
    )
    print("Cells kept by annotation:")
    for annotation, count in kept_by_annotation.items():
        print(f"    {annotation}: {count:,}")

    if filtered.n_obs == 0:
        sys.exit("No cells left after annotation filtering - check annotation names.")

    return filtered


def make_palette(labels: pd.Series, cmap_name: str = "tab20") -> dict:
    unique = sorted(labels.dropna().unique())
    cmap = plt.get_cmap(cmap_name, max(len(unique), 1))
    return {lbl: cmap(i % 20) for i, lbl in enumerate(unique)}


def make_annotation_palette(labels: pd.Series) -> dict:
    """Use the fixed annotation colors, falling back to gray for unexpected labels."""
    unique = sorted(labels.dropna().unique())
    return {lbl: ANNOTATION_COLOR_DICT.get(lbl, "#8D8D8D") for lbl in unique}


def save_figure(fig: plt.Figure, out_path: Path, dpi=None):
    """Save a figure with high-resolution raster data and editable SVG text."""
    out_path = Path(out_path)
    if dpi is None:
        dpi = SVG_DPI if out_path.suffix.lower() == ".svg" else PNG_DPI
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"  Saved -> {out_path}")


def scatter_umap(adata, color_col, palette, title, out_path, pt_size=2.0):
    fig, ax = plt.subplots(figsize=(11, 8), dpi=SVG_DPI)
    coords = adata.obsm["X_umap"]
    labels = adata.obs[color_col]
    for lbl, color in palette.items():
        mask = labels == lbl
        if mask.any():
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                c=[color], s=pt_size, linewidths=0, alpha=0.7, rasterized=True
            )
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.set_aspect("equal", "box")
    handles = [
        mpatches.Patch(color=c, label=str(l))
        for l, c in list(palette.items())[:50]
    ]
    ax.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
        fontsize=6, frameon=False, ncol=max(1, len(handles) // 25)
    )
    fig.tight_layout()
    save_figure(fig, out_path)
    plt.close(fig)


def _marker_values(adata: ad.AnnData, marker: str, layer: str = "marker_expression"):
    """Return a 1D expression vector for one marker."""
    marker_idx = adata.var_names.get_loc(marker)
    matrix = adata.layers[layer] if layer in adata.layers else adata.X
    values = matrix[:, marker_idx]
    if hasattr(values, "toarray"):
        values = values.toarray()
    return np.asarray(values).ravel()


def scatter_marker_umap(
    adata: ad.AnnData,
    marker: str,
    out_path: Path,
    layer: str = "marker_expression",
    cmap: str = MARKER_CMAP,
    pt_size: float = 2.0,
):
    coords = adata.obsm["X_umap"]
    values = _marker_values(adata, marker, layer=layer)
    vmax = np.nanpercentile(values, 99)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = None

    fig, ax = plt.subplots(figsize=(9, 8), dpi=SVG_DPI)
    sc_plot = ax.scatter(
        coords[:, 0], coords[:, 1],
        c=values, s=pt_size, linewidths=0, alpha=0.8,
        cmap=cmap, vmin=0, vmax=vmax, rasterized=True
    )
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title(f"UMAP - {marker} expression", fontsize=13)
    ax.set_aspect("equal", "box")
    fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04, label=marker)
    fig.tight_layout()
    save_figure(fig, out_path)
    plt.close(fig)


def plot_marker_umap_grid(
    adata: ad.AnnData,
    markers: list[str],
    out_path: Path,
    layer: str = "marker_expression",
    cmap: str = MARKER_CMAP,
):
    n_cols = 5
    n_rows = int(np.ceil(len(markers) / n_cols))
    coords = adata.obsm["X_umap"]
    fig = plt.figure(
        figsize=(4.2 * n_cols + 0.8, 3.7 * n_rows), dpi=SVG_DPI
    )
    from matplotlib.gridspec import GridSpec

    gs = GridSpec(
        n_rows,
        n_cols + 1,
        figure=fig,
        width_ratios=[1] * n_cols + [0.08],
        wspace=0.45,
        hspace=0.45,
    )

    # Shared vmax across all markers
    all_values = np.concatenate([_marker_values(adata, m, layer=layer) for m in markers])
    shared_vmax = float(np.nanpercentile(all_values, 99))
    if not np.isfinite(shared_vmax) or shared_vmax <= 0:
        shared_vmax = None

    last_sc = None
    plot_axes = []
    for i, marker in enumerate(markers):
        row, col = divmod(i, n_cols)
        ax = fig.add_subplot(gs[row, col])
        plot_axes.append(ax)
        values = _marker_values(adata, marker, layer=layer)
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=values, s=1.0, linewidths=0, alpha=0.8,
            cmap=cmap, vmin=0, vmax=shared_vmax, rasterized=True
        )
        ax.set_title(marker, fontsize=11)
        ax.set_xlabel("UMAP 1", fontsize=9)
        ax.set_ylabel("UMAP 2", fontsize=9)
        ax.set_aspect("equal", "box")
        last_sc = sc_plot

    # Hide unused cells
    for i in range(len(markers), n_rows * n_cols):
        row, col = divmod(i, n_cols)
        fig.add_subplot(gs[row, col]).axis("off")

    # Single colorbar in its own column, adjacent to the Cytokeratin plot when present.
    if last_sc is not None:
        if "Cytokeratin" in markers:
            cbar_row = markers.index("Cytokeratin") // n_cols
        else:
            cbar_row = (len(markers) - 1) // n_cols
        cbar_ax = fig.add_subplot(gs[cbar_row, n_cols])
        cbar = fig.colorbar(last_sc, cax=cbar_ax)
        cbar.set_label("Expression", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.06, right=0.96, top=0.92, bottom=0.08)
    save_figure(fig, out_path)
    plt.close(fig)


def run_umap_outputs(
    adata_in: ad.AnnData,
    out_dir: Path,
    analysis_label: str,
    filter_annotations: bool = False,
):
    """Run subsampling, UMAP, and plot exports for one analysis mode."""
    out_dir.mkdir(parents=True, exist_ok=True)
    combined = adata_in.copy()

    print(f"\n=== Running UMAP: {analysis_label} ===")
    if filter_annotations:
        combined = filter_definite_annotations(combined, ANNOTATIONS_TO_KEEP)
    else:
        print(f"Keeping all cells: {combined.n_obs:,} cells")

    if SUBSAMPLE_N and SUBSAMPLE_N < combined.n_obs:
        print(f"Subsampling to {SUBSAMPLE_N:,} cells (stratified by slide) ...")
        combined = stratified_subsample(combined, SUBSAMPLE_N, key="slide_id")
        print(f"After subsample: {combined.n_obs:,} cells")

    combined.layers["marker_expression"] = combined.X.copy()

    print("Scaling ...")
    sc.pp.scale(combined, max_value=10)

    print(f"Neighbours (k={N_NEIGHBORS}) ...")
    sc.pp.neighbors(
        combined,
        n_neighbors=N_NEIGHBORS,
        use_rep="X",
        random_state=0,
    )

    print("UMAP ...")
    sc.tl.umap(combined, min_dist=UMAP_MIN_DIST, random_state=0)

    h5ad_out = out_dir / "combined_umap.h5ad"
    combined.write_h5ad(h5ad_out)
    print(f"Saved combined AnnData -> {h5ad_out}")

    print("Generating plots ...")
    cluster_pal = make_palette(combined.obs["cluster_id"], "tab20")
    if filter_annotations:
        annot_pal = make_annotation_palette(combined.obs["annotation"])
    else:
        annot_pal = make_palette(combined.obs["annotation"], "tab20")
    slide_pal = make_palette(combined.obs["slide_id"], "tab20")

    scatter_umap(
        combined, "cluster_id", cluster_pal,
        f"UMAP - cluster ID  (n={combined.n_obs:,})",
        out_dir / "umap_by_cluster.svg"
    )
    scatter_umap(
        combined, "annotation", annot_pal,
        f"UMAP - cell-type annotation  (n={combined.n_obs:,})",
        out_dir / "umap_by_annotation.svg"
    )
    scatter_umap(
        combined, "slide_id", slide_pal,
        f"UMAP - slide ID (batch check)  (n={combined.n_obs:,})",
        out_dir / "umap_by_slide.svg"
    )

    markers_present = [m for m in MARKERS_ALL if m in combined.var_names]
    if markers_present:
        print("Generating marker expression grid ...")
        plot_marker_umap_grid(
            combined,
            markers_present,
            out_dir / "umap_marker_expression_grid.svg",
        )

    print(f"Done with {analysis_label}. Outputs in: {out_dir}")


def main():
    h5ad_dir = Path(H5AD_DIR)
    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
     # --- DEBUG: print what files actually exist ---
    print(f"Looking for h5ad files in: {h5ad_dir}")
    found = sorted(h5ad_dir.glob("*.h5ad"))
    if found:
        print(f"Found {len(found)} h5ad file(s):")
        for f in found:
            print(f"  {f.name}")
    else:
        print("  [!] No .h5ad files found at all - check H5AD_DIR path")

    only_slide = ONLY_SLIDE.strip()
    expected = h5ad_dir / f"{only_slide}_adata.h5ad"
    print(f"Expected file: {expected}")
    print(f"Exists: {expected.exists()}")
    # --- END DEBUG ---

    only_slide = ONLY_SLIDE.strip()
    if only_slide:
        print(f"Filtering to one slide: {only_slide}")

    # 1. Load each slide
    adatas = []
    for idx in range(1, N_SLIDES + 1):
        slide_tag = f"slide_{idx:02d}"
        if only_slide and slide_tag != only_slide:
            continue

        h5ad_path = h5ad_dir / f"{slide_tag}_adata.h5ad"

        if not h5ad_path.exists():
            print(f"[skip] not found: {h5ad_path}")
            continue

        print(f"\nLoading {slide_tag} ...")
        adata = ad.read_h5ad(h5ad_path)

        present = [m for m in MARKERS_ALL if m in adata.var_names]
        missing = [m for m in MARKERS_ALL if m not in adata.var_names]
        if missing:
            print(f"    markers absent (skipped): {missing}")
        if not present:
            print(f"    [skip] no markers in {slide_tag}")
            continue
        adata = adata[:, present].copy()
        adata.obs["slide_id"] = slide_tag

        csv_df = load_cluster_csv(slide_tag)
        adata = attach_metadata(adata, csv_df)

        matched = (adata.obs["cluster_id"] != "unknown").sum()
        print(f"    Cells: {adata.n_obs:,}  |  matched to CSV: {matched:,}")
        adatas.append(adata)

    if not adatas:
        sys.exit("No slides loaded - check H5AD_DIR, CSV_ROOT, and ONLY_SLIDE.")

    # 2. Concatenate
    print(f"\nConcatenating {len(adatas)} slides ...")
    combined = ad.concat(
        adatas,
        join="inner",
        merge="same",
        label="slide_id",
        keys=[a.obs["slide_id"].iloc[0] for a in adatas],
    )
    combined.obs_names_make_unique()
    print(f"Combined: {combined.n_obs:,} cells x {combined.n_vars} markers")

    run_umap_outputs(
        combined,
        out_dir / "all_cells",
        analysis_label="all cells",
        filter_annotations=False,
    )
    run_umap_outputs(
        combined,
        out_dir / "filtered_annotations",
        analysis_label="filtered definite annotations",
        filter_annotations=True,
    )

    print(f"\nDone. All outputs in: {out_dir}")

if __name__ == "__main__":
    main()
