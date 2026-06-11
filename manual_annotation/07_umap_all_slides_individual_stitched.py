"""
07_umap_all_slides_individual_stitched.py

"""

# =============================================================================
# CONFIG - edit these if paths or column names differ
# =============================================================================

import os
from pathlib import Path

# Base directory holding the Annotation `input/` and `output/` folders.
# Defaults to this script's folder so the repository carries no absolute paths.
# Override with the REGISTERED_REPORT_DIR environment variable when the data lives
# elsewhere, e.g.  REGISTERED_REPORT_DIR=/path/to/data python 07_umap_all_slides_individual_stitched.py
BASE_DIR  = Path(os.environ.get("REGISTERED_REPORT_DIR", Path(__file__).resolve().parent))
H5AD_DIR  = str(BASE_DIR / "input" / "h5ad")
CSV_ROOT  = str(BASE_DIR / "output")
OUT_DIR   = str(BASE_DIR / "output" / "umap")

N_SLIDES       = 24
ONLY_SLIDE     = ""          # e.g. "slide_07"; leave blank to keep all slides
SUBSAMPLE_N    = 0     # total cells to keep (stratified per slide); 0 = keep all
N_NEIGHBORS    = 30
UMAP_MIN_DIST  = 0.5
USE_HARMONY    = False        # batch-correct PCA by slide_id before neighbours/UMAP
MARKER_CMAP    = "viridis"   # colormap for marker-expression UMAPs
MIN_SLIDE_FONT_SIZE = 36      # source font size; scaled down when stitched onto A4

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
# =============================================================================

import glob
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
matplotlib.rcParams["svg.fonttype"] = "none"
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


def make_palette(labels: pd.Series, cmap_name: str = "tab20") -> dict:
    unique = sorted(labels.dropna().unique())
    cmap = plt.get_cmap(cmap_name, max(len(unique), 1))
    return {lbl: cmap(i % 20) for i, lbl in enumerate(unique)}


def enforce_min_font_size(fig: plt.Figure, min_size: float = MIN_SLIDE_FONT_SIZE):
    """Raise all Matplotlib text in a figure to at least min_size."""
    for text in fig.findobj(match=plt.Text):
        if text.get_fontsize() < min_size:
            text.set_fontsize(min_size)
    for ax in fig.axes:
        ax.tick_params(axis="both", labelsize=min_size)


def scatter_umap(adata, color_col, palette, title, out_path, pt_size=2.0):
    fig, ax = plt.subplots(figsize=(11, 8))
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
    enforce_min_font_size(fig)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_path}")


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

    fig, ax = plt.subplots(figsize=(9, 8))
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
    enforce_min_font_size(fig)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_path}")


def plot_marker_umap_grid(
    adata: ad.AnnData,
    markers: list[str],
    out_path: Path,
    layer: str = "marker_expression",
    cmap: str = MARKER_CMAP,
):
    n_cols = 4
    n_rows = int(np.ceil(len(markers) / n_cols))
    coords = adata.obsm["X_umap"]
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.7 * n_rows))
    axes = np.asarray(axes).ravel()

    for ax, marker in zip(axes, markers):
        values = _marker_values(adata, marker, layer=layer)
        vmax = np.nanpercentile(values, 99)
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = None
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=values, s=1.0, linewidths=0, alpha=0.8,
            cmap=cmap, vmin=0, vmax=vmax, rasterized=True
        )
        ax.set_title(marker, fontsize=11)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal", "box")
        fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes[len(markers):]:
        ax.axis("off")

    enforce_min_font_size(fig)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_path}")


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

    # 3. Filter to selected cell-type annotations
    #if ANNOTATIONS_TO_KEEP:
       # before_n = combined.n_obs
       # keep_mask = combined.obs["annotation"].isin(ANNOTATIONS_TO_KEEP)
       # combined = combined[keep_mask].copy()
       # print(
       #     "Filtered annotations: "
       #     f"{combined.n_obs:,} / {before_n:,} cells kept "
       #     f"({', '.join(ANNOTATIONS_TO_KEEP)})"
       # )
       # if combined.n_obs == 0:
       #     sys.exit("No cells left after annotation filtering - check annotation names.")

    # 4. Stratified subsample
    if SUBSAMPLE_N and SUBSAMPLE_N < combined.n_obs:
        print(f"Subsampling to {SUBSAMPLE_N:,} cells (stratified by slide) ...")
        combined = stratified_subsample(combined, SUBSAMPLE_N, key="slide_id")
        print(f"After subsample: {combined.n_obs:,} cells")

    # 5. Scale -> PCA -> Harmony -> neighbours -> UMAP
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

    # 5. Save AnnData
    h5ad_out = out_dir / "combined_umap.h5ad"
    combined.write_h5ad(h5ad_out)
    print(f"\nSaved combined AnnData -> {h5ad_out}")

    # 6. UMAP plots
    print("\nGenerating plots ...")
    cluster_pal = make_palette(combined.obs["cluster_id"], "tab20")
    annot_pal   = make_palette(combined.obs["annotation"], "tab20")
    slide_pal   = make_palette(combined.obs["slide_id"],   "tab20")

    scatter_umap(
        combined, "cluster_id", cluster_pal,
        f"UMAP - cluster ID  (n={combined.n_obs:,})",
        out_dir / "umap_by_cluster.png"
    )
    scatter_umap(
        combined, "slide_id", slide_pal,
        f"UMAP - slide ID (batch check)  (n={combined.n_obs:,})",
        out_dir / "umap_by_slide.png"
    )

    # 7. Stitched figure: annotation UMAP (left, tall) + marker grid (right)
    markers_present = [m for m in MARKERS_ALL if m in combined.var_names]
    if markers_present:
        print("\nGenerating stitched annotation + marker-expression figure ...")

        n_cols      = 4
        n_rows      = int(np.ceil(len(markers_present) / n_cols))
        cell_w      = 4.2          # width of each marker subplot (inches)
        cell_h      = 3.7          # height of each marker subplot (inches)
        grid_w      = cell_w * n_cols
        grid_h      = cell_h * n_rows
        annot_w     = grid_w * 0.6   # annotation panel width
        total_w     = annot_w + grid_w
        total_h     = grid_h * 2      # annotation panel is twice as tall

        fig = plt.figure(figsize=(total_w, total_h))

        # Left: annotation UMAP spanning full figure height
        ax_annot = fig.add_axes(
            [0, 0, annot_w / total_w, 1.0]   # [left, bottom, width, height] in figure fraction
        )
        coords = combined.obsm["X_umap"]
        labels = combined.obs["annotation"]
        for lbl, color in annot_pal.items():
            mask = labels == lbl
            if mask.any():
                ax_annot.scatter(
                    coords[mask, 0], coords[mask, 1],
                    c=[color], s=2.0, linewidths=0, alpha=0.7, rasterized=True
                )
        ax_annot.set_xlabel("UMAP 1", fontsize=11)
        ax_annot.set_ylabel("UMAP 2", fontsize=11)
        ax_annot.set_title(f"Cell-type annotation  (n={combined.n_obs:,})", fontsize=13)
        ax_annot.set_aspect("equal", "box")
        handles = [
            mpatches.Patch(color=c, label=str(l))
            for l, c in annot_pal.items()
        ]
        ax_annot.legend(
            handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
            fontsize=8, frameon=False,
        )

        # Right: marker expression grid
        # Build a GridSpec in the right portion of the figure
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(
            n_rows, n_cols,
            figure=fig,
            left=annot_w / total_w + 0.09,
            right=1.0,
            top=1.0 - 0.02,
            bottom=0.0 + 0.02,
            wspace=0.35,
            hspace=0.25,
        )

        # Compute a single shared vmax across all markers (99th percentile)
        all_values = np.concatenate([
            _marker_values(combined, m) for m in markers_present
        ])
        shared_vmax = float(np.nanpercentile(all_values, 99))
        if not np.isfinite(shared_vmax) or shared_vmax <= 0:
            shared_vmax = None

        last_sc = None
        for i, marker in enumerate(markers_present):
            row, col = divmod(i, n_cols)
            ax = fig.add_subplot(gs[row, col])
            values = _marker_values(combined, marker)
            sc_plot = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=values, s=1.0, linewidths=0, alpha=0.8,
                cmap=MARKER_CMAP, vmin=0, vmax=shared_vmax, rasterized=True
            )
            ax.set_title(marker, fontsize=10)
            ax.set_xlabel("UMAP 1", fontsize=7)
            ax.set_ylabel("UMAP 2", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_aspect("equal", "box")
            last_sc = sc_plot

        # Hide unused grid cells
        for i in range(len(markers_present), n_rows * n_cols):
            row, col = divmod(i, n_cols)
            fig.add_subplot(gs[row, col]).axis("off")

        # Single shared colorbar in the bottom-right empty cell or below last row
        if last_sc is not None:
            n_empty = n_rows * n_cols - len(markers_present)
            if n_empty > 0:
                br_row = n_rows - 1
                br_col = n_cols - 1
                cbar_host = fig.add_subplot(gs[br_row, br_col])
                cbar_host.axis("off")
                pos = cbar_host.get_position()
                cbar_ax = fig.add_axes([
                    pos.x0 + pos.width * 0.2,
                    pos.y0 + pos.height * 0.15,
                    pos.width * 0.25,
                    pos.height * 0.6,
                ])
            else:
                pos_last = fig.axes[-1].get_position()
                cbar_ax = fig.add_axes([
                    pos_last.x1 + 0.01,
                    pos_last.y0,
                    0.015,
                    pos_last.height,
                ])
            fig.colorbar(last_sc, cax=cbar_ax, label="Expression")

        stitched_out = out_dir / "umap_annotation_and_markers.png"
        fig.savefig(stitched_out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved -> {stitched_out}")

    print(f"\nDone. All outputs in: {out_dir}")

def run_single_slide(slide_tag: str, out_dir_base: Path) -> "Path | None":
    """Run the full pipeline for one slide and return the stitched image path."""
    import io
    from contextlib import redirect_stdout

    h5ad_dir = Path(H5AD_DIR)
    out_dir = out_dir_base / slide_tag
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = h5ad_dir / f"{slide_tag}_adata.h5ad"
    if not h5ad_path.exists():
        print(f"[skip] not found: {h5ad_path}")
        return None

    print(f"\n{'='*50}")
    print(f"Processing {slide_tag} ...")
    print(f"{'='*50}")

    adata = ad.read_h5ad(h5ad_path)

    present = [m for m in MARKERS_ALL if m in adata.var_names]
    missing = [m for m in MARKERS_ALL if m not in adata.var_names]
    if missing:
        print(f"  markers absent: {missing}")
    if not present:
        print(f"  [skip] no markers in {slide_tag}")
        return None

    adata = adata[:, present].copy()
    adata.obs["slide_id"] = slide_tag

    csv_df = load_cluster_csv(slide_tag)
    adata = attach_metadata(adata, csv_df)

    matched = (adata.obs["cluster_id"] != "unknown").sum()
    print(f"  Cells: {adata.n_obs:,}  |  matched to CSV: {matched:,}")

    # Subsample
    if SUBSAMPLE_N and SUBSAMPLE_N < adata.n_obs:
        print(f"  Subsampling to {SUBSAMPLE_N:,} cells ...")
        adata = stratified_subsample(adata, SUBSAMPLE_N, key="slide_id")

    # Scale -> neighbours -> UMAP
    adata.layers["marker_expression"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=0)
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=0)

    # Save h5ad
    adata.write_h5ad(out_dir / f"{slide_tag}_umap.h5ad")

    # Build stitched figure
    markers_present = [m for m in MARKERS_ALL if m in adata.var_names]
    annot_pal = make_palette(adata.obs["annotation"], "tab20")

    n_cols  = 4
    n_rows  = int(np.ceil(len(markers_present) / n_cols))
    cell_w  = 8.0
    cell_h  = 7.2
    grid_w  = cell_w * n_cols
    grid_h  = cell_h * n_rows
    annot_w = grid_w * 0.65
    gap_w   = 6.0
    total_w = annot_w + gap_w + grid_w
    total_h = grid_h * 2.2

    fig = plt.figure(figsize=(total_w, total_h))

    # Annotation UMAP (left)
    ax_annot = fig.add_axes([0.02, 0.08, annot_w / total_w - 0.03, 0.84])
    coords = adata.obsm["X_umap"]
    labels = adata.obs["annotation"]
    for lbl, color in annot_pal.items():
        mask = labels == lbl
        if mask.any():
            ax_annot.scatter(
                coords[mask, 0], coords[mask, 1],
                c=[color], s=2.0, linewidths=0, alpha=0.7, rasterized=True
            )
    ax_annot.set_xlabel("UMAP 1", fontsize=11)
    ax_annot.set_ylabel("UMAP 2", fontsize=11)
    ax_annot.set_title(f"{slide_tag} — Cell-type annotation  (n={adata.n_obs:,})", fontsize=13)
    ax_annot.tick_params(labelsize=MIN_SLIDE_FONT_SIZE)
    ax_annot.set_aspect("equal", "box")
    handles = [mpatches.Patch(color=c, label=str(l)) for l, c in annot_pal.items()]
    ax_annot.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.98),
        fontsize=MIN_SLIDE_FONT_SIZE,
        frameon=False,
        borderaxespad=0,
        labelspacing=0.6,
        handlelength=1.0,
    )

    # Marker grid (right)
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(
        n_rows, n_cols, figure=fig,
        left=(annot_w + gap_w) / total_w,
        right=0.98,
        top=0.90,
        bottom=0.10,
        wspace=0.85,
        hspace=0.95,
    )

    all_values = np.concatenate([_marker_values(adata, m) for m in markers_present])
    shared_vmax = float(np.nanpercentile(all_values, 99))
    if not np.isfinite(shared_vmax) or shared_vmax <= 0:
        shared_vmax = None

    last_sc = None
    for i, marker in enumerate(markers_present):
        row, col = divmod(i, n_cols)
        ax = fig.add_subplot(gs[row, col])
        values = _marker_values(adata, marker)
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=values, s=1.0, linewidths=0, alpha=0.8,
            cmap=MARKER_CMAP, vmin=0, vmax=shared_vmax, rasterized=True
        )
        ax.set_title(marker, fontsize=MIN_SLIDE_FONT_SIZE)
        ax.set_xlabel("UMAP 1", fontsize=MIN_SLIDE_FONT_SIZE)
        ax.set_ylabel("UMAP 2", fontsize=MIN_SLIDE_FONT_SIZE)
        ax.tick_params(labelsize=MIN_SLIDE_FONT_SIZE)
        ax.set_aspect("equal", "box")
        last_sc = sc_plot

    for i in range(len(markers_present), n_rows * n_cols):
        row, col = divmod(i, n_cols)
        fig.add_subplot(gs[row, col]).axis("off")

    if last_sc is not None:
        n_empty = n_rows * n_cols - len(markers_present)
        if n_empty > 0:
            cbar_host = fig.add_subplot(gs[n_rows - 1, n_cols - 1])
            cbar_host.axis("off")
            pos = cbar_host.get_position()
            cbar_ax = fig.add_axes([
                pos.x0 + pos.width * 0.2, pos.y0 + pos.height * 0.15,
                pos.width * 0.25, pos.height * 0.6,
            ])
        else:
            pos_last = fig.axes[-1].get_position()
            cbar_ax = fig.add_axes([
                pos_last.x1 + 0.01, pos_last.y0, 0.015, pos_last.height,
            ])
        cbar = fig.colorbar(last_sc, cax=cbar_ax, label="Expression")
        cbar.ax.tick_params(labelsize=MIN_SLIDE_FONT_SIZE)
        cbar.set_label("Expression", fontsize=MIN_SLIDE_FONT_SIZE)

    out_path = out_dir / f"{slide_tag}_umap_annotation_and_markers.svg"
    enforce_min_font_size(fig)
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {out_path}")
    return out_path


def stitch_slide_groups_svg(image_paths: list, out_dir: Path):
    """Create A4 portrait SVG pages with four slide UMAP SVGs stacked vertically."""
    import re
    import svgutils.transform as sg

    def size_to_pt(size):
        # Handles strings like "1200pt", "10in", "800px", "210mm"
        match = re.match(r"([0-9.]+)([a-zA-Z]*)", str(size))
        if not match:
            raise ValueError(f"Could not parse SVG size: {size}")
        value = float(match.group(1))
        unit = match.group(2)

        if unit == "pt":
            return value
        if unit == "in":
            return value * 72
        if unit == "mm":
            return value * 72 / 25.4
        if unit == "cm":
            return value * 72 / 2.54
        if unit == "px" or unit == "":
            return value * 72 / 96
        return value

    def set_text_min_font_size(element, min_size=5):
        """Raise explicit SVG font-size attributes below min_size to min_size."""
        for node in element.root.iter():
            style = node.attrib.get("style", "")
            if "font-size:" in style:
                style = re.sub(
                    r"font-size:\s*([0-9.]+)px",
                    lambda m: f"font-size: {max(float(m.group(1)) * 72 / 96, min_size):g}pt",
                    style,
                )
                style = re.sub(
                    r"font-size:\s*([0-9.]+)pt",
                    lambda m: f"font-size: {max(float(m.group(1)), min_size):g}pt",
                    style,
                )
                node.attrib["style"] = style

            font_size = node.attrib.get("font-size")
            if font_size:
                match = re.match(r"([0-9.]+)([a-zA-Z]*)", font_size)
                if match:
                    value = float(match.group(1))
                    unit = match.group(2) or "pt"
                    value_pt = value * 72 / 96 if unit == "px" else value
                    if value_pt < min_size:
                        node.attrib["font-size"] = f"{min_size}pt"

    groups = [
        ("slides_01_04_A4.svg", image_paths[0:4]),
        ("slides_05_08_A4.svg", image_paths[4:8]),
        ("slides_09_12_A4.svg", image_paths[8:12]),
        ("slides_13_16_A4.svg", image_paths[12:16]),
        ("slides_17_20_A4.svg", image_paths[16:20]),
        ("slides_21_24_A4.svg", image_paths[20:24]),
    ]

    # A4 portrait in points: 210 x 297 mm.
    page_w = 210 * 72 / 25.4
    page_h = 297 * 72 / 25.4
    margin = 18
    header_h = 18
    gap = 8
    n_grid_rows = 4
    min_font_size = 11

    panel_w = page_w - 2 * margin
    panel_h = (page_h - 2 * margin - (n_grid_rows * header_h) - ((n_grid_rows - 1) * gap)) / n_grid_rows

    for filename, group in groups:
        loaded = []
        for tag, path in group:
            if path is not None and Path(path).exists():
                svg = sg.fromfile(str(path))
                w, h = svg.get_size()
                set_text_min_font_size(svg, min_size=min_font_size)
                loaded.append((tag, path, svg, size_to_pt(w), size_to_pt(h)))
            else:
                loaded.append((tag, path, None, panel_w, panel_h))

        fig = sg.SVGFigure(f"{page_w}pt", f"{page_h}pt")
        fig.root.set("viewBox", f"0 0 {page_w} {page_h}")
        elements = []

        for idx, (tag, path, svg, w, h) in enumerate(loaded):
            x = margin
            y = margin + idx * (panel_h + header_h + gap)

            # Editable slide label
            label = sg.TextElement(
                x,
                y + min_font_size + 1,
                tag,
                size=min_font_size,
                font="Liberation Sans",
                weight="bold",
            )
            elements.append(label)

            if svg is not None:
                root = svg.getroot()
                scale = min(panel_w / w, panel_h / h)
                scaled_w = w * scale
                plot_x = x + (panel_w - scaled_w) / 2
                root.moveto(plot_x, y + header_h, scale_x=scale, scale_y=scale)
                elements.append(root)
            else:
                missing = sg.TextElement(
                    x + panel_w / 2,
                    y + header_h + panel_h / 2,
                    "missing",
                    size=min_font_size,
                    font="Liberation Sans",
                )
                elements.append(missing)

        fig.append(elements)

        out_path = out_dir / filename
        fig.save(str(out_path))
        print(f"Saved editable SVG grid -> {out_path}")

if __name__ == "__main__":
    base_out = Path(OUT_DIR)

    slide_image_paths = []
    for idx in range(1, N_SLIDES + 1):
        slide_tag = f"slide_{idx:02d}"
        out_path = run_single_slide(slide_tag, base_out)
        slide_image_paths.append((slide_tag, out_path))

    stitch_slide_groups_svg(slide_image_paths, base_out)
