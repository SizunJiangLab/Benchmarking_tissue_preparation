import json
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import importlib

from clustering import subcluster
importlib.reload(subcluster)
from clustering.subcluster import ClusteringResult, ClusteringResultManager

# -------------------- Paths & IO --------------------
annotation_json_path = Path("/registered_report/cluster_annotations.json")
output_root = Path("/registered_report/output")
input_dir   = Path("registered_report/input/h5ad")
output_root.mkdir(parents=True, exist_ok=True)

# -------------------- Slide order (conditions -> slide IDs) --------------------
custom_order = [
    "BIDMC_13","BIDMC_21","BIDMC_5","BIDMC_15","BIDMC_23","BIDMC_17",
    "BIDMC_24","BIDMC_8","BIDMC_7","BIDMC_9","BIDMC_22","BIDMC_16",
    "BIDMC_14","BIDMC_19","BIDMC_1","BIDMC_2","BIDMC_6","BIDMC_20",
    "BIDMC_3","BIDMC_18","BIDMC_4","BIDMC_12","BIDMC_11","BIDMC_10"
]
slide_to_condition = {f"slide_{int(b.split('_')[1]):02d}": b for b in custom_order}
condition_to_slide = {v: k for k, v in slide_to_condition.items()}

# -------------------- Plot colors --------------------
POS_COLOR = "#FF8C00"  # orange (≥ 0)
NEG_COLOR = "#D9D9D9"  # light gray (< 0 or NaN)

# -------------------- Target checks (per Jia Le) --------------------
# (title, target cell type label, candidate marker names, palette key)
PAIRS = [
    ("CD3 in T cells",          "T cells",      ["CD3"],              "T cells"),
    ("CD20 in B cells",         "B cells",      ["CD20"],             "B cells"),
    ("PAX5 in B cells",         "B cells",      ["PAX5", "Pax5"],     "B cells"),
    ("CD20 in T cells",          "T cells",      ["CD20"],              "T cells"),
    ("CD3 in B cells",         "B cells",      ["CD3"],             "B cells"),
]

# ---- T-family handling ----
INCLUDE_TREGS_IN_T_FAMILY = True
T_CELL_FAMILY = (
    "T cells", "CD8+ T cells", "CD4+ T cells"
) + (("Tregs",) if INCLUDE_TREGS_IN_T_FAMILY else ())

def members_for(target_ct: str) -> tuple[str, ...]:
    """Return the set of annotation labels that count as 'in-class' for a target."""
    if target_ct == "T cells":
        # For CD3 we will use the full T-family (T cells + CD8+ + CD8- [+/- Tregs]).
        return T_CELL_FAMILY
    return (target_ct,)

def _find_var(adata: ad.AnnData, candidates) -> str | None:
    """Return the first adata.var that matches any candidate (case-insensitive)."""
    lower_map = {v.lower(): v for v in adata.var_names}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None

def _dense1(col) -> np.ndarray:
    """Return dense 1D np.array from AnnData column (handles sparse)."""
    if hasattr(col, "A1"):  # scipy sparse
        return col.A1
    return np.asarray(col).reshape(-1)

# -------------------- Load & normalize annotations --------------------
with open(annotation_json_path) as f:
    annotations_all = json.load(f)

# normalize annotation strings, drop "Artifacts"
for slide, cluster_dict in annotations_all.items():
    for k, v in list(cluster_dict.items()):
        if v == "Artifacts":
            cluster_dict[k] = ""
        elif isinstance(v, str):
            cluster_dict[k] = v.strip()

# -------------------- Compute per-slide metrics --------------------
# Metric:
#   log2_ratio = log2( mean_in / mean_out )
# where "out" excludes:
#   - "Other" (case-insensitive)
#   - labels containing " and " (case-insensitive), i.e. mixed clusters like
#       "Neutrophils and Other", "Epithelial and B cells", "T cells and X".
#
# For CD8 in T cells:
#   mean_in  = MFI(CD8) in CD8+ T cells
#   mean_out = MFI(CD8) in Others
#             (excluding Other + mixed + the broad "T cells" cluster)

records = []

for slide_key in annotations_all.keys():
    h5ad_path   = input_dir / f"{slide_key}_adata.h5ad"
    out_dir     = output_root / f"{slide_key}_adata_k=200"
    geojson_dir = out_dir / "geojson"

    if not h5ad_path.exists() or not geojson_dir.is_dir():
        continue
    cdirs = [d for d in geojson_dir.iterdir() if d.is_dir()]
    if not cdirs:
        continue
    clustering_id = max(cdirs, key=lambda d: d.stat().st_mtime).name

    # load data + annotations aligned to obs
    _ = ClusteringResult.pop(clustering_id=clustering_id, output_dir=out_dir)
    adata = ad.read_h5ad(h5ad_path)
    manager = ClusteringResultManager(output_dir=out_dir, unit_ids=adata.obs.index)
    ann = manager.summary_df["annotation"].reindex(adata.obs.index).fillna("")

    labels = ann.values.astype(str)
    valid_mask = (labels != "")

    # define "Other" and "mixed" (labels containing " and " case-insensitive)
    labels_lower = np.char.lower(labels)
    unknown_mask = (labels == "Other") | (labels_lower == "other")

    # Base exclusion for the denominator
    exclude_from_out_base = unknown_mask

    # Optional extra: any label that looks like "mixed t cells"
    mixed_t_mask = np.char.find(labels_lower, "mixed t") >= 0

    for (title, target_ct, marker_names, _color_key) in PAIRS:
        var = _find_var(adata, marker_names)
        if var is None:
            records.append({
                "slide": slide_key, "pair": title,
                "log2_ratio": np.nan,
                "mean_in": np.nan, "mean_out": np.nan,
                "n_in": 0, "n_out": 0
            })
            continue

        x_all = _dense1(adata[:, var].X)

        # ---- define in-class labels & denominator exclusion ----
        if title == "CD8 in T cells":
            # SPECIAL CASE FOR CD8:
            # Numerator = CD8+ T cells only
            in_labels = ("CD8+ T cells",)

            # Denominator excludes:
            #  - Other
            #  - mixed clusters containing " and "
            #  - any "mixed t" clusters
            #  - the broad "T cells" cluster
            exclude_from_out = (
                exclude_from_out_base |
                mixed_t_mask |
                (labels == "T cells")
            )
        else:
            # Default for all other markers
            in_labels = members_for(target_ct)
            exclude_from_out = exclude_from_out_base

        mask_in  = np.isin(labels, in_labels) & valid_mask
        mask_out = (~np.isin(labels, in_labels)) & valid_mask

        n_in  = int(mask_in.sum())
        n_out = int(mask_out.sum())

        if n_in == 0 or n_out == 0:
            records.append({
                "slide": slide_key, "pair": title,
                "log2_ratio": np.nan,
                "mean_in": np.nan, "mean_out": np.nan,
                "n_in": n_in, "n_out": n_out
            })
            continue

        xin  = x_all[mask_in]
        xout = x_all[mask_out]

        mean_in  = float(np.mean(xin))
        mean_out = float(np.mean(xout))

        if mean_in > 0 and mean_out > 0:
            log2_ratio = float(np.log2(mean_in / mean_out))
        else:
            log2_ratio = np.nan

        records.append({
            "slide": slide_key, "pair": title,
            "log2_ratio": log2_ratio,
            "mean_in": mean_in, "mean_out": mean_out,
            "n_in": n_in, "n_out": n_out
        })

# -------------------- Tabulate results --------------------
df = pd.DataFrame.from_records(records)

# Save long and wide tables
df.to_csv(output_root / "marker_specificity_log2_metrics_long.csv", index=False)

wide_ratio = df.pivot_table(index="slide", columns="pair", values="log2_ratio", aggfunc="mean")

# Order slides (keep only those present)
ordered_slide_ids = [
    condition_to_slide[c] for c in custom_order
    if c in condition_to_slide and condition_to_slide[c] in wide_ratio.index
]
wide_ratio = wide_ratio.loc[ordered_slide_ids]

wide_ratio.to_csv(output_root / "marker_specificity_log2_ratio_wide.csv")

# -------------------- Plotting --------------------
def _plot_matrix_single_row(wide: pd.DataFrame, fig_title: str, filename_png: str, filename_svg: str):
    """
    One row with len(PAIRS) columns; bars per slide colored by sign (≥0 orange, <0 gray).
    """
    cols = [t for (t, *_rest) in PAIRS if t in wide.columns]
    if not cols:
        return

    ncol = len(cols)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2*ncol, 6.2), sharey=True)
    if ncol == 1:
        axes = np.array([axes])

    def _ok(v):
        return isinstance(v, (int, float, np.floating)) and np.isfinite(v)

    def _plot_col(ax, title):
        vals = wide[title]
        colors = [POS_COLOR if (_ok(v) and v >= 0) else NEG_COLOR for v in vals.values]
        ax.bar(range(len(vals)), vals.values, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_title(title, fontsize=11)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(
            [slide_to_condition.get(s, s) for s in vals.index],
            rotation=90, fontsize=8
        )
        ax.axhline(0, color="black", linewidth=0.8)

    for j, title in enumerate(cols):
        _plot_col(axes[j], title)

    axes[0].set_ylabel("log2(MFI_in / MFI_out)")

    handles = [
        Patch(facecolor=POS_COLOR, edgecolor="black", label="≥ 0 (enriched)"),
        Patch(facecolor=NEG_COLOR, edgecolor="black", label="< 0 (depleted/undef)"),
    ]
    axes[0].legend(handles=handles, loc="upper left", frameon=False, fontsize=9)

    fig.suptitle(fig_title, y=1.02, fontsize=13)
    plt.tight_layout()
    # Save SVG + PNG
    plt.savefig(output_root / filename_svg, dpi=300, bbox_inches="tight")
    plt.savefig(output_root / filename_png, dpi=300, bbox_inches="tight")
    plt.close()

_plot_matrix_single_row(
    wide_ratio,
    fig_title="log2(MFI_in / MFI_out) per slide (others exclude Other, mixed & broad T cells for CD8)",
    filename_png="marker_specificity_log2_ratio_per_slide.png",
    filename_svg="marker_specificity_log2_ratio_per_slide.svg",
)

print("Saved:")
print(output_root / "marker_specificity_log2_metrics_long.csv")
print(output_root / "marker_specificity_log2_ratio_wide.csv")
print(output_root / "marker_specificity_log2_ratio_per_slide.png")
print(output_root / "marker_specificity_log2_ratio_per_slide.svg")
