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

annotation_json_path = Path(
    "/registered_report/cluster_annotations.json"
)
output_root = Path("/registered_report/output")
input_dir   = Path("/registered_report/input/h5ad")
output_root.mkdir(parents=True, exist_ok=True)

# Slide order (conditions -> slide IDs)
custom_order = [
    "BIDMC_13","BIDMC_21","BIDMC_5","BIDMC_15","BIDMC_23","BIDMC_17",
    "BIDMC_24","BIDMC_8","BIDMC_7","BIDMC_9","BIDMC_22","BIDMC_16",
    "BIDMC_14","BIDMC_19","BIDMC_1","BIDMC_2","BIDMC_6","BIDMC_20",
    "BIDMC_3","BIDMC_18","BIDMC_4","BIDMC_12","BIDMC_11","BIDMC_10"
]
slide_to_condition = {f"slide_{int(b.split('_')[1]):02d}": b for b in custom_order}
condition_to_slide = {v: k for k, v in slide_to_condition.items()}

POS_COLOR = "#FF8C00"  # orange
NEG_COLOR = "#D9D9D9"  # light gray

# Target checks (original + cross-specificity)
# (title, target cell type label, candidate marker names, palette key)
ORIGINAL_PAIRS = [
    ("CD20 in B cells",                 "B cells",                 ["CD20"],                "B cells"),
    ("PAX5 in B cells",                 "B cells",                 ["PAX5","Pax5"],         "B cells"),
    ("CD3 in T cells",                  "T cells",                 ["CD3"],                 "T cells"),
    ("CD8 in CD8+ T cells",             "CD8+ T cells",            ["CD8"],                 "CD8+ T cells"),
    ("CD15 in neutrophils",             "Neutrophils and Unknown", ["CD15"],                "Neutrophils and Unknown"),
]
CROSS_SPEC_PAIRS = [
    ("CD20 in T cells",                 "T cells",                 ["CD20"],                "T cells"),
    ("CD3 in B cells",                  "B cells",                 ["CD3"],                 "B cells"),
    ("PAX5 in T cells",                 "T cells",                 ["PAX5","Pax5"],         "T cells"),
    ("CD8 in B cells",                  "B cells",                 ["CD8"],                 "B cells"),
    ("FoxP3 in Neutrophils",            "Neutrophils and Unknown", ["FOXP3","FoxP3"],       "Neutrophils and Unknown"),
    ("CD15 in CD8+ T cells",            "CD8+ T cells",            ["CD15"],                "CD8+ T cells"),
]
PAIRS = ORIGINAL_PAIRS + CROSS_SPEC_PAIRS  # full ordered list

# ---- T-family handling ----
INCLUDE_TREGS_IN_T_FAMILY = True
T_CELL_FAMILY = (
    "T cells", "CD8+ T cells", "CD8- T cells"
) + (("Tregs",) if INCLUDE_TREGS_IN_T_FAMILY else ())

def members_for(target_ct: str) -> tuple[str, ...]:
    """Return the set of annotation labels that count as 'in-class' for a target."""
    if target_ct == "T cells":
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

with open(annotation_json_path) as f:
    annotations_all = json.load(f)

# normalize annotation strings, drop "Artifacts"
for slide, cluster_dict in annotations_all.items():
    for k, v in list(cluster_dict.items()):
        if v == "Artifacts":
            cluster_dict[k] = ""
        elif isinstance(v, str):
            cluster_dict[k] = v.strip()

#per-slide % change vs cohort baseline (in-class mean vs cohort in-class mean)
in_means = {}  # {(slide_key, pair_title): mean_in}
baseline_acc = {title: {"sum": 0.0, "n": 0} for (title, *_r) in PAIRS}

for slide_key in annotations_all.keys():
    h5ad_path  = input_dir / f"{slide_key}_adata.h5ad"
    out_dir    = output_root / f"{slide_key}_adata_k=200"
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

    for (title, target_ct, marker_names, _color_key) in PAIRS:
        var = _find_var(adata, marker_names)
        if var is None:
            in_means[(slide_key, title)] = np.nan
            continue

        x = _dense1(adata[:, var].X)
        in_labels = members_for(target_ct)
        mask_in = np.isin(ann.values, in_labels)  # per-cell-type selection only

        if mask_in.sum() == 0:
            in_means[(slide_key, title)] = np.nan
            continue

        xin = x[mask_in]
        mean_in = float(np.mean(xin))
        n_in    = int(mask_in.sum())

        in_means[(slide_key, title)] = mean_in
        baseline_acc[title]["sum"] += float(np.sum(xin))
        baseline_acc[title]["n"]   += n_in

# Cohort baselines (cell-count weighted)
baseline = {title: (acc["sum"]/acc["n"] if acc["n"] > 0 else np.nan)
            for title, acc in baseline_acc.items()}

# Pass 2: % change vs cohort baseline
records = []
for (slide_key, title), mean_in in in_means.items():
    b = baseline[title]
    if np.isnan(mean_in) or np.isnan(b) or b == 0:
        val = np.nan
    else:
        val = (mean_in - b) / b * 100.0
    records.append({"slide": slide_key, "pair": title, "value": val})

df1 = pd.DataFrame.from_records(records)
wide1 = df1.pivot_table(index="slide", columns="pair", values="value", aggfunc="mean")

# Order slides
ordered_slide_ids = [
    condition_to_slide[c] for c in custom_order
    if c in condition_to_slide and condition_to_slide[c] in wide1.index
]
wide1 = wide1.loc[ordered_slide_ids]

nrow = 2
ncol = max(len(ORIGINAL_PAIRS), len(CROSS_SPEC_PAIRS))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.2*ncol, 6.2*nrow), sharey=True)

# Ensure axes is 2D array
if nrow == 1:
    axes = np.array([axes])
elif ncol == 1:
    axes = axes.reshape(nrow, 1)


def _plot_col(ax, title):
    if title not in wide1.columns:
        ax.set_visible(False)
        return
    vals = wide1[title]
    colors = [POS_COLOR if (isinstance(v, (int, float)) and v >= 0) else NEG_COLOR for v in vals.values]
    ax.bar(range(len(vals)), vals.values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_title(f"{title}\nvs cohort baseline", fontsize=10)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels([slide_to_condition.get(s, s) for s in vals.index], rotation=90, fontsize=7)
    ax.axhline(0, color="black", linewidth=0.8)

# Row 1: originals
for j, (title, *_rest) in enumerate(ORIGINAL_PAIRS):
    ax = axes[0, j]
    _plot_col(ax, title)

# Row 2: cross-specificity
for j, (title, *_rest) in enumerate(CROSS_SPEC_PAIRS):
    ax = axes[1, j]
    _plot_col(ax, title)

# Hide any unused axes
for r in range(nrow):
    for c in range(ncol):
        if r == 0 and c >= len(ORIGINAL_PAIRS):
            axes[r, c].set_visible(False)
        if r == 1 and c >= len(CROSS_SPEC_PAIRS):
            axes[r, c].set_visible(False)

# Y-label on first visible axis
axes[0, 0].set_ylabel("Δ vs cohort mean (%)")

# Legend
handles = [Patch(facecolor=POS_COLOR, edgecolor="black", label="≥ baseline"),
           Patch(facecolor=NEG_COLOR, edgecolor="black", label="< baseline")]
axes[0, 0].legend(handles=handles, loc="upper left", frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(output_root / "marker_specificity_vs_cohort_baseline_per_slide.svg",
            dpi=300, bbox_inches="tight")
plt.savefig(output_root / "marker_specificity_vs_cohort_baseline_per_slide.png",
            dpi=300, bbox_inches="tight")
plt.close()

print("Saved Plot 1 (no whitelist, 2-row layout):",
      output_root / "marker_specificity_vs_cohort_baseline_per_slide.png")
