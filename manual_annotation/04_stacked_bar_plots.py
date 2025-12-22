import json
from pathlib import Path
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import importlib
from clustering import subcluster
importlib.reload(subcluster)
from clustering.subcluster import ClusteringResult, ClusteringResultManager


CELLTYPE_COLOR = {
    # T-cell family
    #"T cells": "#E41A1C",          # Set1 red
    "Macrophages": "#A65628",      # Set1 brown
    "Epithelial": "#984EA3",       # Set1 purple
    "B cells": "#377EB8",          # Set1 blue
    "Tregs": "#FEB24C",            # orange
    "CD4+ T cells": "#FCBBA1",     # lighter salmon
    "CD8+ T cells": "#FB6A4A",     # custom salmon

    # Other lineages (Set1)
    #"Epithelial and B cells": "#4DAF4A",  # Set1 green
    #"Neutrophils and Other": "#F781BF", # Set1 pink

    "Other": "#525252",          # dark gray
}


MIXED_COLOR = "#BDBDBD"  # light gray for Mixed

annotation_json_path = Path("/registered_report/cluster_annotations.json")
output_root = Path("/registered_report/output")
input_dir = Path("/registered_report/input/h5ad")

# Slides visual order (conditions)
custom_order = [
    "BIDMC_13", "BIDMC_21", "BIDMC_5", "BIDMC_15", "BIDMC_23", "BIDMC_17",
    "BIDMC_24", "BIDMC_8", "BIDMC_7", "BIDMC_9", "BIDMC_22", "BIDMC_16",
    "BIDMC_14", "BIDMC_19", "BIDMC_1", "BIDMC_2", "BIDMC_6", "BIDMC_20",
    "BIDMC_3", "BIDMC_18", "BIDMC_4", "BIDMC_12", "BIDMC_11", "BIDMC_10"
]
slide_to_condition = {f"slide_{int(b.split('_')[1]):02d}": b for b in custom_order}
condition_to_slide = {v: k for k, v in slide_to_condition.items()}


def reduce_annotation_for_category(label: str) -> str:
    """
    Keep known labels from CELLTYPE_COLOR, convert unknown-but-nonempty to Mixed,
    preserve explicit 'Other'.
    """
    if label == "Other":
        return "Other"
    if label in CELLTYPE_COLOR:
        return label
    return "Mixed"

# --- Remove Artifacts from plotting ---
with open(annotation_json_path) as f:
    annotations_all = json.load(f)

for slide, cluster_dict in annotations_all.items():
    for k in list(cluster_dict.keys()):
        if cluster_dict[k] == "Artifacts":
            annotations_all[slide][k] = ""

per_slide_annotation_counts: dict[str, dict[str, int]] = {}

for slide_key, annotations in annotations_all.items():
    slide_filename = f"{slide_key}_adata.h5ad"
    h5ad_path = input_dir / slide_filename
    output_dir = output_root / f"{slide_key}_adata_k=200"
    geojson_dir = output_dir / "geojson"
    if not h5ad_path.exists() or not geojson_dir.is_dir():
        continue

    # pick newest clustering id
    clustering_dirs = [d for d in geojson_dir.iterdir() if d.is_dir()]
    if not clustering_dirs:
        continue
    clustering_id = max(clustering_dirs, key=lambda d: d.stat().st_mtime).name

    # load clustering result & apply annotations
    clustering_result = ClusteringResult.pop(clustering_id=clustering_id, output_dir=output_dir)
    df = clustering_result.cluster_df
    normalized_annotations = {str(k).strip(): v for k, v in annotations.items()}
    clustering_result.add_annotation(normalized_annotations)
    clustering_result.save(output_dir)

    # manager for summary
    adata = ad.read_h5ad(h5ad_path)
    manager = ClusteringResultManager(output_dir=output_dir, unit_ids=adata.obs.index)

    # count non-empty annotations
    counts = manager.summary_df.annotation[manager.summary_df.annotation != ""].value_counts()
    per_slide_annotation_counts[slide_key] = {k: int(v) for k, v in counts.items()}

# --- Reduce to (actual labels) + Mixed + Other ---
reduced_annotation_counts: dict[str, dict[str, int]] = {}
for slide, counts in per_slide_annotation_counts.items():
    agg: dict[str, int] = {}
    for label, cnt in counts.items():
        r = reduce_annotation_for_category(label)
        if r == "":
            continue
        agg[r] = agg.get(r, 0) + int(cnt)
    reduced_annotation_counts[slide] = agg

reduced_df = pd.DataFrame(reduced_annotation_counts).T.fillna(0)
reduced_df_frac = reduced_df.div(reduced_df.sum(axis=1), axis=0).fillna(0)

# --- Slide order by condition (keep only present) ---
ordered_slide_ids = [
    condition_to_slide[c] for c in custom_order
    if c in condition_to_slide and condition_to_slide[c] in reduced_df_frac.index
]
reduced_df_frac = reduced_df_frac.loc[ordered_slide_ids]

# --- Label order ---
# 1) Core cell types (everything except the tail group and Other)
core_order = [lab for lab in CELLTYPE_COLOR.keys()
              if lab in reduced_df_frac.columns
              and lab not in {"Other", "Neutrophils and Other", "Epithelial and B cells"}]

# 2) Tail group (place these next to Mixed/Other on the right)
tail_like = [lab for lab in ["Epithelial and B cells", "Neutrophils and Other"] if lab in reduced_df_frac.columns]

# 3) True tail
tail1 = [lab for lab in ["Mixed"] if lab in reduced_df_frac.columns]
tail2 = [lab for lab in ["Other"] if lab in reduced_df_frac.columns]

label_list = tail2 + tail1 + core_order
reduced_df_frac = reduced_df_frac.reindex(columns=label_list, fill_value=0)

# -------- CSV export: counts + fractions for Mixed and Other --------
counts_df = reduced_df.reindex(index=ordered_slide_ids).fillna(0)
counts_df = counts_df.astype(int)

mixed_counts = counts_df.get("Mixed", pd.Series(0, index=counts_df.index))
unknown_counts = counts_df.get("Other", pd.Series(0, index=counts_df.index))

# Total annotated cells per slide (after removing Artifacts)
total_counts = counts_df.sum(axis=1)

# Fractions from the same data used for the plot
mixed_frac = (mixed_counts / total_counts).fillna(0)
unknown_frac = (unknown_counts / total_counts).fillna(0)

export_df = pd.DataFrame({
    "Condition": [slide_to_condition.get(slide_id, slide_id) for slide_id in counts_df.index],
    "TotalCells": total_counts.astype(int),
    "MixedCount": mixed_counts.astype(int),
    "OtherCount": unknown_counts.astype(int),
    "MixedFraction": mixed_frac,
    "OtherFraction": unknown_frac,
})

csv_path = output_root / "mixed_unknown_counts_and_fractions_by_slide.csv"
export_df.to_csv(csv_path, index_label="slide_id")
print(f"Saved CSV with counts and fractions: {csv_path}")


# --- Colors aligned to label_list ---
label_color_map = {
    **{lab: mcolors.to_rgb(CELLTYPE_COLOR[lab]) for lab in CELLTYPE_COLOR if lab in label_list and lab != "Other"},
    "Mixed": mcolors.to_rgb(MIXED_COLOR),
    "Other": mcolors.to_rgb(CELLTYPE_COLOR["Other"]),
}
color_list = [label_color_map[l] for l in label_list]

# --- Plot ---
ax = reduced_df_frac.plot(
    kind="bar",
    stacked=True,
    figsize=(18, 8),
    color=color_list,
    edgecolor="black",
    linewidth=0.4,
)
plt.ylabel("Fraction of annotated cells")
plt.xlabel("Conditions ranked from best to worst (left to right) based on CV values")  # requested x label
plt.title("Annotated cell types per condition")
plt.xticks(rotation=90)
plt.ylim(0, 1.0)
plt.legend(
    title="Annotation",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize="small",
)
plt.tight_layout()
plt.savefig(output_root / "annotation_mixed_vs_single_vs_unknown.svg", bbox_inches="tight")
plt.close()

# --- Plot ---
ax = reduced_df_frac.plot(
    kind="bar",
    stacked=True,
    figsize=(18, 8),
    color=color_list,
    edgecolor="black",
    linewidth=0.4,
)

# Reverse the legend order
handles, labels = ax.get_legend_handles_labels()

plt.ylabel("Fraction of annotated cells")
plt.xlabel("Conditions ranked from best to worst (left to right) based on CV values")
plt.title("Annotated cell types per condition")
plt.xticks(rotation=90)
plt.ylim(0, 1.0)

# Apply reversed handles and labels
plt.legend(
    handles[::-1], 
    labels[::-1],
    title="Annotation",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize="small",
)

plt.tight_layout()
plt.savefig(output_root / "annotation_mixed_vs_single_vs_unknown_reversed.svg", bbox_inches="tight")
plt.close()
