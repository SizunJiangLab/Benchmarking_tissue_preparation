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
    "T cells": "#E41A1C",          # Set1 red
    "CD8+ T cells": "#FB6A4A",     # custom salmon
    "CD8- T cells": "#FCBBA1",     # lighter salmon
    "Tregs": "#FEB24C",            # orange

    # Other lineages (Set1)
    "B cells": "#377EB8",          # Set1 blue
    "Epithelial": "#984EA3",       # Set1 purple
    "Epithelial and B cells": "#4DAF4A",  # Set1 green
    "Macrophages": "#A65628",      # Set1 brown
    "Neutrophils and Unknown": "#F781BF", # Set1 pink

    "Unknown": "#525252",          # dark gray
}

MIXED_COLOR = "#bdbdbd"  # light gray for Mixed

annotation_json_path = Path(
    "/registered_report/cluster_annotations.json"
)
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
    if label == "Unknown":
        return "Unknown"
    if label in CELLTYPE_COLOR:
        return label
    return "Mixed"

#remove Artifacts from plotting
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

# reduce to actual labels + mixed + unknown
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

# slide order by condition, keeping only present
ordered_slide_ids = [
    condition_to_slide[c] for c in custom_order
    if c in condition_to_slide and condition_to_slide[c] in reduced_df_frac.index
]
reduced_df_frac = reduced_df_frac.loc[ordered_slide_ids]

# label order
defined_order = [lab for lab in CELLTYPE_COLOR.keys()
                 if lab in reduced_df_frac.columns and lab != "Unknown"]

tail = [lab for lab in ["Mixed", "Unknown"] if lab in reduced_df_frac.columns]
label_list = defined_order + tail

reduced_df_frac = reduced_df_frac.reindex(columns=label_list, fill_value=0)
-
label_color_map = {
    **{lab: mcolors.to_rgb(CELLTYPE_COLOR[lab]) for lab in defined_order},
    "Mixed": mcolors.to_rgb(MIXED_COLOR),
    "Unknown": mcolors.to_rgb(CELLTYPE_COLOR["Unknown"]),
}
color_list = [label_color_map[l] for l in label_list]

ax = reduced_df_frac.plot(
    kind="bar",
    stacked=True,
    figsize=(18, 8),
    color=color_list,
    edgecolor="black",
    linewidth=0.4,
)
plt.ylabel("Fraction of annotated cells")
plt.title("Annotation Breakdown: Actual Labels vs Mixed vs Unknown (per slide)")
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
