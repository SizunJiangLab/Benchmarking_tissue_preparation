import os
import sys
import argparse
import glob
import numpy as np
import tifffile
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from tqdm import tqdm
import pandas as pd
import matplotlib.colors as mcolors

CELLTYPE_COLOR = {
    # T-cell family
    #"T cells": "#E41A1C",          # Set1 red
    "CD8+ T cells": "#FB6A4A",     # custom salmon
    "CD4+ T cells": "#FCBBA1",     # lighter salmon
    "Tregs": "#FEB24C",            # orange

    # Other lineages (Set1)
    "B cells": "#377EB8",          # Set1 blue
    "Epithelial": "#984EA3",       # Set1 purple
    #"Epithelial and B cells": "#4DAF4A",  # Set1 green
    "Macrophages": "#A65628",      # Set1 brown
    #"Neutrophils and Other": "#F781BF", # Set1 pink

    "Other": "#525252",          # dark gray
}

MIXED_COLOR = "#bdbdbd"  # default/fallback color

# include Mixed explicitly so it shows in the legend
PALETTE = {**CELLTYPE_COLOR, "Mixed": MIXED_COLOR}

celltype_colors_rgb = {ct: mcolors.to_rgb(hx) for ct, hx in PALETTE.items()}
legend_elements = [Patch(facecolor=col, label=ct, edgecolor="black")
                   for ct, col in celltype_colors_rgb.items()]

BASE = "/registered_report"

def pad_slide(slide_num: int) -> str:
    if not (1 <= slide_num <= 24):
        raise ValueError(f"Slide must be in 1..24, got {slide_num}")
    return f"{slide_num:02d}"

def clustering_csv(slide_str: str) -> str:
    """Return the first CSV found for this slide's clustering results (zero-padded path)."""
    pattern = os.path.join(BASE, f"output/slide_{slide_str}_adata_k=200/clustering_results/*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No CSV found matching: {pattern}")
    return files[0]

def output_dir(slide_str: str) -> str:
    d = os.path.join(BASE, f"output/slide_{slide_str}_adata_k=200")
    os.makedirs(d, exist_ok=True)
    return d

def fov_paths(slide_str: str, fov: int):
    """
    MESMER filenames use *un-padded* slide numbers.
    Convert slide_str ('07') -> '7' for overlay/mask file names.
    """
    slide_unpadded = str(int(slide_str))
    outline = os.path.join(BASE, f"MESMER_overlay/overlay_{slide_unpadded}_FOV{fov}.tiff")
    seg     = os.path.join(BASE, f"MESMER_mask/segmentation_mask_{slide_unpadded}_FOV{fov}.tiff")
    return outline, seg

def verify_all_files_present(slide_str: str) -> bool:
    """Check that both FOV1 and FOV2 overlay+seg files exist. If any missing, print and return False."""
    missing = []
    for fov in (1, 2):
        outline, seg = fov_paths(slide_str, fov)
        if not os.path.exists(outline):
            missing.append(outline)
        if not os.path.exists(seg):
            missing.append(seg)
    if missing:
        print(f"For slide {slide_str}, the following required file(s) are missing:")
        for m in missing:
            print(f"   - {m}")
        print("Skipping this slide.")
        return False
    return True

def process_single_core(slide_str: str, core_id: str, annotation_df: pd.DataFrame, out_dir: str):
    """Render a colored segmentation for a single core and save PNG."""
    try:
        fov = 1 if core_id.endswith("_1") else 2
        outline_path, seg_path = fov_paths(slide_str, fov)

        seg = tifffile.imread(seg_path).astype(np.int64)
        outline = tifffile.imread(outline_path)[..., 0]
        outline_mask = (outline == 1)

        core_ann = annotation_df[annotation_df["core"] == core_id].copy()
        core_ann["cell_label"] = pd.to_numeric(core_ann["cell_label"], errors="coerce").astype("Int64")
        core_ann = core_ann.dropna(subset=["cell_label"]).astype({"cell_label": np.int64})

        label_to_type = core_ann.set_index("cell_label")["annotation"].to_dict()

        max_label = int(seg.max()) if seg.size else 0
        color_map = np.zeros((max_label + 1, 3), dtype=np.float32)

        present_labels = np.unique(seg[seg > 0])
        for lab in present_labels:
            ct = label_to_type.get(int(lab), None)

            # normalize/clean the annotation string
            if isinstance(ct, str):
                ct = ct.strip()
            # default to Mixed if missing/empty or not in our palette
            if not ct or ct not in celltype_colors_rgb:
                ct = "Mixed"

            color_map[lab] = celltype_colors_rgb.get(ct, mcolors.to_rgb(MIXED_COLOR))

        colored = color_map[seg]
        colored[outline_mask] = (0, 0, 0)  # white outline

        plt.figure(figsize=(20, 20), facecolor="white")
        plt.imshow(colored)
        plt.axis("off")
        plt.legend(
            handles=legend_elements,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            frameon=False,
            fontsize=14,
            labelcolor="white",
            prop={"weight": "bold"},
        )
        out_png = os.path.join(out_dir, f"{core_id}.png")
        out_svg = os.path.join(out_dir, f"{core_id}.svg")
        plt.savefig(out_png, bbox_inches="tight", dpi=200, facecolor="black", pad_inches=0)
        plt.savefig(out_svg, bbox_inches="tight", facecolor="black", pad_inches=0)
        plt.close()
        print(f"Saved: {out_png}")

    except Exception as e:
        print(f"Failed for {core_id}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Render phenotype maps per slide.")
    parser.add_argument("--slide", type=int, default=None, help="Slide number 1–24 (zero-padded internally).")
    args = parser.parse_args()

    # Resolve slide number (arg > SLURM_ARRAY_TASK_ID)
    slide_num = args.slide
    if slide_num is None:
        env_id = os.getenv("SLURM_ARRAY_TASK_ID")
        if env_id is None:
            print("Provide --slide or run via SLURM array with SLURM_ARRAY_TASK_ID.")
            sys.exit(1)
        slide_num = int(env_id)

    slide_str = pad_slide(slide_num)

    # Check inputs first; if anything missing, exit early
    if not verify_all_files_present(slide_str):
        return

    try:
        csv_path = clustering_csv(slide_str)
    except FileNotFoundError as e:
        print(f"For slide {slide_str}, {e}")
        return

    print(f"Using annotations: {csv_path}")
    out_dir = output_dir(slide_str)

    # Load and filter annotations for this slide only (both cores)
    # Support both '07_1_' and '7_1_' unit_ids
    ann = pd.read_csv(csv_path, usecols=["unit_ids", "cluster_ids", "annotation"])
    slide_unpadded = str(int(slide_str))
    prefixes = (
        f"{slide_str}_1_", f"{slide_str}_2_",
        f"{slide_unpadded}_1_", f"{slide_unpadded}_2_"
    )
    ann = ann[ann["unit_ids"].str.startswith(prefixes)].copy()

    if ann.empty:
        print(f"No matching unit_ids for slide {slide_str} "
              f"(expected prefixes {', '.join(repr(p) for p in prefixes)}).")
        return

    parts = ann["unit_ids"].str.split("_")
    ann["core"] = parts.str[0] + "_" + parts.str[1]
    ann["cell_label"] = parts.str[-1]

    candidate_cores = set(ann["core"].unique())
    cores_to_run = []
    for fov in (1, 2):
        padded = f"{slide_str}_{fov}"
        unpadded = f"{slide_unpadded}_{fov}"
        if padded in candidate_cores:
            cores_to_run.append(padded)
        elif unpadded in candidate_cores:
            cores_to_run.append(unpadded)
        else:
            print(f"No annotations for core {padded} or {unpadded}; skipping.")

    for core_id in tqdm(cores_to_run, desc=f"Rendering slide {slide_str}"):
        process_single_core(slide_str, core_id, ann, out_dir)

    print(f"\nDone with slide {slide_str}. Outputs in: {out_dir}")

if __name__ == "__main__":
    main()
