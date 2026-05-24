#!/usr/bin/env python3
"""
annotate_affinity_histogram.py
==============================
Generates annotated affinity histogram — all structures in one color,
annotations directly on the chart, no legend.

Usage:
  python annotate_affinity_histogram.py --centroid_csv centroid_results.csv

Outputs:
  affinity_histogram_annotated.png
  affinity_histogram_data.csv
"""

import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


DEFAULT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "centroid_results.csv")
DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "affinity_histogram_annotated.png")


def read_csv(path):
    with open(path, "r") as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--centroid_csv", default=DEFAULT_CSV)
    parser.add_argument("--output", default=DEFAULT_OUT)
    parser.add_argument("--metric", default="pIC50",
                        choices=["pIC50", "affinity_prob"])
    parser.add_argument("--bins", type=int, default=50)
    args = parser.parse_args()

    if not os.path.isfile(args.centroid_csv):
        print(f"ERROR: {args.centroid_csv} not found")
        return

    rows = read_csv(args.centroid_csv)
    print(f"Rows: {len(rows)}")

    # Extract values
    all_vals = []
    struct_data = []
    for r in rows:
        try:
            v = float(r[args.metric])
            all_vals.append(v)
            struct_data.append(r["struct_id"])
        except (ValueError, KeyError, TypeError):
            pass

    all_vals = np.array(all_vals)
    print(f"Valid {args.metric}: {len(all_vals)}")

    mu = np.mean(all_vals)
    sigma = np.std(all_vals)
    median = np.median(all_vals)
    q25 = np.percentile(all_vals, 25)
    q75 = np.percentile(all_vals, 75)
    vmin = np.min(all_vals)
    vmax = np.max(all_vals)

    # ── Pocket/far classification (899 centroid + LEU14) ──
    # Compute distance from mean centroid for each structure
    cx_list, cy_list, cz_list = [], [], []
    for r in rows:
        try:
            cx_list.append(float(r["ligand_centroid_x"]))
            cy_list.append(float(r["ligand_centroid_y"]))
            cz_list.append(float(r["ligand_centroid_z"]))
        except:
            pass
    import math
    from collections import namedtuple
    if cx_list:
        mc_x, mc_y, mc_z = np.mean(cx_list), np.mean(cy_list), np.mean(cz_list)
        struct_dists = []
        for r in rows:
            try:
                cx = float(r["ligand_centroid_x"])
                cy = float(r["ligand_centroid_y"])
                cz = float(r["ligand_centroid_z"])
                d = math.sqrt((cx-mc_x)**2 + (cy-mc_y)**2 + (cz-mc_z)**2)
                leu14 = r.get("has_leu14_contact", "").strip().lower() == "true"
                pic50 = float(r[args.metric]) if r.get(args.metric) else None
                struct_dists.append((d, leu14, pic50))
            except:
                pass
        struct_dists.sort(key=lambda x: x[0])
        n_pocket = 899
        pocket_vals = []
        far_vals = []
        for i, (d, leu14, pic50) in enumerate(struct_dists):
            if pic50 is None:
                continue
            if i < n_pocket:
                pocket_vals.append(pic50)
            elif leu14:
                pocket_vals.append(pic50)  # LEU14 recovered → pocket
            else:
                far_vals.append(pic50)
        pocket_vals = np.array(pocket_vals)
        far_vals = np.array(far_vals)
        pocket_min = np.min(pocket_vals) if len(pocket_vals) > 0 else None
        pocket_max = np.max(pocket_vals) if len(pocket_vals) > 0 else None
        pocket_mean = np.mean(pocket_vals) if len(pocket_vals) > 0 else None
        far_min = np.min(far_vals) if len(far_vals) > 0 else None
        far_max = np.max(far_vals) if len(far_vals) > 0 else None
        print(f"Pocket: n={len(pocket_vals)} range=[{pocket_min:.4f},{pocket_max:.4f}]")
        print(f"Far: n={len(far_vals)} range=[{far_min:.4f},{far_max:.4f}]")
    else:
        pocket_min = pocket_max = pocket_mean = None
        far_min = far_max = None
        pocket_vals = far_vals = np.array([])

    # Per-ligand groups
    def extract_ligand(sid):
        parts = sid.split("_", 1)
        return parts[1] if len(parts) > 1 else sid
    lig_groups = defaultdict(list)
    for r in rows:
        sid = r["struct_id"]
        try:
            v = float(r[args.metric])
            lig_groups[extract_ligand(sid)].append(v)
        except:
            pass

    # ── Create figure ──
    fig, ax = plt.subplots(figsize=(14, 7))
    bins = np.linspace(vmin, vmax, args.bins)

    # Single-color histogram
    ax.hist(all_vals, bins=bins, density=True,
            color="#4682B4", alpha=0.65, edgecolor="white", linewidth=0.3)

    # Gaussian fit
    fit_x = np.linspace(vmin, vmax, 300)
    fit_y = stats.norm.pdf(fit_x, mu, sigma)
    ax.plot(fit_x, fit_y, color="#D32F2F", linewidth=2.0, linestyle="--", alpha=0.8)

    # NO legend

    # ── Annotations on the chart ──

    # 1) Stats box (top-left)
    if pocket_min is not None:
        stats_text = (
            f"All (n={len(all_vals)}):\n"
            f"  Mean={mu:.4f}  Std={sigma:.4f}\n"
            f"  Median={median:.4f}\n\n"
            f"Pocket (n={len(pocket_vals)}):\n"
            f"  Mean={pocket_mean:.4f}\n"
            f"  Range: [{pocket_min:.4f}, {pocket_max:.4f}]\n\n"
            f"Far (n={len(far_vals)}):\n"
            f"  Range: [{far_min:.4f}, {far_max:.4f}]"
        )
    else:
        stats_text = (
            f"n = {len(all_vals)}\n"
            f"Mean = {mu:.4f}\n"
            f"Std = {sigma:.4f}\n"
            f"Median = {median:.4f}"
        )
    ax.text(0.02, 0.97, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="#888888"))

    # 2) Range annotation (top-right)
    range_text = f"Range: [{vmin:.4f}, {vmax:.4f}]"
    ax.text(0.98, 0.97, range_text, transform=ax.transAxes,
            fontsize=10, verticalalignment="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF9C4", alpha=0.85, edgecolor="#FBC02D"))

    # 3) Per-ligand mean annotations
    colors = {"6vcb_3_20": "#E53935", "6vcb_3_24": "#1E88E5", "6vcb_3_29": "#43A047"}
    # Use axes fraction for positioning (independent of y-scale)
    for i, lig in enumerate(sorted(lig_groups)):
        vals = lig_groups[lig]
        if len(vals) < 2:
            continue
        lm = np.mean(vals)
        ls = np.std(vals)
        c = colors.get(lig, "#666666")
        # Vertical line at mean
        ax.axvline(x=lm, color=c, linestyle=":", linewidth=1.2, alpha=0.6)
        # Text with mean value
        ax.annotate(f"{lig}: {lm:.4f}",
                    xy=(0.98, 0.88 - i * 0.06),
                    xycoords="axes fraction",
                    fontsize=9, color=c, fontweight="bold",
                    ha="right", va="top",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8, edgecolor=c))

    # 4) Pocket range vertical lines
    if pocket_min is not None:
        ax.axvline(x=pocket_min, color="#2E7D32", linestyle="--", linewidth=1.5, alpha=0.7)
        ax.axvline(x=pocket_max, color="#2E7D32", linestyle="--", linewidth=1.5, alpha=0.7)

    # 5) Threshold markers on x-axis
    thresholds = [0.0, 0.3, 0.5, 0.7, 1.0, 1.3]
    for th in thresholds:
        n_below = int(np.sum(all_vals < th))
        pct = n_below / len(all_vals) * 100
        ax.axvline(x=th, color="#9E9E9E", linestyle="--", linewidth=0.5, alpha=0.3)

    # 5) Binding site center annotation
    # Read centroid center from data (mean of all centroids)
    cx_vals, cy_vals, cz_vals = [], [], []
    for r in rows:
        try:
            cx_vals.append(float(r["ligand_centroid_x"]))
            cy_vals.append(float(r["ligand_centroid_y"]))
            cz_vals.append(float(r["ligand_centroid_z"]))
        except:
            pass
    if cx_vals:
        center_text = (f"Centroid center: ({np.mean(cx_vals):.2f}, "
                       f"{np.mean(cy_vals):.2f}, {np.mean(cz_vals):.2f})")
    else:
        center_text = ""
    ax.text(0.5, -0.14, f"Boltz {args.metric} distribution (n={len(all_vals)})",
            transform=ax.transAxes, fontsize=9, color="#888888",
            ha="center", va="top")

    # ── Axes ──
    ax.set_xlabel(args.metric, fontsize=13)
    ax.set_ylabel("Density", fontsize=13)
    ax.set_title("Boltz Predicted Affinity Distribution", fontsize=15, fontweight="bold")
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Saved: {args.output}")

    # ── Data CSV ──
    data_csv = args.output.replace(".png", "_data.csv")
    with open(data_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["struct_id", args.metric])
        for r in rows:
            try:
                v = float(r[args.metric])
                writer.writerow([r["struct_id"], f"{v:.6f}"])
            except:
                pass
    print(f"Saved: {data_csv}")


if __name__ == "__main__":
    main()
