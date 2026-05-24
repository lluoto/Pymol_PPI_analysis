#!/usr/bin/env python3
"""
boltz_centroid_analysis.py
==========================
Performs three analyses on Boltz-predicted CIF structures:

  [1] Far-centroid cluster analysis
  [2] LEU14 proximity cutoff
  [3] Combined pocket criterion

Usage:
  python boltz_centroid_analysis.py --preds_dir /path/to/predictions \
      --affinity_csv /path/to/affinity_data.csv \
      --out_dir /path/to/output
"""

import argparse
import csv
import math
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed


# ─── Configuration ───────────────────────────────────────────────────────────

QW7_CENTROID = (123.93, 128.79, 119.90)

FAR_CENTROIDS = {
    "36_6vcb_3_20_ligand":  (118.6, 129.2, 118.2),
    "231_6vcb_3_20_ligand": (115.9, 136.3, 118.0),
    "390_6vcb_3_29_ligand": (117.2, 129.6, 123.4),
    "463_6vcb_3_24_ligand": (118.7, 129.6, 117.6),
    "695_6vcb_3_24_ligand": (116.4, 129.5, 120.3),
}

LEU14_CHAIN = "B"
LEU14_RESNAME = "LEU"
LEU14_RESSEQ_CANDIDATES = ("14",)

DIST_CUTOFF = 3.0

DEFAULT_PREDS_DIR = "/media/cuixi/data01/lluoto/boltz_results_input/predictions"
DEFAULT_AFFINITY_CSV = "/media/cuixi/data01/lluoto/boltz_results_input/affinity_data.csv"
DEFAULT_OUT_DIR = "/media/cuixi/data01/lluoto/boltz_results_input/centroid_analysis"


# ─── CIF parsing ─────────────────────────────────────────────────────────────

def parse_cif_simple(filepath):
    """
    Lightweight mmCIF parser — extracts atom_site records.
    Returns (atoms, ligand_resnames, ligand_asym_ids).
    """
    with open(filepath, "r") as f:
        content = f.read()

    lines = content.split("\n")

    atoms = []
    in_loop = False
    atom_fields = []
    ligand_resnames = set()
    loop_section = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if not line or line.startswith("#"):
            i += 1
            continue

        if line.startswith("loop_"):
            in_loop = True
            atom_fields = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                atom_fields.append(lines[i].strip())
                i += 1
            # Determine section from first field
            if atom_fields:
                prefix = atom_fields[0].split(".")[0]
                loop_section = {
                    "_entity": "entity",
                    "_entity_poly": "entity_poly",
                    "_struct_asym": "struct_asym",
                    "_chem_comp": "chem_comp",
                    "_atom_site": "atom_site",
                }.get(prefix, None)
            # Read data rows
            data_rows = []
            while i < len(lines):
                stripped = lines[i].strip()
                if not stripped or stripped.startswith("loop_") or stripped.startswith("_"):
                    break
                data_rows.append(stripped)
                i += 1

            # Parse based on section
            if loop_section == "chem_comp":
                for row in data_rows:
                    vals = row.split()
                    if len(vals) >= 2:
                        comp_id = vals[0]
                        comp_type = vals[1] if len(vals) > 1 else ""
                        if comp_type.lower() in ("non-polymer", "ligand", "other"):
                            ligand_resnames.add(comp_id)

            elif loop_section == "atom_site":
                field_keys = [f.split(".")[-1] for f in atom_fields]
                for row in data_rows:
                    vals = row.split()
                    if len(vals) < 8:
                        continue
                    row_dict = dict(zip(field_keys, vals))
                    try:
                        atom = {
                            "group_pdb": row_dict.get("group_PDB", ""),
                            "type_symbol": row_dict.get("type_symbol", ""),
                            "label_atom_id": row_dict.get("label_atom_id", ""),
                            "label_comp_id": row_dict.get("label_comp_id", ""),
                            "label_asym_id": row_dict.get("label_asym_id", ""),
                            "label_entity_id": row_dict.get("label_entity_id", ""),
                            "label_seq_id": row_dict.get("label_seq_id", ""),
                            "Cartn_x": float(row_dict.get("Cartn_x", 0)),
                            "Cartn_y": float(row_dict.get("Cartn_y", 0)),
                            "Cartn_z": float(row_dict.get("Cartn_z", 0)),
                        }
                        atoms.append(atom)
                    except (ValueError, KeyError):
                        continue

        else:
            i += 1

    # Identify ligand residue names
    if not ligand_resnames:
        het_resnames = set()
        for at in atoms:
            if at["group_pdb"] == "HETATM":
                het_resnames.add(at["label_comp_id"])
        solvent = {"HOH", "WAT", "NA", "CL", "MG", "CA", "ZN", "K", "SO4",
                   "PO4", "GOL", "EDO", "DMS", "ACT", "FMT", "IPA", "PEG"}
        ligand_resnames = het_resnames - solvent

    # Identify ligand asym IDs
    ligand_asym_ids = set()
    for at in atoms:
        if at["label_comp_id"] in ligand_resnames:
            ligand_asym_ids.add(at["label_asym_id"])

    # Fallback: Boltz typically puts ligand in chain C
    if not ligand_asym_ids:
        all_asym = set(at["label_asym_id"] for at in atoms)
        if "C" in all_asym:
            ligand_asym_ids = {"C"}
        else:
            for aid in sorted(all_asym):
                if aid not in ("A", "B"):
                    ligand_asym_ids.add(aid)

    return atoms, ligand_resnames, ligand_asym_ids


def compute_centroid(atoms, asym_ids):
    coords = [(a["Cartn_x"], a["Cartn_y"], a["Cartn_z"])
              for a in atoms if a["label_asym_id"] in asym_ids]
    if not coords:
        return None
    n = len(coords)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)


def distance_3d(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def find_residue_atoms(atoms, asym_id, resname, resseq_candidates):
    result = []
    for at in atoms:
        seq_id = at.get("label_seq_id", "").strip()
        if (at["label_asym_id"] == asym_id and
            at["label_comp_id"] == resname and
            seq_id in resseq_candidates):
            result.append(at)
    return result


def min_distance_between_sets(atoms_a, atoms_b):
    min_d = float("inf")
    for a in atoms_a:
        ax, ay, az = a["Cartn_x"], a["Cartn_y"], a["Cartn_z"]
        for b in atoms_b:
            d = math.sqrt((ax - b["Cartn_x"])**2 +
                          (ay - b["Cartn_y"])**2 +
                          (az - b["Cartn_z"])**2)
            if d < min_d:
                min_d = d
                if min_d < 1e-6:
                    return 0.0
    return min_d


def analyze_cif(args_tuple):
    """Analyze a single CIF file. Returns a result dict."""
    cif_path, struct_id = args_tuple

    result = {
        "struct_id": struct_id,
        "ligand_centroid": None,
        "dist_to_qw7": None,
        "dist_to_far": {},
        "leu14_min_dist": None,
        "has_leu14_contact": False,
        "in_qw7_pocket": False,
        "errors": [],
    }

    if not os.path.isfile(cif_path):
        result["errors"].append(f"File not found: {cif_path}")
        return result

    try:
        atoms, ligand_resnames, ligand_asym_ids = parse_cif_simple(cif_path)
    except Exception as e:
        result["errors"].append(f"Parse error: {e}")
        return result

    if not atoms:
        result["errors"].append("No atoms parsed")
        return result

    # Ligand centroid
    centroid = compute_centroid(atoms, ligand_asym_ids)
    if centroid is None:
        result["errors"].append("No ligand atoms found for centroid")
        return result
    result["ligand_centroid"] = centroid

    # Distance to QW7
    d = distance_3d(centroid, QW7_CENTROID)
    result["dist_to_qw7"] = d
    result["in_qw7_pocket"] = d <= DIST_CUTOFF

    # Distances to far centroids
    for name, fc in FAR_CENTROIDS.items():
        result["dist_to_far"][name] = distance_3d(centroid, fc)

    # LEU14 proximity
    leu14_atoms = find_residue_atoms(
        atoms, LEU14_CHAIN, LEU14_RESNAME, LEU14_RESSEQ_CANDIDATES
    )

    if not leu14_atoms:
        # Try any chain
        for at in atoms:
            seq_id = at.get("label_seq_id", "").strip()
            if at["label_comp_id"] == LEU14_RESNAME and seq_id in LEU14_RESSEQ_CANDIDATES:
                leu14_atoms.append(at)

    if leu14_atoms:
        ligand_atoms = [a for a in atoms if a["label_asym_id"] in ligand_asym_ids]
        if not ligand_atoms:
            ligand_atoms = [a for a in atoms if a["group_pdb"] == "HETATM"]
        if not ligand_atoms:
            leu14_id_set = set(id(a) for a in leu14_atoms)
            ligand_atoms = [a for a in atoms if id(a) not in leu14_id_set]

        if ligand_atoms:
            min_d = min_distance_between_sets(leu14_atoms, ligand_atoms)
            result["leu14_min_dist"] = min_d
            result["has_leu14_contact"] = min_d <= DIST_CUTOFF

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preds_dir", default=DEFAULT_PREDS_DIR)
    parser.add_argument("--affinity_csv", default=DEFAULT_AFFINITY_CSV)
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--max_workers", type=int, default=20)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    # Step 1: Scan predictions directory
    print(f"[1] Scanning: {args.preds_dir}")
    dirs = sorted(os.listdir(args.preds_dir))

    cif_map = {}
    for dname in dirs:
        dpath = os.path.join(args.preds_dir, dname)
        if not os.path.isdir(dpath):
            continue
        # Directory: {seed}_{ligand}_ligand
        # CIF inside: {seed}_{ligand}_ligand_model_0.cif
        cif_candidates = [
            os.path.join(dpath, f"{dname}_model_0.cif"),
        ]
        for cpath in cif_candidates:
            if os.path.isfile(cpath):
                # struct_id = everything before "_ligand"
                base = dname.rsplit("_ligand", 1)[0]
                cif_map[base] = cpath
                break

    print(f"  Found {len(cif_map)} CIF files")

    # Step 2: Match with affinity CSV
    print(f"[2] Reading affinity CSV: {args.affinity_csv}")
    affinity_map = {}

    with open(args.affinity_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["source"] != "boltz":
                continue
            ligand = row["ligand"]
            seed = row["seed"]
            sid = f"{seed}_{ligand}"
            affinity_map[sid] = row

    print(f"  Boltz entries in CSV: {len(affinity_map)}")

    common = sorted(set(cif_map.keys()) & set(affinity_map.keys()))
    missing_cif = sorted(set(affinity_map.keys()) - set(cif_map.keys()))
    missing_csv = sorted(set(cif_map.keys()) - set(affinity_map.keys()))

    print(f"  Matched: {len(common)}")
    if missing_cif:
        print(f"  WARNING: {len(missing_cif)} in CSV but no CIF found")
        for m in missing_cif[:10]:
            print(f"    Missing: {m}")

    to_analyze = [(cif_map[sid], sid) for sid in common]

    if args.limit:
        to_analyze = to_analyze[:args.limit]
        print(f"  Limited to {args.limit}")

    # Step 3: Parallel analysis
    total = len(to_analyze)
    print(f"\n[3] Analyzing {total} structures ({args.max_workers} workers)...")
    all_results = []

    if total == 0:
        print("  ERROR: No structures to analyze. Check --preds_dir and --affinity_csv paths.")
        sys.exit(1)

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        fut_map = {executor.submit(analyze_cif, t): t for t in to_analyze}
        done = 0
        for fut in as_completed(fut_map):
            done += 1
            if done % 200 == 0:
                print(f"  Progress: {done}/{total}")
            try:
                all_results.append(fut.result())
            except Exception as e:
                entry = fut_map[fut]
                all_results.append({
                    "struct_id": entry[1],
                    "ligand_centroid": None,
                    "dist_to_qw7": None,
                    "dist_to_far": {},
                    "leu14_min_dist": None,
                    "has_leu14_contact": False,
                    "in_qw7_pocket": False,
                    "errors": [str(e)],
                })

    print(f"  Done. {len(all_results)} results.")

    # Step 4: centroid_results.csv
    csv_path = os.path.join(args.out_dir, "centroid_results.csv")
    print(f"\n[4] Writing {csv_path}...")

    far_names = sorted(FAR_CENTROIDS.keys())
    csv_fields = [
        "struct_id", "ligand", "seed", "pIC50", "affinity_prob",
        "ligand_centroid_x", "ligand_centroid_y", "ligand_centroid_z",
        "dist_to_qw7", "in_qw7_pocket",
        "leu14_min_dist", "has_leu14_contact", "errors",
    ] + [f"dist_to_{fn}" for fn in far_names]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for r in sorted(all_results, key=lambda x: x["struct_id"]):
            sid = r["struct_id"]
            aff_row = affinity_map.get(sid, {})
            row = {
                "struct_id": sid,
                "ligand": aff_row.get("ligand", ""),
                "seed": aff_row.get("seed", ""),
                "pIC50": aff_row.get("pIC50", ""),
                "affinity_prob": aff_row.get("affinity_prob", ""),
                "ligand_centroid_x": f"{r['ligand_centroid'][0]:.4f}" if r["ligand_centroid"] else "",
                "ligand_centroid_y": f"{r['ligand_centroid'][1]:.4f}" if r["ligand_centroid"] else "",
                "ligand_centroid_z": f"{r['ligand_centroid'][2]:.4f}" if r["ligand_centroid"] else "",
                "dist_to_qw7": f"{r['dist_to_qw7']:.4f}" if r["dist_to_qw7"] is not None else "",
                "in_qw7_pocket": r["in_qw7_pocket"],
                "leu14_min_dist": f"{r['leu14_min_dist']:.4f}" if r["leu14_min_dist"] is not None else "",
                "has_leu14_contact": r["has_leu14_contact"],
                "errors": "; ".join(r.get("errors", [])),
            }
            for fn in far_names:
                val = r["dist_to_far"].get(fn)
                row[f"dist_to_{fn}"] = f"{val:.4f}" if val is not None else ""
            writer.writerow(row)

    # Step 5: Far-centroid clusters
    print(f"\n[5] Far-centroid cluster analysis...")
    far_stats = []
    for name, fc in FAR_CENTROIDS.items():
        count = sum(1 for r in all_results if r["dist_to_far"].get(name, float("inf")) <= DIST_CUTOFF)
        far_stats.append((name, fc, count))
        print(f"  {name}: {count} structures")

    far_path = os.path.join(args.out_dir, "far_cluster_counts.txt")
    with open(far_path, "w") as f:
        f.write(f"Far-Centroid Cluster Analysis (cutoff = {DIST_CUTOFF}A)\n")
        f.write("=" * 60 + "\n\n")
        for name, fc, count in far_stats:
            members = sorted([
                r["struct_id"] for r in all_results
                if r["dist_to_far"].get(name, float("inf")) <= DIST_CUTOFF
            ])
            f.write(f"{name}:\n")
            f.write(f"  Centroid: ({fc[0]:.2f}, {fc[1]:.2f}, {fc[2]:.2f})\n")
            f.write(f"  Within {DIST_CUTOFF}A: {count}\n")
            if members:
                f.write(f"  Members:\n")
                for m in members:
                    f.write(f"    {m}\n")
            f.write("\n")

    # Step 6: LEU14 proximity
    print(f"\n[6] LEU14 proximity analysis...")
    leu14_contact = sum(1 for r in all_results if r["has_leu14_contact"])
    qw7_pocket = sum(1 for r in all_results if r["in_qw7_pocket"])
    both = sum(1 for r in all_results if r["has_leu14_contact"] and r["in_qw7_pocket"])
    recovered = sorted([
        r["struct_id"] for r in all_results
        if r["has_leu14_contact"] and not r["in_qw7_pocket"]
    ])
    combined = sum(1 for r in all_results if r["in_qw7_pocket"] or r["has_leu14_contact"])

    leu14_path = os.path.join(args.out_dir, "leu14_proximity_stats.txt")
    with open(leu14_path, "w") as f:
        f.write(f"LEU14 Proximity Analysis (cutoff = {DIST_CUTOFF}A)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total structures: {len(all_results)}\n\n")
        f.write(f"Centroid within {DIST_CUTOFF}A of QW7: {qw7_pocket}\n")
        f.write(f"Any atom within {DIST_CUTOFF}A of LEU14 (B14): {leu14_contact}\n")
        f.write(f"  - Both QW7 + LEU14: {both}\n")
        f.write(f"  - LEU14 recovered: {len(recovered)}\n")
        f.write(f"\nCombined pocket (QW7 OR LEU14): {combined}\n\n")
        f.write(f"Recovered structures:\n")
        for s in recovered:
            f.write(f"  {s}\n")

    # Step 7: Combined pocket list
    pocket_path = os.path.join(args.out_dir, "combined_pocket_list.txt")
    with open(pocket_path, "w") as f:
        f.write(f"Combined pocket (centroid within {DIST_CUTOFF}A of QW7 OR any atom within {DIST_CUTOFF}A of LEU14)\n")
        f.write("=" * 70 + "\n\n")
        in_pocket = [r for r in all_results if r["in_qw7_pocket"] or r["has_leu14_contact"]]
        f.write(f"Total in pocket: {len(in_pocket)} / {len(all_results)}\n\n")
        for r in sorted(in_pocket, key=lambda x: x["struct_id"]):
            reasons = []
            if r["in_qw7_pocket"]:
                reasons.append(f"QW7 cent {r['dist_to_qw7']:.2f}A")
            if r["has_leu14_contact"]:
                reasons.append(f"LEU14 {r['leu14_min_dist']:.2f}A")
            f.write(f"  {r['struct_id']}  [{', '.join(reasons)}]\n")

    # Summary
    elapsed = time.time() - t0
    summary_path = os.path.join(args.out_dir, "analysis_summary.txt")
    with open(summary_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("BOLTZ CENTROID ANALYSIS - SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total: {len(all_results)}\n")
        f.write(f"Time: {elapsed:.1f}s\n\n")
        f.write(f"QW7: {qw7_pocket} within {DIST_CUTOFF}A\n\n")
        f.write(f"Far clusters:\n")
        for name, fc, count in far_stats:
            f.write(f"  {name}: {count}\n")
        f.write(f"\nLEU14 contacts: {leu14_contact}\n")
        f.write(f"Recovered: {len(recovered)}\n")
        f.write(f"Combined: {combined}\n")

    print(f"\n[Done] {args.out_dir}/")
    print(f"  centroid_results.csv")
    print(f"  far_cluster_counts.txt")
    print(f"  leu14_proximity_stats.txt")
    print(f"  combined_pocket_list.txt")
    print(f"  analysis_summary.txt")
    print(f"  Elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
