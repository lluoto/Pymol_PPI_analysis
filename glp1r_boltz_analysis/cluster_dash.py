import json, os, subprocess, re, math
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

BASE = "/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/interactions_3plus"
OUTPUT = "/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/3plus_clusters_dash"
RECEPTOR = "/media/cuixi/data01/lluoto/6vcb_receptor.pdb"
PYMOL = "/media/cuixi/data01/lluoto/PyMOL-2.6.2_appveyor1945-Linux-x86_64-py311/pymol/pymol"
PLIP_JSON = "/media/cuixi/data01/lluoto/vina_plip_results.json"

os.makedirs(OUTPUT, exist_ok=True)
LIG_DIR = os.path.join(OUTPUT, "ligands")
os.makedirs(LIG_DIR, exist_ok=True)

# Residue -> ProLIF interaction type (matching prolif_pymol_vis.py)
def classify_plip_type(resname):
    name = resname.upper()
    if name in ("LYS", "ARG", "HIS"):
        return "PiCation"
    elif name in ("ASP", "GLU"):
        return "SaltBridge"
    elif name in ("TYR", "TRP", "PHE"):
        return "PiStacking"
    elif name in ("SER", "THR", "ASN", "GLN", "CYS", "MET"):
        return "HydrogenBond"
    elif name == "CL":
        return "HalogenBond"
    else:
        return "Hydrophobic"

# MOLAND colors from prolif_pymol_vis.py (exact match)
MOLAND_FIRST = {
    "HydrogenBond": [0.204, 0.561, 0.655],
    "Hydrophobic":  [0.69, 0.67, 0.72],
    "PiStacking":   [0.812, 0.686, 0.831],
    "EdgeToFace":   [0.6, 0.4, 0.8],
    "FaceToFace":   [0.7, 0.5, 0.9],
    "CationPi":     [0.97, 0.79, 0.49],
    "PiCation":     [0.97, 0.79, 0.49],
    "SaltBridge":   [0.929, 0.759, 0.662],
    "HalogenBond":  [0.251, 0.718, 0.678],
    "Metal":        [1.0, 1.0, 0.0],
    "VdWContact":   [0.89, 0.94, 0.98],
    "WaterBridge":  [0.631, 0.769, 0.992],
}

# Load PLIP for interaction type lookup
print("Loading PLIP results...")
plip_data = {}
if os.path.exists(PLIP_JSON):
    with open(PLIP_JSON) as f:
        plip_all = json.load(f)
    # Index by seed number for quick lookup (key format: {prefix}_{ligand}_{seed}_out)
    for sid, data in plip_all.items():
        m = re.search(r'_(\d+)_out$', sid)
        if m:
            seed = m.group(1)
            plip_data[seed] = data
    print("  Loaded %d PLIP entries" % len(plip_data))
else:
    print("  WARN: PLIP JSON not found at %s" % PLIP_JSON)

def parse_centroid(path):
    coords = []
    for line in open(path):
        if line.startswith("HETATM") and len(line) > 22 and line[21:22] == "L":
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except:
                pass
    if not coords:
        return None
    return (sum(c[0]/len(coords) for c in coords), sum(c[1]/len(coords) for c in coords), sum(c[2]/len(coords) for c in coords))

def extract_ligand(src, dst):
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith("HETATM") and len(line) > 22 and line[21:22] == "L":
                fo.write(line)
    return n

def dist3(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

# Read entries
entries = []
for fn in sorted(os.listdir(BASE)):
    if not fn.endswith(".pdb"):
        continue
    path = os.path.join(BASE, fn)
    c = parse_centroid(path)
    if c is None:
        continue
    m = re.match(r'(\S+)_pic50_([\d.]+)_seed_(\d+)\.pdb', fn)
    pic50 = float(m.group(2)) if m else 0.0
    seed = m.group(3) if m else "0"
    lfn = fn.replace(".pdb", "_lig.pdb")
    lp = os.path.join(LIG_DIR, lfn)
    extract_ligand(path, lp)
    entries.append({"file": fn, "path": path, "lpath": lp, "centroid": c, "pic50": pic50, "seed": seed,
                    "int_json": os.path.join(BASE, fn.replace(".pdb", "_interactions.json"))})

print("Total: %d" % len(entries))

# K-means
X = np.array([e["centroid"] for e in entries])
best_k, best_sil = 2, -1
for k in range(2, min(15, len(entries)//5+1)):
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    sil = silhouette_score(X, km.fit_predict(X))
    if sil > best_sil:
        best_sil, best_k = sil, k
print("Optimal k=%d sil=%.4f" % (best_k, best_sil))

km = KMeans(n_clusters=best_k, random_state=42, n_init=10)
labels = km.fit_predict(X)

cmap = {}
for i, lab in enumerate(labels):
    cmap.setdefault(int(lab), []).append(i)

COLORS = ["marine","red","green","orange","purple","hotpink"]

# Build cluster info (same as before)
cluster_info = []
for label, members in sorted(cmap.items()):
    idx = list(members)
    med = min(idx, key=lambda i: sum(dist3(entries[i]["centroid"], entries[j]["centroid"]) for j in idx))
    pics = [entries[j]["pic50"] for j in idx]
    cid = label + 1
    ci = {"id": cid, "size": len(idx), "medoid_idx": med,
          "medoid_file": entries[med]["file"], "medoid_path": entries[med]["path"],
          "medoid_lpath": entries[med]["lpath"],
          "medoid_seed": entries[med]["seed"],
          "medoid_int_json": entries[med]["int_json"],
          "medoid_pic50": entries[med]["pic50"],
          "pic50_min": min(pics), "pic50_max": max(pics), "pic50_mean": sum(pics)/len(pics)}
    cluster_info.append(ci)
    print("C%d: n=%d medoid=%s seed=%s pic50=%.3f" % (cid, len(idx), entries[med]["file"], entries[med]["seed"], entries[med]["pic50"]))

# Helper: get PLIP interaction type for a given residue pair
def get_plip_type(details_list, chain, resnum, ligand_atom_name):
    """Try to classify interaction from PLIP details based on residue."""
    for d in details_list:
        if (d.get("chain") == chain and str(d.get("res_num")) == str(resnum)):
            return d.get("type", "hydrophobic")
    return "hydrophobic"

# --- Generate PML with interaction dashes, NO visible ligand ---
for ci in cluster_info:
    cid = ci["id"]
    col = COLORS[cid % len(COLORS)]
    pml_c = os.path.join(OUTPUT, "cluster_%d.pml" % cid)
    pse_c = os.path.join(OUTPUT, "cluster_%d.pse" % cid)

    # Load interaction JSON for medoid
    int_data = []
    if os.path.exists(ci["medoid_int_json"]):
        with open(ci["medoid_int_json"]) as f:
            int_all = json.load(f)
        int_data = int_all.get("interactions", [])

    # Chain mapping: interaction JSON (A/B) -> PDB (R/P)
    CHAIN_MAP = {"A": "R", "B": "P"}
    with open(pml_c, "w") as f:
        f.write("from pymol import cmd\n")
        f.write("cmd.bg_color('white')\n")
        f.write("cmd.set('antialias', 1)\n")
        f.write("cmd.set('light_count', 10)\n")
        f.write("cmd.set('spec_count', 1)\n")
        f.write("cmd.set('shininess', 10)\n")
        f.write("cmd.set('specular', 0.6)\n")
        f.write("cmd.set('ambient', 0)\n")
        f.write("cmd.set('direct', 0.45)\n")
        f.write("cmd.set('reflect', 0.75)\n")
        f.write("cmd.set('ray_shadow_decay_factor', 0.1)\n")
        f.write("cmd.set('ray_shadow_decay_range', 2)\n")
        f.write("cmd.set('cartoon_fancy_helices', 1)\n")
        f.write("cmd.set('cartoon_fancy_sheets', 1)\n")
        f.write("cmd.set('cartoon_rect_width', 0.3)\n")
        f.write("cmd.set('cartoon_rect_length', 2)\n")
        f.write("cmd.set('cartoon_loop_radius', 0.15)\n")
        f.write("cmd.set('dash_gap', 0.3)\n")
        f.write("cmd.set('dash_length', 0.3)\n")
        f.write("cmd.set('dash_radius', 0.06)\n")
        # Element rainbow color scheme (CPK-like)
        f.write("cmd.set_color('elem_N', [0.0, 0.0, 1.0])\n")
        f.write("cmd.set_color('elem_O', [1.0, 0.0, 0.0])\n")
        f.write("cmd.set_color('elem_S', [1.0, 1.0, 0.0])\n")
        f.write("cmd.set_color('elem_Cl', [0.0, 1.0, 0.0])\n")
        f.write("cmd.set_color('elem_F', [0.0, 1.0, 1.0])\n")
        f.write("cmd.set_color('elem_P', [1.0, 0.6, 0.0])\n")
        f.write("cmd.set_color('protein_gray', [0.89, 0.94, 0.98])\n")
        f.write("cmd.set_color('ligand_green', [0.83, 0.89, 0.72])\n\n")

        # Color definitions for ProLIF types (from prolif_pymol_vis.py)
        f.write("cmd.set_color('hydrogenbond', [0.204, 0.561, 0.655])\n")
        f.write("cmd.set_color('hydrophobic', [0.69, 0.67, 0.72])\n")
        f.write("cmd.set_color('pistacking', [0.812, 0.686, 0.831])\n")
        f.write("cmd.set_color('pication', [0.97, 0.79, 0.49])\n")
        f.write("cmd.set_color('saltbridge', [0.929, 0.759, 0.662])\n")
        f.write("cmd.set_color('halogenbond', [0.251, 0.718, 0.678])\n")
        f.write("cmd.set_color('elem_N', [0,0,1])\n")
        f.write("cmd.set_color('elem_O', [1,0,0])\n")
        f.write("cmd.set_color('elem_S', [1,1,0])\n\n")

        # Receptor — load from medoid PDB, remove ligand chain L
        f.write("# Receptor (from medoid PDB, chain R + chain P)\n")
        f.write("cmd.load('%s', 'receptor')\n" % ci["medoid_path"])
        f.write("cmd.remove('receptor and chain L')\n")
        f.write("cmd.remove('receptor and not chain R and not chain P')\n")
        f.write("cmd.color('protein_gray', 'receptor')\n")
        f.write("cmd.show('cartoon', 'receptor')\n")
        f.write("cmd.set('transparency', 0.3, 'receptor')\n\n")

        # Load ligand for distance targets
        f.write("# Ligand (for dash targets)\n")
        f.write("cmd.load('%s', 'ligand')\n" % ci["medoid_lpath"])
        f.write("cmd.show('sticks', 'ligand')\n\n")

        f.write("# --- Interaction dashes (colored by ProLIF type) ---\n")
        f.write("# PiCation=orange  PiStacking=purple  SaltBridge=tan  HydrogenBond=teal  Hydrophobic=gray\n")
        dash_count = 0
        for res in int_data:
            if not res.get("interacting"):
                continue
            chain = res.get("chain", "A")
            resid = res.get("resid", "")
            resname = res.get("resname", "")
            pairs = res.get("pairs", [])

            # Classify type from residue name
            plip_type = classify_plip_type(resname)
            color_rgb = MOLAND_FIRST.get(plip_type, [0.69, 0.67, 0.72])

            # Use reference color name (lowercase, no _first suffix)
            col_name = plip_type.lower()

            # Use first pair as representative dash (closest contact)
            if pairs:
                pdb_chain = CHAIN_MAP.get(chain, chain)
                prot_atom = pairs[0][0]
                lig_atom = pairs[0][1]
                dist_val = pairs[0][2]
                dash_name = "d_%s_%s%s" % (plip_type, resname, resid)
                f.write("# %s %s%s - %s  %.2fA (%s)\n" % (resname, chain, resid, lig_atom, dist_val, plip_type))
                f.write("cmd.distance('%s', " % dash_name)
                f.write("'(receptor and chain %s and resi %s and name %s)', " % (pdb_chain, resid, prot_atom))
                f.write("'(ligand and name %s)')\n" % lig_atom)
                f.write("cmd.color('%s', '%s')\n" % (col_name, dash_name))
                # Style the dash
                f.write("cmd.hide('labels', '%s')\n" % dash_name)
                f.write("cmd.set('dash_gap', 0.3, '%s')\n" % dash_name)
                f.write("cmd.set('dash_length', 0.3, '%s')\n" % dash_name)
                f.write("cmd.set('dash_radius', 0.04, '%s')\n" % dash_name)

        # Show interacting residue sidechains (color by interaction type, multi-chain)
        f.write("\n# Interacting residue sidechains (colored by type)\n")
        for res in int_data:
            if not res.get("interacting"):
                continue
            chain = res.get("chain", "A")
            resid = res.get("resid", "")
            resname = res.get("resname", "")
            ptype = classify_plip_type(resname)
            cname = ptype.lower()
            pdb_chain = CHAIN_MAP.get(chain, chain)
            f.write("cmd.select('sc_%s_%s', 'receptor and chain %s and resi %s')\n" % (chain, resid, pdb_chain, resid))
            f.write("cmd.show('sticks', 'sc_%s_%s')\n" % (chain, resid))
            f.write("cmd.color('%s', 'sc_%s_%s')\n" % (cname, chain, resid))

        # Final element color override (applied last so N/O always visible)
        f.write("\n# --- Element color override (applied last) ---\n")
        f.write("cmd.color('elem_N', 'all and element N')\n")
        f.write("cmd.color('elem_O', 'all and element O')\n")
        f.write("cmd.color('elem_S', 'all and element S')\n")
        f.write("cmd.color('elem_Cl', 'all and element Cl')\n")
        f.write("cmd.color('elem_F', 'all and element F')\n\n")

        # View — center on interaction residues since ligand is hidden
        resi_a = [r["resid"] for r in int_data if r.get("interacting") and r.get("chain") == "A"]
        resi_b = [r["resid"] for r in int_data if r.get("interacting") and r.get("chain") == "B"]
        view_sel = "+".join(["chain R and resi %s" % r for r in resi_a] + ["chain P and resi %s" % r for r in resi_b])
        if view_sel:
            f.write("cmd.select('contact_center', '%s')\n" % view_sel)
            f.write("cmd.orient('contact_center')\n")
            f.write("cmd.zoom('contact_center', 4)\n")
        else:
            f.write("cmd.orient('ligand')\n")
            f.write("cmd.zoom('ligand', 6)\n")
        f.write("cmd.save('%s')\n" % pse_c)

    print("PyMOL C%d (dash mode)..." % cid, end=" ", flush=True)
    subprocess.run([PYMOL, "-cq", pml_c], capture_output=True, timeout=600)
    print("OK" if os.path.exists(pse_c) else "FAIL")

# Summary
print("\n=== DONE ===")
print("Output: %s" % OUTPUT)
for ci in cluster_info:
    print("  cluster_%d.pse  (dash visualization)" % ci["id"])
