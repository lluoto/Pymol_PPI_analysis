import json, os, re, math

BASE = "/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/interactions_3plus"
LIG_DIR = "/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/3plus_clusters_dash/ligands"
RECEPTOR = "/media/cuixi/data01/lluoto/6vcb_receptor.pdb"
PYMOL = "/media/cuixi/data01/lluoto/PyMOL-2.6.2_appveyor1945-Linux-x86_64-py311/pymol/pymol"
OUTPUT = "/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/3plus_clusters_dash"

# The 5 hot spots (interaction JSON chain:resid)
HOT5 = [("A", "103"), ("A", "106"), ("A", "159"), ("A", "163"), ("B", "14")]

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

def dist3(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

# Scan all interaction JSONs
candidates = []
for fn in sorted(os.listdir(BASE)):
    if not fn.endswith("_interactions.json"):
        continue
    path = os.path.join(BASE, fn)
    with open(path) as f:
        data = json.load(f)
    interactions = data.get("interactions", [])
    
    # Build set of (chain, resid) that are interacting
    active = set()
    for r in interactions:
        if r.get("interacting"):
            active.add((r.get("chain",""), r.get("resid","")))
    
    # Check if all 5 hot spots are present
    all5 = all(h in active for h in HOT5)
    if not all5:
        continue
    
    # Also need the PDB for centroid
    pdb_fn = fn.replace("_interactions.json", ".pdb")
    pdb_path = os.path.join(BASE, pdb_fn)
    if not os.path.exists(pdb_path):
        continue
    
    centroid = parse_centroid(pdb_path)
    if centroid is None:
        continue
    
    m = re.match(r'(\S+)_pic50_([\d.]+)_seed_(\d+)\.pdb', pdb_fn)
    pic50 = float(m.group(2)) if m else 0.0
    seed = m.group(3) if m else "?"
    ligand = m.group(1) if m else "?"
    
    candidates.append({
        "file": pdb_fn, "path": pdb_path,
        "centroid": centroid, "pic50": pic50,
        "seed": seed, "ligand": ligand,
        "n_interactions": len(active),
    })

print("=== Structures with ALL 5 hot spots ===")
print("Hot spots:", HOT5)
print("Found: %d candidates" % len(candidates))
print()

# Sort by pic50 (prefer mid-range as representative)
candidates.sort(key=lambda c: c["pic50"])
for c in candidates:
    print("  %s  ligand=%s  seed=%s  pic50=%.3f" % (
        c["file"], c["ligand"], c["seed"], c["pic50"]))

# Pick medoid of the candidates (closest to their mean centroid)
if candidates:
    mean_cx = sum(c["centroid"][0] for c in candidates) / len(candidates)
    mean_cy = sum(c["centroid"][1] for c in candidates) / len(candidates)
    mean_cz = sum(c["centroid"][2] for c in candidates) / len(candidates)
    mean_c = (mean_cx, mean_cy, mean_cz)
    
    best = min(candidates, key=lambda c: dist3(c["centroid"], mean_c))
    
    print("\nBest representative (closest to mean centroid of all-5 structures):")
    print("  File: %s" % best["file"])
    print("  Ligand: %s" % best["ligand"])
    print("  Seed: %s" % best["seed"])
    print("  pIC50: %.3f" % best["pic50"])
    print("  Centroid: (%.2f, %.2f, %.2f)" % best["centroid"])
    
    # Write the pick info
    with open(os.path.join(OUTPUT, "hot5_representative.json"), "w") as f:
        json.dump({
            "hot5_residues": HOT5,
            "total_candidates": len(candidates),
            "representative": {
                "file": best["file"],
                "ligand": best["ligand"],
                "seed": best["seed"],
                "pic50": round(best["pic50"], 4),
                "centroid": [round(c, 2) for c in best["centroid"]],
            }
        }, f, indent=2)
    
    # Generate PML for this representative
    pml_path = os.path.join(OUTPUT, "hot5_representative.pml")
    pse_path = os.path.join(OUTPUT, "hot5_representative.pse")
    
    # Load interaction data
    int_fn = best["file"].replace(".pdb", "_interactions.json")
    int_path = os.path.join(BASE, int_fn)
    with open(int_path) as f:
        int_data = json.load(f).get("interactions", [])
    
    # Classify function
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
    
    MOLAND_FIRST = {
        "HydrogenBond": [0.204, 0.561, 0.655],
        "Hydrophobic":  [0.69, 0.67, 0.72],
        "PiStacking":   [0.812, 0.686, 0.831],
        "PiCation":     [0.97, 0.79, 0.49],
        "SaltBridge":   [0.929, 0.759, 0.662],
        "HalogenBond":  [0.251, 0.718, 0.678],
    }
    
    CHAIN_MAP = {"A": "R", "B": "P"}
    
    with open(pml_path, "w") as f:
        f.write("from pymol import cmd\n")
        f.write("cmd.bg_color('white')\n")
        f.write("cmd.set('antialias', 1)\n")
        f.write("cmd.set('dash_gap', 0.3)\n")
        f.write("cmd.set('dash_length', 0.3)\n")
        f.write("cmd.set('dash_radius', 0.06)\n\n")
        
        # Colors
        f.write("cmd.set_color('ligand_green', [0.83, 0.89, 0.72])\n")
        f.write("cmd.set_color('protein_gray', [0.89, 0.94, 0.98])\n")
        for itype, rgb in MOLAND_FIRST.items():
            f.write("cmd.set_color('%s', %s)\n" % (itype.lower(), rgb))
        f.write("cmd.set_color('elem_N', [0,0,1])\n")
        f.write("cmd.set_color('elem_O', [1,0,0])\n")
        f.write("cmd.set_color('elem_S', [1,1,0])\n\n")
        
        # Receptor
        f.write("cmd.load('%s', 'receptor')\n" % best["path"])
        f.write("cmd.remove('receptor and chain L')\n")
        f.write("cmd.remove('receptor and not chain R and not chain P')\n")
        f.write("cmd.color('protein_gray', 'receptor')\n")
        f.write("cmd.show('cartoon', 'receptor')\n")
        f.write("cmd.set('transparency', 0.3, 'receptor')\n\n")
        
        # Ligand for dash targets
        lig_fn = best["file"].replace(".pdb", "_lig.pdb")
        lig_path = os.path.join(LIG_DIR, lig_fn)
        f.write("cmd.load('%s', 'ligand')\n" % lig_path)
        
        # Dashes — named by contact type
        f.write("\n# --- Interaction dashes (named by type:residue) ---\n")
        dash_idx = 0
        for r in int_data:
            if not r.get("interacting"):
                continue
            chain = r.get("chain", "A")
            resid = r.get("resid", "")
            resname = r.get("resname", "")
            pairs = r.get("pairs", [])
            if not pairs:
                continue
            
            ptype = classify_plip_type(resname)
            pdb_ch = CHAIN_MAP.get(chain, chain)
            pa, la, dv = pairs[0][0], pairs[0][1], pairs[0][2]
            
            # Dash name includes type
            dash_name = "d_%s_%s%s" % (ptype, chain, resid)
            
            f.write("# %s %s%s - %s %.2fA (%s)\n" % (resname, chain, resid, la, dv, ptype))
            f.write("cmd.distance('%s', '(receptor and chain %s and resi %s and name %s)', '(ligand and name %s)')\n" % (
                dash_name, pdb_ch, resid, pa, la))
            f.write("cmd.color('%s', '%s')\n" % (ptype.lower(), dash_name))
            f.write("cmd.hide('labels', '%s')\n" % dash_name)
            f.write("cmd.set('dash_gap', 0.3, '%s')\n" % dash_name)
            f.write("cmd.set('dash_length', 0.3, '%s')\n" % dash_name)
            f.write("cmd.set('dash_radius', 0.04, '%s')\n" % dash_name)
            
            # Also show sidechain
            f.write("cmd.select('sc_%s%s', 'receptor and chain %s and resi %s')\n" % (chain, resid, pdb_ch, resid))
            f.write("cmd.show('sticks', 'sc_%s%s')\n" % (chain, resid))
            f.write("cmd.color('%s', 'sc_%s%s')\n" % (ptype.lower(), chain, resid))
            dash_idx += 1
        
        # Element override
        f.write("\ncmd.color('elem_N', 'all and element N')\n")
        f.write("cmd.color('elem_O', 'all and element O')\n\n")
        
        # View
        hot5_resis_a = "+".join(["chain R and resi %s" % r for _, r in HOT5 if _ == "A"])
        hot5_resis_b = "+".join(["chain P and resi %s" % r for _, r in HOT5 if _ == "B"])
        f.write("cmd.select('hot5', '%s or %s')\n" % (hot5_resis_a, hot5_resis_b))
        f.write("cmd.orient('hot5')\n")
        f.write("cmd.zoom('hot5', 4)\n\n")
        
        f.write("cmd.save('%s')\n" % pse_path)
        f.write("print('Saved: %s')\n" % pse_path)
    
    print("\nRunning PyMOL...", end=" ", flush=True)
    import subprocess
    r = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, timeout=600)
    if os.path.exists(pse_path):
        print("OK -> hot5_representative.pse")
    else:
        print("FAIL")

print("\nDone.")
