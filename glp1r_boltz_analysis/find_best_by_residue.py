"""Find best-affinity structure for each residue with interactions."""
import json, os, glob

PLIP_JSON = "/media/cuixi/data01/lluoto/boltz_results_input/plip_results.json"
PREDICTIONS = "/media/cuixi/data01/lluoto/boltz_results_input/predictions"

with open(PLIP_JSON) as f:
    plip_data = json.load(f)

# Query: structures where each target residue has interactions
# Target: chain A resi 170 (full seq 202, K)
target_residues = [
    ("A", 170, "GLP-1R 202 Lys"),
]

results = {}
for ch, rn, label in target_residues:
    # Find all structures with this residue interacting
    candidates = []
    for name, info in plip_data.items():
        if "details" not in info: continue
        for det in info["details"]:
            if det.get("chain") == ch and det.get("res_num") == rn:
                # Read affinity JSON
                aff_path = os.path.join(PREDICTIONS, name, f"affinity_{name}.json")
                affinity_val = None
                ligand_type = None
                
                # Determine which affinity_pred_value to use based on structure name
                if "3_20" in name:
                    ligand_type = "3_20"
                    key = "affinity_pred_value"
                elif "3_24" in name:
                    ligand_type = "3_24"
                    key = "affinity_pred_value1"
                elif "3_29" in name:
                    ligand_type = "3_29"
                    key = "affinity_pred_value2"
                else:
                    continue
                
                if os.path.exists(aff_path):
                    with open(aff_path) as af:
                        aff_data = json.load(af)
                    affinity_val = aff_data.get(key)
                
                # Also get interaction count for tie-breaking
                total_int = sum(info.get("interactions", {}).values())
                
                candidates.append((name, affinity_val, ligand_type, total_int))
                break  # One entry per structure
    
    # Sort by affinity (lower = better), then by interaction count (higher = better)
    valid_candidates = [c for c in candidates if c[1] is not None]
    valid_candidates.sort(key=lambda x: (x[1], -x[3]))
    
    print(f"\n=== {ch}{rn} ({label}) ===")
    print(f"Total structures with interactions: {len(candidates)}")
    print(f"With affinity data: {len(valid_candidates)}")
    
    if valid_candidates:
        print(f"\nTop 10 by affinity (lower=better):")
        for name, aff, lig, n_int in valid_candidates[:10]:
            print(f"  {name}: affinity_pred={aff:.4f} ({lig}) interactions={n_int}")
        
        best = valid_candidates[0]
        print(f"\n=== BEST: {best[0]} ===")
        print(f"  affinity_pred={best[1]:.4f} ({best[2]}) interactions={best[3]}")
        results[(ch, rn)] = best
    else:
        print("  NO STRUCTURES FOUND")
        results[(ch, rn)] = None

# Also check A110 and A113 with broader approach
print("\n\n=== CHECKING A110 and A113 with ANY interaction ===")
for ch, rn, label in [("A", 110, "GLP-1R 142"), ("A", 113, "GLP-1R 145")]:
    count = 0
    for name, info in plip_data.items():
        if "details" not in info: continue
        for det in info["details"]:
            if det.get("chain") == ch and det.get("res_num") == rn:
                count += 1
                break
    print(f"  {ch}{rn} ({label}): {count} structures with interactions")

# Generate PSE command for the best one
if ("A", 170) in results and results[("A", 170)]:
    best_name = results[("A", 170)][0]
    print(f"\n\n=== PSE for best structure ===")
    print(f"Structure: {best_name}")
    
    # Generate PyMOL script
    project_dir = f"/media/cuixi/data01/lluoto/boltz_results_input/pse_best_by_residue/A170_{best_name}"
    os.makedirs(project_dir, exist_ok=True)
    pse_path = f"{project_dir}/A170_{best_name}.pse"
    pml_path = f"{project_dir}/A170_{best_name}.pml"
    
    info = plip_data[best_name]
    details = info.get("details", [])
    
    from collections import defaultdict
    type_residues = defaultdict(set)
    type_dashes = defaultdict(list)
    
    for det in details:
        itype = det.get("type", "")
        rnum = det.get("res_num")
        p_idx = det.get("protein_atom_idx")
        l_idx = det.get("ligand_atom_idx")
        if rnum: type_residues[itype].add(rnum)
        if p_idx and l_idx and int(p_idx) > 0 and int(l_idx) > 0:
            type_dashes[itype].append((int(p_idx), int(l_idx)))
    
    cif_dir = os.path.join(PREDICTIONS, best_name)
    cif_files = glob.glob(os.path.join(cif_dir, "*.cif"))
    
    if cif_files:
        with open(pml_path, "w") as f:
            f.write("from pymol import cmd\nimport pymol\n\n")
            f.write(f"cmd.load('{cif_files[0]}')\n\n")
            
            # Color defs (prolif scheme)
            f.write("cmd.set_color('protein_gray', [0.89, 0.94, 0.98])\n")
            f.write("cmd.set_color('ligand_green', [0.83, 0.89, 0.72])\n")
            for itype, rgb in [("hbond", [0.204, 0.561, 0.655]), ("hydrophobic", [0.69, 0.67, 0.72]), ("pistacking", [0.812, 0.686, 0.831]), ("pication", [0.97, 0.79, 0.49]), ("saltbridge", [0.929, 0.759, 0.662]), ("halogen", [0.251, 0.718, 0.678])]:
                f.write(f"cmd.set_color('{itype}_dash', {rgb})\n")
                f.write(f"cmd.set_color('{itype}_stick', {rgb})\n")
            for elem, rgb in [("O", [1,0,0]), ("N", [0,0,1]), ("S", [1,1,0]), ("Cl", [0,1,0]), ("F", [0,1,1])]:
                f.write(f"cmd.set_color('elem_{elem}', {rgb})\n")
            f.write("cmd.bg_color('white')\n\n")
            
            # Base
            f.write("cmd.hide('everything')\n")
            f.write("cmd.show('cartoon', 'polymer.protein')\n")
            f.write("cmd.color('protein_gray', 'polymer.protein')\n")
            f.write("cmd.show('sticks', 'hetatm')\n")
            f.write("cmd.color('ligand_green', 'hetatm')\n\n")
            
            # Highlight target residue A170
            f.write("# Highlight target residue A170 (GLP-1R 202 Lys)\n")
            f.write("cmd.select('target_res', 'chain A and resi 170')\n")
            f.write("cmd.show('sticks', 'target_res')\n")
            f.write("cmd.color('deeppurple', 'target_res')\n")
            f.write("cmd.show('spheres', 'target_res and name CA')\n\n")
            
            # Interacting residues
            for itype, resns in type_residues.items():
                color = f"{itype}_stick"
                res_str = " or ".join([f"resi {r}" for r in sorted(resns)])
                if res_str:
                    sel = f"int_{itype}"
                    f.write(f"cmd.select('{sel}', 'polymer and ({res_str})')\n")
                    f.write(f"cmd.show('sticks', '{sel}')\n")
                    f.write(f"cmd.color('{color}', '{sel}')\n")
            
            # Element colors
            for elem in ["O", "N", "S", "Cl", "F"]:
                f.write(f"cmd.color('elem_{elem}', 'elem {elem}')\n")
            
            # Dashes
            for itype, dash_list in type_dashes.items():
                color = f"{itype}_dash"
                for i, (p, l) in enumerate(dash_list):
                    dn = f"dash_{itype}_{i}"
                    f.write(f"cmd.distance('{dn}', 'index {p}', 'index {l}')\n")
                    f.write(f"cmd.color('{color}', '{dn}')\n")
                    f.write(f"cmd.hide('labels', '{dn}')\n")
            
            f.write("\ncmd.orient()\n")
            f.write(f"cmd.save('{pse_path}')\n")
            f.write("cmd.delete('all')\n")
            f.write("print('Done')\n")
        
        print(f"PML: {pml_path}")
        print(f"PSE: {pse_path}")
        print(f"\nRun: /media/cuixi/data01/lluoto/PyMOL-2.6.2_appveyor1945-Linux-x86_64-py311/pymol/pymol -cq {pml_path}")
    else:
        print("CIF not found!")
