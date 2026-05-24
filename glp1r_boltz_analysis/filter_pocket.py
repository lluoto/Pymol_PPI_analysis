"""Filter boltz_results_input CIFs by pocket distance (30A from chain B resi 10 CA)."""
import gemmi
import os, sys, json, math
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

PREDICTIONS = "/media/cuixi/data01/lluoto/boltz_results_input/predictions"
OUT_JSON = "/media/cuixi/data01/lluoto/boltz_results_input/pocket_pass.json"
CUTOFF = 30.0  # Angstrom

# Chain B residue 10 CA - the reference point
REF_CHAIN = "B"
REF_RESI = "10"
REF_ATOM = "CA"


def check_one_cif(subdir_path):
    """Check if any ligand atom in this CIF is within CUTOFF of the reference point.
    Returns (subdir_name, passes: bool, ligand_chain: str, min_dist: float)."""
    subdir_name = os.path.basename(subdir_path)
    
    # Find the CIF file
    cif_files = [f for f in os.listdir(subdir_path) if f.endswith(".cif")]
    if not cif_files:
        return (subdir_name, False, None, None, "no_cif")
    
    cif_path = os.path.join(subdir_path, cif_files[0])
    
    try:
        doc = gemmi.cif.read_file(cif_path)
    except Exception as e:
        return (subdir_name, False, None, None, f"parse_error: {e}")
    
    block = doc.sole_block()
    
    try:
        loop = block.find_mmcif_category("_atom_site")
    except Exception:
        return (subdir_name, False, None, None, "no_atom_site")
    
    if len(loop) == 0:
        return (subdir_name, False, None, None, "empty_atom_site")
    
    # Find reference atom (chain B, residue 10, CA)
    ref_pos = None
    ligand_atoms = []  # (element, x, y, z, chain)
    
    for row in loop:
        auth_chain = row[15]
        auth_resi = row[7]
        atom_name = row[3]
        comp_name = row[5]
        group = row[0]
        
        # Reference point
        if auth_chain == REF_CHAIN and auth_resi == REF_RESI and atom_name == REF_ATOM:
            try:
                ref_pos = (float(row[10]), float(row[11]), float(row[12]))
            except ValueError:
                pass
        
        # Ligand atoms (non-polymer HETATM)
        if group == "HETATM" and comp_name != "HOH" and comp_name != "WAT":
            try:
                x, y, z = float(row[10]), float(row[11]), float(row[12])
                ligand_atoms.append((row[2], x, y, z, auth_chain))
            except ValueError:
                pass
    
    if ref_pos is None:
        return (subdir_name, False, None, None, "ref_not_found")
    
    if not ligand_atoms:
        return (subdir_name, False, None, None, "no_ligand")
    
    # Calculate minimum distance
    rx, ry, rz = ref_pos
    min_dist = float("inf")
    min_chain = None
    
    for elem, lx, ly, lz, chain in ligand_atoms:
        dx = lx - rx
        dy = ly - ry
        dz = lz - rz
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        if dist < min_dist:
            min_dist = dist
            min_chain = chain
    
    passes = min_dist <= CUTOFF
    return (subdir_name, passes, min_chain, round(min_dist, 2), None)


def main():
    # Get all subdirectories
    entries = sorted([
        os.path.join(PREDICTIONS, d)
        for d in os.listdir(PREDICTIONS)
        if os.path.isdir(os.path.join(PREDICTIONS, d))
    ])
    
    print(f"Found {len(entries)} prediction subdirectories")
    
    results = {}
    
    # Process sequentially (gemmi is fast enough, 3000 * ~0.1s = 5 min)
    for i, entry in enumerate(entries):
        name, passes, lig_chain, min_dist, error = check_one_cif(entry)
        results[name] = {
            "passes": passes,
            "ligand_chain": lig_chain,
            "min_dist_A": min_dist,
            "error": error
        }
        if (i+1) % 100 == 0:
            passed = sum(1 for v in results.values() if v.get("passes"))
            print(f"  [{i+1}/{len(entries)}] passed so far: {passed}")
    
    # Save results
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    
    passed_count = sum(1 for v in results.values() if v.get("passes"))
    failed_count = sum(1 for v in results.values() if not v.get("passes"))
    error_count = sum(1 for v in results.values() if v.get("error"))
    
    print(f"\nResults saved to {OUT_JSON}")
    print(f"PASSED: {passed_count} / {len(entries)}")
    print(f"FAILED: {failed_count} (errors: {error_count})")
    
    # Show some passing examples
    pass_names = sorted([k for k, v in results.items() if v.get("passes")])
    if pass_names:
        print(f"\nFirst 10 passing:")
        for name in pass_names[:10]:
            r = results[name]
            print(f"  {name} - chain {r['ligand_chain']}, dist={r['min_dist_A']}A")


if __name__ == "__main__":
    main()
