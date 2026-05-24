#!/usr/bin/env python3
"""
PLIP-based interaction detection for prolif_pymol_vis.py.

Called as subprocess from PyMOL. Analyzes a PDB file for:
  - Protein-ligand interactions (using PLIP natively)
  - Protein-protein interactions (using OpenBabel-based detection via PLIP's mol)

Usage:
  python3 _plip_detect.py <pdb_file> <receptor_chains> <partner_chains> [--mode ligand|ppi] [--cutoff 5.0]

Output: JSON to stdout with interaction data.
"""
import json, sys, os, math
from collections import defaultdict

# --- OpenBabel atom utility ---
def _atom_name(obatom):
    """Get PDB atom name from an OpenBabel atom."""
    res = obatom.GetResidue()
    if res:
        return res.GetAtomID(obatom).strip()
    return ''

def _residue_info(obatom):
    """Get (chain, resname, resnum, atom_name) from an OpenBabel atom."""
    res = obatom.GetResidue()
    if not res:
        return ('?', '?', '?', '?')
    return (res.GetChain(), res.GetName().strip(), res.GetNum(), res.GetAtomID(obatom).strip())

def _atom_coords(obatom):
    """Get (x, y, z) from an OpenBabel atom."""
    return (obatom.GetX(), obatom.GetY(), obatom.GetZ())

def _dist3(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

# --- Residue classification (matching prolif_pymol_vis.py) ---
def classify_plip_type(resname):
    name = resname.upper().strip()
    if name in ("LYS", "ARG", "HIS", "HIP"):
        return "PiCation"
    elif name in ("ASP", "GLU", "CYM"):
        return "SaltBridge"
    elif name in ("TYR", "TRP", "PHE"):
        return "PiStacking"
    elif name in ("SER", "THR", "ASN", "GLN", "CYS", "MET"):
        return "HydrogenBond"
    else:
        return "Hydrophobic"

def is_hydrophobic(resname):
    return resname.upper().strip() in ("ALA", "VAL", "LEU", "ILE", "PHE", "PRO", "MET")

def is_aromatic(resname):
    return resname.upper().strip() in ("PHE", "TYR", "TRP")

def is_charged_pos(resname):
    return resname.upper().strip() in ("LYS", "ARG", "HIS", "HIP")

def is_charged_neg(resname):
    return resname.upper().strip() in ("GLU", "ASP", "CYM")

def is_hbond_donor(resname):
    return resname.upper().strip() in ("SER", "THR", "ASN", "GLN", "CYS", "MET", "LYS", "ARG", "HIS", "TRP", "TYR")

def is_hbond_acceptor(resname):
    return resname.upper().strip() in ("SER", "THR", "ASN", "GLN", "CYS", "ASP", "GLU", "HIS", "TYR")

# --- Atom-level classification ---
def _is_aromatic_atom(atom_name):
    return atom_name.upper() in ('CG','CD1','CD2','CE1','CE2','CE3','CZ','CZ2','CZ3','CH2')

def _is_hbond_donor_atom(atom_name, resname):
    rn = resname.upper().strip()
    an = atom_name.upper()
    # Standard donor atoms
    if rn == "ARG" and an in ("NE", "NH1", "NH2"): return True
    if rn == "ASN" and an in ("ND2",): return True
    if rn == "GLN" and an in ("NE2",): return True
    if rn == "HIS" and an in ("ND1", "NE2"): return True
    if rn == "LYS" and an == "NZ": return True
    if rn == "SER" and an == "OG": return True
    if rn == "THR" and an == "OG1": return True
    if rn == "TRP" and an == "NE1": return True
    if rn == "TYR" and an == "OH": return True
    if rn == "CYS" and an == "SG": return True
    if rn == "MET" and an == "SD": return True
    return False

def _is_hbond_acceptor_atom(atom_name, resname):
    rn = resname.upper().strip()
    an = atom_name.upper()
    if rn in ("ASP",) and an in ("OD1", "OD2"): return True
    if rn in ("GLU",) and an in ("OE1", "OE2"): return True
    if rn in ("ASN",) and an in ("OD1",): return True
    if rn in ("GLN",) and an in ("OE1",): return True
    if rn in ("HIS",) and an in ("ND1", "NE2"): return True
    if rn in ("SER",) and an == "OG": return True
    if rn in ("THR",) and an == "OG1": return True
    if rn in ("TYR",) and an == "OH": return True
    if rn in ("CYS",) and an == "SG": return True
    return False

# --- Main detection ---

def detect_ppi(mol, receptor_chains, partner_chains, cutoff=5.0):
    """
    Detect PPI between receptor and partner chains using OpenBabel molecule.
    
    Args:
        mol: plip.structure.preparation.PDBComplex (loaded)
        receptor_chains: list of chain IDs (e.g. ['R'])
        partner_chains: list of chain IDs (e.g. ['P'])
        cutoff: max distance in Angstrom
    
    Returns:
        list of interaction dicts
    """
    rc_set = set(receptor_chains)
    pc_set = set(partner_chains)
    
    # Collect atoms by chain
    rec_atoms = []  # list of (obatom, chain, resname, resnum, atom_name, coords)
    par_atoms = []
    
    for idx in sorted(mol.atoms.keys()):
        atom = mol.get_atom(idx)
        if atom is None:
            continue
        oba = atom.OBAtom
        res = oba.GetResidue()
        if res is None:
            continue
        chain = res.GetChain()
        resname = res.GetName().strip()
        resnum = res.GetNum()
        atom_name = res.GetAtomID(oba).strip()
        coords = (oba.GetX(), oba.GetY(), oba.GetZ())
        entry = (oba, chain, resname, resnum, atom_name, coords)
        
        if chain in rc_set:
            rec_atoms.append(entry)
        elif chain in pc_set:
            par_atoms.append(entry)
    
    results = []
    
    # For each pair of atoms (receptor atom, partner atom), check interactions
    for ra in rec_atoms:
        ra_oba, ra_chain, ra_resname, ra_resnum, ra_aname, ra_xyz = ra
        for pa in par_atoms:
            pa_oba, pa_chain, pa_resname, pa_resnum, pa_aname, pa_xyz = pa
            
            d = _dist3(ra_xyz, pa_xyz)
            if d > cutoff:
                continue
            
            # --- Classify interaction type ---
            rn_r = ra_resname.upper().strip()
            rn_p = pa_resname.upper().strip()
            an_r = ra_aname.upper()
            an_p = pa_aname.upper()
            
            iname = None
            
            # 1. Salt bridge: charged positive + charged negative
            if (is_charged_pos(rn_r) and is_charged_neg(rn_p)) or \
               (is_charged_neg(rn_r) and is_charged_pos(rn_p)):
                if d < 4.0:
                    iname = "SaltBridge"
            
            # 2. Pi-cation: aromatic ring + cation
            if not iname:
                if d < 4.5:
                    if (is_aromatic(rn_r) and _is_aromatic_atom(an_r) and is_charged_pos(rn_p)) or \
                       (is_aromatic(rn_p) and _is_aromatic_atom(an_p) and is_charged_pos(rn_r)):
                        iname = "PiCation"
            
            # 3. Pi-stacking: aromatic + aromatic
            if not iname:
                if d < 5.0:
                    if is_aromatic(rn_r) and is_aromatic(rn_p):
                        if _is_aromatic_atom(an_r) and _is_aromatic_atom(an_p):
                            iname = "PiStacking"
            
            # 4. Hydrogen bond: donor + acceptor
            if not iname:
                if d < 3.5:
                    if (_is_hbond_donor_atom(an_r, rn_r) and _is_hbond_acceptor_atom(an_p, rn_p)) or \
                       (_is_hbond_donor_atom(an_p, rn_p) and _is_hbond_acceptor_atom(an_r, rn_r)):
                        iname = "HydrogenBond"
            
            # 5. Hydrophobic: both hydrophobic side chain atoms
            if not iname:
                if d < 4.0:
                    if is_hydrophobic(rn_r) and is_hydrophobic(rn_p):
                        # Only count hydrophobic atoms (carbon, not polar)
                        if ra_oba.GetAtomicNum() == 6 and pa_oba.GetAtomicNum() == 6:
                            iname = "Hydrophobic"
            
            if iname:
                # Determine which is receptor-side and partner-side
                # Receptor and partner could be either atom in the pair
                # We maintain the convention: first atom = receptor chain atom
                if ra_chain in rc_set:
                    results.append({
                        "type": iname,
                        "receptor_chain": ra_chain,
                        "receptor_resname": ra_resname,
                        "receptor_resnum": str(ra_resnum),
                        "receptor_atom": ra_aname,
                        "partner_chain": pa_chain,
                        "partner_resname": pa_resname,
                        "partner_resnum": str(pa_resnum),
                        "partner_atom": pa_aname,
                        "distance": round(d, 2)
                    })
                else:
                    # Swap
                    results.append({
                        "type": iname,
                        "receptor_chain": pa_chain,
                        "receptor_resname": pa_resname,
                        "receptor_resnum": str(pa_resnum),
                        "receptor_atom": pa_aname,
                        "partner_chain": ra_chain,
                        "partner_resname": ra_resname,
                        "partner_resnum": str(ra_resnum),
                        "partner_atom": ra_aname,
                        "distance": round(d, 2)
                    })
    
    return results


def detect_ligand(mol, cutoff=5.0):
    """
    Use PLIP native analysis for protein-ligand interactions.
    Returns interactions between the auto-detected ligand and protein chains.
    """
    try:
        mol.analyze()
    except Exception as e:
        return {"error": str(e)}
    
    results = []
    for lig_id, iset in mol.interaction_sets.items():
        lig_chain = iset.ligand.chain
        lig_hetid = iset.ligand.hetid
        
        for attr_name, output_type in [
            ('hydrophobic_contacts', 'Hydrophobic'),
            ('hbonds_ldon', 'HydrogenBond'),
            ('hbonds_pdon', 'HydrogenBond'),
            ('pistacking', 'PiStacking'),
            ('pication_laro', 'PiCation'),
            ('pication_paro', 'PiCation'),
            ('saltbridge_lneg', 'SaltBridge'),
            ('saltbridge_pneg', 'SaltBridge'),
            ('halogen_bonds', 'HalogenBond'),
        ]:
            items = getattr(iset, attr_name, [])
            for item in items:
                try:
                    # Extract atom info from the interaction object
                    bsa = item.bsatom.OBAtom  # receptor side
                    lga = item.ligatom.OBAtom  # ligand side
                    
                    bsa_res = bsa.GetResidue()
                    
                    entry = {
                        "type": output_type,
                        "receptor_chain": getattr(item, 'reschain', bsa_res.GetChain()),
                        "receptor_resname": getattr(item, 'restype', bsa_res.GetName().strip()),
                        "receptor_resnum": str(getattr(item, 'resnr', bsa_res.GetNum())),
                        "receptor_atom": bsa_res.GetAtomID(bsa).strip(),
                        "partner_chain": lig_chain,
                        "partner_resname": getattr(item, 'restype_l', 'LIG'),
                        "partner_resnum": str(getattr(item, 'resnr_l', '1')),
                        "partner_atom": item.ligatom_orig_idx,
                        "distance": round(float(item.distance), 2),
                        "ligand_mode": True,
                    }
                    results.append(entry)
                except Exception:
                    continue
    
    return results


# --- Main entry point ---
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='PLIP-based PPI/PL detection')
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('receptor_chains', help='Receptor chain IDs (comma-separated, e.g. R,A,B)')
    parser.add_argument('partner_chains', help='Partner/ligand chain IDs (comma-separated, e.g. P,L)')
    parser.add_argument('--mode', choices=['auto', 'ppi', 'ligand'], default='auto',
                        help='Detection mode: ppi=protein-protein, ligand=PLIP native, auto=try both')
    parser.add_argument('--cutoff', type=float, default=5.0, help='Distance cutoff')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.pdb_file):
        print(json.dumps({"error": "PDB file not found: %s" % args.pdb_file}))
        sys.exit(1)
    
    receptor_list = [c.strip() for c in args.receptor_chains.split(',') if c.strip()]
    partner_list = [c.strip() for c in args.partner_chains.split(',') if c.strip()]
    
    # Load structure via PLIP's PDBComplex
    from plip.structure.preparation import PDBComplex
    
    mol = PDBComplex()
    try:
        mol.load_pdb(args.pdb_file)
    except Exception as e:
        print(json.dumps({"error": "Failed to load PDB: %s" % str(e)}))
        sys.exit(1)
    
    output = {
        "pdb": args.pdb_file,
        "receptor_chains": receptor_list,
        "partner_chains": partner_list,
        "mode": args.mode,
        "cutoff": args.cutoff,
    }
    
    results = []
    
    if args.mode in ('auto', 'ppi'):
        ppi_results = detect_ppi(mol, receptor_list, partner_list, args.cutoff)
        results.extend(ppi_results)
    
    if args.mode in ('auto', 'ligand'):
        lig_results = detect_ligand(mol, args.cutoff)
        # For "auto" mode, only include ligand results if no PPI-specific chains found
        if args.mode == 'auto':
            # Check if the auto-detected ligand has at least one interaction
            # and our PPI detection found nothing (possibly because both sides
            # are the same protein complex)
            if not ppi_results:
                results = lig_results
                output['active_mode'] = 'ligand'
            else:
                output['active_mode'] = 'ppi'
        else:
            results = lig_results if args.mode == 'ligand' else results
            output['active_mode'] = args.mode
    
    # Deduplicate: same residues + same type = keep only shortest distance
    seen = {}
    deduped = []
    for r in results:
        key = (r['type'], r['receptor_resnum'], r['receptor_chain'], 
               r['partner_resnum'], r.get('partner_chain', ''))
        old = seen.get(key)
        if old is None or r['distance'] < old['distance']:
            seen[key] = r
    
    output['interactions'] = list(seen.values())
    output['n_interactions'] = len(output['interactions'])
    
    print(json.dumps(output, indent=2))
