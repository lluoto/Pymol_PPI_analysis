#!/usr/bin/env python
# coding: utf-8
'''
    Author: Chenxi Wang (chenxi.wang@salilab.org) Yulin Luo
    Date: 2024-6-21
    Updated: 2025 - Added all ProLIF interaction types with moland color system
    
    The residue contacts in the GLP-1R-Gs interface is calculated by Pymol.
    Now supports ProLIF interaction types: HydrogenBond, Hydrophobic, PiStacking, 
    EdgeToFace, FaceToFace, CationPi, PiCation, SaltBridge, HalogenBond, Metal, VdWContact
'''

from pymol import cmd
from pymol import stored
import time
import csv
import os
import re
import ast
import subprocess
import json
import tempfile

# --- PLIP integration constants ---
# PLIP detection script path (system Python, not PyMOL's Python)
PLIP_SCRIPT = '/tmp/_plip_detect.py'
PLIP_PYTHON = '/usr/bin/python3'
PLIP_TEMP_PDB = '/tmp/prolif_plip_temp.pdb'

cmd.set_color('protein_gray', [0.89, 0.94, 0.98])
cmd.set_color('ligand_green', [0.83, 0.89, 0.72])
cmd.set_color('hydrophobic_gray', [0.69, 0.67, 0.72])
cmd.set_color('hydrogen_bond_blue', [0.204, 0.561, 0.655])
cmd.set_color('pi_cation_orange', [0.97, 0.79, 0.49])
cmd.set_color('pi_stack_purple', [0.812, 0.686, 0.831])
cmd.set_color('salt_bridge_red', [0.929, 0.759, 0.662])
cmd.set_color('H_color', [0.7, 0.7, 0.7])
cmd.set_color('oxygen_color', [1.0, 0.0, 0.0])
cmd.set_color('nitrogen_color', [0.0, 0.0, 1.0])
cmd.set_color('fluorine_color', [0.251, 0.718, 0.678])
cmd.set_color('sulfur_color', [1.0, 1.0, 0.0])

cmd.set('light_count',10)
cmd.set('spec_count',1)
cmd.set('shininess',10)
cmd.set('specular',0.4)
cmd.set('ambient',0)
cmd.set('direct',0.45)
cmd.set('reflect',0.75)
cmd.set('ray_shadow_decay_factor',0.1)
cmd.set('ray_shadow_decay_range',2)
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_fancy_sheets', 1)
cmd.set('cartoon_rect_width', 0.3)
cmd.set('cartoon_rect_length', 2)
cmd.set('cartoon_loop_radius', 0.15)
cmd.set('cartoon_rect_width', 0.3)
cmd.set('specular', 0.6)

cmd.bg_color('white')

start=time.time()

amino_acid_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
_AA_ONE_TO_THREE = {v: k for k, v in amino_acid_map.items()}

hbonds_atoms = {
    'ARG': ['NE', 'NH1', 'NH2'],
    'ASN': ['OD1', 'ND2'],
    'ASP': ['OD1', 'OD2'],
    'CYS': ['SG'],
    'GLN': ['OE1', 'NE2'],
    'GLU': ['OE1', 'OE2'],
    'HIS': ['ND1', 'NE2'],
    'LYS': ['NZ'],
    'MET': ['SD'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TRP': ['NE1'],
    'TYR': ['OH']
}

MOLAND_COLORS = {
    'HydrogenBond': {
        'first': [0.204, 0.561, 0.655],
        'second': [0.404, 0.761, 0.855],
    },
    'Hydrophobic': {
        'first': [0.69, 0.67, 0.72],
        'second': [0.89, 0.87, 0.92],
    },
    'PiStacking': {
        'first': [0.812, 0.686, 0.831],
        'second': [0.906, 0.843, 0.915],
    },
    'EdgeToFace': {
        'first': [0.6, 0.4, 0.8],
        'second': [0.8, 0.6, 0.9],
    },
    'FaceToFace': {
        'first': [0.7, 0.5, 0.9],
        'second': [0.85, 0.75, 0.95],
    },
    'CationPi': {
        'first': [0.97, 0.79, 0.49],
        'second': [0.985, 0.895, 0.745],
    },
    'PiCation': {
        'first': [0.97, 0.79, 0.49],
        'second': [0.985, 0.895, 0.745],
    },
    'SaltBridge': {
        'first': [0.929, 0.759, 0.662],
        'second': [0.965, 0.88, 0.831],
    },
    'HalogenBond': {
        'first': [0.251, 0.718, 0.678],
        'second': [0.525, 0.859, 0.839],
    },
    'Metal': {
        'first': [1.0, 1.0, 0.0],
        'second': [1.0, 1.0, 0.5],
    },
    'VdWContact': {
        'first': [0.89, 0.94, 0.98],
        'second': [0.945, 0.97, 0.99],
    },
    'WaterBridge': {
        'first': [0.631, 0.769, 0.992],
        'second': [0.815, 0.884, 0.996],
    },
}

for inter_type, colors in MOLAND_COLORS.items():
    color_name = f"{inter_type.lower()}_first"
    cmd.set_color(color_name, colors['first'])
    color_name_second = f"{inter_type.lower()}_second"
    cmd.set_color(color_name_second, colors['second'])

PROLIF_INTERACTION_COLORS = {
    'HBDonor': 'hydrogen_bond_first',
    'HBAcceptor': 'hydrogen_bond_first',
    'Hydrophobic': 'hydrophobic_first',
    'PiStacking': 'pi_stack_first',
    'FaceToFace': 'pi_stack_first',
    'EdgeToFace': 'pistacking_first',
    'CationPi': 'pi_cation_first',
    'PiCation': 'pi_cation_first',
    'Cationic': 'salt_bridge_first',
    'Anionic': 'salt_bridge_first',
    'VdWContact': 'vdwcontact_first',
    'XBAcceptor': 'halogenbond_first',
    'XBDonor': 'halogenbond_first',
    'MetalDonor': 'metal_first',
    'MetalAcceptor': 'metal_first',
    'WaterBridge': 'waterbridge_first',
}

def list_electron_interaction(selection, channel, partners, selection2=None, cutoff=4.0, angle=150, mode=1, hb_list_name='hbonds'):
    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    posresList = ["LYS","ARG","HIS","HIP"]
    negresList = ["GLU","ASP","CYM"]
    posatomList = ["NE","NH1","NH2","NZ","ND1","NE2"]
    negatomList = ["SG","OE1","OE2","OD1","OD2"]
    
    selection = selection + " & e. n+o+s"
    
    sb = cmd.find_pairs(selection, selection, mode=mode, cutoff=cutoff, angle=angle)
    sb.sort(key=lambda x:x[0][1])
    
    result_h = []
    result_s = []
    
    for pairs in sb:
        stored.x = []
        stored.y = []
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
        
        if stored.x[0] in channel and stored.y[0] in partners:
            cmd.show('stick', "chain %s and resi %s chain %s and resi %s" % (stored.x[0], stored.x[1], stored.y[0], stored.y[1]))
            cmd.select('should_label', "chain %s and resi %s or should_label" % (stored.x[0], stored.x[1]))
            cmd.select('should_label', "chain %s and resi %s or should_label" % (stored.y[0], stored.y[1]))
            
            if (stored.x[2] in posresList and stored.x[3] in posatomList and stored.y[2] in negresList and stored.y[3] in negatomList):
                d = cmd.distance('saltbridge', "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]), cutoff=4.0)
                result_s.append([stored.x, stored.y, float(d)])
                cmd.color('salt_bridge_red', 'saltbridge')
            elif stored.x[2] in hbonds_atoms.keys() and stored.x[2] != 'MET' and stored.y[2] in hbonds_atoms.keys() and stored.y[3] in hbonds_atoms[stored.y[2]]:
                d = cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]), cutoff=3.5)
                result_h.append([stored.x, stored.y, float(d)])
                cmd.color('hydrogen_bond_blue', hb_list_name)
    
    return result_s, result_h


def list_hydrophobic(selection, channel, partners, selection2=None, cutoff=4.0, angle=180, mode=1, hb_list_name='hydrophobic'):
    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    hydrophobiclist = ["ALA","VAL","LEU","ILE","PHE","PRO","MET"]
    aromaticlist = ["PHE","TYR","TRP"]
    aromaticatoms = ['CG','CZ','CD1','CD2','CE1','CE2']
    
    selection = selection + " & e. c+n"
    hp = cmd.find_pairs(selection, selection, mode=mode, cutoff=cutoff, angle=angle)
    hp.sort(key=lambda x:x[0][1])
    
    result = []
    for pairs in hp:
        stored.x = []
        stored.y = []
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
        
        if (stored.x[2] in hydrophobiclist and stored.y[2] in hydrophobiclist):
            if stored.x[0] in channel and stored.y[0] in partners:
                cmd.show('stick', "chain %s and resi %s chain %s and resi %s" % (stored.x[0], stored.x[1], stored.y[0], stored.y[1]))
                cmd.select('should_label', "chain %s and resi %s or should_label" % (stored.x[0], stored.x[1]))
                cmd.select('should_label', "chain %s and resi %s or should_label" % (stored.y[0], stored.y[1]))
                
                if stored.x[2] in aromaticlist and stored.y[2] in aromaticlist:
                    d = cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]), cutoff=5.0, mode=6)
                    result.append([stored.x, stored.y, float(d)])
                    cmd.color('pi_stack_purple', hb_list_name)
                elif stored.x[2] in aromaticlist and stored.x[3] in aromaticatoms:
                    d = cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]), cutoff=4.0, mode=7)
                    result.append([stored.x, stored.y, float(d)])
                    cmd.color('pi_cation_orange', hb_list_name)
                elif stored.x[2] not in aromaticlist and stored.y[2] not in aromaticlist:
                    d = cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]))
                    result.append([stored.x, stored.y, float(d)])
                    cmd.color('hydrophobic_gray', hb_list_name)
        else:
            continue

    return result


def select_interface_pair(bonds, filepath, various, allpairs, channel, partners, part1, part2):
    R_A = []
    specific_pair = [(subunit, partner) for subunit in channel for partner in partners]

    for item in bonds:
        pairs = (item[0][0], item[1][0])
        if pairs in specific_pair:
            R_A += [item]

    if R_A != []:
        chains = set([f'{x[0][0]},{y[1][0]}' for x in R_A for y in R_A if x[0][0] != y[1][0] and x[0][0].isdigit() == False and y[1][0].isdigit() == False])
        allpairs[various] = {}
        for chain in chains:
            allpairs[various][chain] = []
        for item in R_A:
            pair = f'{item[0][1]}{amino_acid_map[item[0][2]]},{item[1][1]}{amino_acid_map[item[1][2]]}'
            part1.append(f'{item[0][0]},{item[0][1]}')
            part2.append(f'{item[1][0]},{item[1][1]}')
            allpairs[various][f'{item[0][0]},{item[1][0]}'] += [pair]


def _normalize_resn_code(resn):
    resn = str(resn).strip().upper()
    if not resn:
        return None
    if len(resn) == 1:
        return _AA_ONE_TO_THREE.get(resn)
    if len(resn) == 3:
        return resn
    return None


def _parse_residue_token(token):
    token = str(token).strip().strip('"').strip("'")
    if not token:
        return None, None

    # Try the two patterns that work for: GLU214, 214A, 214GLU
    patterns = [
        r'^(?P<resn>[A-Za-z]{1,3})(?P<resi>-?\d+[A-Za-z]?)$',
        r'^(?P<resi>-?\d+[A-Za-z]?)(?P<resn>[A-Za-z]{1,3})$',
    ]
    for pattern in patterns:
        m = re.match(pattern, token)
        if m:
            resn = m.groupdict().get('resn')
            resi = m.groupdict().get('resi')
            if resn and resi:
                return _normalize_resn_code(resn), resi

    # For ambiguous tokens like "195SER" or "454D":
    # try all splits: keep only digits in resi, everything else is resn
    m = re.match(r'^(-?\d+)([A-Za-z]+)$', token)
    if m:
        return _normalize_resn_code(m.group(2)), m.group(1)

    # If nothing matched but there are digits, use them as resi
    m = re.search(r'-?\d+', token)
    if m:
        return None, m.group(0)

    return None, token


def _closest_atom_point(obj_name, selection_expr, fallback_expr=None):
    coords = cmd.get_model(selection_expr).atom
    if not coords and fallback_expr:
        coords = cmd.get_model(fallback_expr).atom
    if not coords:
        return False

    best = coords[0]
    for atom in coords[1:]:
        if atom.index < best.index:
            best = atom

    cmd.delete(obj_name)
    cmd.pseudoatom(obj_name, pos=list(best.coord))
    cmd.hide('everything', obj_name)
    return True


def draw_prolif_contacts(contact_file, selection="all", cutoff=6.0, show_labels=True):
    if not os.path.exists(contact_file):
        print(f"Contact file {contact_file} not found")
        return

    interactions_drawn = {}
    contact_count = {'total': 0}
    label_sel = show_labels if isinstance(show_labels, str) else ('should_label' if show_labels else None)
    if label_sel:
        cmd.select(label_sel, 'none')

    with open(contact_file, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if len(row) == 1 and (not row[0].strip() or row[0].startswith('#')):
                continue

            try:
                contacts = []
                if len(row) >= 9 and row[0].strip() in ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'):
                    lig_chain, lig_resid_raw = row[0].strip(), row[1].strip()
                    lig_name = row[2].strip() if len(row) > 2 and row[2].strip() else '*'
                    prot_chain, prot_resid_raw = row[3].strip(), row[4].strip()
                    prot_name = row[5].strip() if len(row) > 5 and row[5].strip() else '*'
                    interaction = row[6].strip()
                    color = row[7].strip() if len(row) > 7 and row[7].strip() else 'white'
                    distance = float(row[8]) if len(row) > 8 and row[8].strip() else cutoff
                    contacts.append((lig_chain, lig_resid_raw, lig_name, prot_chain, prot_resid_raw, prot_name, interaction, color, distance))
                elif len(row) == 3:
                    interaction = row[0].strip()
                    chain_pair = [c.strip() for c in row[1].split(',', 1)]
                    if len(chain_pair) != 2:
                        continue
                    pair_set = ast.literal_eval(row[2])
                    if isinstance(pair_set, str):
                        pair_set = {pair_set}
                    for pair in sorted(str(pair).strip() for pair in pair_set if str(pair).strip()):
                        pair_bits = [p.strip() for p in pair.split(',', 1)]
                        if len(pair_bits) != 2:
                            continue
                        contacts.append((chain_pair[0], pair_bits[0], '*', chain_pair[1], pair_bits[1], '*', interaction, 'white', cutoff))
                else:
                    continue

                for lig_chain, lig_resid_raw, lig_name, prot_chain, prot_resid_raw, prot_name, interaction, color, distance in contacts:
                    lig_resn, lig_resi = _parse_residue_token(lig_resid_raw)
                    prot_resn, prot_resi = _parse_residue_token(prot_resid_raw)
                    if lig_resi is None or prot_resi is None:
                        print(f"Skipping line with unparseable residue ids: {row}")
                        continue

                    inter_key = f"{interaction}_{lig_chain}_{lig_resi}_{prot_chain}_{prot_resi}"
                    if inter_key in interactions_drawn:
                        continue

                    base_sel1 = f"chain {lig_chain} and resi {lig_resi}"
                    base_sel2 = f"chain {prot_chain} and resi {prot_resi}"
                    if lig_resn:
                        base_sel1 += f" and resn {lig_resn}"
                    if prot_resn:
                        base_sel2 += f" and resn {prot_resn}"

                    stick_sel1 = f"({base_sel1}) and (sidechain or name CA)"
                    stick_sel2 = f"({base_sel2}) and (sidechain or name CA)"

                    point_sel1 = stick_sel1
                    point_sel2 = stick_sel2
                    if lig_name != '*':
                        point_sel1 = f"({base_sel1}) and name {lig_name}"
                    if prot_name != '*':
                        point_sel2 = f"({base_sel2}) and name {prot_name}"

                    point1 = f"pt1_{inter_key}"
                    point2 = f"pt2_{inter_key}"
                    if not _closest_atom_point(point1, point_sel1, stick_sel1):
                        continue
                    if not _closest_atom_point(point2, point_sel2, stick_sel2):
                        continue

                    dist_name = f"contact_{inter_key}"
                    cmd.distance(dist_name, point1, point2, cutoff=distance + 0.5)

                    pymol_color = PROLIF_INTERACTION_COLORS.get(color, color)
                    if pymol_color in PROLIF_INTERACTION_COLORS:
                        pymol_color = PROLIF_INTERACTION_COLORS[color]
                    if pymol_color in cmd.get_color_indices():
                        cmd.color(pymol_color, dist_name)

                    cmd.show('sticks', stick_sel1)
                    cmd.show('sticks', stick_sel2)

                    interactions_drawn[inter_key] = True
                    contact_count['total'] += 1

                    if label_sel:
                        cmd.select(label_sel, f"{base_sel1} or {label_sel}")
                        cmd.select(label_sel, f"{base_sel2} or {label_sel}")

                    if interaction not in contact_count:
                        contact_count[interaction] = 0
                    contact_count[interaction] += 1

            except Exception as e:
                print(f"Error processing line: {','.join(row)}, error: {e}")
                continue

    print(f"Drew {len(interactions_drawn)} unique contacts from {os.path.basename(contact_file)}")
    return interactions_drawn


def get_interaction_color(interaction_type, use_second=False):
    """Get moland color for interaction type"""
    suffix = '_second' if use_second else '_first'
    return PROLIF_INTERACTION_COLORS.get(interaction_type, 'white')



cmd.extend("draw_prolif_contacts", draw_prolif_contacts)
cmd.extend("get_interaction_color", get_interaction_color)


def _parse_chain_arg(chains):
    if chains is None:
        return None
    if isinstance(chains, str):
        return [c.strip() for c in chains.split(",") if c.strip()]
    return [str(c).strip() for c in chains if str(c).strip()]


def _apply_default_contact_style(receptor_chains, partner_chains):
    receptor_sel = "+".join(receptor_chains)
    partner_sel = "+".join(partner_chains)
    if receptor_sel:
        cmd.color("protein_gray", f"chain {receptor_sel}")
    if partner_sel:
        cmd.color("pi_cation_orange", f"chain {partner_sel}")
    cmd.color("nitrogen_color", "element n")
    cmd.color("oxygen_color", "element o")
    cmd.color("sulfur_color", "element s")
    cmd.color("fluorine_color", "element f")
    cmd.set("cartoon_transparency", 0.7)


def analyze_loaded_structure(selection="all", receptor_chains="A,B,C,D", partner_chains="E",
                             output_csv=None, add_hydrogen=True, use_plip=False):
    receptor_chains = _parse_chain_arg(receptor_chains) or []
    partner_chains = _parse_chain_arg(partner_chains) or []
    if not receptor_chains or not partner_chains:
        print("receptor_chains and partner_chains are required")
        return {}

    if cmd.count_atoms(selection) == 0:
        print(f"Selection not found or empty: {selection}")
        return {}

    # --- PLIP mode ---
    if use_plip:
        temp_pdb = _save_temp_pdb(selection)
        if not temp_pdb:
            print("PLIP: Cannot save temp PDB, falling back to PyMOL detection")
            use_plip = False
        else:
            plip_data = _call_plip(temp_pdb, receptor_chains, partner_chains, mode='ppi')
            if plip_data and plip_data.get('interactions'):
                _draw_plip_results(plip_data)
                print("PLIP: Analyzed selection: %s" % selection)
                print("PLIP: Chains: R=%s P=%s" % (receptor_chains, partner_chains))
                print("PLIP: Total interactions: %d" % plip_data['n_interactions'])
                return plip_data
            else:
                print("PLIP: No interactions found or analysis failed, falling back")

    # --- Original PyMOL detection (fallback or use_plip=False) ---
    if add_hydrogen:
        cmd.h_add(selection)

    _apply_default_contact_style(receptor_chains, partner_chains)

    allpairs = {}
    part1 = []
    part2 = []
    cmd.select('should_label', 'none')

    sb, hb = list_electron_interaction(selection, receptor_chains, partner_chains)
    hp = list_hydrophobic(selection, receptor_chains, partner_chains)

    select_interface_pair(hb, "", "hbonds", allpairs, receptor_chains, partner_chains, part1, part2)
    select_interface_pair(sb, "", "saltbridges", allpairs, receptor_chains, partner_chains, part1, part2)
    select_interface_pair(hp, "", "hydrophobic", allpairs, receptor_chains, partner_chains, part1, part2)

    model_name = str(selection).replace(" ", "_")
    cmd.group(f"{model_name}_contact", f"hbonds hydrophobic saltbridge* {model_name}")

    if output_csv:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("a", newline="") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow([model_name])
            for type1 in allpairs:
                for chain in allpairs[type1]:
                    pairs = sorted(set(allpairs[type1][chain]))
                    writer.writerow([type1, chain, "{" + ", ".join(repr(p) for p in pairs) + "}"])
            writer.writerow([])

    summary = {k: sum(len(set(v)) for v in allpairs[k].values()) for k in allpairs}
    print(f"Analyzed selection: {selection}")
    print(f"Receptor chains: {receptor_chains}; partner chains: {partner_chains}")
    print(f"Contact summary: {summary}")
    return allpairs


def load_structures_from_dir(input_dir, pattern="*.pdb", object_prefix=""):
    from glob import glob

    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        print(f"Input directory not found: {input_dir}")
        return []

    loaded = []
    for file_path in sorted(glob(os.path.join(input_dir, pattern))):
        obj_name = object_prefix + Path(file_path).stem
        cmd.load(file_path, obj_name)
        loaded.append(obj_name)

    print(f"Loaded {len(loaded)} structures from {input_dir}")
    return loaded


def draw_loaded_contacts(selection="all", receptor_chains="A,B,C,D", ligand_chains="E", 
                          output_csv=None, use_plip=False):
    return analyze_loaded_structure(selection=selection, receptor_chains=receptor_chains,
                                    partner_chains=ligand_chains, output_csv=output_csv,
                                    use_plip=use_plip)


def prolif_pymol_agent(mode="loaded", input_dir=None, pattern="*.pdb", selection="all",
                       receptor_chains="A,B,C,D", partner_chains="E", output_csv=None,
                       object_prefix="", use_plip=False):
    mode = str(mode).strip().lower()
    if mode == "loaded":
        return analyze_loaded_structure(selection=selection, receptor_chains=receptor_chains,
                                        partner_chains=partner_chains, output_csv=output_csv,
                                        use_plip=use_plip)

    if mode == "dir":
        if not input_dir:
            print("input_dir is required when mode='dir'")
            return {}
        loaded_objects = load_structures_from_dir(input_dir, pattern=pattern, object_prefix=object_prefix)
        results = {}
        for obj_name in loaded_objects:
            results[obj_name] = analyze_loaded_structure(
                selection=obj_name,
                receptor_chains=receptor_chains,
                partner_chains=partner_chains,
                output_csv=output_csv,
                use_plip=use_plip,
            )
        return results

    print("mode must be 'loaded' or 'dir'")
    return {}


cmd.extend("draw_loaded_contacts", draw_loaded_contacts)
cmd.extend("prolif_pymol_agent", prolif_pymol_agent)


# ======================== PLIP Integration ========================

def _save_temp_pdb(selection="all", output_path=PLIP_TEMP_PDB):
    """Save a PyMOL selection to a temp PDB file for PLIP analysis."""
    try:
        cmd.save(output_path, selection)
        return output_path
    except Exception as e:
        print("PLIP: Failed to save temp PDB:", e)
        return None


def _call_plip(pdb_file, receptor_chains, partner_chains, mode='ppi',
               cutoff=5.0, timeout=120):
    """
    Call the PLIP detection script as subprocess (system Python).
    Returns parsed JSON dict.
    """
    if not os.path.exists(PLIP_SCRIPT):
        print("PLIP: Script not found at %s" % PLIP_SCRIPT)
        return None
    if not os.path.exists(pdb_file):
        print("PLIP: PDB not found at %s" % pdb_file)
        return None

    rc = ','.join(receptor_chains) if isinstance(receptor_chains, (list, tuple)) else receptor_chains
    pc = ','.join(partner_chains) if isinstance(partner_chains, (list, tuple)) else partner_chains

    cmd_args = [
        PLIP_PYTHON, PLIP_SCRIPT,
        pdb_file, rc, pc,
        '--mode', mode,
        '--cutoff', str(cutoff),
    ]
    try:
        output = subprocess.check_output(cmd_args, timeout=timeout,
                                         stderr=subprocess.PIPE)
        data = json.loads(output.decode('utf-8'))
        return data
    except subprocess.TimeoutExpired:
        print("PLIP: Analysis timed out (%ds)" % timeout)
    except subprocess.CalledProcessError as e:
        print("PLIP: Error (%d): %s" % (e.returncode, e.stderr.decode()[:200]))
    except json.JSONDecodeError as e:
        print("PLIP: JSON parse error:", e)
    return None


def _draw_plip_results(plip_data, receptor_ns='receptor', partner_ns='partner'):
    """
    Draw interactions from PLIP JSON results in PyMOL.
    Uses existing MOLAND colors and interaction type conventions.

    Args:
        plip_data: parsed JSON dict from _plip_detect.py
        receptor_ns: name for receptor chain object (for labeling)
        partner_ns: name for partner chain object
    """
    # PLIP type → PROLIF color key mapping
    PLIP_TYPE_COLOR_MAP = {
        'HydrogenBond': 'HBDonor',
        'Hydrophobic': 'Hydrophobic',
        'PiStacking': 'PiStacking',
        'PiCation': 'PiCation',
        'SaltBridge': 'Cationic',
        'HalogenBond': 'XBDonor',
        'WaterBridge': 'WaterBridge',
    }

    if not plip_data or 'interactions' not in plip_data:
        print("PLIP: No interactions to draw")
        return

    interactions = plip_data['interactions']
    drawn_types = set()
    contact_count = {}
    label_sel = 'should_label'
    cmd.select(label_sel, 'none')

    # Draw each interaction as a distance dash
    for i, entry in enumerate(interactions):
        try:
            itype = entry['type']
            r_chain = entry['receptor_chain']
            r_resnum = entry['receptor_resnum']
            r_resname = entry['receptor_resname']
            r_atom = entry['receptor_atom']
            p_chain = entry['partner_chain']
            p_resnum = entry['partner_resnum']
            p_resname = entry['partner_resname']
            p_atom = entry['partner_atom']
            dist = entry['distance']

            # Skip dashes that would hit numeric atom indices (from PLIP ligand mode)
            if isinstance(p_atom, int) or (isinstance(p_atom, str) and p_atom.isdigit() and int(p_atom) > 1000):
                continue

            # Build PyMOL selections for both atoms
            sel_r = '(all and chain %s and resi %s and name %s)' % (r_chain, r_resnum, r_atom)
            sel_p = '(all and chain %s and resi %s and name %s)' % (p_chain, p_resnum, p_atom)

            # Generate unique dash name by type + residue identifier
            dash_name = 'd_%s_%s%s_%s%s' % (itype, r_resname, r_resnum, p_resname, p_resnum)

            # Map PLIP type to PROLIF color key
            prolif_key = PLIP_TYPE_COLOR_MAP.get(itype, itype)
            try:
                color_name = PROLIF_INTERACTION_COLORS.get(prolif_key, 'white')
            except:
                color_name = 'white'

            # Create the distance dash
            cmd.distance(dash_name, sel_r, sel_p)
            cmd.color(color_name, dash_name)
            cmd.hide('labels', dash_name)
            cmd.set('dash_gap', 0.3, dash_name)
            cmd.set('dash_length', 0.3, dash_name)
            cmd.set('dash_radius', 0.04, dash_name)

            # Show sidechain sticks for both residues (colored by type)
            sel_r_sc = '(all and chain %s and resi %s)' % (r_chain, r_resnum)
            sel_p_sc = '(all and chain %s and resi %s)' % (p_chain, p_resnum)
            cmd.show('sticks', sel_r_sc)
            cmd.show('sticks', sel_p_sc)
            cmd.color(color_name, sel_r_sc)
            cmd.color(color_name, sel_p_sc)

            # Track for labeling
            cmd.select(label_sel, '%s or %s or %s' % (label_sel, sel_r_sc, sel_p_sc))

            # Type counts
            drawn_types.add(itype)
            contact_count[itype] = contact_count.get(itype, 0) + 1

        except Exception as e:
            print("PLIP: Error drawing interaction %d: %s" % (i, e))

    print("PLIP: Drawn %d interactions (%s)" % (
        len(interactions),
        ', '.join('%s=%d' % (t, contact_count[t]) for t in sorted(contact_count))))


# ======================== Command extensions ========================