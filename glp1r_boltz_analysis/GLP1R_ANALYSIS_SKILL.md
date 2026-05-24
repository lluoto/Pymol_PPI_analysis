# GLP-1R Boltz Structure Analysis — Agent & Skill Definition

## Overview

Complete analysis pipeline for 3001 Boltz-predicted GLP-1R structures from Schr&ouml;dinger 6VCB docking experiment. Includes centroid clustering, affinity histogram, LEU14 proximity filtering, PLIP interaction analysis, and PyMOL PSE generation with contact-type-colored interaction dashes.

## Environment

| Resource | Path/Config |
|---|---|
| Remote Host | `cuixi` (SSH 10.19.25.48, user cuixi) |
| Data Root | `/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/` |
| Receptor PDB | `/media/cuixi/data01/lluoto/6vcb_receptor.pdb` |
| PLIP Results | `/media/cuixi/data01/lluoto/vina_plip_results.json` |
| PyMOL | `/media/cuixi/data01/lluoto/PyMOL-2.6.2_appveyor1945-Linux-x86_64-py311/pymol/pymol` |
| Local Scripts | `C:\Users\lluoto\AppData\Local\Temp\opencode\` |

## Key Variables & Constants

```python
# Centroid (Boltz frame, mean across 3001 structures)
MEAN_CENTROID = (0.1713, 0.1209, -0.0928)

# Pocket classification
POCKET = 913  # 899 (centroid-based) + 14 (LEU14-recovered)
FAR = 2088     # 3001 - 913

# pIC50 ranges
POCKET_PIC50 = [0.0325, 1.6805]
FAR_PIC50 = [-0.9206, 1.5156]

# Chain mapping: Interaction JSON → PDB
CHAIN_MAP = {"A": "R", "B": "P"}  # A=receptor chain R, B=peptide chain P

# Key hot spots (5 residues of interest)
HOT5_RESIDUES = [("A", 103), ("A", 106), ("A", 159), ("A", 163), ("B", 14)]
# Residue names: LEU103, TYR106, ASP159, LYS163, LEU14

# k-means clusters (894 3plus structures)
K = 2  # silhouette = 0.5147
C1 = {"size": 369, "medoid_seed": 925, "ligand": "6vcb_3_20", "pic50": 0.686}
C2 = {"size": 525, "medoid_seed": 420, "ligand": "6vcb_3_24", "pic50": 1.200}

# Hot5 representative (432/894 structures have all 5 hot spots)
HOT5_REP = {"seed": 899, "ligand": "6vcb_3_24", "pic50": 1.244, "centroid": (152.99, 146.47, 127.3)}
```

## ProLIF Color Scheme (from prolif_pymol_vis.py)

```python
MOLAND_COLORS = {
    "HydrogenBond": [0.204, 0.561, 0.655],   # teal
    "Hydrophobic":  [0.69, 0.67, 0.72],       # gray-purple
    "PiStacking":   [0.812, 0.686, 0.831],     # purple
    "PiCation":     [0.97, 0.79, 0.49],        # orange
    "SaltBridge":   [0.929, 0.759, 0.662],      # tan
    "HalogenBond":  [0.251, 0.718, 0.678],      # green-teal
}
```

## Residue → Interaction Type Classification

```python
def classify_plip_type(resname):
    name = resname.upper()
    if name in ("LYS", "ARG", "HIS"): return "PiCation"
    elif name in ("ASP", "GLU"): return "SaltBridge"
    elif name in ("TYR", "TRP", "PHE"): return "PiStacking"
    elif name in ("SER", "THR", "ASN", "GLN", "CYS", "MET"): return "HydrogenBond"
    elif name == "CL": return "HalogenBond"
    else: return "Hydrophobic"
```

## Script Inventory

| Script | Location | Purpose |
|---|---|---|
| `boltz_centroid_analysis.py` | local temp | Process 3001 CIF files for centroid extraction |
| `annotate_affinity_histogram.py` | Schr&ouml;dinger dir | Affinity histogram with pocket/far classification |
| `filter_pocket.py` | local temp | LEU14 proximity-based pocket filtering |
| `cluster_dash.py` | `/tmp/` (remote) | **Primary**: k-means clustering + PyMOL PSE generation with contact-type dash objects |
| `find_5hotspot.py` | `/tmp/` (remote) | Find structures with all 5 key hot spot residues; generate hot5_representative.pse |

## Output Files

| File | Location | Description |
|---|---|---|
| `affinity_histogram_annotated.png` | local | Final histogram (blue bars, red Gaussian fit, text stats) |
| `cluster_1.pse` | remote `3plus_clusters_dash/` | C1 medoid (seed 925), 7 dashes named by type |
| `cluster_2.pse` | remote `3plus_clusters_dash/` | C2 medoid (seed 420), 6 dashes named by type |
| `hot5_representative.pse` | remote `3plus_clusters_dash/` | Hot5 rep (seed 899), 10 dashes named by type |
| `centroid_results.csv` | local | Per-structure centroid data |
| `far_centroid_clusters/` | local | Far-cluster analysis output |

## PyMOL PSE Generation Workflow

The `cluster_dash.py` script performs the full pipeline:

1. **Load PLIP data** — Parse `vina_plip_results.json`, index by seed number
2. **Parse centroids** — Extract ligand centroid from each PDB (chain L HETATM)
3. **k-means clustering** — On centroid coordinates, auto-select k via silhouette score
4. **Medoid selection** — Per cluster, pick structure with minimum sum of intra-cluster distances
5. **Extract ligands** — Write ligand-only PDBs (chain L) to `ligands/` directory
6. **Generate PML** — Write PyMOL command files with:
   - Receptor cartoon (chain R/P, protein_gray, 30% transparency)
   - Ligand sticks (chain L)
   - Interaction dashes: `cmd.distance('d_{Type}_{Resname}{Resid}', ...)`
   - "Contact-type-colored dashes (hydrogenbond/hydrophobic/pistacking/pication/saltbridge/halogenbond)
   - Sidechain sticks colored by interaction type
   - Element color override (N=blue, O=red) applied last
   - Camera view centered on interaction contacts
7. **Render PSE** — Invoke PyMOL in batch mode, save `.pse`

### PML Conventions

- **Dash object naming**: `d_{Type}_{Resname}{Resid}` (e.g., `d_Hydrophobic_LEU103`, `d_PiStacking_TYR106`)
- **Sidechain selection naming**: `sc_{chain}_{resid}` (e.g., `sc_A_103`, `sc_B_14`)
- **Chain mapping**: Interaction JSON chain `A` → PDB chain `R`; JSON chain `B` → PDB chain `P`
- **Labels**: Hidden on dashes (`cmd.hide('labels', ...)`)
- **No `forbidden hide`**: All selections/manipulations remain enabled

## Common Operations

### Re-run cluster PSE generation
```bash
ssh cuixi "python3 /tmp/cluster_dash.py"
```

### Re-run hot5 representative finder
```bash
ssh cuixi "python3 /tmp/find_5hotspot.py"
```

### Copy PSE files locally
```bash
scp "cuixi:/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/3plus_clusters_dash/*.pse" .
```

## Interaction JSON Structure

Per-structure JSON at `{BASE}/{ligand}_pic50_{value}_seed_{seed}_interactions.json`:
```json
{
  "interactions": [
    {"chain": "A", "resid": "103", "resname": "LEU", "interacting": true,
     "pairs": [["CG", "C47", 3.63], ...]},
    ...
  ]
}
```

## Troubleshooting

| Issue | Fix |
|---|---|
| Dashes show distance label instead of type | `cmd.hide('labels', 'd_*')` — already applied |
| Sidechains not visible | Check chain mapping (A→R, B→P); verify PDB has sidechain atoms |
| Element colors overwritten | Ensure `cmd.color('elem_N/O/S', ...)` is the FINAL command in PML |
| Receptor/ligand not loading | Verify `medoid_path` and `medoid_lpath` are valid paths |
| `forbidden hide manipulate` error | Ensure no PyMOL settings file restricts `hide` or `manipulate` |
