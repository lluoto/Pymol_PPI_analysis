# Pymol_PPI_analysis

PyMOL-based protein-protein interaction (PPI) analysis script. Detect and visualize hydrogen bonds, salt bridges, and hydrophobic contacts between receptor and partner chains. Supports both single-structure analysis and high-throughput batch processing.

## Features

- **Analyze loaded structures** — detect contacts in structures already loaded in PyMOL
- **Batch directory mode** — load all PDBs from a directory and analyze each one
- **Three interaction types**: Hydrogen bonds, Salt bridges, Hydrophobic contacts
- **Automatic color styling** with ProLIF/moland color scheme
- **CSV export** of contact residue pairs for downstream analysis
- **Sidechain stick display** for interacting residues

## Dependencies

- [PyMOL](https://pymol.org/) 2.4+ (Schrodinger or open-source)
- Python 3.8+ (bundled with PyMOL)

## Files

| File | Description |
|---|---|---|
| `prolif_pymol_vis.py` | Main PPI analysis script (unified agent entrypoint,final script) |
| `residue_contact_ring.py` | Ring-specific contact analysis variant |
| `findcontactsV2-3.tcl` | VMD Tcl script for trajectory contact/SASA analysis (forked from [FindContacts](https://github.com/TMiclot/FindContacts)) |
| `glp1r_boltz_analysis/` | GLP-1R Boltz structure analysis pipeline — centroid clustering, affinity histogram, PLIP interaction analysis, PyMOL PSE generation with ProLIF-colored dashes |

## Quick Start

### 1. Load the script in PyMOL

```pymol
run prolif_pymol_vis.py
```

### 2. Analyze an already-loaded structure

```pymol
prolif_pymol_agent loaded, selection=all, receptor_chains=A,B,C,D, partner_chains=E
```

### 3. Load structures from a directory and analyze

```pymol
prolif_pymol_agent dir, input_dir=/path/to/pdbs, pattern=*.pdb, receptor_chains=A,B,C,D, partner_chains=E, output_csv=/path/to/contacts.csv
```

### 4. Shortcut command

```pymol
draw_loaded_contacts all, A,B,C,D, E
```

## Demo

### Example: GPCR–G protein complex

```pymol
# Load the script
run prolif_pymol_vis.py

# Load a GPCR-G protein PDB
fetch 7W0N

# Analyze contacts: receptor chains A,B,C,D, partner chain E
prolif_pymol_agent loaded, selection=7W0N, receptor_chains=A,B,C,D, partner_chains=E
```

### Example: Batch analysis of docking results

```pymol
run prolif_pymol_vis.py

# Analyze all PDBs in a docking output directory
prolif_pymol_agent dir, input_dir=/mnt/sdb2/lluoto/final_top, receptor_chains=A,B,C,D, partner_chains=E, output_csv=/mnt/sdb2/lluoto/all_contacts.csv
```

## Input Data Format

### Structure files
- Standard PDB format
- Must have chain identifiers (e.g., `A`, `B`, `C`, `D` for receptor; `E` for ligand/partner)
- Multi-model PDBs are loaded as separate objects

### Chain convention
- **Receptor chains**: Typically `A,B,C,D` (the first 4 chains)
- **Partner/Ligand chain**: Typically `E` (the last chain)

## Output

### In PyMOL
- Colored sticks for interacting sidechains
- Distance measurements drawn between contact atoms
- Auto-generated group objects per interaction type

### CSV export (optional)
When `output_csv` is specified, writes contact pairs in the format:

| interaction_type | chain_pair | residue_pairs |
|---|---|---|
| hbonds | A,E | {"454D,123T", "607E,118R", ...} |
| saltbridges | A,E | {"401R,235E", "448R,35E"} |
| hydrophobic | A,E | {"709L,135I", "720L,130I", ...} |

Residue notation: `{residue_number}{one-letter amino acid code}` (e.g., `454D` = residue 454, ASP).

## Command Reference

### `prolif_pymol_agent`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mode` | str | `loaded` | `"loaded"` or `"dir"` |
| `input_dir` | str | `None` | Directory with PDB files (required for `dir` mode) |
| `pattern` | str | `*.pdb` | File glob pattern for `dir` mode |
| `selection` | str | `all` | PyMOL selection expression |
| `receptor_chains` | str | `A,B,C,D` | Comma-separated receptor chain IDs |
| `partner_chains` | str | `E` | Comma-separated partner chain IDs |
| `output_csv` | str | `None` | Path to output CSV file |
| `object_prefix` | str | `""` | Prefix for object names in `dir` mode |

### `draw_loaded_contacts`

```pymol
draw_loaded_contacts selection, receptor_chains, ligand_chains, output_csv
```

## Color Scheme

| Interaction | Color |
|---|---|
| Hydrogen bonds | Hydrogen bond blue |
| Salt bridges | Salt bridge red |
| Hydrophobic | Hydrophobic gray |
| Aromatic (π-stacking) | Pi stack purple |
| Aromatic–cation | Pi cation orange |

## Related Tools

- [FindContacts](https://github.com/TMiclot/FindContacts) — VMD Tcl scripts for trajectory-based contact and SASA analysis
- [ProLIF](https://github.com/chemosim-lab/ProLIF) — Python library for protein-ligand interaction fingerprints

## Authors

- Chenxi Wang (chenxi.wang@salilab.org)
- Yulin Luo

## License

MIT
