# GLP-1R Boltz Structure Analysis — Project README

## Summary

Analysis of 3001 GLP-1R structures from Boltz-predicted docking experiments (6VCB, Schrödinger). 
Includes centroid clustering, affinity histogram, LEU14 proximity classification, PLIP interaction analysis,
and PyMOL visualizations with contact-type-colored interaction dashes.

## Key Results

| Metric | Value |
|---|---|
| Total structures | 3001 |
| Pocket (centroid + LEU14) | 913 |
| Far (>3Å centroid) | 2088 |
| 3+ interaction structures | 894 |
| Has all 5 hot spots | 432 |
| k-means clusters (894) | k=2, silhouette=0.5147 |
| C1 medoid | seed 925, 6vcb_3_20, pIC50=0.686 |
| C2 medoid | seed 420, 6vcb_3_24, pIC50=1.200 |
| Hot5 representative | seed 899, 6vcb_3_24, pIC50=1.244 |

## 5 Key Hot Spots

| Residue | Chain (JSON→PDB) | Top Type | Contact Freq |
|---|---|---|---|
| LEU 103 | A→R | Hydrophobic | 1028 |
| TYR 106 | A→R | PiStacking | 1485 |
| ASP 159 | A→R | SaltBridge | — |
| LYS 163 | A→R | PiCation | 2316 |
| LEU 14 | B→P | Hydrophobic | 1132 |

## Output Files

### PyMOL Sessions (remote)
| File | Dashes | Medoid |
|---|---|---|
| `3plus_clusters_dash/cluster_1.pse` | 7 | seed 925 |
| `3plus_clusters_dash/cluster_2.pse` | 6 | seed 420 |
| `3plus_clusters_dash/hot5_representative.pse` | 10 | seed 899 |

### Plots & Data (local)
- `affinity_histogram_annotated.png` — Final histogram
- `centroid_results.csv` — Per-structure centroids
- `affinity_histogram_annotated_data.csv` — Histogram bin data

### Scripts
- `C:\Users\lluoto\AppData\Local\Temp\opencode\cluster_dash.py` — Main PSE generator (v2, 303 lines)
- `C:\Users\lluoto\AppData\Local\Temp\opencode\find_5hotspot.py` — Hot5 finder
- `C:\Users\lluoto\Documents\Schrodinger\glp1r_6vcb_dock_2\annotate_affinity_histogram.py`

## Environment

- **Remote**: `ssh cuixi` (10.19.25.48)
- **PyMOL**: `/media/cuixi/data01/lluoto/PyMOL-2.6.2_appveyor1945-Linux-x86_64-py311/pymol/pymol`
- **Receptor**: `/media/cuixi/data01/lluoto/6vcb_receptor.pdb`
- **3plus data**: `/media/cuixi/data01/lluoto/vina_plip_pse/boltz_pse_899/interactions_3plus/`
- **PLIP**: `/media/cuixi/data01/lluoto/vina_plip_results.json`

## Chain Mapping

| Source | Target |
|---|---|
| Interaction JSON chain A | PDB chain R (receptor) |
| Interaction JSON chain B | PDB chain P (peptide) |

## ProLIF Interaction Types

| Residue | Type | Color |
|---|---|---|
| LYS, ARG, HIS | PiCation | Orange |
| ASP, GLU | SaltBridge | Tan |
| TYR, TRP, PHE | PiStacking | Purple |
| SER, THR, ASN, GLN, CYS, MET | HydrogenBond | Teal |
| CL | HalogenBond | Green-teal |
| Other (LEU, VAL, ILE, etc.) | Hydrophobic | Gray-purple |

## Quick Re-run

```bash
# Regenerate cluster PSEs
scp cluster_dash.py cuixi:/tmp/
ssh cuixi "python3 /tmp/cluster_dash.py"
```

## See Also

- `GLP1R_ANALYSIS_SKILL.md` — Full agent/skill definition for reproducibility
