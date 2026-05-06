# PPI_analysis

simple pymol base protein-protein interaction analyze script, could work in high throughtput
precisional denote should utility on plip and others server/software

besides fork others findcontact script for trajectory SASA/contact calcualtion
https://github.com/TMiclot/FindContacts/tree/main

run /home/cuixi/prolif_pymol_vis.py
Analyze already loaded structure:
prolif_pymol_agent loaded, selection=all, receptor_chains=A,B,C,D, partner_chains=E
Load from directory and analyze:
prolif_pymol_agent dir, input_dir=/path/to/pdbs, pattern=*.pdb, receptor_chains=A,B,C,D, partner_chains=E
Shortcut:
draw_loaded_contacts all, A,B,C,D, E
If you want, I can next make chain inference automatic instead of requiring A,B,C,D and E.
