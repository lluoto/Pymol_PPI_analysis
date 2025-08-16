'''
    Author: Chenxi Wang (chenxi.wang@salilab.org) Yulin Luo
    Date: 2024-6-21

    The residue contacts in the GLP-1R-Gs interface is calculated by Pymol.
    
'''

#!/usr/bin/env python
# coding: utf-8

from pymol import cmd
from pymol import stored
from os import listdir
from os.path import isfile, join
import time
from mpi4py import MPI
import csv
start=time.time()
amino_acid_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
def list_hbonds(selection,channel,partners,selection2=None,cutoff=3.5,angle=150,mode=1,hb_list_name='hbonds'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)

    selection = selection + " & e. n+o+s"
    
    if not selection2:
        hb = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
        pass
    else:
        selection2 = selection2 + " & e. n+o"
        hb = cmd.find_pairs(selection,selection2,mode = mode,cutoff = cutoff,angle = angle)
        
    # sort the list for easier reading
    hb.sort(key=lambda x:x[0][1])
    
    result = []
    
    for pairs in hb:
        stored.x = []
        stored.y = []
        cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
        cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
        
        if stored.x[0] in channel and stored.y[0] in partners:
            d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
            cmd.show('stick',"chain %s and resi %s chain %s and resi %s" % (stored.x[0],stored.x[1],stored.y[0],stored.y[1]))
            result.append([stored.x,stored.y,float(d)])
            cmd.color('skyblue',hb_list_name)
    return result
    
    


def list_saltbridge(selection,channel,partners,selection2=None,cutoff=4.0,angle=180,mode=1,hb_list_name='saltbridges'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)

    posresList = ["LYS","ARG","HIS","HIP"]
    negresList = ["GLU","ASP","CYM"]
    posatom='N'
    negatomList = ['S','O']
    
    selection = selection + " & e. n+o+s"
    
    sb = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
    sb.sort(key=lambda x:x[0][1])
    
    result = []
    for pairs in sb:
            stored.x = []
            stored.y = []
            cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name,elem]')
            cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name,elem]')
            if (stored.x[2] in posresList and stored.x[4]==posatom and stored.y[2] in negresList and stored.y[4] in negatomList):
                if stored.x[0] in channel and stored.y[0] in partners:
                    d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
                    cmd.show('stick',"chain %s and resi %s chain %s and resi %s" % (stored.x[0],stored.x[1],stored.y[0],stored.y[1]))
                    result.append([stored.x,stored.y,float(d)])
                    cmd.color('purple',hb_list_name)
            else:
                continue
    

    return result




def list_hydrophobic(selection,channel,partners,selection2=None,cutoff=4.0,angle=180,mode=1,hb_list_name='hydrophobic'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    hydrophobiclist = ["ALA","VAL","LEU","ILE","PHE","PRO","MET"]
    aromaticlist=["PHE","TYR","TRP"]
    aromaticatoms=['CG','CZ','CD1','CD2','CE1','CE2']
    selection = selection + " & e. c+n+o+s"
    hp = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
    hp.sort(key=lambda x:x[0][1])
    result = []
    for pairs in hp:
            stored.x = []
            stored.y = []
            cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
            cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
            
            if (stored.x[2] in hydrophobiclist and stored.y[2] in hydrophobiclist):
                if stored.x[0] in channel and stored.y[0] in partners:
                    cmd.show('stick',"chain %s and resi %s chain %s and resi %s" % (stored.x[0],stored.x[1],stored.y[0],stored.y[1]))
                    '''if stored.x[2] in aromaticlist: # the workflow of pseudoatoms
                        cmd.deselect()
                        cmd.select(f'chain {stored.x[0]} and resi {stored.x[1]} and name CZ+CG+CD1+CD2+CE1+CE2')
                        cmd.pseudoatom(f'{stored.x[1]}{stored.x[2]}','sele')
                        cmd.deselect()
                    if stored.y[2] in aromaticlist:
                        cmd.deselect()
                        cmd.select(f'chain {stored.y[0]} and resi {stored.y[1]} and name CZ+CG+CD1+CD2+CE1+CE2')
                        cmd.pseudoatom(f'{stored.y[1]}{stored.y[2]}','sele')
                        cmd.deselect()
                    if stored.x[2] in aromaticlist and stored.y[2] in aromaticlist:
                        d=cmd.distance(hb_list_name,f'{stored.x[1]}{stored.x[2]}',f'{stored.y[1]}{stored.y[2]}')
                        if d>5:
                            cmd.delete(hb_list_name)
                            continue
                        else:
                            result.append([[stored.x,stored.y,float(d)]])
                            cmd.color('green',hb_list_name)

                    elif stored.x[2] not in aromaticlist and stored.y[2] not in aromaticlist:
                        d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
                        result.append([stored.x,stored.y,float(d)])
                        cmd.color('green',hb_list_name)
                    else:
                        if stored.x[2] in aromaticlist and stored.x[3] in aromaticatoms:
                            d=cmd.distance(hb_list_name,f'{stored.x[1]}{stored.x[2]}',"%s and index %s" % (pairs[1][0],pairs[1][1]))
                            if d>4:
                                cmd.delete(hb_list_name)
                                continue
                            else:
                                result.append([[stored.x,stored.y,float(d)]])
                                cmd.color('green',hb_list_name)
                        elif stored.y[2] in aromaticlist and stored.y[3] in aromaticatoms:
                            d=cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),f'{stored.y[1]}{stored.y[2]}')
                            if d>4:
                                cmd.delete(hb_list_name)
                                continue
                            else:
                                result.append([[stored.x,stored.y,float(d)]])
                                cmd.color('green',hb_list_name)
                        else:'''
                    if stored.x[2] in aromaticlist and stored.y[2] in aromaticlist:
                        d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]),cutoff=5.0,mode=6)
                        result.append([stored.x,stored.y,float(d)])
                        cmd.color('green',hb_list_name)
                    elif stored.x[2] in aromaticlist and stored.x[3] in aromaticatoms:
                        d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]),cutoff=4.0,mode=7)
                        result.append([stored.x,stored.y,float(d)])
                        cmd.color('green',hb_list_name)
                    elif stored.x[2] not in aromaticlist and stored.y[2] not in aromaticlist:    
                        d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
                        result.append([stored.x,stored.y,float(d)])
                        cmd.color('green',hb_list_name)
            else:
                continue

    return result
    
def select_interface_pair(bonds,filepath,various,allpairs,channel,partners,part1,part2):

    R_A= []
    specific_pair=[(subunit,partner) for subunit in channel for partner in partners]

    for item in bonds:
        pairs = (item[0][0],item[1][0])
        # change the index to chain index

        if pairs in specific_pair: R_A += [item]
        #if pairs in [('I','D'),('D','I')]: R_A += [item]


    file = open(filepath, 'w')
    if R_A != []:
        chains=set([f'{x[0][0]},{y[1][0]}' for x in R_A for y in R_A if x[0][0]!=y[1][0] and x[0][0].isdigit()==False and y[1][0].isdigit()==False])
        allpairs[various]={}
        for chain in chains:
            allpairs[various][chain]=[]
        for item in R_A:
            pair=f'{item[0][1]}{amino_acid_map[item[0][2]]},{item[1][1]}{amino_acid_map[item[1][2]]}'
            print('R_A', item[0][0], item[0][1], item[0][2], item[0][3], item[1][0], item[1][1], item[1][2], item[1][3], item[2], sep=' ', file=file)
            part1.append(f'{item[0][0]},{item[0][1]}')
            part2.append(f'{item[1][0]},{item[1][1]}')
            allpairs[various][f'{item[0][0]},{item[1][0]}']+=[pair]
    file.close()



########################### main ###########################

input_pdb_path = '/home/lluoto/Downloads/'
inputfiles = [f for f in listdir(input_pdb_path) if isfile(join(input_pdb_path, f)) if '.pdb' in f]
output_dir = '/home/lluoto/interact/'


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    # Master process
    # Distribute files to workers
    files_per_process = len(inputfiles) // size
    extra = len(inputfiles) % size
    offsets = [0] * size
    counts = [files_per_process] * size
    for i in range(extra):
        counts[i] += 1
    for i in range(1, size):
        offsets[i] = offsets[i - 1] + counts[i - 1]
    scattered_data = [inputfiles[offsets[i]:offsets[i] + counts[i]] for i in range(size)]
else:
    scattered_data = None

# Scatter the data to all processes
local_files = comm.scatter(scattered_data, root=0)
colors = ["red", "blue", "green", "yellow", "orange"]

# Each process works on their part of the list
for pdbname in local_files:
    print(f"Process {rank} working on {pdbname}")
    cmd.load(input_pdb_path + pdbname, pdbname)
    allpairs={}
    part1=[]
    part2=[]
    #channel=['A','B','C','D']
    #partners=['E','F','G','H']
    channel=['A','B','C','D']
    partners=['E','F']
    if 'ard' in pdbname:
        channel=['A']
        partners=['B','C','D']
    cmd.h_add()
    cmd.extend("list_hbonds", list_hbonds)
    hb = list_hbonds(pdbname,channel,partners)
    select_interface_pair(hb,output_dir + 'hbonds/hbonds_' + pdbname[:-4],'hbonds',allpairs,channel,partners,part1,part2)
    print(f'all finished,took {time.time()-start}s')
    cmd.extend("list_saltbridge", list_saltbridge)
    sb = list_saltbridge(pdbname,channel,partners)
    select_interface_pair(sb, output_dir + 'saltbridges/saltbridges_' + pdbname[:-4],'saltbridges',allpairs,channel,partners,part1,part2)
    print(f'all finished,took {time.time()-start}s')
    cmd.extend("list_hydrophobic", list_hydrophobic)
    hp = list_hydrophobic(pdbname,channel,partners)
    select_interface_pair(hp, output_dir + 'hydrophobic/distance/hydrophobic_' + pdbname[:-4],'hydrophobic',allpairs,channel,partners,part1,part2)
    print(f'all finished,took {time.time()-start}s')
    cmd.group(f'{pdbname}_contact', f'hbonds hydrophobic saltbridges {pdbname}')
    cmd.bg_color('white')
    cmd.set('cartoon_transparency',0.7)
    chains = cmd.get_chains()
    for i, chain_c in enumerate(chains):
        color_index = i % len(colors)  
        color = colors[color_index]
        selection = f"chain {chain_c}"
        cmd.color(color, selection)
    #cmd.save(f'{input_pdb_path}' + f'{pdbname[:-4]}' + '.pse')
    with open(f'{output_dir}'+f'test.csv','a') as f:
        csv.writer(f,delimiter=',').writerow([pdbname[:-4]])
        for type1 in allpairs:
            for chain in allpairs[type1]:
                allww=set(allpairs[type1][chain])
                lis=[type1]+[chain]+[set(allpairs[type1][chain])]
                csv.writer(f, delimiter=',').writerow(lis)
        csv.writer(f,delimiter=',').writerow('\n')
    #for part in set(part1):
    #    face1=open(input_pdb_path+pdbname[:-4]+'.face1','a')
    #    print(part.split(',')[0],part.split(',')[1],'_',sep=' ',file=face1)
    #for part in set(part2):
    #    face2=open(input_pdb_path+pdbname[:-4]+'.face2','a')
    #    print(part.split(',')[0],part.split(',')[1],'_',sep=' ',file=face2)

    #cmd.delete('all')
cmd.save(f'{input_pdb_path}' + f'{pdbname[:-4]}' + '.pse')
# Gather results or finalize
if rank == 0:
    print("Master process collecting results...")
    # Additional code to collect and finalize results if necessary
else:
    print(f"Process {rank} finished, took {time.time()-start}s")
print(f'all finished,took {time.time()-start}s')