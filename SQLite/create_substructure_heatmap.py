#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 14:34:43 2018
Creates a similarity heatmap with substructures from class XXX generated with
generate_data_for_vector.py
Command line: python3 create_substructure_heatmap.py /mnt/nexenta/stokm006/class_Akua_str_sub.txt
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import seaborn as sns


def similarity(smiles_data):
    # Generate substructure list
    sub_smiles = smiles_data.split('\n')
    sub_smiles = sub_smiles[2:-1]
   
    sub_smiles_list0 = []
    for sub_smile in sub_smiles:
        sub_smile = sub_smile.split('\t')
        sub_smile = sub_smile[1]
        if sub_smile not in sub_smiles_list0:
            sub_smiles_list0 += [sub_smile]
  
    sub_smiles_list = []
    sub_mol_list = []
    for sub_mol in sub_smiles_list0:
        sm = Chem.MolFromSmiles(sub_mol)
        if sm != None:
            sub_mol_list += [sm]
            sm = Chem.MolToSmiles(sm)
            sub_smiles_list += [sm]
    
    sim_score_list2 = []
    fps = [FingerprintMols.FingerprintMol(x) for x in sub_mol_list]
    for i in range(len(fps)):
        sim_score_list = []
        for x in range(len(fps)):
            sim_score = DataStructs.FingerprintSimilarity(fps[x], fps[i])     
            sim_score_list += [sim_score]
        sim_score_list2 += [sim_score_list]
    
    
    xylabels = []
    for i in range(len(sub_smiles_list)):
        xylabels += [i]
    
    sim_score_array = np.asarray(sim_score_list2)
    
    return sim_score_array, xylabels

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        vector_array, sub_list = similarity(input_file)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(35,30)) 
    ax = sns.heatmap(vector_array, linewidths=.5, xticklabels=True, yticklabels=True, cmap="YlGnBu")
    ax.set_yticklabels(sub_list, fontsize=5)
    ax.set_xticklabels(sub_list, fontsize=5)
    plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    plt.xlabel('Substructure')
    plt.ylabel('Substructure')
    fig.savefig('/mnt/scratch/stokm006/SQLite/Akuammilan_sub_Heatmap.png')
 #   plt.show()

