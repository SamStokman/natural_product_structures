#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 15:08:28 2018
Make a vector array for a class (x: substructure, y:structures). This vector
array can be visualized with a heatmap.(see generate_data_for_vector.py)
Command line: python3 create_vector.py /mnt/nexenta/stokm006/Akua_str_sub.txt
@author: stokm006
"""

from __future__ import print_function
import sqlite3
from rdkit import Chem
from sys import argv
import numpy as np
import matplotlib.pyplot as plt


def vector(data_input):
    """ Uses the data created with generate_data_for_vector.py and creates a
    vector (0/1) and a heatmap.
    
    data_input: txt file from generate_data_for_vector.py
    """
    
    # Generate substructure list
    sub_smiles = data_input.split('\n')
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
         
    # Connect to SQLite database
    conn = sqlite3.connect("/path/to/SQLiteDatabse/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Generate the nr of structures for class XXX
    numrows = c.execute("SELECT count(*) FROM structure WHERE Class = 'XXX';")
    numrows = str(numrows.fetchone()).lstrip("(").rstrip(",)")
    numrows = int(numrows)

    # Select all data for class XXX
    c.execute("SELECT * FROM structure WHERE Class = 'XXX';")
    
    # Generate structure list
    str_mol_list = []
    str_smiles_list = []
    for x in range(0, numrows):
        row = c.fetchone()
        str_smiles_list += [row[3]]
        m = Chem.MolFromSmiles(row[3])
        str_mol_list += [m]
 
    # Substructure matching
    vector_list = []
    for structure in str_mol_list:
        vector2_list = []
        for substructure in sub_mol_list:
            if structure.HasSubstructMatch(substructure) == True:
                vector2_list += [1]
            if structure.HasSubstructMatch(substructure) == False:
                vector2_list += [0]
        vector_list += [vector2_list]
        

    vector_array = np.array(vector_list) # if you print this you get the vector [0,1,etc] overview

    # Create heatmap
    column_labels = sub_smiles_list
    row_labels = str_smiles_list
    fig, axis = plt.subplots()
    heatmap = axis.pcolor(vector_array, cmap=plt.cm.Blues)
    axis.set_yticks(np.arange(vector_array.shape[0])+0.5, minor=False)
    axis.set_xticks(np.arange(vector_array.shape[1])+0.5, minor=False)

    axis.set_yticklabels(row_labels, fontsize=8, minor=False)
    axis.set_xticklabels(column_labels, fontsize=5, minor=False)
    plt.setp(axis.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    plt.xlabel('Substructure')
    plt.ylabel('Structure')
    fig.set_size_inches(340, 4)
    fig.grid(which="minor", color="w", linestyle='-', linewidth=33)
    plt.colorbar(heatmap)
    fig.savefig('/path/to/store/vectorHeatmap.png')
    plt.show()

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        vector(input_file)
