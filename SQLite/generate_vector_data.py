#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:30:28 2018
Generates all substructures per class, and for that class, it creates the 
junction table (structure_has_substructure table).
Command line: python3 generate_vector_data.py
@author: stokm006
"""

from __future__ import print_function
import sqlite3
from rdkit import Chem

def generate_substructures():
    """ Uses the structure table from the NP SQLite database and generates
    substructures data which will be stored in a text file.
    
    """

    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Generate the nr of structures for class XXX
    numrows = c.execute("SELECT count(*) FROM structure WHERE Class = 'XXX';")
    numrows = str(numrows.fetchone()).lstrip("(").rstrip(",)")
    numrows = int(numrows)

    # Select all data for class XXX
    c.execute("SELECT * FROM structure WHERE Class = 'XXX';")
    
    str_mol_list = []
    sub_smiles_list = []
    for x in range(0, numrows):
        row = c.fetchone()
        m = Chem.MolFromSmiles(row[3])
        str_mol_list += [m]
        nr_of_atoms = m.GetNumAtoms()

        substructures_list = []
        for i in range(nr_of_atoms):
            for j in range(nr_of_atoms):
                mol = Chem.FindAtomEnvironmentOfRadiusN(m,i,j)
                substructures_list += [mol]

        smile_list = []
        for mol in substructures_list:
            amap={}
            submol=Chem.PathToSubmol(m, mol, atomMap=amap)
            p = Chem.MolToSmiles(submol, canonical=True)
            # prevent overlapping substructures per structure
            if p != '' and p not in smile_list:
                smile_list += [p]
       
        for sm in smile_list:
            # prevent overlapping substructure for all structures
            if sm not in sub_smiles_list:
                sub_smiles_list += [sm]
    
    sub_mol_list = []
    for sub_struc in sub_smiles_list:
        sm = Chem.MolFromSmiles(sub_struc)
        if sm != None:
            sub_mol_list += [sm]

    with open("/path/to/store/substructures_from_class_XXX.txt", 'w') as db_file:
        db_file.write("Structure has Substructure" + '\n\n')
        for structure in str_mol_list:
            for substructure in sub_mol_list:
                if structure.HasSubstructMatch(substructure) == True:
                    struc = Chem.MolToSmiles(structure)
                    substruc = Chem.MolToSmiles(substructure)
                    db_file.write(struc + '\t' + substruc + '\n')

    # Close the connection
    conn.close()


if __name__ == '__main__':
    generate_substructures()