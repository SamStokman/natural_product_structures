#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 08:58:48 2018
Generates substructures for all structures from the SQLite NP database. The 
data file which will be created contains the class, the structure, the 
substructure and the nr of substructure matches per structure.
Command line: python3 substructures_for_all_structures.py
@author: stokm006
"""

from __future__ import print_function
import sqlite3
from rdkit import Chem

def generate_substructures():
    """ Uses the structure table from the NP SQLite database and generates
    substructures data which will be stored in a text file.
    
    """

    # Connect to SQlite database
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Generate the nr of structures present in the SQLite database
    numrows = c.execute("SELECT count(*) FROM structure;")
    numrows = str(numrows.fetchone()).lstrip("(").rstrip(",)")
    numrows = int(numrows)
    
    # Select all data from the structure table
    c.execute("SELECT * FROM structure;")
    
    str_mol_list = []
    sub_smiles_list = []
    class_name_list = []
    for x in range(0, numrows):
        row = c.fetchone()
        # Select the classes
        class_name = [row[9]]
        class_name = str(class_name).lstrip("['").rstrip("']")
        class_name_list += [class_name]
        # Select the smiles
        m = Chem.MolFromSmiles(row[3])
        str_mol_list += [m]
        nr_of_atoms = m.GetNumAtoms()

        # Generate all possible mol environments per structure
        mol_env_list = []
        for i in range(nr_of_atoms):
            for j in range(nr_of_atoms):
                env = Chem.FindAtomEnvironmentOfRadiusN(m,i,j)
                mol_env_list += [env]

        # Generate all possible substructures based on the mol envs.
        smile_list = []
        for env in mol_env_list:
            amap={}
            submol=Chem.PathToSubmol(m, env, atomMap=amap)
            # Generate the mol of each substructure
            mol = Chem.MolToSmiles(submol, canonical=True)
            # Prevent overlapping substructures per structure (the nr of 
            # generated substructures here is not correct yet)
            if mol != '' and mol not in smile_list:
                smile_list += [mol]
        
        # Add the substructure to the 'all substructures list'
        for smile in smile_list:
            # Prevent overlapping substructure for all structures 
            if smile not in sub_smiles_list:
                sub_smiles_list += [smile]

    # Convert all substructure smiles to mols
    sub_mol_list = []
    for sub_struc in sub_smiles_list:
        sm = Chem.MolFromSmiles(sub_struc)
        if sm != None:
            sub_mol_list += [sm]
        
    # Create textfile 
    with open("/mnt/scratch/stokm006/generate_substructures/structure_substructure_table.txt", 'w') as db_file:
        db_file.write("Structure"+'\t'+"Substructure"+'\t'+\
                      "Nr of substructures matches in structure"+'\n\n')
        for i, structure in enumerate(str_mol_list):
            for substructure in sub_mol_list:
                # Match the correct substructure to the correct structure
                if structure.HasSubstructMatch(substructure) == True:
                    # Generate nr of substructure matches per structure
                    match = structure.GetSubstructMatches(substructure)
                    quantity = len(match)
                    # Convert mols to smiles
                    struc = Chem.MolToSmiles(structure)
                    substruc = Chem.MolToSmiles(substructure)
                    db_file.write(str(class_name_list[i])+'\t'+ str(struc)\
                                  +'\t'+str(substruc) +'\t'+str(quantity)\
                                  +'\n')

    # Close the connection
    conn.close()


if __name__ == '__main__':
    generate_substructures()