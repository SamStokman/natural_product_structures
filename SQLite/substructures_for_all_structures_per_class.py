#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 10:19:39 2018
Generates substructures for all structures per class(!) from the SQLite NP 
database. The data file which will be created contains the class, 
the structure, the substructure and the nr of substructure matches per 
structure.
Command line: python3 substructures_for_all_structures_per_class.py
@author: stokm006
"""

from __future__ import print_function
import sqlite3
from rdkit import Chem

def generate_substructures():
    """ Uses the structure table from the NP SQLite database and generates
    substructures data which will be stored in a text file.
    
    """
    
    # Create textfile 
    with open("/path/to/data/file/structure_substructure_per_class_table.txt", 'w') as db_file:
         db_file.write("Class"+'\t'+"Structure"+'\t'+"Substructure"+'\t'+\
                       "Nr of substructures matches in structure" + '\n\n')

    # Connect to SQlite database
    conn = sqlite3.connect("/path/tpSQLiteDatabase/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Generate the nr of structures present in the SQLite database
    numrows = c.execute("SELECT count(*) FROM structure;")
    numrows = str(numrows.fetchone()).lstrip("(").rstrip(",)")
    numrows = int(numrows)
    
    # Select all data from the structure table
    c.execute("SELECT * FROM structure;")
    
    class_list = []
    class_struc_list = []
    for x in range(0, numrows):
        row = c.fetchone()
        class_name = [row[9]]
        class_name = str(class_name).lstrip("['").rstrip("']")
        # Create a list with unique classes
        if class_name not in class_list:
            class_list += [class_name]
        m = Chem.MolFromSmiles(row[3])
        # Create a list of list with the class and structure
        class_struct_combo = [class_name, m]
        class_struc_list += [class_struct_combo]
    
    # Generate the substructures per class
    for class_name in class_list:
        str_mol_list = []
        for class_struc in class_struc_list:
            # Group all structures for the specific class
            if class_struc[0] == class_name:
                str_mol_list += [class_struc]
        
        # Generate the mols for each structure in the class
        sub_smiles_list = []
        for m in str_mol_list:
            m = m[1]
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
                # Prevent overlapping substructures within the same structure
                # (the nr of generated substructures here is not correct yet)
                if mol != '' and mol not in smile_list:
                    smile_list += [mol]
 
            # Add the substructure to the 'all substructures list' (per class)
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

        # Add (sub)structure data per class to text file 
        with open("/path/to/data/file/structure_substructure_per_class_table.txt", 'a') as db_file:
            for structure in str_mol_list:
                structure = structure[1]
                for substructure in sub_mol_list:
                    # Match the correct substructure to the correct structure
                    if structure.HasSubstructMatch(substructure) == True:
                        # Generate nr of substructure matches per structure
                        match = structure.GetSubstructMatches(substructure)
                        quantity = len(match)
                        # Convert mols to smiles
                        struc = Chem.MolToSmiles(structure)
                        substruc = Chem.MolToSmiles(substructure)
                        db_file.write(str(class_name) +'\t'+ str(struc) +'\t'+ str(substruc) +'\t'+\
                                     str(quantity) +'\n')
 
    # Close the connection
    conn.close()

if __name__ == '__main__':
    generate_substructures()
