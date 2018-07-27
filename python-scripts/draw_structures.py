#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 12:11:51 2018
This script generates substructures from a given structure and draws them in a
png file.
Command line: python3 draw_substuctures.py 
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Draw 


def draw_substruc():
    """ Uses the given structure and generates its substructure and draws the
    structure and all substructures in a png file.
    
    """
    
    # Structure examples
    # Rubradirin
    #’CC1C=C(C(=O)C2=C3C(=CC(=C2O)C)C(=O)C4=C(C3=O)OC(CN4)(C(=O)C(C1OC(=O)C5=C(C(=CC(=N5)C(=O)NC6=C(C7=C(C=C(C=C7)OC)OC6=O)O)OC8CC(C(C(O8)C)OC)(C)[N+](=O)[O-])O)O)C)C’

    # Rifamycin
    #’CC1C=CC=C(C(=O)NC2=CC(=C3C(=C2O)C(=C(C4=C3C(=O)C(O4)(OC=CC(C(C(C(C(C(C1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)C’

    # Polyketomycin
    #’CC1C(CCC(O1)OC2C(=O)C(=C(C3(C2(CC4=C(C5=C(C(=O)C=C(C5=O)OC)C(=C4C3=O)O)C)O)O)O)C(=O)C)OC6CC(C(C(O6)C)OC(=O)C7=C(C=CC(=C7O)C)C)(C)O’

    # Everninomicin
    #’CC1C(C(CC(O1)OC2C(OC3(CC2O)OC4C(OC(C(C4(O3)C)O)OC5C(C(OC(C5OC)C)OC6C(OC(C(C6O)OC)OC7C(C8C(CO7)OC9(O8)C1C(C(CO9)OC(=O)C2=C(C=C(C=C2C)O)O)OCO1)O)COC)O)C)C)OC1CC(C(C(O1)C)OC)(C)[N+](=O)[O-])OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C’

    # Simocyclinone
    #‘CC1C(C(CC(O1)C2=C(C3=C(C=C2)C(=O)C45C(C3O)(O4)C(CC6(C5(C(=O)C=C(C6)C)O)O)O)O)OC(=O)C=CC=CC=CC=CC(=O)NC7=C(C8=C(C(=C(C=C8)O)Cl)OC7=O)O)OC(=O)C’

    
    first_smiles_list =["CC1C=C(C(=O)C2=C3C(=CC(=C2O)C)C(=O)C4=C(C3=O)OC(CN4)(C(=O)C(C1OC(=O)C5=C(C(=CC(=N5)C(=O)NC6=C(C7=C(C=C(C=C7)OC)OC6=O)O)OC8CC(C(C(O8)C)OC)(C)[N+](=O)[O-])O)O)C)C", \
                        "CC1C=CC=C(C(=O)NC2=CC(=C3C(=C2O)C(=C(C4=C3C(=O)C(O4)(OC=CC(C(C(C(C(C(C1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)C", \
                        "CC1C(CCC(O1)OC2C(=O)C(=C(C3(C2(CC4=C(C5=C(C(=O)C=C(C5=O)OC)C(=C4C3=O)O)C)O)O)O)C(=O)C)OC6CC(C(C(O6)C)OC(=O)C7=C(C=CC(=C7O)C)C)(C)O", \
                        "CC1C(C(CC(O1)OC2C(OC3(CC2O)OC4C(OC(C(C4(O3)C)O)OC5C(C(OC(C5OC)C)OC6C(OC(C(C6O)OC)OC7C(C8C(CO7)OC9(O8)C1C(C(CO9)OC(=O)C2=C(C=C(C=C2C)O)O)OCO1)O)COC)O)C)C)OC1CC(C(C(O1)C)OC)(C)[N+](=O)[O-])OC(=O)C1=C(C(=C(C(=C1OC)Cl)O)Cl)C", \
                        "CC1C(C(CC(O1)C2=C(C3=C(C=C2)C(=O)C45C(C3O)(O4)C(CC6(C5(C(=O)C=C(C6)C)O)O)O)O)OC(=O)C=CC=CC=CC=CC(=O)NC7=C(C8=C(C(=C(C=C8)O)Cl)OC7=O)O)OC(=O)C"]

    structure = first_smiles_list[0] # Select the structure
    structure = Chem.MolFromSmiles(structure)

    m = Chem.MolFromSmiles(first_smiles_list[0]) # Select the structure
    nr_of_atoms = m.GetNumAtoms()
    
    # Generate all possible mol environments per structure
    substructures_list = []
    for i in range(nr_of_atoms):
        for j in range(nr_of_atoms):
            env = Chem.FindAtomEnvironmentOfRadiusN(m,i,j)
            substructures_list += [env]

    # Generate all possible substructures based on the mol envs.
    smile_list = []
    for env in substructures_list:
        amap={}
        submol=Chem.PathToSubmol(m, env, atomMap=amap)
        mol = Chem.MolToSmiles(submol, canonical=True)
        # Prevent overlapping substructures 
        if mol != '' and mol not in smile_list:
            smile_list += [mol]
    
    # Add the substructure to the 'all substructures list'
    sub_list = [structure]
    for smile in smile_list:
        x = Chem.MolFromSmiles(smile)
        if x != None:
             sub_list += [x]   
             
    # Draw and save (sub)structures
    multiple_molecules = Draw.MolsToGridImage(sub_list,molsPerRow=5,subImgSize=(200,100))
    multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/Simocyclinone_substructures.png')


if __name__ == '__main__':
    draw_substruc()