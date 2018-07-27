#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 16:28:26 2018
This script generates all substructures (and its frequency) from a structure
and creates a histogram with that data (x: substructure frequency, 
y: substructure smile.)
command line: python3 histogram_substructure_frequency_per_structure.py
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
import matplotlib.pyplot as plt

def create_histogram():
    """ Uses the given structure and generates its substructure and the 
    substructures frequency. This data is visualized with a histogram 
    (png file).
    
    """

    first_smiles_list =['C=C(C)C1CCC(C)=CCCc2coc(c2)CC2(C)OC2C1',\
                        'c1nccc2n1ccc2','OCC=CC(=O)O','OC1C2C1CC2']
    structure = first_smiles_list[0] # Select structure
    structure = Chem.MolFromSmiles(structure)
    
    for smile in first_smiles_list:
        m = Chem.MolFromSmiles(smile)
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
            if mol != '' and mol not in smile_list:
                smile_list += [mol]
     
    # Add the substructure to the 'all substructures list'
    sub_list = []
    for smile in smile_list:
        x = Chem.MolFromSmiles(smile)
        if x != None:
             sub_list += [x]         


    nr_of_matches = 0
    sub_dict = {}
    for substructure in sub_list:
        match = structure.GetSubstructMatches(substructure)
        nr_of_matches += len(match)
        mol = Chem.MolToSmiles(substructure)
        sub_dict[mol] = len(match)

    # Create and save histogram
    fig, ax = plt.subplots(figsize=(10,5)) 
    plt.bar(list(sub_dict.keys()), sub_dict.values(), color='b')
    plt.xticks(fontsize=7, rotation=90)
    xlabel = plt.xlabel('Substructure smile')
    plt.ylabel('Substructure frequency')
    plt.title("Substructure frequency for structure XXX")
    fig.savefig('/path/to/histogram.png', bbox_extra_artists=[xlabel], bbox_inches='tight')


if __name__ == '__main__':
    create_histogram()