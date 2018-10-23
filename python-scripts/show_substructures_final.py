#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 09:11:22 2018
This script shows the highlighted substructure, as input the molblock input and
output are needed, also you have to specify which structure you want to
visualize (int)

Command line: python3 show_substructures.py pyra_structures.txt filtered_pyra_substructures.txt 3
Command line: python3 show_substructures.py boro_structures.txt filtered_boro_substructures.txt 3
Command line: python3 show_substructures.py 2_flav_structures.txt filtered_2flav_cycle_1_substructures.txt 3 n

@author: stokm006
"""
from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw    
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG = True


def visualize_substructures(input_data, result_data, structure_nr, boolean):

    input_data = input_data.split('\n')
    input_data = input_data[structure_nr:structure_nr+2]
    input_list = []
    for row in input_data:
        row = row.split('\t')
        input_list += [row[0]]
    all_lines = result_data.split('\n')
    all_lines = all_lines[structure_nr:structure_nr+2]
    all_info_input_list = []
    all_info_output_list = []
    

    for x, line in enumerate(all_lines):
        line = line.split('\t')
        smile = line[0]
        smile = smile.split(".")
        nr_of_substructures = len(smile)
        input_molecule = []
        input_molecule = nr_of_substructures * [input_list[x]]
        input_mol_list = []
        for m in input_molecule:
            m = Chem.MolFromSmiles(m)
            input_mol_list += [m]
        all_info_input_list += [input_mol_list]
 
        output_list = []
        for i in range(len(smile)):
            smile_str = smile[i]
            output_list += [smile_str]
        mol_result_list = []
        for smile in output_list:
            m2 = Chem.MolFromSmiles(smile, sanitize=False)
            mol_result_list += [m2]
        all_info_output_list += [mol_result_list]
    

        
    # just show the substructures    
    if boolean == 'n':

        im = Draw.MolsToGridImage(all_info_output_list[0], 
                                  useSVG=False,
                                  molsPerRow=4,
                                  subImgSize = (300,200))
        im.show()
  
    # shows the (highlighted) substructures in structure (does not work for every substructure)  
    if boolean == 'y':
        [Chem.GetSSSR(mol) for mol in all_info_output_list[0]]
        [AllChem.Compute2DCoords(mol) for mol in all_info_output_list[0]]
     
        mol_input = all_info_input_list[0][0] 

        
        highlight_lists = []
        for mol in all_info_output_list[0]: 
            AllChem.GenerateDepictionMatching2DStructure(mol_input, mol, acceptFailure = True)
            highlight_lists += [mol_input.GetSubstructMatch(mol)]

        input_molecule = []
        for i in range(len(highlight_lists)): 
            input_molecule += [mol_input]

        im = Draw.MolsToGridImage(input_molecule, 
                                  highlightAtomLists = highlight_lists,
                                  useSVG=False,
                                  molsPerRow=2,
                                  subImgSize = (600,400))
        im.show()    
    
    # shows all (highlighted) substructures in structure (does not work for every substructure)  
    if boolean == 'a':
        [Chem.GetSSSR(mol) for mol in all_info_output_list[0]]
        [AllChem.Compute2DCoords(mol) for mol in all_info_output_list[0]]
     
        mol_input = all_info_input_list[0][0] 

        
        highlight_lists = []
        for mol in all_info_output_list[0]: 
            AllChem.GenerateDepictionMatching2DStructure(mol_input, mol, acceptFailure = True)
            if len(mol_input.GetSubstructMatches(mol)) == 1:
                highlight_lists += [mol_input.GetSubstructMatch(mol)]
            if len(mol_input.GetSubstructMatches(mol)) > 1:
                matches = mol_input.GetSubstructMatches(mol)
                for match in matches:
                    highlight_lists += [match]
                               
        input_molecule = []
        for i in range(len(highlight_lists)): 
            input_molecule += [mol_input]

        im = Draw.MolsToGridImage(input_molecule, 
                                  highlightAtomLists = highlight_lists,
                                  useSVG=False,
                                  molsPerRow=2,
                                  subImgSize = (600,400))
        im.show()


if __name__ == "__main__":
    
    structure_nr = int(argv[3])
    boolean = str(argv[4])
    
    with open(argv[1]) as file_object:
        input_file = file_object.read()
    
    with open(argv[2]) as file_object:
        result_file = file_object.read()
        visualize_substructures(input_file, result_file, structure_nr, boolean)
