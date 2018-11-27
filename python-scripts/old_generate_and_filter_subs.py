#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 11:00:54 2018
Script that first generated all possible substructures out of the input
structures and then filters the substructures on mol weigth, presence of hetero
atoms, nr of rings in the substructure, abundancy and nr of the bonds which
are connected to the structure.

Command line: python3 old2_gener_and_filter_subs.py simple_structures.txt n n
@author: stokm006
"""


from __future__ import print_function
from sys import argv
import itertools as it
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw 
from rdkit import RDLogger

def generate_substructures(input_file, create_files):
    """ takes all text from input file containing the structures' smile string
    and identifier. Returns structure info list and a dictionary with all 
    possibles substructure per structure.
  
    input_file: structure txt file
    create_files: if 'y', files with (sub)structure png and txt files will be 
    created
    """
    official_subs_dict = {}
    # Create a structure list
    all_lines = input_file.split('\n')
    structure_smile_list = []
    structure_mol_list = []
    structure_combo_list = []
    for line in all_lines[0:5]:
  #  for line in all_lines[:-1]:
        line = line.split('\t')
        structure_id = line[1]
        structure_mol = Chem.MolFromSmiles(line[0])
        structure_smile = Chem.MolToSmiles(structure_mol)
        structure_smile_list += [structure_smile]
        structure_mol_list += [structure_mol]
        structure_combo_list += [[structure_smile, structure_mol, structure_id]]
    
    # TO DRAW ALL INPUT STRUCTURES
    if create_files == 'y':
       multiple_molecules = Draw.MolsToGridImage(structure_mol_list, molsPerRow=1,subImgSize=(500,500), legends = structure_smile_list)
       multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/official/non_perm_structures.png')       
       
    # Generate the mols for each structure in the class
    draw_list = []
    draw_legend_list = []
    for i, structure_info in enumerate(structure_combo_list):
        valid_sub_list = []
        valid_sub_mol_list = []
        structure_smile = structure_info[0]
        structure_mol = structure_info[1]
        structure_id = structure_info[2]
        
        nr_of_atoms = structure_mol.GetNumAtoms()

        # Generate all possible mol environments per structure
        mol_env_list = []
        for j in range(nr_of_atoms):
            for k in range(nr_of_atoms):
                env = Chem.FindAtomEnvironmentOfRadiusN(structure_mol,j,k)
                mol_env_list += [env]

        # Generate all possible substructures based on the mol envs
        for env in mol_env_list:
            submol=Chem.PathToSubmol(structure_mol, env)
            # Generate the mol of each substructure
            sub_smile = Chem.MolToSmiles(submol)
            submol = Chem.MolFromSmiles(sub_smile)
            if sub_smile != '' and sub_smile != structure_smile:
                lg = RDLogger.logger()
                lg.setLevel(RDLogger.CRITICAL)
                try:
                    Chem.SanitizeMol(submol)
                    if sub_smile not in valid_sub_list and structure_mol.HasSubstructMatch(submol) == True:
                        valid_sub_list += [sub_smile]
                        valid_sub_mol_list += [submol]
                except:
                    pass
        # Write each substructure per structure in a dictionary and also generate the draw_list

        for i, valid_substructure in enumerate(valid_sub_list):
            if valid_substructure not in draw_list:
                draw_list += [valid_sub_mol_list[i]]
                draw_legend_list += [valid_substructure]
            if structure_id in official_subs_dict:
                official_subs_dict[structure_id].append(valid_substructure)
            if structure_id not in official_subs_dict:
                official_subs_dict[structure_id] = [valid_substructure]
        if structure_id not in official_subs_dict:
            official_subs_dict[structure_id] = ['<NA>']

    # To write all possible substructures in MOLBLOCKS format
    if create_files == 'y':
        official_subs_dict_sorted = sorted(official_subs_dict)
        with open("non_perm_substructures.txt", 'w') as db_file:
            for name in official_subs_dict_sorted:
                for key in official_subs_dict.keys():
                    if key == name:
                        for value in official_subs_dict[key]:
                            db_file.write(value + ".")
                        db_file.write('\t' + key + '\n')
        # To draw all_substructures       
        multiple_molecules = Draw.MolsToGridImage(draw_list, molsPerRow=8,subImgSize=(150,100), legends = draw_legend_list)
        multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/official/non_perm_all_substructures.png')       
    
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('All possible substructures')
    nr_of_subs = 0
    for key, value in official_subs_dict.items():
        # print (key,'\t',  value)
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
  
    return structure_combo_list, official_subs_dict


def filter1(substructure_dict):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on mol weigth.
  
    substructure_dict: dict with substructures and structure identifier
    """
    filter1_dict = {}
    for key, values in substructure_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                mol_weigth = Descriptors.ExactMolWt(sub_mol)
                if mol_weigth >= 0:    # adjust this value
                    if structure_id in filter1_dict:
                            filter1_dict[structure_id].append(smile)
                    if structure_id not in filter1_dict:
                            filter1_dict[structure_id] = [smile]
        if structure_id not in filter1_dict:
            filter1_dict[structure_id] = ['<NA>']
    
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 1')
    nr_of_subs = 0
    for key, value in filter1_dict.items():
      for val in value:
          nr_of_subs += 1
    print (nr_of_subs)    
    return filter1_dict
        
def filter2(filter1_dict):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the presence of
    heteroatoms.
  
    filter1_dict: dict with substructures and structure identifier
    """
    
    filter2_dict = {}
    for key, values in filter1_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                nr_of_hetero_atoms = Lipinski.NumHeteroatoms(sub_mol)
                if nr_of_hetero_atoms >= 0:    # adjust this value
                    if structure_id in filter2_dict:
                        filter2_dict[structure_id].append(smile)
                    if structure_id not in filter2_dict:
                        filter2_dict[structure_id] = [smile]
        if structure_id not in filter2_dict:
            filter2_dict[structure_id] = ['<NA>']
  
  
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 2')
    nr_of_subs = 0
    for key, value in filter2_dict.items():
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
    return filter2_dict

def filter3(filter2_dict):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the presence of
    cyclic structures.
  
    filter2_dict: dict with substructures and structure identifier
    """
    
    filter3_dict = {}
    for key, values in filter2_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                nr_ring_matches = Chem.GetSSSR(sub_mol)
                if nr_ring_matches >= 0:    # adjust this value 
                    if structure_id in filter3_dict:
                            filter3_dict[structure_id].append(smile)
                    if structure_id not in filter3_dict:
                            filter3_dict[structure_id] = [smile]
        if structure_id not in filter3_dict:
            filter3_dict[structure_id] = ['<NA>']
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 3')
    nr_of_subs = 0
    for key, value in filter3_dict.items():
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
    return filter3_dict

def filter4(structure_info, filter3_dict):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the abundancy
    (the percentage of structures that contains the substructure)
    
    structure_info: list of lists with structure_smile, structrure_mol and 
    structure_id.
    filter3_dict: dict with substructures and structure identifier
    """
    filter4_dict = {}
    
    structure_mol_list = []
    for structure_combo in structure_info:
        structure_smile = structure_combo[0]
        structure_mol = Chem.MolFromSmiles(structure_smile)
        structure_mol_list += [structure_mol]
        
    quant_pre_dict = {}
    for key, values in filter3_dict.items():
        structure_id = key
        quant_pre_dict[structure_id] = []
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                for structure in structure_mol_list:
                    if structure.HasSubstructMatch(sub_mol) == True:
                        substruc = Chem.MolToSmiles(sub_mol)
                        if substruc not in quant_pre_dict[key]:
                            quant_pre_dict[key].append(substruc)
                            
    quantity_substructure_list = []
    for values in quant_pre_dict.values():
        for value in values:
            quantity_substructure_list += [value]
    
    # The same smiles are grouped together in the list
    if quantity_substructure_list == []:
        quantity_substructure_list = ['<NA>']
    sorted_list = sorted(quantity_substructure_list) 

    # put the unique substructure smile and the nr of structures it occurs in in dict
    x = 0
    
    quantity_substructure_dict = {}
    substruct_smile = sorted_list[0]
    quantity_substructure_dict[substruct_smile] = [1]    
    for line in sorted_list:
        if substruct_smile != line:        
            substruct_smile = line
            x = 0
        if substruct_smile == line:
            x += 1
            quantity_substructure_dict[substruct_smile] = [x]
            substruct_smile = line

    abundancy_dict = {}
    for key, value in quantity_substructure_dict.items():
        for val in value:
            abundancy = val/len(structure_mol_list)*100    
        abundancy_dict[key] = [value, abundancy]

    abundancy_list = []
    for key, value in abundancy_dict.items(): 
        if value[1] >= 0 and value[1] <= 100:    # adjust these values
            abundancy_list += [key]

    for key, values in filter3_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                if smile in abundancy_list:
                    if structure_id in filter4_dict:
                            filter4_dict[structure_id].append(smile)
                    if structure_id not in filter4_dict:
                            filter4_dict[structure_id] = [smile]
        if structure_id not in filter4_dict:
            filter4_dict[structure_id] = ['<NA>']
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 4')
    nr_of_subs = 0
    for key, value in filter4_dict.items():
      #  print (key, '\t', len(value))
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
    return filter4_dict

def filter5(structure_info, filter4_dict, draw_substructures):
    """ takes all text from dictionary with the substructures and creates a 
    dictionary with substructures which are filtered based on the
    substructures' nr of bonds which are connected to structure. Also generates
    a txt file with all generated and filtered substructures.
    
    structure_info: list of lists with structure_smile, structrure_mol and 
    structure_id.
    filter4_dict: dict with substructures and structure identifier
    draw_substructures: if 'y', the substructures will be drawn in a png file.
    """
   
    draw_list = []
    draw_legend_list = []
    filter5_dict = {}
    for structure_combo in structure_info:
        ref_dict = {}
        empty_dict = {}
        structure_smile = structure_combo[0]
        structure_id = structure_combo[2]
        input_mol = Chem.MolFromSmiles(structure_smile)
        whole_struct = input_mol.GetSubstructMatch(input_mol)

        
        # determine the nr of bonds per atom in the whole structure (reference)
        for i in range(len(whole_struct)):
            atom =  input_mol.GetAtomWithIdx(i)
            nr_of_neighbors = len([y.GetAtomicNum() for y in atom.GetNeighbors()])
            ref_dict[i] = nr_of_neighbors
            empty_dict[i] = 0
        
        # determine the nr of bonds per atom in the substructure + dummy atoms
        for key, values in filter4_dict.items():     
            if key == structure_id:
                for smile in values:
                    if smile != '<NA>':
                        nr_of_bonds_list = []
                        sub_mol = Chem.MolFromSmiles(smile)
                        if input_mol.HasSubstructMatch(sub_mol) == True:
                            matches = input_mol.GetSubstructMatches(sub_mol)
                            for match in matches:
                                submol_dict = empty_dict
                                for j, num in enumerate(match):
                                    
                                    if num in submol_dict:
                                        atom =  sub_mol.GetAtomWithIdx(j)         
                                        nr_of_neighbors = len([z.GetAtomicNum() for z in atom.GetNeighbors()])
                                        submol_dict[num] = nr_of_neighbors           

                                # determine the nr of bonds per atom in the substructure        
                                difference_dict = {}
                                for key, value in ref_dict.items():
                                    difference = int(value) - int(submol_dict[key])
                                    difference_dict[key] = difference
                                nr_of_bonds = 0                               
                                
                                for num in match:
                                    if num in difference_dict and difference_dict[num] != 0:
                                        nr_of_bonds += int(difference_dict[num])
                                nr_of_bonds_list += [nr_of_bonds]
                            min_nr_of_bonds = min(nr_of_bonds_list)
                            if min_nr_of_bonds <= 10000:    # adjust this value
                                if draw_substructures == 'y':
                                    if smile not in draw_list:
                                        draw_mol = Chem.MolFromSmiles(smile)
                                        draw_list += [draw_mol]
                                    draw_legend_list += [smile]
                                if structure_id in filter5_dict:
                                    filter5_dict[structure_id].append(smile)
                                if structure_id not in filter5_dict:
                                    filter5_dict[structure_id] = [smile]
        if structure_id in filter5_dict:
            filter5_dict[structure_id].append(structure_smile)
        
        if structure_id not in filter5_dict:
            filter5_dict[structure_id] = [structure_smile]
        
   
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 5')
    nr_of_subs = 0
    for key, value in filter5_dict.items():
       # print (key, '\t', len(value))
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
    # To write all filtered substructures in MOLBLOCKS format
    filter5_dict_sorted = sorted(filter5_dict)
 #   with open("final_nonpermutation_e10_substructures.txt", 'w') as db_file:
 #       for name in filter5_dict_sorted:
 #           for key in filter5_dict.keys():
 #               if key == name:
 #                   for value in filter5_dict[key]:
 #                      db_file.write(value + ".")
 #                   db_file.write('\t' + key + '\n')
                    
    if draw_substructures == 'y':
        # To draw all_substructures       
        multiple_molecules = Draw.MolsToGridImage(draw_list, molsPerRow=5,subImgSize=(250,150), legends = draw_legend_list)
        multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/official/nonpermutation_e10_substructures.png')       

 
if __name__ == '__main__':
    
    create_structure_files = str(argv[2])
    draw_substructures = str(argv[3])
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        structure_info, substructure_dict = generate_substructures(input_file, create_structure_files)
        
        # filter 1: mol weight
        filter1_dict = filter1(substructure_dict)
        
        # filter 2: hetero atoms
        filter2_dict = filter2(filter1_dict)
        
        # filter 3: nr. of rings
        filter3_dict = filter3(filter2_dict)
        
        # filter 4: abundancy
        filter4_dict = filter4(structure_info, filter3_dict)
        
        # filter 5: connected bonds
        filter5_dict = filter5(structure_info, filter4_dict, draw_substructures)
        
        
        
