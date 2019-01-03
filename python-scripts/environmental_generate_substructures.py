#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 11:00:54 2018
Script that first generated all possible substructures using the 
environmental method and then filters the substructures base on mol weight, 
presence of hetero atoms, nr of rings in the substructure, abundancy and 
nr of the bonds which are connected to the rest structure.

Command line: python3 environmental_generate_substructures.py -i [input file] -m [minimum mol weight] -h [nr of hetero atoms] -r [nr of rings] -a [abundancy] -b [nr of bonds]
@author: stokm006
"""


from __future__ import print_function
import sys, getopt
import itertools as it
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw 
from rdkit import RDLogger

def generate_substructures(input_file):
    """ takes all text from input file containing the structures' smile string
    and identifier. Returns structure info list and a dictionary with all 
    possibles substructure per structure.
  
    input_file: structure txt file
    """
    
    official_subs_dict = {}
   
    with open(input_file) as file_object:
        input_file = file_object.read()  
    
    # Create a structure list
    all_lines = input_file.split('\n')
    structure_smile_list = []
    structure_mol_list = []
    structure_combo_list = []
  #  for line in all_lines[0:5]:
    for line in all_lines[:-1]:
        line = line.split('\t')
        structure_id = line[1]
        structure_mol = Chem.MolFromSmiles(line[0])
        structure_smile = Chem.MolToSmiles(structure_mol)
        structure_smile_list += [structure_smile]
        structure_mol_list += [structure_mol]
        structure_combo_list += [[structure_smile, structure_mol, structure_id]]
    


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


    official_subs_dict_sorted = sorted(official_subs_dict)
    with open("all_test_substructures.txt", 'w') as db_file:
        for name in official_subs_dict_sorted:
            for key in official_subs_dict.keys():
                if key == name:
                    value_string = ''
                    for value in official_subs_dict[key]:
                        value_string += value + "."
                    value_string = value_string[:-1]
                    db_file.write(value_string + '\t' + key + '\n')
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('All possible substructures')
    nr_of_subs = 0
    for key, value in official_subs_dict.items():
        for val in value:
            nr_of_subs += 1
    print (nr_of_subs)
  
    return structure_combo_list, official_subs_dict


def filter1(substructure_dict, input_mol_weight):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on mol weight.
  
    substructure_dict: dict with substructures and structure identifier
    input_mol_weight: int, mol weight input paramater value
    """
    filter1_dict = {}
    for key, values in substructure_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                mol_weight = Descriptors.ExactMolWt(sub_mol)
                if mol_weight >= input_mol_weight:
                    if structure_id in filter1_dict:
                            filter1_dict[structure_id].append(smile)
                    if structure_id not in filter1_dict:
                            filter1_dict[structure_id] = [smile]
        if structure_id not in filter1_dict:
            filter1_dict[structure_id] = ['<NA>']
    
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 1 [mol weight]')
    nr_of_subs = 0
    for key, value in filter1_dict.items():
      for val in value:
          if val != '<NA>':
              nr_of_subs += 1
    print (nr_of_subs)    
    return filter1_dict
        
def filter2(filter1_dict, input_nr_of_hetero_atoms):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the presence of
    heteroatoms.
  
    filter1_dict: dict with substructures and structure identifier
    input_nr_of_hetero_atoms: integer, number of hetero atoms input paramater value
    """
    
    filter2_dict = {}
    for key, values in filter1_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                nr_of_hetero_atoms = Lipinski.NumHeteroatoms(sub_mol)
                if nr_of_hetero_atoms >= input_nr_of_hetero_atoms:
                    if structure_id in filter2_dict:
                        filter2_dict[structure_id].append(smile)
                    if structure_id not in filter2_dict:
                        filter2_dict[structure_id] = [smile]
        if structure_id not in filter2_dict:
            filter2_dict[structure_id] = ['<NA>']
  
  
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 2 [hetero atoms]')
    nr_of_subs = 0
    for key, value in filter2_dict.items():
        for val in value:
            if val != '<NA>':
                nr_of_subs += 1
    print (nr_of_subs)
    return filter2_dict

def filter3(filter2_dict, input_nr_of_rings):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the presence of
    cyclic structures.
  
    filter2_dict: dict with substructures and structure identifier
    input_nr_of_rings: int, number of rings input paramater value
    """
    
    filter3_dict = {}
    for key, values in filter2_dict.items():
        structure_id = key
        for smile in values:
            if smile != '<NA>':
                sub_mol = Chem.MolFromSmiles(smile)
                nr_ring_matches = Chem.GetSSSR(sub_mol)
                if nr_ring_matches >= input_nr_of_rings:
                    if structure_id in filter3_dict:
                            filter3_dict[structure_id].append(smile)
                    if structure_id not in filter3_dict:
                            filter3_dict[structure_id] = [smile]
        if structure_id not in filter3_dict:
            filter3_dict[structure_id] = ['<NA>']
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 3 [nr of rings]')
    nr_of_subs = 0
    for key, value in filter3_dict.items():
        for val in value:
            if val != '<NA>':
                nr_of_subs += 1
    print (nr_of_subs)
    return filter3_dict

def filter4(structure_info, filter3_dict, input_abundancy):
    """ takes all text from dictionary with the substructures and returns a 
    dictionary with substructures which are filtered based on the abundancy
    (the percentage of structures that contains the substructure)
    
    structure_info: list of lists with structure_smile, structrure_mol and 
    structure_id.
    filter3_dict: dict with substructures and structure identifier
    input_abundancy: str, abundancy input paramater value
    """
    filter4_dict = {}
    
    min_abundancy, max_abundancy = input_abundancy.split('-')
    min_abundancy = int(min_abundancy)
    max_abundancy = int(max_abundancy)
    
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
        if value[1] >= min_abundancy and value[1] <= max_abundancy:
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
    print ('nr of substructures after filter 4 [abundancy]')
    nr_of_subs = 0
    for key, value in filter4_dict.items():
        for val in value:
            if val != '<NA>':
                nr_of_subs += 1
    print (nr_of_subs)
    return filter4_dict

def filter5(structure_info, filter4_dict, input_nr_of_bonds, parameter_info):
    """ takes all text from dictionary with the substructures and creates a 
    dictionary with substructures which are filtered based on the
    substructures' nr of bonds which are connected to structure. Also generates
    a txt file with all generated and filtered substructures.
    
    structure_info: list of lists with structure_smile, structrure_mol and 
    structure_id.
    filter4_dict: dict with substructures and structure identifier
    input_nr_of_bonds: int, number of connected bonds input paramater value
    parameter_info: str with all input parameter information
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
                            if min_nr_of_bonds <= input_nr_of_bonds:
                                if structure_id in filter5_dict:
                                    filter5_dict[structure_id].append(smile)
                                if structure_id not in filter5_dict:
                                    filter5_dict[structure_id] = [smile]
     #   if structure_id in filter5_dict:
     #       filter5_dict[structure_id].append(structure_smile)
        
        if structure_id not in filter5_dict:
            filter5_dict[structure_id] = ['<NA>']
        
   
    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print ('nr of substructures after filter 5 [nr of connected bonds]')
    nr_of_subs = 0
    for key, value in filter5_dict.items():
        for val in value:
            if val != '<NA>':
                nr_of_subs += 1
    print (nr_of_subs)
    
    # To write all filtered substructures in MOLBLOCKS format
    filter5_dict_sorted = sorted(filter5_dict)
    with open("final_test_output.txt", 'w') as db_file:
        
        db_file.write('Environmental method;\t' + parameter_info + '\n')
        for name in filter5_dict_sorted:
            for key in filter5_dict.keys():
                if key == name:
                    value_string = ''
                    for value in filter5_dict[key]:
                        value_string += value + "."
                    value_string = value_string[:-1]
                    db_file.write(value_string + '\t' + key + '\n')
                    
                    
def main(argv):
    """
    Main function that calls all other functions with the corresponding
    argument parameter values.
    
    argv: list with all input arguments
    """
    try:
        opts, args = getopt.getopt(argv,"i:m:h:r:a:b:")
        if len(opts) != 6:
           print ('You did not enter the right argument')
           sys.exit(2)
        opt_list = []
        for opt, arg in opts:
           if opt not in opt_list:
               opt_list += [opt]
        if '-i' not in opt_list or '-m' not in opt_list or '-h' not in opt_list or '-r' not in opt_list or '-a' not in opt_list or '-b' not in opt_list:
           print ('You did not enter the right arguments')
           sys.exit(2)      
    except getopt.GetoptError:
        print ('You did not enter the right arguments')
        sys.exit(2)
      
    for opt, arg in opts:
       if opt in ['-i']:
          input_file = arg
       if opt in ['-m']:
          mol_weight = int(arg)
       elif opt in ['-h']:
          nr_of_hetero_atoms = int(arg)
       elif opt in ['-r']:
           nr_of_rings = int(arg)
       elif opt in ['-a']:
           abundancy = arg
       elif opt in ['-b']:
           nr_of_bonds = int(arg)
           
    parameter_info = 'mol weight >= {}'.format(mol_weight) + '\t hetero atoms >= {}'.format(nr_of_hetero_atoms) + '\t rings >= {}'.format(nr_of_rings) + '\t abundany: {}'.format(abundancy) + '\t bonds <= {}'.format(nr_of_bonds)
    
    structure_info, substructure_dict = generate_substructures(input_file)
        
    # filter 1: mol weight
    filter1_dict = filter1(substructure_dict, mol_weight)
        
    # filter 2: hetero atoms
    filter2_dict = filter2(filter1_dict, nr_of_hetero_atoms)
       
    # filter 3: nr. of rings
    filter3_dict = filter3(filter2_dict, nr_of_rings)
        
    # filter 4: abundancy
    filter4_dict = filter4(structure_info, filter3_dict, abundancy)
        
    # filter 5: nr. of connected bonds
    filter5_dict = filter5(structure_info, filter4_dict, nr_of_bonds, parameter_info)
        
if __name__ == '__main__':
    
   main(sys.argv[1:])


        
        
