#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 09:12:53 2018
Filter that checks several rules
Command line: python3 generate_substructures.py pyra_structures.txt
Command line: python3 generate_substructures.py boro_structures.txt
Command line: python3 generate_substructures.py indo_structures.txt
Command line: python3 generate_substructures.py 2_flav_structures.txt
Command line: python3 generate_substructures.py flavo2_structures.txt
@author: stokm006
"""


from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw 
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

def generate_substructures(input_file):
    
    # Create a structure list
    all_lines = input_file.split('\n')
    structure_list = []
    structure_mol_list = []
    structure_combo_list = []
    for line in all_lines[:-1]:
        line = line.split('\t')
        structure_list += [line[0]]
        mol = Chem.MolFromSmiles(line[0])
        structure_mol_list += [mol]
        structure_combo_list += [[mol, line[1]]]

    sorted_list = sorted(structure_combo_list, key = lambda x:(x[1]))

  #  print (structure_combo_list)

  #  multiple_molecules = Draw.MolsToGridImage(structure_mol_list, molsPerRow=10,subImgSize=(250,150), legends = structure_list)
  #  multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/indo_names_structures.png')       
       
      

    # Generate the mols for each structure in the class
    substructure_list = []
    for mol in structure_mol_list:
        nr_of_atoms = mol.GetNumAtoms()
        

        # Generate all possible mol environments per structure
        mol_env_list = []
        for i in range(nr_of_atoms):
            for j in range(nr_of_atoms):
                env = Chem.FindAtomEnvironmentOfRadiusN(mol,i,j)
                mol_env_list += [env]


        # Generate all possible substructures based on the mol envs.
        substructure_list_temp = []
        for env in mol_env_list:
            amap={}
            submol=Chem.PathToSubmol(mol, env, atomMap=amap)
            # Generate the mol of each substructure
            mol1 = Chem.MolToSmiles(submol, canonical=True)

            # Prevent overlapping substructures within the same structure
            # (the nr of generated substructures here is not correct yet)
            if mol1 != '' and mol1 not in substructure_list_temp:
                substructure_list_temp += [mol1]
        # Add the substructure to the 'all substructures list'
        for smile in substructure_list_temp:
            # Prevent overlapping substructure for all structures 
            if smile not in substructure_list:
                substructure_list += [smile]


    # Convert all substructure smiles to mols
    substructure_mol_list = []
    for sub_struc in substructure_list:
        sm = Chem.MolFromSmiles(sub_struc, sanitize=False)
        if sm != None:
            substructure_mol_list += [sm]
 
       
    all_subs_dict= {}
    for structure_comb in structure_combo_list:
        
        for substructure in substructure_mol_list:
            # Match the correct substructure to the correct structure
            structure = structure_comb[0]
            name = structure_comb[1]
            if structure.HasSubstructMatch(substructure) == True:
                substruc = Chem.MolToSmiles(substructure)
            #    struc = Chem.MolToSmiles(structure)
                if name in all_subs_dict:
                    all_subs_dict[name].append(substruc)
                else:
                    all_subs_dict[name] = [substruc]

 #   all_subs_dict_sorted = sorted(all_subs_dict)

   # werkt 
 #   with open("all_flavo_substructures.txt", 'w') as db_file:
 #      for key in all_subs_dict_sorted:         
 #         for i in range(0, len(all_subs_dict[key])-1, 1):
 #             db_file.write(all_subs_dict[key][i-1] + ".")
 #         db_file.write(all_subs_dict[key][len(all_subs_dict[key])-1])
 #         db_file.write('\t' + key + '\n')
   
  #  print (len(substructure_mol_list))
  #  multiple_molecules = Draw.MolsToGridImage(substructure_mol_list, molsPerRow=25,subImgSize=(50,50))#, legends = substructure_list)
  #  multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/indo7_all_substructures.png')       
     
  
    # Filter; Oxygens, nitrogens and halogens
    # And mol weigth, NOTE: hydrogen ignored

    chlorine_check_list2 = []
    chlorine_check = substructure_list[0][1]
    bromine_check_list2 = []
    bromine_check = substructure_list[0][1]
    for smile in substructure_list:
        chlorine_check_list = []
        bromine_check_list = []
        for i in range(0, len(smile)-1, 1):
            chlorine_check = smile[i+1]
            chlorine_check_list +=[chlorine_check]
            bromine_check = smile[i+1]
            bromine_check_list +=[bromine_check]
        chlorine_check_list.insert(-1, 'X')
        chlorine_check_list2 += [chlorine_check_list]
        bromine_check_list.insert(-1, 'X')
        bromine_check_list2 += [bromine_check_list]

    filter1_list = []
    for i in range(len(substructure_list)):
        mw = 0
        c_count = 0
        n_count = 0
        o_count = 0
        Cl_count = 0
        F_count = 0
        Br_count = 0
        I_count = 0
        At_count = 0
        Zr_count = 0
        B_count = 0
        
        if substructure_list[i][-1] == 'l': # for Cl at the end of the smile
            c_count += -1
            
        for j in range(len(substructure_list[i])):
            if substructure_list[i][j] == 'c':
                c_count += 1
                mw += 12.0107
            if substructure_list[i][j] == 'C' and chlorine_check_list2[i][j] != 'l': # Ignore Cl atoms
                c_count += 1
                mw += 12.0107
            if substructure_list[i][j] == 'n':
                n_count += 1
                mw += 14.0067
            if substructure_list[i][j] == 'N':
                n_count += 1
                mw += 14.0067
            if substructure_list[i][j] == 'o':
                o_count += 1
                mw += 15.999
            if substructure_list[i][j] == 'O':
                o_count += 1
                mw += 15.999
            if substructure_list[i][j] == 'I':
                I_count += 1   
                mw += 126.905
            if substructure_list[i][j] == 'F':
                F_count += 1
                mw += 18.998
            if substructure_list[i][j] == 'B' and bromine_check_list2[i][j] != 'r': # Ignore Br atoms
                B_count += 1
                mw += 10.811

        for k in range(0, len(substructure_list[i]), 2):
            two_atoms = substructure_list[i][k:k+2]       
            if two_atoms == 'Cl':
                Cl_count += 1
                mw += 35.453
            if two_atoms == 'Br':
                Br_count += 1
                mw += 79.904
            if two_atoms == 'At':
                At_count += 1
                mw += 210.000
            if two_atoms == 'Zr':
                Zr_count += 1
                mw += 91.224

            
        
        for k in range(1, len(substructure_list[i]), 2):
            two_atoms = substructure_list[i][k:k+2] 
            if two_atoms == 'Cl':
                Cl_count += 1
                mw += 35.453
            if two_atoms == 'Br':
                Br_count += 1
                mw += 79.904
            if two_atoms == 'At':
                At_count += 1
                mw += 210.000
            if two_atoms == 'Zr':
                Zr_count += 1
                mw += 91.224

   #     print (substructure_list[i])
   #     print ('nr. of Cs   ',c_count)
   #     print ('nr. of Cls   ',Cl_count)
   #    print ('nr. of Ns   ',n_count)
   #     print ('nr. of Os   ',o_count)
   #     print ('nr. of Brs   ',Br_count)
   #     print ('nr. of Bs   ',B_count)
   #     print ('nr. of Zrs   ',Zr_count)
   #     print ('nr. of Fls   ',F_count)
   #     print ('mol weigth  ', mw)
            
      
        if mw > 40:
            filter1_list += [substructure_list[i]]
            
        if substructure_list[i] not in filter1_list:
            if n_count >= 1 or o_count >= 1 or Cl_count >= 1 or F_count >= 1 \
            or Br_count >= 1 or I_count >= 1 or At_count >= 1 or B_count >= 1 \
            or Zr_count >= 1:
                filter1_list += [substructure_list[i]]


    # Filter2: for abundancy (substructure match)
    # and ring

    substructure_mol_list = []
    for sub_struc in filter1_list:
        sm = Chem.MolFromSmiles(sub_struc, sanitize=False)
        substructure_mol_list += [sm]
 

    
    quanitity_substructure_list = []
    for structure in structure_mol_list:
        for substructure in substructure_mol_list:
            # Match the correct substructure to the correct structure
            if structure.HasSubstructMatch(substructure) == True:
                substruc = Chem.MolToSmiles(substructure)
                quanitity_substructure_list += [substruc] 
   
    
    # The same smiles are grouped together in the list
    sorted_list = sorted(quanitity_substructure_list) 

    # put the unique substructure smart and the nr of structures it occurs in in dict
    x = 0
    
    quanitity_substructure_dict = {}
    substruct_smart = sorted_list[0]

    quanitity_substructure_dict[substruct_smart] = [1]    
    for line in sorted_list:
        if substruct_smart != line:        
            substruct_smart = line
            x = 0
        if substruct_smart == line:
            x += 1
            quanitity_substructure_dict[substruct_smart] = [x]
            substruct_smart = line

    # abundancy = percentage of structures where the substructure occurs in
    abundancy_dict = {}
    for key, value in quanitity_substructure_dict.items():
        for val in value:
            abundancy = val/len(structure_list)*100    
        abundancy_dict[key] = [value, abundancy]

    filter2_list = []
    filter2_mol_list = []
    for key, value in abundancy_dict.items():
        if value[1] >= 20 and value[1] <= 80:    # adjust these values
            key_mol = Chem.MolFromSmarts(key)
            nr_ring_matches = Chem.GetSSSR(key_mol)
            if nr_ring_matches >= 1:
                filter2_list += [key]
                filter2_mol_list += [key_mol]
     
               
    filtered_subs_dict= {}
    for structure_comb in structure_combo_list:
        for substructure in filter2_mol_list:
            # Match the correct substructure to the correct structure
            structure = structure_comb[0]
            name = structure_comb[1]
            if structure.HasSubstructMatch(substructure) == True:
                substruc = Chem.MolToSmiles(substructure)
            #    struc = Chem.MolToSmiles(structure)
                if name in filtered_subs_dict:
                    filtered_subs_dict[name].append(substruc)
                else:
                    filtered_subs_dict[name] = [substruc]
   
    for structure_comb in structure_combo_list:
        name = structure_comb[1]
        if name not in filtered_subs_dict:
            filtered_subs_dict[name] = ['<NA>']
            
        
    filtered_subs_sorted = sorted(filtered_subs_dict)

   # werkt 
    with open("filtered_2flav_cycle_1_substructures.txt", 'w') as db_file:
       for key in filtered_subs_sorted:       
          for i in range(0, len(filtered_subs_dict[key])-1, 1):
              db_file.write(filtered_subs_dict[key][i] + ".")
          db_file.write(filtered_subs_dict[key][len(filtered_subs_dict[key])-1])
          db_file.write('\t' + key + '\n')

   
    # To draw molecules    
    legends_list = []
    draw_mol_list = []
    for smart in filter2_list:
        mol_smart = Chem.MolFromSmarts(smart)
        smile = Chem.MolToSmiles(mol_smart)
        legends_list += [smile]
        mol_smile = Chem.MolFromSmiles(smile, sanitize=False)
        draw_mol_list += [mol_smile]

    

  #  print (legends_list)
    
  #  multiple_molecules = Draw.MolsToGridImage(draw_mol_list, molsPerRow=5,subImgSize=(200,150), legends = legends_list)
  #  multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/filtered_indo2_substructures.png')       

if __name__ == '__main__':
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        generate_substructures(input_file)
        
# =============================================================================
#         
#                 # nr_of_matches: multiple substructures per structures meegerekend
#                 nr_of_matches = len(match)                
#                 if nr_of_matches == 1:
#                    struc = Chem.MolToSmiles(structure)
#                    quanitity_structure_list += [struc]
#                    substruc = Chem.MolToSmiles(substructure)
#                    quanitity_substructure_list += [substruc] 
#                 
#                 if nr_of_matches >= 2:
#                    struc = Chem.MolToSmiles(structure)
#                    quanitity_structure_list += [struc] * nr_of_matches
#                    substruc = Chem.MolToSmiles(substructure)
#                    quanitity_substructure_list += [substruc] * nr_of_matches
#         

 #   multiple_molecules = Draw.MolsToGridImage(substructure_mol_list, molsPerRow=5,subImgSize=(200,150), legends = substructure_list)
 #   multiple_molecules.save('/mnt/scratch/stokm006/generate_substructures/filter/filtered_substructures.png')       
       
        
  
# =============================================================================
#     similarity_list = []
#     i_list = []
#     j_list = []
#     fps = [FingerprintMols.FingerprintMol(xx) for xx in filter2_mol_list]
#     for i in range(len(filter2_mol_list)):
#         for j in range(len(filter2_mol_list)):         
#             if DataStructs.FingerprintSimilarity(fps[i],fps[j]) >= 0.9:
#                 print (DataStructs.FingerprintSimilarity(fps[i],fps[j]))
#                 similarity_list += [DataStructs.FingerprintSimilarity(fps[i],fps[j])]
#                 if filter2_mol_list[i] not in i_list:
#                     i_list += [filter2_mol_list[i]]
#                 if filter2_mol_list[j] not in j_list:
#                     j_list += [filter2_mol_list[j]]
#         
#     for value in similarity_list[0:10]:
#         print (value)
#    
#     for smart in i_list[0:10]:
#         sm = Chem.MolToSmarts(smart)
#         print (sm)
# 
#     for smart in j_list[0:10]:
#         sm = Chem.MolToSmarts(smart)
#         print (sm)
#    
# =============================================================================
    
