#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:37:51 2018
This script generates the mol_mass, inchi, inchi_key1, inchi_key2, mol_formula,
from a given smile and creates a text file which can be used as input to add
new structures to the NP sqlite database.

The name of the input txt file (minus “.txt”) will be used as the source 
name later on in the sqlite database.

The input file needs to be a text file with the smile and source identifier 
which are tab-seperated, then each pair has to be seperated by a new line.
Of cource this script can be costumized

Command line: python3 generate_new_structure_data.py dummydb_input.txt
@author: stokm006
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdinchi
from sys import argv

def generate_data(input_file):
    """ takes all text from the input structure data file and returns a list of
    lists with all generated data needed for the sqlite database.
  
    input_file: input structure txt file
    """
       
    mol_mass_list = []
    inchi_list = []
    SMILES_list = []
    identifier_list = []
    inchi_key1_list = [] 
    inchi_key2_list = [] 
    mol_formula_list = []
    NA_list = []
    
    pre_SMILES_list = []
    identifier_list = []
    all_lines = input_file.split('\n')
    if all_lines[-1] == '':
        all_lines = all_lines[:-1]
    for line in all_lines:
        line = line.split('\t')

        #Convert to mol and remove invalid structures 
        smile_string = ''
        id_string = ''
        m = line[0]
        id_name = line[1]
        mol = Chem.MolFromSmiles(m)
        if mol != None:
            smile_string += m
            id_string += id_name
        pre_SMILES_list += [smile_string]
        
        #Source identifiers
        identifier_list += [id_string]
        
    pre_inchi_list = []
    for smile in pre_SMILES_list:
        #Generate mol
        m = Chem.MolFromSmiles(smile)
        #SMILES, canonical
        sm = Chem.MolToSmiles(m)
        SMILES_list += [sm]
        #Monoisotopic mass
        mol_weigth = Descriptors.ExactMolWt(m)
        mol_mass_list += [mol_weigth]
        #Mol Forumula
        mol_formula = rdMolDescriptors.CalcMolFormula(m)
        mol_formula_list += [mol_formula]
        # InChI         
        inchi = rdinchi.MolToInchi(m)
        pre_inchi_list += [inchi[0]]  
    
   
    # InChIKey1 and InChIKey2
    for inchi in pre_inchi_list:
        if not str(inchi).startswith('InCh'):
            inchi = 'NA'
        inchi_list += [inchi]
    
    pre_inchi_key_list =[]
    for inchi2 in inchi_list: 
        if inchi2 == 'NA':
            inchi_key = "NA-NA"
            pre_inchi_key_list += [inchi_key]
        if inchi2 != 'NA':
            inchi_key = rdinchi.InchiToInchiKey(inchi2)
            pre_inchi_key_list += [inchi_key]
    
    for inchi_key in pre_inchi_key_list:
        inchi_key = inchi_key.split('-')
        inchi_key2 = inchi_key[1]
        inchi_key2_list += [inchi_key2]
        inchi_key1 = inchi_key[0]
        inchi_key1_list += [inchi_key1]

    # NA list    
    nr_of_structures = len(SMILES_list)
    NA_list += ['NA'] * nr_of_structures

    overall_list = [mol_mass_list]+[inchi_list]+[SMILES_list]+\
    [identifier_list]+[inchi_key2_list]+[inchi_key1_list]+[mol_formula_list]+\
    [NA_list]+[NA_list]+[NA_list]+[NA_list]
    
    return  overall_list  
     
def write_CLASS_db(file_name, data):
    """ takes list of lists with data from generate_data() and writes a txt
    file with all data which can be used to add new structures to the sqlite 
    database.
  
    file_name: name of the input txt file
    data: all data parsed and created in parse_data()
    """
 
    file_name = file_name[:-4] + "_with_generated_data.txt"
    
    # Create a txt file with all structure data
    with open(file_name, 'w') as db_file:
        data2 = list(map(list, zip(*data))) # to 'rotate' the list of lists
        for line in data2:
            if not str(line[0]).startswith('0'):
                db_file.write(str(line[0])+ '\t')
                db_file.write(str(line[1])+ '\t')
                db_file.write(str(line[2])+ '\t')
                db_file.write(str(line[3])+ '\t')
                db_file.write(str(line[4])+ '\t')
                db_file.write(str(line[5])+ '\t')
                db_file.write(str(line[6])+ '\t')
                db_file.write(str(line[7])+ '\t')
                db_file.write(str(line[8])+ '\t')
                db_file.write(str(line[9])+ '\t')
                db_file.write(str(line[10])+ '\n')


if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        txt_file_name = argv[1]
        data_list = generate_data(input_file)
        write_CLASS_db(txt_file_name, data_list)
        