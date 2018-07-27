#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:37:51 2018
Norine CLASS parser which creates a CLASS format database for the Norine data.
Command line: python3 norine_CLASS_parser.py NORINE.txt
@author: stokm006
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdinchi
from sys import argv

def parse_data(input_file):
    """ takes all text from norine database file and returns a list of lists 
    with all CLASS data and an attribute list.
  
    input_file: norine database txt file
    """
    attribute_names = ['MonoisotopicMass', 'InChI', 'SMILES', 'Identifier',\
    'InChIKey2', 'InChIKey1', 'MolecularFormula', 'kingdom_name',\
    'superclass_name', 'class_name', 'subclass_name']
        
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
    all_lines = all_lines[2:]
    for line in all_lines:
        line = line.split('\t')
        
        #Convert to mol and remove invalid structures 
        smile_string = ''
        id_string = ''
        m = line[2]
        id_name = line[0]
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
        #SMILES
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
    
    return  attribute_names, overall_list  
    
def write_CLASS_db(file_name, attribute_names, data):
    """ takes list of lists with data and attributes list from farse_data() and
    writes a CLASS format database.
  
    file_name: the file name of the txt file which will be written
    attribute_names: attribute names created in parse_data()
    data: all data parsed and created in parse_data()
    """
    
    with open(file_name, 'w') as db_file:
        for i in range(len(attribute_names)-1):
            db_file.write(str(attribute_names[i]) + '\t')
        db_file.write(str(attribute_names[10]))
        db_file.write('\n')
        
        data2 = list(map(list, zip(*data)))
        
        for line in (data2):
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
                db_file.write(str(line[10]))
                db_file.write('\n')


if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        att_names, data_list = parse_data(input_file)
        write_CLASS_db("norineCLASS.txt", att_names, data_list)
        
