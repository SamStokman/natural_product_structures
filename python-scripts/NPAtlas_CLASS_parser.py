#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 13:41:03 2018
NPAtlas CLASS parser which creates a CLASS format database for the NPAtlas data.
this script parses the NPAtlas_DB into the uniform CLASS one.
Command line: python3 parse_NPAtlas.py NPAtlas_DB_last_version.tsv
@author: stokm006
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from sys import argv

def parse_data(input_file):
    """ takes all text from NPAtlas database file and returns a list of lists 
    with all CLASS data and an attribute list.
  
    input_file: NPAtlas database txt file
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
    
    all_lines = input_file.split('\n')
    all_lines = all_lines[1:-1]
    for line in  all_lines:
        line = line.split('\t')
        
        #SMILE
        m = line[1]
        m = Chem.MolFromSmiles(m)
        sm = Chem.MolToSmiles(m)
        SMILES_list += [sm]
        #Monoisotopic mass
        mol_weigth = Descriptors.ExactMolWt(m)
        mol_mass_list += [mol_weigth]
        #Source identifiers
        identifier_list += [line[0]]
        #Inchi
        inchi_list += [line[2]]
        #InchiKeys
        inchi_key = line[3].split('-')
        inchi_key2 = inchi_key[1]
        inchi_key2_list += [inchi_key2]
        inchi_key1 = inchi_key[0]
        inchi_key1_list += [inchi_key1]
        #Mol Forumula
        mol_formula = rdMolDescriptors.CalcMolFormula(m)
        mol_formula_list += [mol_formula]

    # NA list    
    nr_of_structures = len(SMILES_list)
    NA_list += ['NA'] * nr_of_structures


    overall_list = [mol_mass_list]+[inchi_list]+[SMILES_list]+\
    [identifier_list]+[inchi_key2_list]+[inchi_key1_list]+[mol_formula_list]+\
    [NA_list]+[NA_list]+[NA_list]+[NA_list]
    
    return  attribute_names, overall_list  
    
def write_CLASS_db(file_name, attribute_names, data):
    """ takes list of lists with data and attributes list from farse_file() and
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
        write_CLASS_db("NPAtlasCLASS.txt", att_names, data_list)
        