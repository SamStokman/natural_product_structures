# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 11:17:39 2018
nanpdb CLASS parser, structures 78 and 2297 are removed from the orginal db, as
they could not be converted with their SMILE. This script creates a CLASS 
format database for the nanpdb data.
Command line: python3 nanpdb_CLASS_parser.py nanpdbDATA.txt
@author: stokm006
"""

from __future__ import print_function
from sys import argv
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdinchi
from rdkit.Chem import rdMolDescriptors

def parse_file(input_file, db_name):
    """ takes all text from nanpdb database file and returns a list of lists 
    with NPs which is easy to use.
  
    input_file: nanpdb database txt file
    db_name: database name
    """
    
    all_lines = input_file.split('\n')
    all_lines = all_lines[:-1]
    all_info_list = []
    for line in all_lines:
        line = line.split('\t')
        info_per_row_list = []
        for value in line:
            my_string = ""
            if len(value) == 0:
                value = "NA"
            my_string += value
            info_per_row_list += [my_string]
        info_per_row_list += [db_name]    
        all_info_list += [info_per_row_list]
    
        
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
    
    for line in all_info_list:
        # generate molecules
        m = Chem.MolFromSmiles(line[0])
        
        # MonoisotopicMass
        mol_mass = str(Descriptors.ExactMolWt(m))[:-5]
        mol_mass_list += [mol_mass]
        
        # InChI         
        inchi = rdinchi.MolToInchi(m)
        inchi_list += [inchi[0]]        
        
        # SMILES
        SMILES_list += [line[0]]
        
        # Identifier
        identifier_list += [line[1]]
       
        # MolecularFormula
        mol_formula = rdMolDescriptors.CalcMolFormula(m)
        mol_formula_list += [mol_formula]
        
        # NA list    
    nr_of_structures = len(all_info_list)
    NA_list += ['NA'] * nr_of_structures  
           
            # InChIKey
    inchi_key_list = [] 
    inchi_key_list2 = []    
    for inchi in inchi_list:
        inchi_key = rdinchi.InchiToInchiKey(inchi)
        inchi_key_list2 += [inchi_key]
    inchi_key_list += inchi_key_list2
    
        # InChiKey1 and InChiKey2
    for inchikey in inchi_key_list:
        inchikey = inchikey.split('-')
        inchikey1 = inchikey[0]
        inchikey2 = inchikey[1]
        inchi_key1_list += [inchikey1]
        inchi_key2_list += [inchikey2]
        
        
    overall_list = [mol_mass_list]+[inchi_list]+[SMILES_list]+\
    [identifier_list]+[inchi_key2_list]+[inchi_key1_list]+[mol_formula_list]+\
    [NA_list]+[NA_list]+[NA_list]+[NA_list]
    
    return  attribute_names, overall_list  
    
def write_CLASS_db(file_name, attribute_names, data):
    """ takes list of lists with data and attributes list from farse_file() and
    writes a ClASS format database.
  
    file_name: the file name of the txt file which will be written
    attribute_names: attribute names created in parse_file()
    data: all data parsed and created in parse_file()
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
        db_name = argv[1]
        db_name = db_name[:-4]
        input_file = file_object.read()
        att_names, data_list = parse_file(input_file, db_name)
        write_CLASS_db("nanpdbCLASS.txt", att_names, data_list)
        
