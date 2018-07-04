# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 08:07:13 2018
Streptomedb CLASS parser (structure with smile 'O=C1NC(=O)c2c1c1c3ccccc3cnc1c
1c2c2ccccc2[n-]1.CCOC(=O)C1=C(C)[CH-]C(=C1)C.C$O.[Ru+2] is removed, as it was
not recognized). This script creates a CLASS format database for
the streptomedb data.

Command line: python3 streptomedb_CLASS_parser.py streptomedb-2try.sdf
@author: stokm006
"""

from __future__ import print_function
from sys import argv
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdinchi    
from rdkit.Chem import rdMolDescriptors


def make_dict(input_file):
    """ takes all text from streptomedb file and returns a dictionary with
    NPs which is easy to use
    
    input_file: streptomedb sdf file
    """    
    my_string = ""
    all_lines = input_file.split('\n')
    for lines in all_lines:
        if not lines.startswith ('$$$$'):
            my_string += lines.strip('\n')
        if lines.startswith ('$$$$'):
            my_string += ('|')
    my_string = my_string.split('|')
    my_string = my_string[:-1]
    
    # remove 'IC50 >' to allow splitting on '>' later on
    all_info_list = []
    for line in my_string:
        line = line.split('END>')
        line = line[1].split('IC50 >')
        list_1 = []
        if len(line) != 1:
            list_2 = ""
            for i in range (len(line)):
                list_2 += line[i]
            list_1 += [list_2]
        if len(line) == 1:
            list_1 += line
        all_info_list += list_1    

    all_info_list2 = []
    for info_line in all_info_list:
        info_line = info_line.split('>')
        list_3 = []
        for value in info_line:
            my_string = ""
            value = value.strip('\'"')
            value = value.lstrip('  <')
            if len(value) == 0:
                 value = "NA"
            my_string += value
            list_3 += [my_string]
        all_info_list2 += [list_3]
           
    column_name_list = ['canonical_smiles', 'name', 'compound_id', 'molwt', \
    'HBD', 'HBA', 'rotatable bonds', 'pubchem cid', 'organism', 'pmids', \
    'activities', 'synthesizing routes']
            
    
    column_list1 = []    
    column_list2 = []
    column_list3 = []
    column_list4 = []
    column_list5 = []
    column_list6 = []
    column_list7 = []
    column_list8 = []
    column_list9 = []
    column_list10 = []
    column_list11 = []
    column_list12 = []
    
    strep_dict = {}
    for lists in all_info_list2:
        if len(lists) == 24:
            column_list1 += [lists[1]]        
            column_list2 += [lists[3]]
            column_list3 += [lists[5]]
            column_list4 += [lists[7]]        
            column_list5 += [lists[9]]
            column_list6 += [lists[11]]
            column_list7 += [lists[13]]        
            column_list8 += [lists[15]]
            column_list9 += [lists[17]]
            column_list10 += [lists[19]]
            column_list11 += [lists[21]]
            column_list12 += [lists[22]]
            
        if len(lists) == 22: # one structure has no synthesizing routes info
            column_list1 += [lists[1]]        
            column_list2 += [lists[3]]
            column_list3 += [lists[5]]
            column_list4 += [lists[7]]        
            column_list5 += [lists[9]]
            column_list6 += [lists[11]]
            column_list7 += [lists[13]]        
            column_list8 += [lists[15]]
            column_list9 += [lists[17]]
            column_list10 += [lists[19]]
            column_list11 += [lists[21]]
            column_list12 += ['NA']
     
        strep_dict[column_name_list[0]] = column_list1
        strep_dict[column_name_list[1]] = column_list2
        strep_dict[column_name_list[2]] = column_list3
        strep_dict[column_name_list[3]] = column_list4
        strep_dict[column_name_list[4]] = column_list5
        strep_dict[column_name_list[5]] = column_list6
        strep_dict[column_name_list[6]] = column_list7
        strep_dict[column_name_list[7]] = column_list8
        strep_dict[column_name_list[8]] = column_list9
        strep_dict[column_name_list[9]] = column_list10
        strep_dict[column_name_list[10]] = column_list11
        strep_dict[column_name_list[11]] = column_list12
    return (strep_dict)


def create_CLASS_data(data_dict):
    """ Generates CLASS data for the strepto data present in the strep_dict.
    
    input_file: streptodb dictionary
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
    
    # Identifier
    identifier_list = data_dict['compound_id']
    
    # SMILES
    SMILES_list = data_dict['canonical_smiles']
    
    for SMILE in SMILES_list:
        # generate molecules
        m = Chem.MolFromSmiles(SMILE)
    
        # MonoisotopicMass
        mol_mass = str(Descriptors.ExactMolWt(m))[:-3]
        mol_mass_list += [mol_mass]

        # InChI         
        inchi = rdinchi.MolToInchi(m)
        inchi_list += [inchi[0]]   
        
        # MolecularFormula
        mol_formula = rdMolDescriptors.CalcMolFormula(m)
        mol_formula_list += [mol_formula]

    # NA list    
    nr_of_structures = len(data_dict['canonical_smiles'])
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
    """ takes list of lists with data and attributes list from create_CLASS_\
    data and writes a ClASS format database.
  
    file_name: the file name of the txt file which will be written
    attribute_names: attribute names created in create_CLASS_data()
    data: all data parsed and created in create_CLASS_data()
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
        my_dictionary = make_dict(input_file)
        att_names, data_list = create_CLASS_data(my_dictionary)
        write_CLASS_db("streptomedbCLASS.txt", att_names, data_list)
    
