# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 10:46:23 2018
Takes the information from easy NP database and creates all tables for the 
structure database which can be used in mySQL workbench.
Command line: python3 create_tables_NPdata_sql.py Test_Structure_Database_File
Command line: python3 create_tables_NPdata_sql.py Structure_Database_File
@author: stokm006
"""

from __future__ import print_function
from sys import argv
from rdkit import Chem

def use_data(NP_data):
    """ takes all data from the 'easy' NP database and returns a sorted list
    with the NP data.
    
    NP_data: data from the 'easy' NP database
    """

    all_lines = NP_data.split('\n')
    all_lines = all_lines[1:-1]

    succes_list = []
    fail_list = []

    all_info = []
    all_lines = all_lines
    for line in all_lines:
        line = line.split('\t')
        if len(line) == 12: # for the empty lines in CLASS.txt files
            m = Chem.MolFromSmiles(line[2])
            m = str(m)
            if m.startswith('<'):              
                all_info = [line[0], line[1], line[2], line[3], line[4], line[5], \
                line[6], line[7], line[8], line[9], line[10], line[11]]
                succes_list += [all_info]
            if not m.startswith('<'):
                all_info = [line[0], line[1], line[2], line[3], line[4], line[5], \
                line[6], line[7], line[8], line[9], line[10], line[11]]
                fail_list += [all_info]

    all_info2 = []
    for line in succes_list:
        m = Chem.MolFromSmiles(line[2])
        n = Chem.MolToSmiles(m)
        all_info = [line[0], line[1], n, line[3], line[4], line[5], line[6],\
        line[7], line[8], line[9], line[10], line[11]]
        all_info2 += [all_info]
    
    sorted_list = sorted(all_info2, key = lambda x:(x[2]))
    
    return (sorted_list)
    
def create_structure_table(all_info_list):
    """ takes all NP data and finds unique structures and adds structure iden-
    tifiers, also the Structure_table will be created in a txt file.
    
    all_info_list: the sorted data from use_data()
    """
    
    smile = 'SMILE'
    structure_dict = {}
    structure_dict2 = {}
    for i, line in enumerate(all_info_list):
        if smile != line[2]:   
            structure_dict['NP_ID_{}'.format(len(structure_dict)+1)] = \
            [line[0], line[1], line[2], line[4], line[5], line[6], line[7],\
            line[8], line[9], line[10], line[11]]
            structure_dict2['NP_ID_{}'.format(len(structure_dict2)+1)] = \
            [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7],\
            line[8], line[9], line[10], line[11]]
            smile = line[2]
        if smile == line[2]:
            smile = line[2]
    
    # create Structure_table txt file         
    with open("/mnt/nexenta/stokm006/structure_data.txt", 'w') as db_file:
        for key, value in structure_dict.items():
            db_file.write(str(key) + '\t')
            for i in range(len(value)-1):
                db_file.write(str(value[i]) + '\t')
            db_file.write(str(value[10]) + '\n')
    
    return structure_dict2
 
def create_structure_database_table(smile_list, all_info_dict):
    """ takes all NP data and correlates the right structure identifier to the 
    right structure. This function creates the structure_database_table.
    
    smile_list: the sorted data from use_data()
    all_info_dict: the data dict created in create_structure_table()
    """
    
    # Isolate the structure and the corresponding structure identifier    
    all_info_list2 = []
    for key, value in all_info_dict.items():
        all_info_list2 += [[key, value[2]]] 
    
    smile_list2 = smile_list
    for i in range(len(all_info_list2)):
        for j, line in enumerate(smile_list):
            if line[2] == all_info_list2[i][1]:
                smile_list2[j].append(all_info_list2[i][0])
    
    # Find unique structures and put them together will all info in a dict
    y = 0
    smile = smile_list2[0][2]
    unique_smile_list = [smile_list2[0][2]]
    
    smile_dict = {}
    for line in smile_list2:
        if smile != line[2]:      
            smile = line[2]
            unique_smile_list += [smile]
            y = 0
        if smile == line[2]:
            y += 1
            smile_dict[smile] = [y]
            smile = line[2]

    smile_dict2 = smile_dict
    for line in smile_list2:
        for key in smile_dict:
             if key == line[2]:
                 smile_dict2[line[2]].append(line[12])
                 smile_dict2[line[2]].append(line[3])
                 smile_dict2[line[2]].append(line[11])               
    
    # Add the structure identifier to database identifier and database in a 
    # list
    my_list1 = []
    my_list2 = []
    for values in smile_dict2.values():
        values = values[1:]
        for value in values:
            my_list1 += [value]
    
    # Make a list of lists
    for n in range (0, len(my_list1), 3):
        my_list2 += [[my_list1[n], my_list1[n+1], my_list1[n+2]]]
    
    # Add the data to a dict and add the structure_database identifier as key.
    structure_database_dict = {}
    for line in my_list2:
        structure_database_dict['NP_DB_ID_{}'.format(len(structure_database_dict)+1)]\
        = line
    
    # create Structure_Database_table txt file  
    with open("/mnt/nexenta/stokm006/structure_database_data.txt", 'w') as db_file:
        for key, value in structure_database_dict.items():
            db_file.write(str(key) + '\t')
            for i in range(len(value)-1):
                db_file.write(str(value[i]) + '\t')
            db_file.write(str(value[2]) + '\n')
  
    return structure_database_dict
    
def create_database_table(structure_database_dict):
    """ takes the structure_database_dict and determined the nr of structures
    per database, with that information, the database_table is created. 

    strcuture_database_dict: the data dict from 
    create_Structure_Database_table()
    """
    
    # Get the databases and put them in a list
    database_list = []
    for value in structure_database_dict.values():
        database_list += [value[2]]
    
    # Checks for overlap in db and write the db and the nr of structure into
    # dict
    database_dict = {}
    quantity_list = []
    database_list2 = []
    for i in database_list:
        if i in database_dict: 
            continue 
        else:
            database_dict[i] = database_list.count(i)
            quantity_list.append(database_list.count(i))
            database_list2.append(i)
    
    # create Database_table txt file  
    with open("/mnt/nexenta/stokm006/database_data.txt", 'w') as db_file:
        for i in range(len(quantity_list)):
            db_file.write(str(database_list2[i]) + '\t')
            db_file.write(str(quantity_list[i]) + '\n')

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        all_info_list = use_data(input_file)
        all_info_dict = create_structure_table(all_info_list)
        structure_database_dict = create_structure_database_table(all_info_list, all_info_dict)
        create_database_table(structure_database_dict)