#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 09:08:53 2018
Script to add new structures to the NP SQlite database. The input data txt file
(example_database_with_generated_data.txt) has to be generated first with
generate_new_structure_data.py
Command line: python3 add_new_structures.py dummy_data_CLASS.txt
@author: stokm006
"""

from sys import argv
import sqlite3


def parse_data(input_file, db_name):
    """ takes all text from pre generated database file and returns an easy to 
    use list of lists with all structure data
    
    input_file: pre generated data txt file from generate_new_structure_data.py
    """

    all_lines = input_file.split('\n')
    if all_lines[-1] == '':
        all_lines = all_lines[:-1]
    
    all_info_list = []
    for line in all_lines:
        line = line.split('\t')
        line.insert(len(line), db_name)
        all_info_list += [line]
    
    all_info_list = sorted(all_info_list, key = lambda x:(x[2]))
        
    return all_info_list    

def generate_structure_table_data(data):
    """ Checks whether the structure is already present in the structure table.
    If already present; no new structure_id will be created
    If not present; a new structure_id will be created
    A dictionary with all structure data will be returned
    
    input_file: structure data list created with parse_data()
    """
    
    # Connect to sqlite database
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Retrieve all smiles and structure_id's from sqlite database
    c.execute('SELECT * from structure')
    smiles_sqlite_list = []
    structure_id_sqlite_list = []
    for row in c:
        smiles_sqlite_list += [row[3]]
        structure_id_sqlite_list += [row[0]]
    
    # For structure_id's: delete "NP_" and only keep the nr
    structure_id_list = []
    for structure_id in structure_id_sqlite_list:
        structure_id_list += [int(structure_id[6:])] #for updated database this must be 3
    
    # Get the highest structure_id nr    
    structure_id_nr = max(structure_id_list)
    
    # Dictionary that temporarily stores the structure table data
    new_structure_dict = {}  
    for line in data:
        # If structure already in database, retrieve stucture_id:
        if line[2] in smiles_sqlite_list:
            get_structure_id = 'SELECT structure_id from structure where smile = "%s";'% (line[2])
            c.execute(get_structure_id)
            rows = c.fetchall()
            for row in rows:
                row = str(row)
                row = row.lstrip("( '").rstrip("' ,)")
                structure_id_old = row
                new_structure_dict[structure_id_old] = [line[0], line[1],\
                                   line[2], line[3], line[4], line[5],\
                                   line[6], line[7], line[8], line[9],\
                                   line[10], line[11]]

        # If structure not yet in database, create new structue_id
        if line[2] not in smiles_sqlite_list:
            structure_id_nr = structure_id_nr+1
            new_structure_id = 'NP_ID_{}'.format(structure_id_nr) #for updated database this must be "NP_{}"
            new_structure_dict[new_structure_id] = [line[0], line[1], line[2],\
                               line[3], line[4], line[5], line[6], line[7],\
                               line[8], line[9], line[10], line[11]]

    # Close the connection
    conn.close()

    return new_structure_dict, structure_id_nr


def structure_table(structure_dict, structure_id_nr):
    """ adds the new structures to the structure table if the structure is not
    yet present in the structure table.
    
    input_file: new_structure_dict from generate_structure_table_data()
    """
    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor() 
 
    sql_list = [] 
    for key, value in structure_dict.items():
        key = int(key.lstrip('NP_ID_'))
        
        # If the structure is not present in the database yet, a statement will 
        # be created with all data (including new structure_id)
        if key > structure_id_nr:
            value9 = value[9].replace("'", "")
            value10 = value[10].replace("'", "")
            sql_string = "INSERT INTO structure VALUES ('%s', '%s', '%s',\
            '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"% (key, value[0],\
             value[1], value[2], value[4], value[5], value[6], value[7],\
             value[8], value9, value10)
            sql_list += [sql_string]

        
    # Add new structure data to structure table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
   
    # Commit changes
    conn.commit()

    # Close the connection
    conn.close()


def generate_structure_has_source_table_data(new_structure_dict):
    """ generates the new structures-data source combination
    if source_id is already in the db nothing happens
    if source_id is unknown a new structures-data source combination is created
    returns a dict
    
    source_id = data source identifier given by the external database where the
    new structure data is retrieved from (e.g. Super Natural II).
    structure_source_id = a number that represents a structure and data source 
    combination which will be used as primary key in the 
    structure_has_data_source table.
        
    input_file: new_structure_dict from generate_structure_table_data()
    """

    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Retrieve all structure_source_id's and source_id's from sqlite database
    c.execute('SELECT * from structure_has_data_source')
    structure_source_id_sqlite_list = []
    source_id_list = []
    
    for row in c:
        structure_source_id_sqlite_list += [row[0]]
        source_id_list += [row[2]]
    
    # Get the highest structure_source_id nr       
    structure_source_id_list = []
    for structure_source_id in structure_source_id_sqlite_list:
        structure_source_id_list += [int(structure_source_id[9:])] #for updated database this must be 7
    structure_source_id_nr = max(structure_source_id_list)
    
    # Dictionary that temporarily stores the structure_has_data_source table data
    new_structure_source_dict = {}
    # If source_id not yet in sqlite database, a new structure_source_id is 
    # created
    for key, value in new_structure_dict.items(): 
        if value[3] not in source_id_list:
            structure_source_id_nr = structure_source_id_nr+1
            new_structure_source_id = 'NP_BD_ID_{}'.format(structure_source_id_nr) # for updated database NP_SRC_
            new_structure_source_dict[new_structure_source_id] = [key, value[3],\
                                      value[11]]
    return new_structure_source_dict


def structure_source_table(structure_source_dict): 
    """ adds the new structures to the structure_has_data_source table if the 
    source_id is not yet present in the structure source table. 
        
    input_file: new_structure_source_dict from generate_structure_has_source_table_data()
    """
   
    # Statements with all structure_has_data_source data will be created
    sql_list = [] 
    for key, value in structure_source_dict.items():
        sql_string = "INSERT INTO structure_has_data_source VALUES\
        ('%s','%s', '%s', '%s');"% (key, value[0], value[1], value[2])
        sql_list += [sql_string] 
        
  
    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Add new data to structure_has_data_source table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
   
    # Commit changes
    conn.commit()

    # Close the connection
    conn.close()

def generate_data_source_table_data():
    """ updates the data source table, adds new database name if needed and its
    number of structures, or adds the new nr of structures to the old one
        
    """
    
     # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Retrieve all source_names from structure_has_data_source table    
    c.execute('SELECT * from structure_has_data_source')
    source_name_list = []
    for row in c:
        source_name_list += [row[3]]
    
    source_name_list_sorted = sorted(source_name_list)
    source_name = source_name_list_sorted[0]

    # Count all structures per data source and create a dictionary with that
    # information
    x = 0
    source_dict= {}
    source_dict[source_name] = [1]    
    for line in source_name_list_sorted:
        if source_name != line:        
            source_name = line
            x = 0
        if source_name == line:
            x += 1
            source_dict[source_name] = x
            source_name = line
    return source_dict

    
def source_table(source_dict):
    """ adds the new database name and nr of structures or updates the nr of 
    stuctures in the data_source table.
        
    input_file: new_source_dict from generate_data_source_table_data()
    """
     
    # If source_name is already in the data_source table, then the row that 
    # contains that data needs to be deleted first. The statements for deletion
    # are stored in the delete list.
    delete_list = [] 
    #  List for statements for all data_source data
    sql_list = [] 
    for key, value in source_dict.items():
        delete_string = "DELETE from data_source where source_name = ('%s');"% (key)
        delete_list +=  [delete_string]
        sql_string = "INSERT INTO data_source VALUES ('%s', %s);"% (key, str(value))
        sql_list += [sql_string]


    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
    c = conn.cursor()
    
    # Remove data from and add data to data_source table
    for i in range(len(sql_list)):
      c.execute(delete_list[i]) 
      c.execute(sql_list[i])
   
    # Commit changes
    conn.commit()

    # Close the connection
    conn.close()
    

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        db_name = argv[1]
        db_name = db_name[:-24]
        input_file = file_object.read()
        data = parse_data(input_file, db_name)
        # For structure table
        structure_dict, structure_id_nr = generate_structure_table_data(data)
        structure_table(structure_dict, structure_id_nr)
        # For structure_has_data_source_table
        structure_source_dict = generate_structure_has_source_table_data(structure_dict)
        structure_source_table(structure_source_dict)
        # For data_source table
        source_dict = generate_data_source_table_data()
        source_table(source_dict)
  
         