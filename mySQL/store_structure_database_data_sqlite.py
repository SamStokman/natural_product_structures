# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 08:54:51 2018
Stores the structure database data in de sqlite database file
Command line: python store_structure_database_data_sqlite.py structure_database_data.txt
Command line: python store_structure_database_data_sqlite.py /mnt/nexenta/stokm006/official_structure_database_data.txt
@author: stokm006
"""

import sqlite3
from sys import argv

def use_sqlite(data):
    all_lines = data.split('\n')
    all_lines = all_lines[:-1]
    all_info_list = []
    for line in all_lines:
        line = line.split('\t')
        all_info_list += [line]

    sql_list = [] 
    for line in all_info_list:
        sql_string = "INSERT INTO structure_has_data_source VALUES\
        ('%s','%s', '%s', '%s');"% (line[0], line[1], line[2], line[3])
        sql_list += [sql_string] 
  
    # Connect
    conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")

    c = conn.cursor()
     
    # Deletes the table, which allows creation of a newer version
    c.execute("DROP TABLE IF EXISTS structure_has_data_source")
    
    # Add table with atrribute names
    c.execute('''CREATE TABLE structure_has_data_source
              (structure_source_id text PRIMARY KEY, structure_id text,\
              source_id text, source_name text,
              FOREIGN KEY(structure_id) REFERENCES structure(structure_id),
              FOREIGN KEY(source_id) REFERENCES data_source(source_id))''')

    # Add data to table
    for i in range(len(sql_list)):
        c.execute(sql_list[i])
 
    # Commit your changes
    conn.commit()

    # Close the connection
    conn.close()

if __name__ == '__main__':
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        use_sqlite(input_file)