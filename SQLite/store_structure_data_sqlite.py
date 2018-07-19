# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:25:40 2018
Stores the structure data in the sqlite database file.
Command line: python store_structure_data_sqlite.py structure_data.txt
Command line: python store_structure_data_sqlite.py /mnt/nexenta/stokm006/official_structure_data.txt
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
         line9 = line[9].replace("'", "")
         line10 = line[10].replace("'", "")
         sql_string = "INSERT INTO structure VALUES ('%s', '%s', '%s', '%s',\
         '%s', '%s', '%s', '%s', '%s', '%s', '%s');"% (line[0], line[1], \
         line[2], line[3], line[4], line[5], line[6], line[7], line[8], line9, line10)
         sql_list += [sql_string] 

     conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")
     c = conn.cursor()
     
     # Deletes the table, which allows creation of a newer version
     c.execute("DROP TABLE IF EXISTS structure")
    
     # Add table with atrribute names
     c.execute('''CREATE TABLE structure 
              (structure_id text PRIMARY KEY, monoisotopic_mass num, inchi text,\
              smile text, inchi_key2 text, inchi_key1 text,\
              molecular_formula text, kingdom text,\
              superclass text, class text, subclass text)''')

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