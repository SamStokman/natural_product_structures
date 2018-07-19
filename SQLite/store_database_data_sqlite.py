# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 16:41:57 2018
Stores the database data in de sqlite database file
Command line: python store_database_data_sqlite.py database_data.txt
Command line: python store_database_data_sqlite.py /mnt/nexenta/stokm006/official_database_data.txt
@author: stokm006
"""


import sqlite3
from sys import argv

def create_database_data(data):
     all_lines = data.split('\n')
     all_lines = all_lines[:-1]
     all_info_list = []
     for line in all_lines:
         line = line.split('\t')
         all_info_list += [line]
         

      # Generate database data
     sql_list = [] 
     for line in all_info_list:
         line = line
         sql_string = "INSERT INTO data_source VALUES ('%s', %s);"% (line[0], line[1])
         sql_list += [sql_string]
 

     # Connect
     conn = sqlite3.connect("/mnt/nexenta/stokm006/Natural_Product_Structure.sqlite")

     c = conn.cursor()
     
     # Deletes the table, which allows creation of a newer version
     c.execute("DROP TABLE IF EXISTS data_source")
    
     # Add table with atrribute names
     c.execute('''CREATE TABLE data_source 
              (source_name text PRIMARY KEY, nr_of_structures num)''')
    
#     c.execute("INSERT INTO data_source VALUES ('sam', 2)")
#     c.execute("INSERT INTO data_source VALUES ('yolo', 1)")
              
     
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
        create_database_data(input_file)