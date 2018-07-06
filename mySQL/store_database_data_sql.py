# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 16:41:57 2018
Stores the database data in de structure database on the sql server.
Command line: python store_database_data_sql.py database_data.txt
@author: stokm006
"""

import MySQLdb
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
        sql_string = "INSERT INTO `dbstokm006structuredb`.`Database` (`Database_Name`, `Nr_Of_Structures`)\
        VALUES ('%s', %s);"% (line[0], line[1])
        sql_list += [sql_string]

    # Connect
    db = MySQLdb.connect(host = "dmcourse.cgdm9vqizsxq.eu-central-1.rds.amazonaws.com", 
                     user="stokm006", 
                     passwd="7b24d294f52c1dd2667202841f30a76ce259a967", 
                     db="dbstokm006structuredb")

    cursor = db.cursor()
    
    # Add table with atrribute names
    cursor.execute("CREATE TABLE IF NOT EXISTS `dbstokm006structuredb`.`Database`\
    (`Database_Name` VARCHAR(45) NOT NULL,\
  `Nr_Of_Structures` INT(11) NOT NULL,\
  PRIMARY KEY (`Database_Name`),\
  UNIQUE INDEX `Database_UNIQUE` (`Database_Name` ASC))\
  ENGINE = InnoDB\
  DEFAULT CHARACTER SET = latin1;")

    # Add data to table
    for i in range(len(sql_list)):
        cursor.execute(sql_list[i])

    # Commit your changes
    db.commit()

    # Close the connection
    db.close()

if __name__ == '__main__':
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        create_database_data(input_file)
