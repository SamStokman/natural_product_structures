# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:25:40 2018
Stores the structure data in de structure database on the sql server.
Command line: python store_structure_data_sql.py structure_data.txt
Command line: python store_structure_data_sql.py /mnt/nexenta/stokm006/structure_data.txt
@author: stokm006
"""

import MySQLdb
from sys import argv


def try_sql(data):
    all_lines = data.split('\n')
    all_lines = all_lines[:-1]
    all_info_list = []
    for line in all_lines:
        line = line.split('\t')
        all_info_list += [line]
     
    sql_list = [] 
    for line in all_info_list:
        sql_string = "INSERT INTO `dbstokm006structuredb`.`Structure` (`Structure_ID`, `MonoisotopicMass`, `InChI`, `SMILES`, `InChIKey2`, `InChIKey1`, `MolecularFormula`, `Kingdom`, `Superclass`, `Class`, `Subclass`)\
        VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s');"\
        % (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10])
        sql_list += [sql_string] 

    # Connect
    db = MySQLdb.connect(host = "dmcourse.cgdm9vqizsxq.eu-central-1.rds.amazonaws.com", 
                     user="stokm006", 
                     passwd="7b24d294f52c1dd2667202841f30a76ce259a967", 
                     db="dbstokm006structuredb")

    cursor = db.cursor()
    
    # Add table with atrribute names
    cursor.execute("CREATE TABLE IF NOT EXISTS `dbstokm006structuredb`.`Structure` (\
  `Structure_ID` VARCHAR(20) NOT NULL,\
  `MonoisotopicMass` VARCHAR(45) NOT NULL,\
  `InChI` VARCHAR(500) NOT NULL,\
  `SMILES` VARCHAR(500) NOT NULL,\
  `InChIKey2` VARCHAR(10) NOT NULL,\
  `InChIKey1` VARCHAR(14) NOT NULL,\
  `MolecularFormula` VARCHAR(45) NOT NULL,\
  `Kingdom` VARCHAR(45) NOT NULL,\
  `Superclass` VARCHAR(45) NOT NULL,\
  `Class` VARCHAR(45) NOT NULL,\
  `Subclass` VARCHAR(200) NOT NULL,\
  PRIMARY KEY (`Structure_ID`),\
  UNIQUE INDEX `Structure_ID_UNIQUE` (`Structure_ID` ASC))\
  ENGINE = InnoDB\
  DEFAULT CHARACTER SET = latin1;")
  
    # Add data to table
    for i in range(len(sql_list)):
        cursor.execute(sql_list[i])
 

    # Commit changes
    db.commit()

    # Close the connection
    db.close()

if __name__ == '__main__':
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        try_sql(input_file)
