# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 08:54:51 2018
Stores the structure database data in de structure database on the sql server.
Command line: python store_structure_database_data_sql.py structure_database_data.txt
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
        sql_string = "INSERT INTO `dbstokm006structuredb`.`Structure_Database`\
        (`Structure_Database_ID`, `Structure_ID`, `Database_Identifier`, \
        `Database_Name`) VALUES ('%s','%s', '%s', '%s');"\
        % (line[0], line[1], line[2], line[3])
        sql_list += [sql_string] 
  

    # Connect
    db = MySQLdb.connect(host = "dmcourse.cgdm9vqizsxq.eu-central-1.rds.amazonaws.com", 
                     user="stokm006", 
                     passwd="7b24d294f52c1dd2667202841f30a76ce259a967", 
                     db="dbstokm006structuredb")

    cursor = db.cursor()
    
    # Add table with atrribute names
    cursor.execute("CREATE TABLE IF NOT EXISTS `dbstokm006structuredb`.\
    `Structure_Database` (`Structure_Database_ID` VARCHAR(45) NOT NULL,\
  `Structure_ID` VARCHAR(45) NOT NULL,\
  `Database_Identifier` VARCHAR(45) NOT NULL,\
  `Database_Name` VARCHAR(45) NOT NULL,\
  PRIMARY KEY (`Structure_Database_ID`),\
  UNIQUE INDEX `Structure_Database_ID_UNIQUE` (`Structure_Database_ID` ASC),\
  INDEX `Structure_ID_idx` (`Structure_ID` ASC),\
  INDEX `Database_Name_idx` (`Database_Name` ASC),\
  CONSTRAINT `Structure_ID`\
  FOREIGN KEY (`Structure_ID`)\
  REFERENCES `dbstokm006structuredb`.`Structure` (`Structure_ID`)\
  ON DELETE NO ACTION\
  ON UPDATE NO ACTION,\
  CONSTRAINT `Database_Name`\
  FOREIGN KEY (`Database_Name`)\
  REFERENCES `dbstokm006structuredb`.`Database` (`Database_Name`)\
  ON DELETE NO ACTION\
  ON UPDATE NO ACTION)\
  ENGINE = InnoDB;")

    
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
