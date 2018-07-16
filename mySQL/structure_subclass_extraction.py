#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 13:09:57 2018
Get structures based on subclass
Command line: python structure_subclass_extraction.py
@author: stokm006
"""


import MySQLdb

def try_sql():

    # Connect
    db = MySQLdb.connect(host = "dmcourse.cgdm9vqizsxq.eu-central-1.rds.amazonaws.com", 
                     user="stokm006", 
                     passwd="7b24d294f52c1dd2667202841f30a76ce259a967", 
                     db="dbstokm006structuredb")

    cursor = db.cursor()
    
      
    cursor.execute("SELECT SMILES FROM dbstokm006structuredb.Structure\
                    WHERE Subclass = 'Halomethanes';")
    
                    
    numrows = cursor.rowcount


    with open("Structure_subclasss_TEST.txt", 'w') as db_file:
        for x in range(0, numrows):
            row = cursor.fetchone()
            row = str(row)
            row = row.lstrip('(\'').rstrip('\',)')
            db_file.write(row + '\n')
 
    # Close the connection
    db.close()

if __name__ == '__main__':
    try_sql()