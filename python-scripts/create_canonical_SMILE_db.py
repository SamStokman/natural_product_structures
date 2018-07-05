# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 10:56:27 2018
Create a CLASS database with the canonical SMILE of each structure, based on
these SMILES, structures are grouped. This script can be used for the combined
NP database (with the nine NP dbs in it).
command line: python3 create_canonical_SMILE_db.py Structure_Database_File.txt
@author: stokm006
"""

from __future__ import print_function
from sys import argv
from rdkit import Chem

def use_data(NP_data):
    """ takes all data from the NP database and writes a sorted structure txt
    file, the structures are sorted based on their canonical smile.
    
    NP_data: data from the NP database
    """

    all_lines = NP_data.split('\n')
    attribute_list = [all_lines[0]]
    all_lines = all_lines[1:-1]
    succes_list = []

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

    all_info2 = []
    for line in succes_list:
        m = Chem.MolFromSmiles(line[2])
        n = Chem.MolToSmiles(m)
        all_info = [line[0], line[1], n, line[3], line[4], line[5], line[6],\
        line[7], line[8], line[9], line[10], line[11]]
        all_info2 += [all_info]

    # The same smiles are grouped together in the list
    sorted_list = sorted(all_info2, key = lambda x:(x[2]))
    

    with open("Canonical_db_file.txt", 'w') as db_file:
        for line in attribute_list:
            db_file.write(str(line) + '\n')
        for line in sorted_list:
            for i in range (len(line)-1):
                db_file.write(str(line[i]) + '\t')
            db_file.write(str(line[11]) + '\n')


if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        use_data(input_file)
