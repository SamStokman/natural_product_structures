# -*- coding: utf-8 -*-

"""
Created on Thu Jun 28 11:12:48 2018
This script can be used to get information about the database in
Structure_Database_File.txt from create_structure_db.py
Command line: python3 parse_structure_db.py Structure_Database_File.txt
@author: stokm006
"""

from __future__ import print_function
from sys import argv
from rdkit import Chem

def use_data(NP_data):

    all_lines = NP_data.split('\n')
    all_lines = all_lines[1:-1]

    x = 0
    y = 0

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
                x += 1
                all_info = [line[0], line[1], line[2], line[3], line[4], line[5], \
                line[6], line[7], line[8], line[9], line[10], line[11]]
                succes_list += [all_info]
            if not m.startswith('<'):
                y += 1
                all_info = [line[0], line[1], line[2], line[3], line[4], line[5], \
                line[6], line[7], line[8], line[9], line[10], line[11]]
                fail_list += [all_info]

    # whether the smile is recognized by RDkit
    print ('nr of recognized structures: ', x)
    print ('nr of unrecognized structures :', y)
    
    all_info2 = []
    for line in succes_list:
        m = Chem.MolFromSmiles(line[2])
        n = Chem.MolToSmiles(m)
        all_info = [line[0], line[1], n, line[3], line[4], line[5], line[6],\
        line[7], line[8], line[9], line[10], line[11]]
        all_info2 += [all_info]
    
    # the overlapping smiles are grouped together in the list
    sorted_list = sorted(all_info2, key = lambda x:(x[2]))

    # put the unique smile and the nr of dbs it occurs in in dict
    x = 1
    
    smile = sorted_list[0][2]
    unique_smile_list = [sorted_list[0][2]]
    
    smile_dict = {}
    for line in sorted_list:
        if smile != line[2]:          
            smile = line[2]
            x = 0
        if smile == line[2]:
            x += 1
            smile_dict[smile] = [x]
            smile = line[2]
        
    # returns a dict with smile as key and its quantity, identiefier and 
    # database as value
    smile_dict2 = smile_dict
    for line in sorted_list:
        for key in smile_dict:
             if key == line[2]:
                 smile_dict2[line[2]].append(line[3])
                 smile_dict2[line[2]].append(line[11])
 
    # only print the structures which where found in multiple dbs
    a=b=c=d=e=f=g=h=i=j=0
    for key, value in smile_dict.items():
        if value[0] == 1:
            a += 1
        if value[0] == 2:
            b += 1
        if value[0] == 3:
            c += 1
        if value[0] == 4:
            d += 1
        if value[0] == 5:
            e += 1
        if value[0] == 6:
            f += 1
        if value[0] == 7:
            g += 1
        if value[0] == 8:
            h += 1
        if value[0] == 9:
            i += 1
        if value[0] > 9:
            j += 1
            
    print ('Nr of unique SMILES: ', len(smile_dict))    
    print ('Nr of structures that occurs once: ', a)
    print ('Nr of structures that occurs twice: ', b)
    print ('Nr of structures that occurs three times: ', c)
    print ('Nr of structures that occurs four times: ', d)
    print ('Nr of structures that occurs five times: ', e)
    print ('Nr of structures that occurs six times: ', f)
    print ('Nr of structures that occurs seven times: ', g)
    print ('Nr of structures that occurs eight times: ', h)
    print ('Nr of structures that occurs nine times: ', i)
    print ('Nr of structures that occurs more than nine times: ', j)
    
    
if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        use_data(input_file)
