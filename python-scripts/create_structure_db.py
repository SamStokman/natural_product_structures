# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 16:27:08 2018
Makes a text file with all data from the original CLASS NP databases which is
tab-separatable). All NP databases can be added one by one.
Can be used for the databases: chebi, drugdb, gnps, hmdb, nanpdb, nubbe, 
streptomedb, supernatural, ymdb
Command line: python3 create_structure_db.py CLASSDATABASEFILE.txt
@author: stokm006
"""

from sys import argv
import os.path


def parse_file(input_file, db_name):
    """ takes all text from CLASS database file and returns a list of lists 
    with NPs which is easy to use
    
    input_file: CLASS database txt file
    """
    
    all_lines = input_file.split('\n')
    all_info_list = []
    attribute_list = []
    for line in all_lines:
        line = line.split('\t')
        info_per_row_string = ""
        for value in line:
            my_string = ""
            value = value.strip('\'"')
            if len(value) == 0:
                value = "NA"
            my_string += value
            info_per_row_string += my_string+'\t'
        all_info_list += [info_per_row_string+db_name]
        attribute_list += [info_per_row_string]
    all_info_list2 = [all_info_list[1:]]
    attribute_list2 = [attribute_list[0]]
    
    # remove empty lines, for nubbe db   
    if db_name.startswith('nubb'):
        all_info_list3 = [] 
        for info_list in all_info_list2:
            for line in info_list:
                if len(line) > 75:
                    all_info_list3 += [line]
                    all_info_list2 = [all_info_list3]
    
    return all_info_list2, attribute_list2


    
def write_CLASS_txtfile(input_file, data, attribute_names):
    """ takes all text from info list and writes an easy to use CLASS database.
    
    input_file_name: name of txt file that will be created
    data: database list created with parse_file()
    attribute_names: list with attribute names
    """
    
    if not os.path.isfile(input_file):
        with open(input_file, 'w') as db_file:
            for line in attribute_names:
                db_file.write(str(line) +'Database_name'+ '\n')        
        
    with open(input_file, 'a') as db_file:
        for line in data:
            for i in range (len(line)):
                db_file.write(str(line[i]) + '\n')

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        db_name = argv[1]
        db_name = db_name[:-4]
        input_file = file_object.read()      
        parsed_data, attr_list = parse_file(input_file, db_name)
        write_CLASS_txtfile("Structure_Database_FriWedd.txt", parsed_data, attr_list)




