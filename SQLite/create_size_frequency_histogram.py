#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 14:10:10 2018
Creates a histogram with substructures from class XXX generated with
generate_data_for_vector.py. (x: molecule size, y: substructure frequency)
Command line: python3 create_size_frequency_histogram.py /mnt/nexenta/stokm006/class_Akua_str_sub.txt
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
import matplotlib.pyplot as plt
from sys import argv

def create_histogram(data_input):
    """ Uses the data created with generate_data_for_vector.py and creates a
    histogram (substructure frequency vs. molecule size).
    
    data_input: txt file from generate_data_for_vector.py
    """
    all_lines = data_input.split('\n')
    all_lines = all_lines[2:-1]
    
    atom_list = []
    for line in all_lines:
        line = line.split('\t')
        m = Chem.MolFromSmiles(line[1])
        atom_list += [m.GetNumAtoms()]
        
    atom_list = sorted(atom_list)

    atom_dict = {}  
    nr_of_atoms = atom_list[0]
    x = 0
    for nr_of_atoms in atom_list:
        if nr_of_atoms not in atom_dict.keys():
            x = 0
            atom_dict[nr_of_atoms] = 1
        if nr_of_atoms in atom_dict.keys():
            x += 1
            atom_dict[nr_of_atoms] = x

    # Create histogram
    fig, ax = plt.subplots(figsize=(10,10)) 
    plt.bar(list(atom_dict.keys()), atom_dict.values(), color='r')
    plt.xlabel('Nr. of atoms in substructure (H atoms excluded)')
    plt.ylabel('Substructure frequency')
    plt.title("Substructure frequency for the class 'XXX'")
    plt.grid(True)
    fig.savefig('/path/to/save/frequency_histogram.png')
  #  plt.show()

if __name__ == "__main__":
    with open(argv[1]) as file_object:
        input_file = file_object.read()
        create_histogram(input_file)