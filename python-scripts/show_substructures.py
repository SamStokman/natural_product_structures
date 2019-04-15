"""
Created on Mon Sep 10 09:11:22 2018
This script shows the highlighted substructure, which requires the structure
input and substructure output.

Usage:
Command line: python3 show_substructures.py -i [input file] -o [output_file] -s [structure nr] -m [method of substructure showing]
-s: int, 0 - nr. of input structures, ascending order.
-m: n = show the substructures.
    y = show the unique substructure with the whole structure as background, for each substructure, one match is shown.
    a = show the substructures with the whole structure as backgroud, all substructure matches are shown.


@author: stokm006
"""

import sys, getopt
from rdkit import Chem
from rdkit.Chem import Draw    
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG = True


def visualize_substructures(input_file, result_data, structure_nr, method):
    """
    This function visualizes substructures for quick and easy substructure
    evaluation.
        
    input_file: input text file with structures' SMILES and structure_id's
    result_data: output text file with substructures' SMILES and original
    structure's id
    structure_nr: the structure's number which substructures should be shown
    method: the approach to show the substructures
    """
   
    with open(input_file) as file_object:
        input_file = file_object.read()  
        
    input_data = input_file.split('\n')
    input_data = input_data[structure_nr:structure_nr+2]
    input_list = []
    for row in input_data:
        row = row.split('\t')
        input_list += [row[0]]
    
    with open(result_data) as file_object:
        result_data = file_object.read()    
   
    all_lines = result_data.split('\n')
    all_lines = all_lines[structure_nr:structure_nr+2]
    all_info_input_list = []
    all_info_output_list = []
    



    for x, line in enumerate(all_lines):
        line = line.split('\t')
        smile = line[0]
        smile = smile.split(".")
        nr_of_substructures = len(smile)
        input_molecule = []
        input_molecule = nr_of_substructures * [input_list[x]]
        input_mol_list = []
        for m in input_molecule:
            m = Chem.MolFromSmiles(m)
            input_mol_list += [m]
        all_info_input_list += [input_mol_list]
 
        output_list = []
        for i in range(len(smile)):
            smile_str = smile[i]
            output_list += [smile_str]
        mol_result_list = []
        for smile in output_list:
            m2 = Chem.MolFromSmiles(smile, sanitize=False)
            mol_result_list += [m2]
        all_info_output_list += [mol_result_list]


    if all_info_output_list[0][0] is None: # if no substructures are present, the whole structure is shown
        all_info_output_list = all_info_input_list

    # just show the substructures    
    if method == 'n':

        im = Draw.MolsToGridImage(all_info_output_list[0], 
                                  useSVG=False,
                                  molsPerRow=4,
                                  subImgSize = (300,200))
        im.show()
  
    # shows the (highlighted) substructures in structure (does not work for every substructure)  
    if method == 'y':
        [Chem.GetSSSR(mol) for mol in all_info_output_list[0]]
        [AllChem.Compute2DCoords(mol) for mol in all_info_output_list[0]]
     
        mol_input = all_info_input_list[0][0] 

        
        highlight_lists = []
        for mol in all_info_output_list[0]: 
            AllChem.GenerateDepictionMatching2DStructure(mol_input, mol, acceptFailure = True)
            highlight_lists += [mol_input.GetSubstructMatch(mol)]

        input_molecule = []
        for i in range(len(highlight_lists)): 
            input_molecule += [mol_input]

        im = Draw.MolsToGridImage(input_molecule, 
                                  highlightAtomLists = highlight_lists,
                                  useSVG=False,
                                  molsPerRow=2,
                                  subImgSize = (600,400))
        im.show()    
    
    # shows all (highlighted) substructures in structure (does not work for every substructure)  
    if method == 'a':
        [Chem.GetSSSR(mol) for mol in all_info_output_list[0]]
        [AllChem.Compute2DCoords(mol) for mol in all_info_output_list[0]]
     
        mol_input = all_info_input_list[0][0] 

        
        highlight_lists = []
        for mol in all_info_output_list[0]: 
            AllChem.GenerateDepictionMatching2DStructure(mol_input, mol, acceptFailure = True)
            if len(mol_input.GetSubstructMatches(mol)) == 1:
                highlight_lists += [mol_input.GetSubstructMatch(mol)]
            if len(mol_input.GetSubstructMatches(mol)) > 1:
                matches = mol_input.GetSubstructMatches(mol)
                for match in matches:
                    highlight_lists += [match]
                               
        input_molecule = []
        for i in range(len(highlight_lists)): 
            input_molecule += [mol_input]

        im = Draw.MolsToGridImage(input_molecule, 
                                  highlightAtomLists = highlight_lists,
                                  useSVG=False,
                                  molsPerRow=3,
                                  subImgSize = (350,250))
        im.show()




def main(argv):
    """
    Main function that calls all other functions with the corresponding
    argument parameter values.
    
    argv: list with all input arguments
    """
    

    try:
        opts, args = getopt.getopt(argv,"i:o:s:m:")
        if len(opts) != 4:
           print ('You did not enter the correct argument')
           sys.exit(2)
        opt_list = []
        for opt, arg in opts:
           if opt not in opt_list:
               opt_list += [opt]
        if '-i' not in opt_list or '-o' not in opt_list or '-s' not in opt_list or '-m' not in opt_list:
           print ('You did not enter the correct arguments')
           sys.exit(2)      
    except getopt.GetoptError:
        print ('You did not enter the correct arguments')
        sys.exit(2)
      
    for opt, arg in opts:
       if opt in ['-i']:
          input_file = arg
       if opt in ['-o']:
          output_file = arg
       elif opt in ['-s']:
          structure_nr = int(arg)
       elif opt in ['-m']:
           method = arg

    visualize_substructures(input_file, output_file, structure_nr, method)

if __name__ == '__main__':
    
   main(sys.argv[1:])
