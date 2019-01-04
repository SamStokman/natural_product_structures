"""
04-01-'19
Script that parses the molBLOCKS substructure output, the substructures
SMILES are converted to the canonical RDkit SMILES and the option 
can be selected to add the whole structure as a substructure.

Command line: python3 process_molBLOCKS_output.py -i [input structures file] -s [input substructure file] -a [add orginal structure y/n] -o [output file name]
"""

import sys, getopt
from rdkit import Chem


def parse_substructures(input_structure_data, input_substructure_data, add_whole_structure, output_file_name):
    """
    Takes the (sub)structures SMILES and identifier and creates a output
    file with their canonical RDkit SMILES and identifier. Also the 
    original structure can be added as substructure. 
    
    input_structure_data: text file with input structure SMILES
    input_substructure_data: text file with the output substructure
    SMILES from molBLOCKS
    add_whole_structure: [y/n] add the whole structure as substructure
    output_file_name: the name of the output text file
    """

    with open(input_structure_data) as file_object:
        input_file = file_object.read()
  
    input_data = input_file.split('\n')
    input_list = []
    for row in input_data[0:-1]:
        row = row.split('\t')
        mol = Chem.MolFromSmiles(row[0])
        SMILES = Chem.MolToSmiles(mol)
        input_list += [SMILES]

    with open(input_substructure_data) as file_object:
        input_subs_file = file_object.read()

    input_subs_data = input_subs_file.split('\n')
    subs_list = []
    structure_id_list = []
    for row in input_subs_data[0:-1]:
        row = row.split('\t')
        if row[0] != '<NA>':
            mol = Chem.MolFromSmiles(row[0])
            SMILES = Chem.MolToSmiles(mol)
            subs_list += [SMILES]
        if row[0] == '<NA>':
            subs_list += [row[0]]
        structure_id_list += [row[1]]

    combo_list= []
    for i in range(len(input_list)):
        combo = [subs_list[i], input_list[i]]
        combo_list += [combo]
    
    sub_out_list = []
    for combo in combo_list:
        if combo[0] == '<NA>':
            if add_whole_structure == 'y':
                subs_out = combo[1]
            if add_whole_structure == 'n':
                subs_out = '<NA>'
                
        if combo[0] != '<NA>':
            if add_whole_structure == 'y':
                subs_out = combo[0]+'.'+combo[1]
            if add_whole_structure == 'n':
                subs_out = combo[0]
        sub_out_list += [subs_out]
    
    with open(output_file_name, 'w') as db_file:
        for i in range(len(input_list)):
           db_file.write(sub_out_list[i]+'\t'+structure_id_list[i]+'\n')

def main(argv):
    """
    Main function that calls all other functions with the corresponding
    argument parameter values.
    
    argv: list with all input arguments
    """
    try:
        opts, args = getopt.getopt(argv,"i:s:a:o:")
        if len(opts) != 4:
           print ('You did not enter the right arguments')
           sys.exit(2)
        opt_list = []
        for opt, arg in opts:
           if opt not in opt_list:
               opt_list += [opt]
        if '-i' not in opt_list or '-s' not in opt_list or '-o' not in opt_list or '-a' not in opt_list:
           print ('You did not enter the right arguments')
           sys.exit(2)      
    except getopt.GetoptError:
        print ('You did not enter the right arguments')
        sys.exit(2)
      
    for opt, arg in opts:
       if opt in ['-i']:
          input_structure_file = arg
       elif opt in ['-s']:
          input_substructure_file = arg
       elif opt in ['-a']:
          add_whole_structure = arg
       elif opt in ['-o']:
          output_file = arg

    parse_substructures(input_structure_file, input_substructure_file, add_whole_structure, output_file)

if __name__ == '__main__':
    
   main(sys.argv[1:])
