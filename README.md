# natural_product_structures
The repository natural_product_structures stores data and python scripts for natural product (sub)structures.



## Creating the Natural Product Structure Database
Used databases:
- GNPS (CLASS)
- Super Natural II (CLASS)         
- NuBBe (CLASS)
- Human Metabolome (CLASS)
- Yeast Metabolome (CLASS)
- ChEBI (CLASS)
- DrugBank (CLASS)
- NANPDB (converted to CLASS, see python-scripts/nanpdb_CLASS_parser.py)
- Streptomedb2 (converted to CLASS, see python-scripts/streptomedb_CLASS_parser.py))

These CLASS databases are merged together in a file (Data/Structure_Database_File) with the script python-scripts/create_structure_db.py. The file is tab-separated and contains overlapping structures (193 944 KB). 

The script python-scripts/parse_structure_db.py gives information about the number of (un)recognized and unique smiles in the Structure_Database_File. Also the number of structures that occurs multiple times is generated. 

In total, 477.349 structures are recognized by their SMILE. 592 structures are not recognized, 348 of those structures are empty lines which originate from NuBBe.

All SMILES from the recognized structures are made uniform and converted into their canonical SMILE. Based on these canonical SMILES, the number of unique SMILES is determined, which is 312.938.

The Structure_Database_File is also used to create the tables for the SQL database (see Data/DatabaseDesign). The tables are created by the python-scripts/create_tables_NPdata_sql.py. 



## RDkit
The directory RDkit is used to store several python scripts that test the most important functionalities of RDkit. Also a RDkit docx file is added with some furter explanation about the installation and the functionalities. 
