# natural_product_structures
The repository natural_product_structures stores data and python scripts for natural product (sub)structures.



## Creating the Natural Product Structure Database
Used databases:
- GNPS (CLASS)
- Super Natural II (CLASS)         
- NuBBE (CLASS)
- Human Metabolome (CLASS)
- Yeast Metabolome (CLASS)
- ChEBI (CLASS)
- DrugBank (CLASS)
- NANPDB (converted to CLASS, see python-scripts/nanpdb_CLASS_parser.py)
- Streptomedb2 (converted to CLASS, see python-scripts/streptomedb_CLASS_parser.py)
- NP Atlas (converted to CLASS, see python-scripts/NPAtlas_CLASS_parser.py)
- Norine (converted to CLASS, see python-scripts/norine_CLASS_parser.py)

These CLASS databases are merged together in a file (Data/Structure_Database_File) with the script python-scripts/create_structure_db.py. The file is tab-separated and contains overlapping structures (193.944 KB). (new version is 200.061 KB, not uploaded yet)

HAS TO BE UPDATED*
The script python-scripts/get_db_information.py gives information about the number of (un)recognized and unique smiles in the Structure_Database_File. Also the number of structures that occurs multiple times is generated. In total, 477.349* structures are recognized by their SMILE and 244* structures are not recognized. The number of unique SMILES is 312.938*.

All SMILES from the recognized structures are made uniform and converted into their canonical SMILE. The structures are ordered based on these canonical SMILE and stored in Data/Canonical_db_file.txt, the script python-script/create_canonical_SMILE_db.py is used to create this file.

The Canonical_db_file.txt file is used to create the tables for the sqlite database (see Data/DatabaseDesign). The tables are created by the python-scripts/sqlite_data_table_creator.py. These tables are then converted (with the scripts in the SQLite dir) into sqlite tables and added to the Natural_Product_Structure.sqlite database (Dropbox).

## RDkit
The directory RDkit is used to store several python scripts that test the most important functionalities of RDkit. Also a RDkit docx file is added with some furter explanation about the installation and the functionalities. 
