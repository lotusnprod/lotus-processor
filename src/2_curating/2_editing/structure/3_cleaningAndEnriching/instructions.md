# Standardization of structural inputs

This is done using modules of the MolVS package, a molecule validation and standardization tool, written in Python. (<https://molvs.readthedocs.io/en/latest/>)

## Running the script

### Requirements

To do: make requirements files (mainly needs rdkit, molvs). Note that the molvs **input**.py file needs to be mannually modified to launch the function we call in the script (see script)

## fragment.py

in /Users/USERNAME/opt/anaconda3/lib/python3.7/site-packages/molvs/ needs to be replaced with the version present in this folder

### What the scripts does

Given a tab delimited file as input having a SMILES and an InChI column it will firstgenerate a ROmol object for SMILES and InChI.
It will then:

- validate the ROmol and apend a new column with the validation log (allowing to see what are the problems).
- standardize the ROmol object
- fetch the largest fragment
- uncharge the molecule

Once these steps are realized the following fields are generated from the sanitized output:

- SMILES
- InChI
- InChIKey
- Short InChI key
- Molecular Formula
- Exact Mass
- Protonated Exact Mass
- Deprotonated Exact Mass

### Command line

In your appropriate conda env run

- General

`python opennpdb_sanitizer_script.py /path/to/input_file.tsv /path/to/output_file.tsv smiles_column_header inchi_column_header`

- Example

`nohup python opennpdb_sanitizer_script.py ../outputs/tables/open_NP_db.tsv ../outputs/tables/open_NP_db_sanitized.tsv smiles inchi`

NB the nohup command allows to launch a command via ssh and let it run in the background.
