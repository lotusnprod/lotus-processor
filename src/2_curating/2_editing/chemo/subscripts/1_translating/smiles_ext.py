#!/usr/bin/env python
# coding: utf-8

# # RDKIT to generate InChIs

# loading packages
import sys
from rdkit import Chem
import pandas as pd
import gzip
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import numpy as np

# here we deine input, outputs and other variables as sys arguments
try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    smiles_column_header = sys.argv[3]

    print('Parsing gziped tab separated file'
          + input_file_path 
          + 'with column: '
          + smiles_column_header
          + 'as SMILES column.'
          + 'Proceeding to the the ROMol object conversion and returning the InChI in file :'
          + ouput_file_path)
except:
    print('Please add input and output file path as first and second argument and SMILES column header as third argument.')


# loading data
myZip = gzip.open('../data/interim/tables_min/0_original/smiles.tsv.zip')

df = pd.read_csv(
	myZip,
	sep = '\t') 

# df.head()
# df.info()

# keeping non-null entries
df_i = df[df['structureOriginalSmiles'].notnull()]

# replacing unwanted characters
df_i['structureOriginalSmiles'].replace(regex = True,
										inplace = True,
										to_replace = r'"',
										value = r''
										)


# generating ROMOL
df_i['ROMol'] = df_i['structureOriginalSmiles'].map(Chem.MolFromSmiles)

df = df_i

# removing entries for which no ROMol object could be generated
df = df[~df['ROMol'].isnull()]

# generating the InChI
df['inchi'] = df['ROMol'].map(Chem.MolToInchi)

# renaming
df['structureTranslatedSmiles'] = df['inchi']

# dropping old columns
df = df.drop(['inchi', 'ROMol'], axis=1)

# exporting
df.to_csv(
    "../data/interim/tables_min/1_translated/smiles.tsv.zip",
    sep = '\t',
    index = False,
    compression = 'gzip'
)

