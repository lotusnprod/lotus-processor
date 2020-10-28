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
myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
	sep = '\t') 

if len(df) == 1:
    df[smiles_column_header] = '[Pu]'
    print('your dataframe is empty, plutonium loaded')
else:
  print('your dataframe is not empty :)')

# df.head()
# df.info()

# keeping non-null entries
df = df[df[smiles_column_header].notnull()]

# replacing unwanted characters
df[smiles_column_header].replace(regex = True,
										inplace = True,
										to_replace = r'"',
										value = r''
										)

##


# generating ROMOL
df['ROMol'] = df[smiles_column_header].map(Chem.MolFromSmiles)

# removing entries for which no ROMol object could be generated
df = df[~df['ROMol'].isnull()]

# generating the InChI
df['inchi'] = df['ROMol'].map(Chem.MolToInchi)

# renaming
# naming not OK checkz with Adriano 

df['structureTranslated_smiles'] = df['inchi']

# dropping old columns
df = df.drop(['inchi', 'ROMol'], axis=1)

# exporting

import os
import errno

filename = ouput_file_path
if not os.path.exists(os.path.dirname(filename)):
    try:
        os.makedirs(os.path.dirname(filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
            
df.to_csv(
    ouput_file_path,
    sep = '\t',
    index = False,
    compression = 'gzip'
)

