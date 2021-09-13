#!/usr/bin/env python
# coding: utf-8

# # RDKIT to generate InChIs

# loading packages
import errno
import gzip
import os
import sys

import pandas as pd
from rdkit import Chem

# here we deine input, outputs and other variables as sys arguments
try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    smiles_column_header = sys.argv[3]

    print('Parsing gziped tab separated file '
          + '\n'
          + input_file_path
          + 'with column: '
          + '\n'
          + smiles_column_header
          + ' as SMILES column. \n'
          + ' Proceeding to the the ROMol object conversion and returning the InChI in file : '
          + '\n'
          + ouput_file_path)
except:
    print(
        'Please add input and output file path as first and second argument and SMILES column header as third argument. \n')

# loading data
myZip = gzip.open(input_file_path)

df = pd.read_csv(
    myZip,
    sep='\t')

if len(df) == 1:
    df[smiles_column_header] = '[Pu]'
    print('your dataframe is empty, plutonium loaded \n')
else:
    print('your dataframe is not empty :) \n')

# df.head()
# df.info()

# keeping non-null entries
df = df[df[smiles_column_header].notnull()]

# replacing unwanted characters
df[smiles_column_header].replace(regex=True,
                                 inplace=True,
                                 to_replace=r'"',
                                 value=r''
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
if not os.path.exists(os.path.dirname(ouput_file_path)):
    try:
        os.makedirs(os.path.dirname(ouput_file_path))
    except OSError as exc:  # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

df.to_csv(
    ouput_file_path,
    sep='\t',
    index=False,
    compression='gzip'
)
