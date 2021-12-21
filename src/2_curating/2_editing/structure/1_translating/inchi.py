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
    inchi_column_header = sys.argv[3]

    print('Parsing gziped tab separated file '
          + '\n'
          + input_file_path
          + 'with column: '
          + '\n'
          + inchi_column_header
          + ' as InChI column. \n'
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
    sep='\t',
    encoding='utf-8',
    on_bad_lines='error'
    )

if (len(df) == 1) and (df.empty):
    df[inchi_column_header] = 'InChI=1S/Pu'
    print('your dataframe is empty, plutonium loaded \n')
else:
    print('your dataframe is not empty :) \n')

# df.head()
# df.info()

# keeping non-null entries
df = df[df[inchi_column_header].astype(str).str.startswith('InChI')]

# replacing unwanted characters
df[inchi_column_header].replace(regex=True,
                                 inplace=True,
                                 to_replace=r'"',
                                 value=r''
                                 )

##


# generating ROMOL
df['ROMol'] = df[inchi_column_header].map(Chem.MolFromInchi)

# removing entries for which no ROMol object could be generated
df = df[~df['ROMol'].isnull()]

# generating the InChI
df['smiles'] = df['ROMol'].map(Chem.MolToSmiles)

# renaming
df['structureTranslated_inchi'] = df['smiles']

# dropping old columns
df = df.drop(['smiles', 'ROMol'], axis=1)

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
    compression='gzip',
    encoding='utf-8'
)
