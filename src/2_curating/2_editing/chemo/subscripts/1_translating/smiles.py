#!/usr/bin/env python
# coding: utf-8

# # RDKIT to generate InChIs

# loading packages

from rdkit import Chem
import pandas as pd
import gzip
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import numpy as np

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

