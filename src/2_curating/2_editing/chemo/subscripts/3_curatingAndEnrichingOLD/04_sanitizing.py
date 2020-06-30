#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

## importing packages 
from rdkit import Chem
import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from molvs import Standardizer
from molvs import LargestFragmentChooser #requires to modify __init__.py in /Users/USERNAME/opt/anaconda3/lib/python3.7/site-packages/molvs/ accordingly
from molvs import FragmentRemover #requires to modify __init__.py in /Users/USERNAME/opt/anaconda3/lib/python3.7/site-packages/molvs/ accordingly
from molvs import Validator
from molvs import Uncharger ## add this to init.py see above from .fragment import LargestFragmentChooser from .charge import Uncharger
import sys
import gzip
import yaml

try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    inchi_column_header = sys.argv[3]

    print('Parsing tab separated file'
          + input_file_path 
          + 'with column: '
          + inchi_column_header
          + 'as InChI column.'
          + 'Proceeding to the validation, standardization, fragment choosing and uncharging of the ROMol object and returning the sanitized outputs in file :'
          + ouput_file_path)
except:
    print('Please add input and output file path as first and second argument and SMILES and InChI column header as third and forth argument.')


# This load the associated params.yaml file were we can specify some variables of the script, here (head_lenght)
# The pipeline can then just be rerunned using dvc repro

with open("params.yaml", 'r') as fd:
    params = yaml.safe_load(fd)

head_lenght = params['head_lenght']

## Loading the df with inchi columns 
myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
	sep = '\t')

## eventually filter display some info, comment according to your needs
#df = df[df['originaldb'] == 'tcm']
df = df.head(head_lenght)
df.columns
df.info()

# defining the validator log output format
fmt = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'

# save the Standardizer and LargestFragmentChooser classes as variables
validator = Validator(log_format = fmt)
s = Standardizer()
#lf = LargestFragmentChooser()
lf = FragmentRemover()
uc = Uncharger()

#we make sure we have inchi for each row
df = df[~df[inchi_column_header].isnull()]

# we generate ROMol object from smiles and or inchi
df['ROMol'] = df[inchi_column_header].map(Chem.MolFromInchi)

# we eventually remove rows were no ROMol pobject was generated
df = df[~df['ROMol'].isnull()]

## and now apply the validation, standardization, fragment chooser and uncharging scripts as new columns.
# Note that these are sequentially applied
df['validatorLog'] = df['ROMol'].apply(validator.validate)
df['ROMolSanitized'] = df['ROMol'].apply(s.standardize)
#df['ROMolSanitizedLargestFragment'] = df['ROMolSanitized'].apply(lf.choose)
df['ROMolSanitizedLargestFragment'] = df['ROMolSanitized'].apply(lf.remove)
df['ROMolSanitizedLargestFragmentUncharged'] = df['ROMolSanitizedLargestFragment'].apply(uc.uncharge)

# outputting smiles, inchi, molecular formula, exact mass and protonated and deprotonated exactmasses from the latest object of the above scripts
df['smilesSanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].map(Chem.MolToSmiles)
print(df['ROMolSanitized'])
df['inchiSanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].map(Chem.MolToInchi)
df['inchikeySanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].map(Chem.MolToInchiKey)
df['shortikSanitized'] = df['inchikeySanitized'].str.split("-", n = 1, expand = True)[0]

df['formulaSanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].apply(rdMolDescriptors.CalcMolFormula)

df['exactmassSanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].apply(Descriptors.ExactMolWt)

df['xlogpSanitized'] = df['ROMolSanitizedLargestFragmentUncharged'].apply(Chem.Crippen.MolLogP)

# outputing final df inofs
df.info()

# dropping some irrelevant columns prior to export
colstodrop = ['ROMol',
              'ROMolSanitized',
              'ROMolSanitizedLargestFragment',
              'ROMolSanitizedLargestFragmentUncharged'
              ]

# colstodrop = ['ROMol',
#               'ROMol_largest_fragment',
#               'ROMol_largest_fragment_sanitized',
#               'ROMol_largest_fragment_sanitized_uncharged'
#               ]

df = df.drop(colstodrop, axis = 1)

# exporting df
df.to_csv(
	ouput_file_path, 
	sep = '\t', 
	index = False,
	compression = 'gzip'
	)
