#!/usr/bin/env python

'''
File: /Users/pma/Dropbox/Research_UNIGE/git_repos/chemo_sanitizer/src/chemosanitizer.py
Project: /Users/pma/Dropbox/Research_UNIGE/git_repos/chemo_sanitizer/src
Created Date: Monday June 8th 2020
Author: PMA
-----
Last Modified: Monday June 8th 2020 5:24:23 pm
Modified By: the developer formerly known as PMA at <you@you.you>
-----
Copyright (c) 2020 Your Company
-----
HISTORY:
'''

# importing packages
import errno
import os

from chemosanitizer_functions import *

# defining the command line arguments
try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    smiles_column_header = sys.argv[3]

    print('Parsing tab separated file'
          + input_file_path
          + 'with column: '
          + smiles_column_header
          + 'as SMILES column.'
          + 'Proceeding to the stereocenters counting of the ROMol object and returning the counts outputs in file :'
          + ouput_file_path)
except:
    print(
        '''Please add input and output file path as first and second argument, SMILES column header as third argument.
        Example :
        python stereocounter.py ~/translatedStructureRdkit.tsv ./test.tsv smilesSanitized''')

# Loading the df with inchi columns
# df = pd.read_csv(input_file_path, sep='\t')
# input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/1_translated/structure/unique.tsv.zip'
myZip = gzip.open(input_file_path)

df = pd.read_csv(
    myZip,
    sep='\t')

if len(df) == 1:
    df['structureTranslated'] = 'InChI=1S/Pu'
    print('your dataframe is empty, plutonium loaded')
else:
    print('your dataframe is not empty :)')

# eventually filter display some info, comment according to your needs
# df = df[df['originaldb'] == 'tcm']
# df = df.head(50000)
df.columns
df.info()

# df = df[df[inchi_column_header].astype(str).str.startswith('InChI')]

# df.columns
# df.info()

# timer is started
start_time = time.time()

df['ROMol'] = df[smiles_column_header].map(Chem.MolFromSmiles)

df = df[~df['ROMol'].isnull()]

# this one "Returns the number of unspecified atomic stereocenters"
df['count_unspecified_atomic_stereocenters'] = df['ROMol'].map(
    Chem.rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters)

# this one "Returns the total number of atomic stereocenters (specified and unspecified)"
df['count_atomic_stereocenters'] = df['ROMol'].map(
    Chem.rdMolDescriptors.CalcNumAtomStereoCenters)

# We export this final df, after dropping the ROMol object

df.drop('ROMol', axis=1, inplace=True)

# timer is stopped
print(" Above command executed in --- %s seconds ---" %
      (time.time() - start_time))

df.info()

# exporting df
filename = ouput_file_path
if not os.path.exists(os.path.dirname(filename)):
    try:
        os.makedirs(os.path.dirname(filename))
    except OSError as exc:  # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

# df.to_csv(ouput_file_path, sep='\t', index=False)
df.to_csv(
    ouput_file_path,
    sep='\t',
    index=False,
    compression='gzip'
)
