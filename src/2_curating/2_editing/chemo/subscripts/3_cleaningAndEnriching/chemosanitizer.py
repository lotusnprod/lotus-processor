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


import chemosanitizer_functions
from chemosanitizer_functions import *

# defining the command line arguments

try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    inchi_column_header = sys.argv[3]
    cpus = sys.argv[4]

    print('Parsing tab separated file'
          + input_file_path
          + 'with column: '
          + inchi_column_header
          + 'as InChI column.'
          + 'Parralelized on '
          + cpus
          + ' cores.'
          + 'Proceeding to the validation, standardization, fragment choosing and uncharging of the ROMol object and returning the sanitized outputs in file :'
          + ouput_file_path)
except:
    print(
        '''Please add input and output file path as first and second argument, InChI column header as third argument and finally the number of cpus you want to use.
        Example :
        python chemosanitizer.py ~/translatedStructureRdkit.tsv ./test.tsv structureTranslated 6''')

# Loading the df with inchi columns
#df = pd.read_csv(input_file_path, sep='\t')

myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
	sep = '\t')

# eventually filter display some info, comment according to your needs
# df = df[df['originaldb'] == 'tcm']
#df = df.head(10000)
df.columns
df.info()


# timer is started
start_time = time.time()


# here we define the multiprocessing wrapper for the function. Beware to set the number of running tasks according to your cpu number

if __name__ == "__main__":
    # with multiprocessing.Pool(multiprocessing.cpu_count() - 2 ) as pool:
    with multiprocessing.Pool(int(cpus)) as pool:
        
        # # we generate ROMol object from smiles and or inchi
        df['ROMol'] = pool.map(MolFromInchi_fun, df[inchi_column_header])
        # # we eventually remove rows were no ROMol pobject was generated
        df = df[~df['ROMol'].isnull()]
        # # and now apply the validation, standardization, fragment chooser and uncharging scripts as new columns.
        # # Note that these are sequentially applied
        df['validatorLog'] = pool.map(validator_fun, df['ROMol'])
        df['ROMolSanitized'] = pool.map(standardizor_fun, df['ROMol'])
        df.drop('ROMol', axis=1, inplace=True)
        df['ROMolSanitizedLargestFragment'] = pool.map(fragremover_fun, df['ROMolSanitized'])
        df.drop('ROMolSanitized', axis=1, inplace=True)
        df['ROMolSanitizedLargestFragmentUncharged'] = pool.map(uncharger_fun, df['ROMolSanitizedLargestFragment'])
        df.drop('ROMolSanitizedLargestFragment', axis=1, inplace=True)
        # # outputting smiles, inchi, molecular formula, exact mass and protonated and deprotonated exactmasses from the latest object of the above scripts
        df['smilesSanitized'] = pool.map(MolToSmiles_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        # for the inchi and IK since some specific structures are raising issues we use the ***_fun_safe functions (see associated chemosanitizer_function.py)
        # df['inchi_sanitized'] = pool.map(MolToInchi_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        df['inchiSanitized'] = pool.starmap(MolToInchi_fun_safe, zip(df['smilesSanitized'], df['ROMolSanitizedLargestFragmentUncharged']))
        #df['inchikeySanitized'] = pool.map(MolToIK_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        df['inchikeySanitized'] = pool.starmap(MolToIK_fun_safe, zip(df['smilesSanitized'], df['ROMolSanitizedLargestFragmentUncharged']))
        df['shortikSanitized'] = df['inchikeySanitized'].str.split("-", n=1, expand=True)[0]
        df['formulaSanitized'] = pool.map(MolToMF_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        df['exactmassSanitized'] = pool.map(MolToEmass_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        df['xlogpSanitized'] = pool.map(MolToLogP_fun, df['ROMolSanitizedLargestFragmentUncharged'])
        
        pool.close()
        pool.join()

# timer is stopped
print(" Above command executed in --- %s seconds ---" %
      (time.time() - start_time))




# # dropping some irrelevant columns prior to export


# colstodrop = ['ROMol',
#               'ROMolSanitized',
#               'ROMolSanitizedLargestFragment',
#               'ROMolSanitizedLargestFragmentUncharged'
#               ]

colstodrop = ['ROMolSanitizedLargestFragmentUncharged']

df.drop(colstodrop, axis=1, inplace=True)

df.info()
# exporting df

#df.to_csv(ouput_file_path, sep='\t', index=False)

df.to_csv(
	ouput_file_path, 
	sep = '\t', 
	index = False,
	compression = 'gzip'
	)

