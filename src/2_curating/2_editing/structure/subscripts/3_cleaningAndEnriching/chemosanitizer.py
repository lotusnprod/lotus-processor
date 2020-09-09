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
# df = pd.read_csv(input_file_path, sep='\t')
# input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/1_translated/structure/unique.tsv.zip'
myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
	sep = '\t')

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

df = df[df[inchi_column_header].astype(str).str.startswith('InChI')]

df.columns
df.info()

# the full df is splitted and each subdf are treated sequentially as df > 900000 rows retruned errors 
# (parralel treatment of these subdf should improve performance)
n = 20000  # chunk row size
list_df = [df[i:i+n] for i in range(0,df.shape[0],n)]

# timer is started
start_time = time.time()

for i in range(0, len(list_df)):

    # here we define the multiprocessing wrapper for the function. Beware to set the number of running tasks according to your cpu number

    if __name__ == "__main__":
        # with multiprocessing.Pool(multiprocessing.cpu_count() - 2 ) as pool:
        with multiprocessing.Pool(int(cpus)) as pool:

            # list_df[i].to_csv(
            #     "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/2_cleaned/structure/outfirst38%i.csv" % i , 
            #     sep = '\t', 
            #     index = False,
            #     compression = 'gzip'
            #     )

            # # we generate ROMol object from smiles and or inchi
            list_df[i]['ROMol'] = pool.map(MolFromInchi_fun, list_df[i][inchi_column_header])
            # # we eventually remove rows were no ROMol pobject was generated
            list_df[i] = list_df[i][~list_df[i]['ROMol'].isnull()]
            # # and now apply the validation, standardization, fragment chooser and uncharging scripts as new columns.
            # # Note that these are sequentially applied
            list_df[i]['validatorLog'] = pool.map(validator_fun, list_df[i]['ROMol'])
            list_df[i]['ROMolSanitized'] = pool.map(standardizor_fun, list_df[i]['ROMol'])
            list_df[i].drop('ROMol', axis=1, inplace=True)
            list_df[i]['ROMolSanitizedLargestFragment'] = pool.map(fragremover_fun, list_df[i]['ROMolSanitized'])
            list_df[i].drop('ROMolSanitized', axis=1, inplace=True)
            list_df[i]['ROMolSanitizedLargestFragmentUncharged'] = pool.map(uncharger_fun, list_df[i]['ROMolSanitizedLargestFragment'])
            list_df[i].drop('ROMolSanitizedLargestFragment', axis=1, inplace=True)
            # # outputting smiles, inchi, molecular formula, exact mass and protonated and deprotonated exactmasses from the latest object of the above scripts
            list_df[i]['smilesSanitized'] = pool.map(MolToSmiles_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            # for the inchi and IK since some specific structures are raising issues we use the ***_fun_safe functions (see associated chemosanitizer_function.py)
            # list_df[i]['inchi_sanitized'] = pool.map(MolToInchi_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            list_df[i]['inchiSanitized'] = pool.starmap(MolToInchi_fun_safe, zip(list_df[i]['smilesSanitized'], list_df[i]['ROMolSanitizedLargestFragmentUncharged']))
            #list_df[i]['inchikeySanitized'] = pool.map(MolToIK_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            list_df[i]['inchikeySanitized'] = pool.starmap(MolToIK_fun_safe, zip(list_df[i]['smilesSanitized'], list_df[i]['ROMolSanitizedLargestFragmentUncharged']))
            list_df[i]['shortikSanitized'] = list_df[i]['inchikeySanitized'].str.split("-", n=1, expand=True)[0]
            list_df[i]['formulaSanitized'] = pool.map(MolToMF_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            list_df[i]['exactmassSanitized'] = pool.map(MolToEmass_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            list_df[i]['xlogpSanitized'] = pool.map(MolToLogP_fun, list_df[i]['ROMolSanitizedLargestFragmentUncharged'])
            list_df[i].drop('ROMolSanitizedLargestFragmentUncharged', axis=1, inplace=True)
            
            pool.close()
            pool.join()

            # list_df[i].to_csv(
            #     "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/2_cleaned/structure/oouthoupla_%i.csv" % i , 
            #     sep = '\t', 
            #     index = False,
            #     compression = 'gzip'
            #     )

# timer is stopped
print(" Above command executed in --- %s seconds ---" %
      (time.time() - start_time))

# we merge the previously obtained df
df = pd.concat(list_df)

# # dropping some irrelevant columns prior to export
# colstodrop = ['ROMol',
#               'ROMolSanitized',
#               'ROMolSanitizedLargestFragment',
#               'ROMolSanitizedLargestFragmentUncharged'
#               ]
# colstodrop = ['ROMolSanitizedLargestFragmentUncharged']
# df.drop(colstodrop, axis=1, inplace=True)

df.info()

# exporting df
import os
import errno
filename = ouput_file_path
if not os.path.exists(os.path.dirname(filename)):
    try:
        os.makedirs(os.path.dirname(filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
            
# df.to_csv(ouput_file_path, sep='\t', index=False)
df.to_csv(
	ouput_file_path, 
	sep = '\t', 
	index = False,
	compression = 'gzip'
	)
