#!/usr/bin/env python
# coding: utf-8

'''
Author: PMA
Contributor: JB, AR
'''
import errno
import numpy as np
import os
import sys
from chemosanitizer_functions import *
from tqdm import tqdm

# defining the command line arguments
try:
    input_file_path = sys.argv[1]
    ouput_file_path = sys.argv[2]
    smiles_column_header = sys.argv[3]
    cpus = int(sys.argv[4])

    print(f'Parsing tab separated file {input_file_path}'
          + f' with column: {smiles_column_header}'
          + f' as SMILES column. Paralelized on {cpus} cores.'
          + ' Proceeding to the validation, standardization, fragment choosing and uncharging of the ROMol object and returning the sanitized outputs in file :'
          + ouput_file_path)
except:
    print(
        '''Please add input and output file path as first and second argument, SMILES column header as third argument and finally the number of cpus you want to use.
        Example :
        python chemosanitizer.py ~/translatedStructureRdkit.tsv ./test.tsv structureTranslated 6''')
    sys.exit(1)

if __name__ == "__main__":
    myZip = gzip.open(input_file_path)
    df = pd.read_csv(
        myZip,
        sep='\t',
        encoding='utf-8',
        on_bad_lines='error'
    )

    if (len(df) == 1) and df.empty:
        df['structureTranslated'] = '[Pu]'
        print('your dataframe is empty, plutonium loaded')
    else:
        print('your dataframe is not empty :)')

    df = df[df[smiles_column_header].notnull()]
    chunk_size = len(df) // cpus + 1
    df_chunks = [df.iloc[i:i+chunk_size] for i in range(0, len(df), chunk_size)]
    f = CleaningFunc(smiles_column_header).f
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        processed_list = list(tqdm(executor.map(f, df_chunks), total=len(df_chunks)))

    processed_df = pd.concat(processed_list, ignore_index=True)

    output_path = os.path.dirname(ouput_file_path)

    if output_path != '' and not os.path.exists(output_path):
        try:
            os.makedirs(os.path.dirname(ouput_file_path))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    processed_df.to_csv(
        ouput_file_path,
        sep='\t',
        index=False,
        compression='gzip',
        encoding='utf-8'
    )
