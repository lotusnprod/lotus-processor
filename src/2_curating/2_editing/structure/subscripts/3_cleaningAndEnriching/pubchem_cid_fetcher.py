#import pubchempy as pcp
from pubchempy import Compound, get_compounds


## generic modules
import time
import multiprocessing
import pandas as pd
import sys
import gzip



# c = pcp.Compound.from_cid(5090)


# print(c.molecular_formula)
# print(c.molecular_weight)
# print(c.isomeric_smiles)
# print(c.xlogp)
# print(c.iupac_name)
# print(c.synonyms)

get_compounds('C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1', 'smiles')


input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables_min/3_curated/table.tsv.gz'

myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
	sep = '\t')

# eventually filter display some info, comment according to your needs
# df = df[df['originaldb'] == 'tcm']
df.columns
df.info()

df = df[~df['structureCleanedSmiles'].isnull()]
df.drop_duplicates('structureCleanedSmiles', inplace = True)
df = df.head(100)
df.info()



def SMILES_to_PCID_fun(smiles):
    m = get_compounds(smiles, 'smiles')
    if m:
        return m
    return None


# the full df is splitted and each subdf are treated sequentially as df > 900000 rows retruned errors 
# (parralel treatment of these subdf should improve performance)

# n = 50  #chunk row size
# list_df = [df[i:i+n] for i in range(0,df.shape[0],n)]

# timer is started
start_time = time.time()
cpus = 5


if __name__ == "__main__":
    # with multiprocessing.Pool(multiprocessing.cpu_count() - 2 ) as pool:
    with multiprocessing.Pool(int(cpus)) as pool:

        # # we generate ROMol object from smiles and or inchi
        df['PCID'] = pool.map(SMILES_to_PCID_fun, df['structureCleanedSmiles'])
        
        pool.close()
        pool.join()

# timer is stopped
print(" Above command executed in --- %s seconds ---" %
      (time.time() - start_time))



df['PCID'] = map(SMILES_to_PCID_fun, df['structureCleanedSmiles'])