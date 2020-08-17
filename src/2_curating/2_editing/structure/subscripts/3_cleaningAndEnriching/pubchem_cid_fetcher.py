#import pubchempy as pcp
from pubchempy import Compound, get_compounds


## generic modules
import time
import multiprocessing
import pandas as pd
import sys
import gzip

from functools import reduce



# c = pcp.Compound.from_cid(5090)


# print(c.molecular_formula)
# print(c.molecular_weight)
# print(c.isomeric_smiles)
# print(c.xlogp)
# print(c.iupac_name)
# print(c.synonyms)

get_compounds('C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1', 'smiles')


input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables_min/3_curated/table.tsv.gz'
input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/structureMetadata.tsv.gz'

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
#df = df.head(1000)
df.info()




def SMILES_to_PCID_fun(smiles):
    m = get_compounds(smiles, 'smiles')
    if m:
        return m
    return None


# the full df is splitted and each subdf are treated sequentially as df > 900000 rows retruned errors 
# (parralel treatment of these subdf should improve performance)


# timer is started
start_time = time.time()
cpus = 5


n = 500000  #chunk row size
list_df = [df[i:i+n] for i in range(0,df.shape[0],n)]

# timer is started
start_time = time.time()


for i in range(0, len(list_df)):

    if __name__ == "__main__":
        # with multiprocessing.Pool(multiprocessing.cpu_count() - 2 ) as pool:
        with multiprocessing.Pool(int(cpus)) as pool:

            # # we generate ROMol object from smiles and or inchi
            list_df[i]['PCID'] = pool.map(SMILES_to_PCID_fun, list_df[i]['structureCleanedSmiles'])
            
            pool.close()
            pool.join()

# timer is stopped
print(" Above command executed in --- %s seconds ---" %
    (time.time() - start_time))


df = pd.concat(list_df)

df['PCID']


df['PCID'] = map(SMILES_to_PCID_fun, df['structureCleanedSmiles'])



# lets try to work on chunks 

myZip = gzip.open(input_file_path)

chunks = pd.read_csv(
    myZip, chunksize=1000, usecols=[
        "structureCleanedSmiles",
        "organismCleaned"
    ], sep = "\t"
)

# 2. Map. For each chunk, calculate the per-street counts:
def get_counts(chunk):
    by_party = chunk.groupby("structureCleanedSmiles")
    street = by_party["organismCleaned"]
    return street.value_counts()

processed_chunks = map(get_counts, chunks)


def get_PCID(chunk):
    by_party = chunk.groupby("structureCleanedSmiles")
    street = by_party["organismCleaned"]
    return street.value_counts()
    
processed_chunks = map(get_counts, chunks)

# 3. Reduce. Combine the per-chunk voter counts:
def add(previous_result, new_result):
    return previous_result.add(new_result, fill_value=0)
result = reduce(add, processed_chunks)

# 4. Post-process.
result.sort_values(ascending=False, inplace=True)

print(result)


### preparing input for https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

df = df['structureCleanedSmiles']


ouput_file_path0 = "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/smiles_0.gz"
ouput_file_path1 = "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/smiles_1.gz"

list_df[0].to_csv(
	ouput_file_path0, 
	sep = '\t', 
	index = False,
    header = False,
	compression = 'gzip'
	)

list_df[1].to_csv(
	ouput_file_path1, 
	sep = '\t', 
	index = False,
    header = False,
	compression = 'gzip'
	)



input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/structureMetadata.tsv.gz'

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
df = df['structureCleanedSmiles']
#df = df.head(1000)
df.info()

n = 100000  #chunk row size
list_df = [df[i:i+n] for i in range(0,df.shape[0],n)]


for i in range(0, len(list_df)):

    list_df[i].to_csv(
	"/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/smiles_%s.txt" % i, 
	sep = '\t', 
	index = False,
    header = False
	)


## checking the ratio of pcided smiles 

input_file_path = '/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/pcided_smiles.txt.gz'

myZip = gzip.open(input_file_path)

df = pd.read_csv(
	myZip,
    header = None,
	sep = '\t')


df.info()
pd.set_option('display.width', 1000)

df[0][793409]