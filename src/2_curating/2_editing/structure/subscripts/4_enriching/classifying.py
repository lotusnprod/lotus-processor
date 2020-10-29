# npClassifier (adapted from https://github.com/mwang87/NP-Classifier/blob/master/notebooks/BatchClassificationNotebook.ipynb

# importing packages
import pandas as pd
import grequests
import urllib.parse
from tqdm import tqdm
import sys
import gzip

# defining command line arguments
try:
    input_file_path = sys.argv[1]
    smiles_column_header = sys.argv[2]
    ouput_file_path = sys.argv[3]

    print('Parsing tab separated file'
          + input_file_path
          + 'with column: '
          + smiles_column_header
          + 'as smiles column.'
          + 'Proceeding to classification and returning the classified results in file :'
          + ouput_file_path)
except:
    print(
        '''Please add input file path as first argument, smiles column header as second and output file path as third.
        Example :
        python 2_curating/2_editing/structure/subscripts/4_enriching/classify.py ../data/interim/tables/2_cleaned/structure/cleaned.tsv.gz structureCleanedSmiles ../data/interim/tables/2_cleaned/structure/classified.tsv.gz''')

# importing file
myZip = gzip.open(input_file_path)

dfFull = pd.read_csv(
    myZip,
    sep='\t')

# keeping smiles column only
dfFull.rename(columns={smiles_column_header: 'smiles'}, inplace=True)
dfSmiles = dfFull[['smiles']]
df = dfSmiles.drop_duplicates()
df = df.reset_index(drop=True)

# server url
SERVER_URL = "https://npclassifier.ucsd.edu"

# getting classification
all_urls = []

for entry in tqdm(df.to_dict(orient="records")):
    smiles = str(entry["smiles"])
    # if len(smiles) > 5:
    request_url = "{}/classify?smiles={}".format(
        SERVER_URL, urllib.parse.quote(smiles))
    all_urls.append(request_url)

rs = (grequests.get(u) for u in all_urls)
responses = grequests.map(rs, size=20)
all_responses_list = [response.json() for response in responses]

# getting dataframe
classes = pd.json_normalize(all_responses_list)

# concatenating with initial file
finalDf = pd.concat([df, classes], axis=1)

# exporting
finalDf.to_csv(
    ouput_file_path,
    sep='\t',
    index=False,
    compression='gzip'
)
