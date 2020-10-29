# %%
# ## generic modules
import re
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from functools import reduce
import gzip


# loading the input dataframes

# %%
opennpdb_table_path = gzip.open(
    "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/interim/tables/3_curated/table.tsv.gz"
)

wikidata_inchilist_path = (
    "/home/EPGL.UNIGE.LOCAL/allardp/opennaturalproductsdb/data/external/query.tsv"
)


# %%
df_onpdb_table = pd.read_csv(
    opennpdb_table_path, sep="\t", error_bad_lines=False, low_memory=False
)
df_wd_inchi = pd.read_csv(
    wikidata_inchilist_path, sep="\t", error_bad_lines=False, low_memory=False
)
# df_GNPS_output = pd.read_csv(GNPS_output_path,
#                              sep='\t', error_bad_lines=False, usecols=['cluster index', 'componentindex'],
#                              low_memory=False)

# %%
# We start by adding a sik column to the wd table

df_wd_inchi["shortik"] = df_wd_inchi["inchikey"].str.split(
    "-", n=1, expand=True)[0]
df_wd_inchi.info()

# %%
# Here we drop duplicated inchi in the wd table

df_wd_inchi.drop_duplicates("inchikey", inplace=True)
df_wd_inchi.info()

# %%
# Here we drop duplicated inchi in the openpdb table

df_onpdb_table.drop_duplicates("inchikeySanitized", inplace=True)
df_onpdb_table.info()

# %%
# Now we will append the wd compound field to the opnnpdb table after merging on the ik column

df_onpdb_table_wded = pd.merge(
    df_onpdb_table,
    df_wd_inchi,
    left_on="inchikeySanitized",
    right_on="inchikey",
    how="left",
)
df_onpdb_table_wded.info()

# %%


df = px.data.gapminder().query("year == 2007").query("continent == 'Europe'")
df.loc[
    df["pop"] < 2.0e6, "country"
] = "Other countries"  # Represent only large countries
fig = px.pie(
    df, values="pop", names="country", title="Population of European continent"
)
fig.show()

# %%
